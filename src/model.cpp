/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
//#include <Eigen/Dense>

#include <iostream>
#include <stdint.h>

#include "eigen_wrap.h"

#include "libQBD/inc/libQBD.hpp"

#include "main.h"

#include <Rcpp.h>

#include <RcppEigen.h>

#if defined(_OPENMP)
//#include <omp.h>
#include <thread>
#include <cpuid.h>
#endif

using namespace libQBD;
using namespace Eigen;
using namespace std;
using namespace Rcpp;

//using Eigen::VectorXd;

//' @name Rcpp_ModelTransient
//' @aliases ModelTransient
//' @aliases Rcpp_ModelTransient-class
//' @rdname ModelTransient-class
//' @title ModelTransient Class
//' @description Internal representation for multi-server job model with random speed scaling for transient analysis.
//' @slot h - Returns current step size.
//' @slot mean_clients.
//' \itemize{
//' \item Parameter: max_time - maximum time for which you want to calculate mean customers in system;
//' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
//' \item Returns: mean clients in system.
//' }
//' @slot mean_queue
//' \itemize{
//' \item Parameter: max_time - maximum time for which you want to calculate mean queue length;
//' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
//' \item Returns: mean queue length
//' }
//' @slot distribution 
//' \itemize{
//' \item Parameter: max_time - maximum level for which you want to calculate distribution;
//' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
//' \item Returns: list of lists of vectors containing transient distribution.
//' }
class ModelTransient
{
    private:
    TaylorSeriesTransient<double> td;
    
    std::vector<Eigen::VectorX<double>> *queue_size_vector;
    vector<VectorX<double>> pi_0;

    vector<VectorX<double>> &cast_pi0(const Rcpp::List &pi_0)
    {
        this->pi_0.clear();
        for(auto it = pi_0.begin(); it != pi_0.end(); it++){
            this->pi_0.push_back(*it);
        }
        return this->pi_0;
    }

    public:
    ModelTransient(const QBD<double> &proc, uint8_t order, double step, std::vector<Eigen::VectorX<double>> *q)
    {
        this->queue_size_vector = q;
        td.bind(proc, order, step);
    }

    double h()
    {
        return td.get_step();
    }
    
    vector<double> get_mean_clients(double max_time, Rcpp::List pi_0)
    {
        return td.get_mean_clients(max_time, cast_pi0(pi_0));
    }

    vector<double> get_mean_queue(double max_time, Rcpp::List pi_0)
    {
        return td.get_mean_queue(*queue_size_vector, max_time, cast_pi0(pi_0));
    }

    Rcpp::List get_distribution(double max_time, Rcpp::List pi_0)
    {
        Rcpp::List ret;
        vector<vector<VectorX<double>>> tmp = td.get_dist(max_time, cast_pi0(pi_0));
        for(auto it = tmp.begin(); it != tmp.end(); it++){
            Rcpp::List dft;
            for(auto itt = it->begin(); itt != it->end(); itt++){
                dft.push_back(*itt);
            }
            ret.push_back(dft);
        }
        return ret;
    }
};

//' @name Rcpp_Model
//' @aliases Model
//' @aliases Rcpp_Model-class
//' @rdname Model-class
//' @title Model Class
//' @description Internal representation for multi-server job model with random speed scaling.
//' @slot rho Returns rho computated by using Neuts ergodicity criteria.
//' @slot mean_clients Returns mean clients in system.
//' @slot mean_queue Returns mean queue length.
//' @slot sum_pi_from_c_to_inf Returns sum of distribution vectors from level N to infinity level.
//' @slot pi_0_c Returns distribution for first N levels.
//' @slot distribution 
//' \itemize{
//' \item Parameter: max_level - maximum level for which you want to calculate distribution;
//' \item Returns: list of vectors containing stationary distribution.
//' }
//' @slot level_description 
//' \itemize{
//' \item Parameter: level -  level for which you want to get a description;
//' \item Returns: human-readable description of level in form "m|s_1,s_2,..." where "m" is a speed mode and "s_i" is
//' a number of serviced customers of "i" class.
//' }
//' @slot level_busy_servers 
//' \itemize{
//' \item Parameter: level -  level for which you want to get a description;
//' \item Returns: vector of numbers of busy servers for each state at specified level.
//' }
//' @slot level_serviced_clients 
//' \itemize{
//' \item Parameter: level -  level for which you want to get a description;
//' \item Returns: vector of numbers of serviced servers for each state at specified level.
//' }
//' @slot transient_analysis 
//' \itemize{
//' \item Parameter: order - numerical order of Taylor series method
//' \item Parameter: stepsize - If positive, then  step length of the algorithm. If negative - the multiplier of the 
//' maximum step of the algorithm.
//' \item Returns: internal representation for model transient analysis (ModelTransient class).
//' }
class Model
{
    private:
    QBD<double> process;

    unsigned int c;
    StationaryDistribution<double>  sd;
    vector<vector<state>> states;

    void init(double lambda, unsigned int c, const server_dist &dist, const vector<double> &f, const MatrixXd &P_a, const MatrixXd &P_d)
    {
        // Initialize OpenMP
        #ifdef _OPENMP
        #if defined(__x86_64__) || defined(__i686__)
        unsigned int eax, ebx, ecx, edx;
        __get_cpuid(1, &eax, &ebx, &ecx, &edx);
        bool hyper_th =  (edx & (1 << 28)) > 0;
        #else
        bool hyper_th = false;
        #endif
        unsigned int cores_log = thread::hardware_concurrency();
        unsigned int cores_ph = hyper_th ? cores_log >> 1 : cores_log;
        Eigen::setNbThreads(cores_ph);
        #endif
        // Initilize model
        this->c = c;
        // Prepare states descriptions:
        for(unsigned int k = 0; k <= (c + 1); k++){
            states.push_back(generate_all_states(P_d.rows(), c, k, dist));
        }
        auto cutted_dists = prepare_dist(dist,c);
        process.add_A_plus(create_forward_matrix(states[0], states[1], lambda, P_a, dist, 0));
        for(unsigned int k = 1; k <= c; k++){
            process.add_A_plus(create_forward_matrix(states[k], states[k+1], lambda, P_a, dist, k));
            process.add_A_minus(create_backward_matrix(states[k], states[k-1], P_d, f, cutted_dists, k, c));
        }
        process.add_A_minus(create_backward_matrix(states[c + 1], states[c], P_d, f, cutted_dists, c + 1, c));
        process.auto_A_0();
        sd.bind(process);
    }

    bool is_queue_size_vector_comp = false;
    std::vector<Eigen::VectorX<double>> queue_size_vector;
    bool is_mean_queue_comp = false;
    double mean_queue = 0;
    void computate_queue_size_vec(void)
    {
        if(is_queue_size_vector_comp != true){
            for(unsigned int k = 0; k <= c; k++){
                Eigen::VectorX<double> v = Eigen::VectorX<double>::Constant(states[k].size(), 1, 0);
                for(auto it = states[k].begin(); it != states[k].end(); it++)
                {
                    v(it-states[k].begin(), 0) = k - it->apps();
                }
                queue_size_vector.push_back(v);
            }
            is_queue_size_vector_comp = true;
        }
    }

    bool is_pi_0_c_comp = false;
    Rcpp::List pi_0_c;

    public:

    Model(double lambda, unsigned int c, NumericMatrix classes, vector<double> f, MatrixXd P_a, MatrixXd P_d)
    {
        unsigned int len = 0;
        for(int k = 0; k < classes.cols(); k++){
            if(classes(1,k) != 0.0){
                len++;
            }
        }
        server_dist dist(len);
        unsigned int pos = 0;
        for(int k = 0; k < classes.cols(); k++){
            if(classes(1,k) != 0.0){
                dist.mu[pos] = classes(2,k);
                dist.prob[pos] = classes(1,k);
                dist.serv_count[pos] = classes(0,k);
                pos++;
            }
        }
        this->init(lambda, c, dist, f, P_a, P_d);
    }

    double get_rho()
    {
        return sd.get_rho();
    }

    double get_mean_clients()
    {
        return sd.get_mean_clients();
    }

    VectorXd dist_sum_from_c_to_inf()
    {
        return sd.get_sum_from_c_to_inf();
    }

    Rcpp::List get_pi_0_c()
    {
        if(!is_pi_0_c_comp){
            auto dist = sd.get_pi_0_c();
            for(auto it = dist.begin(); it != dist.end(); it++){
                pi_0_c.push_back(*it);
            }
            is_pi_0_c_comp = true;
        }
        return pi_0_c;
    }

    double get_mean_queue()
    {
        if(!is_mean_queue_comp){
            computate_queue_size_vec();
            mean_queue = sd.get_mean_queue(queue_size_vector);
            is_mean_queue_comp = true;
        }
        return mean_queue;
    }

    Rcpp::List get_distribution(unsigned int max_level)
    {
        Rcpp::List ret;
        auto dist = sd.get_dist(max_level);
        for(auto it = dist.begin(); it != dist.end(); it++){
            ret.push_back(*it);
        }
        return ret;
    }

    Rcpp::NumericVector level_busy_servers(unsigned int level)
    {
        Rcpp::NumericVector ret;
        if(level >= states.size()){
            //ret.reserve(states.back().size());
            for(auto it = states.back().begin(); it != states.back().end(); it++){
                ret.push_back(it->busy());
            }
        } else{
            //ret.reserve(states[level].size());
            for(auto it = states[level].begin(); it != states[level].end(); it++){
                ret.push_back(it->busy());
            }
        }
        return ret;
    }
    
    Rcpp::NumericVector level_serviced_clients(unsigned int level)
    {
        Rcpp::NumericVector ret;
        if(level >= states.size()){
            //ret.reserve(states.back().size());
            for(auto it = states.back().begin(); it != states.back().end(); it++){
                ret.push_back(it->apps());
            }
        } else{
            //ret.reserve(states[level].size());
            for(auto it = states[level].begin(); it != states[level].end(); it++){
                ret.push_back(it->apps());
            }
        }
        return ret;
    }

    Rcpp::StringVector level_description(unsigned int level)
    {
        Rcpp::StringVector ret;
        if(level >= states.size()){
            //ret.reserve(states.back().size());
            for(auto it = states.back().begin(); it != states.back().end(); it++){
                ret.push_back((string)(*it));
            }
        } else{
            //ret.reserve(states[level].size());
            for(auto it = states[level].begin(); it != states[level].end(); it++){
                ret.push_back((string)(*it));
            }
        }
        return ret;
    }

    ModelTransient transient_analysis(uint8_t order, double step)
    {
        return ModelTransient(this->process, order, step, &queue_size_vector);
    }
};

RCPP_EXPOSED_CLASS(ModelTransient)
RCPP_EXPOSED_CLASS(Model)

RCPP_MODULE(master){
    using namespace Rcpp ;

    class_<ModelTransient>("ModelTransient")

    .property("h", &ModelTransient::h, "Returns current step size.")

    .method("mean_clients", &ModelTransient::get_mean_clients, "Returns vector of mean customers in system.")
    .method("mean_queue", &ModelTransient::get_mean_queue, "Returns vector of mean queue length.")
    .method("distribution", &ModelTransient::get_distribution, "Returns list of transient distributions.")
    ;

    class_<Model>( "Model" )

    .constructor<double, unsigned int, NumericMatrix, vector<double>, MatrixXd, MatrixXd>()

    .property("rho", &Model::get_rho, "Returns rho computated by using Neuts ergodicity criteria.")
    .property("mean_clients", &Model::get_mean_clients, "Returns mean clients in system.")
    .property("mean_queue", &Model::get_mean_queue, "Returns mean queue length.")
    .property("sum_pi_from_c_to_inf", &Model::dist_sum_from_c_to_inf, "Returns sum of distribution vectors from level N to infinity level.")
    .property("pi_0_c", &Model::get_pi_0_c, "Returns distribution for first N levels.")

    .method("distribution", &Model::get_distribution, "Computate stationary distribution.")
    .method("level_description", &Model::level_description, "Returns human-readable level description.")
    .method("level_busy_servers", &Model::level_busy_servers, "Returns busy servers for certain level.")
    .method("level_serviced_clients", &Model::level_serviced_clients, "Returns number of serviced customers for certain level.")
    .method("transient_analysis", &Model::transient_analysis, "Prepare for transient analysis")
    ;
}