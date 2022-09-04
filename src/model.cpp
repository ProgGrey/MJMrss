/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#include "libQBD/inc/libQBD.hpp"
#include <iostream>
#include <stdint.h>

#include "main.hpp"

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

using Eigen::VectorXd;

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
    
    vector<double> get_mean_clients(double max_time, const Rcpp::List &pi_0)
    {
        return td.get_mean_clients(max_time, cast_pi0(pi_0));
    }

    vector<double> get_mean_queue(double max_time, const Rcpp::List &pi_0)
    {
        return td.get_mean_queue(*queue_size_vector, max_time, cast_pi0(pi_0));
    }

    Rcpp::List get_distribution(double max_time, const Rcpp::List &pi_0)
    {
        Rcpp::List ret;
        vector<vector<VectorX<double>>> tmp = td.get_dist(max_time, cast_pi0(pi_0));
        for(auto it = tmp.begin(); it != tmp.end(); it++){
            Rcpp::List dft;
            for(auto itt = it->begin(); itt != it->end(); it++){
                dft.push_back(*it);
            }
            ret.push_back(dft);
        }
        return ret;
    }
};

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
    double mean_queue;
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
    Model()
    {
        unsigned int c = 10;
        double lambda = 0.1;
        vector<double> f;
        f.push_back(1.0);
        f.push_back(2.2);

        MatrixXd P_a{{0.8, 0.2},
                    {0.0, 1.0}};
        MatrixXd P_d{{1.0, 0.0},
                    {0.2, 0.8}};

        server_dist dist(c);
        for(unsigned int k = 0; k < c; k++){
            dist.prob[k] = 1.0/c;
            dist.serv_count[k] = k + 1;
            dist.mu[k] = (double)(k + 1);
        }
        this->init(lambda, c, dist, f, P_a, P_d);
    }

    Model(double lambda, unsigned int c, const NumericMatrix &classes, vector<double> f, MatrixXd P_a, MatrixXd P_d)
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

    .property("h", &ModelTransient::h)

    .method("mean_clients", &ModelTransient::get_mean_clients)
    .method("mean_queue", &ModelTransient::get_mean_queue)
    .method("distribution", &ModelTransient::get_distribution)
    ;

    class_<Model>( "Model" )

    .constructor()
    .constructor<double, unsigned int, const NumericMatrix &, vector<double>, MatrixXd, MatrixXd>()

    .property("rho", &Model::get_rho)
    .property("mean_clients", &Model::get_mean_clients)
    .property("mean_queue", &Model::get_mean_queue)
    .property("sum_pi_from_c_to_inf", &Model::dist_sum_from_c_to_inf)
    .property("pi_0_c", &Model::get_pi_0_c)

    .method("distribution", &Model::get_distribution)
    .method("level_description", &Model::level_description)
    .method("level_busy_servers", &Model::level_busy_servers)
    .method("level_serviced_clients", &Model::level_serviced_clients)
    .method("transient_analysis", &Model::transient_analysis)
    ;
}