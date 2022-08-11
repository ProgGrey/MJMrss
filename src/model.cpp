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

// [[Rcpp::depends(RcppEigen)]]

/*
useDynLib(MJMrss, .registration=TRUE)
import(methods, Rcpp, RcppEigen)
exportPattern("^[[:alpha:]]+")
//*/
//////' @exportPattern "^[[:alpha:]]+"

using namespace libQBD;
using namespace Eigen;
using namespace std;
using namespace Rcpp;

using Eigen::VectorXd;

class Model
{
    private:
    QBD<double> process;

    unsigned int c;
    StationaryDistribution<double>  sd;
    vector<vector<state>> states;

    void init(double lambda, unsigned int c, const server_dist &dist, const vector<double> &f, const MatrixXd &P_a, const MatrixXd &P_d)
    {
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

    bool is_mean_queue_comp = false;
    double mean_queue;

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

    Model(double lambda, unsigned int c, const NumericVector &mu, const NumericVector &clients_server_distribution, vector<double> f, MatrixXd P_a, MatrixXd P_d)
    {
        unsigned int len = 0;
        for(auto it = clients_server_distribution.begin(); it != clients_server_distribution.end(); it++){
            if(*it != 0.0){
                len++;
            }
        }
        server_dist dist(len);
        unsigned int pos = 0;
        for(unsigned int k = 0; k < clients_server_distribution.size(); k++){
            if(clients_server_distribution[k] != 0.0){
                dist.mu[pos] = mu[k];
                dist.prob[pos] = clients_server_distribution[k];
                dist.serv_count[pos] = k + 1;
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
            std::vector<Eigen::VectorX<double>> queue_size_vector;
            for(unsigned int k = 0; k <= c; k++){
                Eigen::VectorX<double> v = Eigen::VectorX<double>::Constant(states[k].size(), 1, 0);
                for(auto it = states[k].begin(); it != states[k].end(); it++)
                {
                    v(it-states[k].begin(), 0) = k - it->apps();
                }
                queue_size_vector.push_back(v);
            }
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
};


RCPP_MODULE(master){
    using namespace Rcpp ;

    class_<Model>( "Model" )

    .constructor()
    .constructor<double, unsigned int, const NumericVector&, const NumericVector&, vector<double>, MatrixXd, MatrixXd>()

    .property("rho", &Model::get_rho)
    .property("mean_clients", &Model::get_mean_clients)
    .property("mean_queue", &Model::get_mean_queue)
    .property("sum_pi_from_c_to_inf", &Model::dist_sum_from_c_to_inf)
    .property("pi_0_c", &Model::get_pi_0_c)

    .method("distribution", &Model::get_distribution)
    .method("level_description", &Model::level_description)
    .method("level_busy_servers", &Model::level_busy_servers)
    ;
}