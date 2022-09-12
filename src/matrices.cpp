/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#include "matrices.h"
#include "state.h"

#include <iostream>

using namespace Eigen;
using namespace std;


vector<server_dist> prepare_dist(server_dist dist, uint16_t c)
{
    vector<server_dist> ret;
    for(unsigned int k = 0; k <= c; k++){
        ret.push_back(dist.cut_prob(k));
    }
    return ret;
}

inline unsigned int find_index(const state &what, const vector<state> &where)
{
    auto it = lower_bound(where.begin(), where.end(), what);
    return (unsigned int)(it - where.begin());
}

MatrixXd create_forward_matrix(const vector<state> &from, const vector<state> &to, double lambda, const MatrixXd &P_a, const server_dist &dist, uint16_t level)
{
    uint16_t mode_count = P_a.rows();
    MatrixXd A_plus = MatrixXd::Constant(from.size(), to.size(), 0.0);
    for(auto it = from.begin(); it != from.end(); it++){
        unsigned int fr = it - from.begin();
        // All posible arrivals of apps
        for(unsigned int k = 0; k < dist.length; k++){
            state what = *it;
            if(what.apps() == level){
                // Addition of new app is posible only if app in queue is not blocked this
                what.add_app(dist.serv_count[k]);
            }
            // All posible mode switching
            for(unsigned int j = 0; j < mode_count; j++){
                if(P_a(it->m, j) != 0.0){
                    what.m = j;
                    A_plus(fr, find_index(what, to)) += P_a(it->m, j) * dist.prob[k] * lambda;
                }
            }
        }
    }
    return A_plus;
}

void review_all_arrivals(const state &what, const server_dist &state_dist, const server_dist &orig_dist, double rate, uint16_t level, uint16_t c, unsigned int fr, const vector<state> &to, MatrixXd &A_minus)
{
    if((what.busy() == c) || (what.apps() == level)){
        A_minus(fr, find_index(what, to)) += rate;
    } else if(what.apps() < level){
        for(unsigned int k = 0; k < state_dist.length; k++){
            state copy(what);
            copy.s[state_dist.serv_count[k] - 1]++;
            if(copy.busy() <= c){
                review_all_arrivals(copy, orig_dist, orig_dist, rate * state_dist.prob[k], level, c, fr, to, A_minus);
            }else{
                A_minus(fr, find_index(what, to)) += rate * state_dist.prob[k];
            }
            //*/
        }
    }
}

MatrixXd create_backward_matrix(const vector<state> &from, const vector<state> &to, const MatrixXd &P_d, const vector<double> &velocity, const vector<server_dist> &dists, uint16_t level, uint16_t c)
{
    uint16_t mode_count = P_d.rows();
    MatrixXd A_minus = MatrixXd::Constant(from.size(), to.size(), 0.0);
    // For all states
    for(auto it = from.begin(); it != from.end(); it++){
        unsigned int fr = it - from.begin();
        server_dist state_dist = dists[c - it->busy()];
        // for all switching
        for(unsigned int k = 0; k < mode_count; k++){
            if(P_d(it->m, k) != 0){
                // For all apps, that can exit
                for(unsigned int j = 0; j < it->s_len; j++){
                    if(it->s[j] != 0){
                        state what = *it;
                        what.s[j]--;
                        what.m = k;
                        double rate = P_d(it->m, k) * it->s[j] * dists.front().get_mu(j+1) * velocity[it->m];
                        review_all_arrivals(what, state_dist, dists.front(), rate, level - 1, c, fr, to, A_minus);
                    }
                }
            }
        }
    }
    return A_minus;
}