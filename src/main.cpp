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

using namespace libQBD;
using namespace Eigen;
using namespace std;

class Cluster
{
    private:
    QBD<double> process;

    public:
    StationaryDistribution<double>  sd;
    vector<vector<state>> states;

    Cluster(double lambda, unsigned int c, const server_dist &dist, const vector<double> &f, const MatrixXd &P_a, const MatrixXd &P_d)
    {
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
        //*/
    }
};

int main()
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
        //dist.mu[k] = (double)(k + 1)*pow(10,k);
        dist.mu[k] = (double)(k + 1);
    }
    //dist.mu[0] = 1;
    //dist.mu[1] = 2;
    //dist.mu[2] = 100;

    //cout << (string)dist << endl;
    //prepare_dist(dist,c);
    //create_forward_matrix(lambda, P_a, dist, 5, c);
    //create_backward_matrix(P_d, f, prepare_dist(dist,c), 1, c);
    Cluster *model = new Cluster(lambda, c, dist, f, P_a, P_d);
    cout.precision(22);
    cout << model->sd.get_mean_clients() << endl;
    delete model;
    return 0;
}