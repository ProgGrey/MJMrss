/*
 * Copyright (c) 2022-2023 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#pragma once
#ifndef __MATRICES_HPP__
#define __MATRICES_HPP__

#include <vector>
#include <cstdint>
#include "eigen_wrap.h"

#include "dist.h"
#include "state.h"

std::vector<server_dist> prepare_dist(server_dist dist, uint16_t c);
Eigen::MatrixXd create_forward_matrix(const std::vector<state> &from, const std::vector<state> &to, double lambda, const Eigen::MatrixXd &P_a, const server_dist &dist, uint16_t level);
Eigen::MatrixXd create_backward_matrix(const std::vector<state> &from, const std::vector<state> &to, const Eigen::MatrixXd &P_d, const std::vector<double> &velocity, const std::vector<server_dist> &dists, uint16_t level, uint16_t c);


#endif