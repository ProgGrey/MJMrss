/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#pragma once
#ifndef __DIST_HPP__
#define __DIST_HPP__

#include <cstdint>
#include <cstdlib>
#include <string>

class server_dist
{
    public:
    uint16_t* serv_count = NULL;
    double* prob = NULL;
    double* mu = NULL;
    uint16_t length = 0;

    explicit server_dist(uint16_t length);
    server_dist(const server_dist &right);
    ~server_dist();

    server_dist cut_prob(uint16_t free_servers);
    double get_mu(uint16_t serv) const;

    inline uint16_t max_serv() const
    {
        return serv_count[length - 1];
    }

    server_dist& operator= (const server_dist &right);

    operator std::string() const;
};
#endif