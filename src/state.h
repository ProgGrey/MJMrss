/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#pragma once
#ifndef __STATE_HPP__
#define __STATE_HPP__

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "dist.h"

class state
{
    public:
    uint8_t m;
    uint16_t *s = NULL;
    uint16_t s_len = 0;
    
    explicit state(uint16_t len);
    state(const state &parent);
    ~state();

    uint32_t busy() const;
    uint32_t apps() const;

    void add_app(uint16_t serv);

    state& operator= (const state &right);

    operator std::string() const;
};

bool operator< (const state &a, const state &b);
inline bool operator>(const state &a, const state &b)
{
    return b < a;
}
inline bool operator<=(const state &a, const state &b)
{
    return !(b < a);
}
inline bool operator>=(const state &a, const state &b)
{
    return !(a < b);
}


bool operator== (const state &a, const state &b);
inline bool operator!= (const state &a, const state &b)
{
    return !(a == b);
}


std::vector<state> generate_all_states(uint8_t mode_count, uint16_t max_servers, uint16_t level, const server_dist &dist);

#endif