/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#include "state.hpp"
#include <sstream>
#include <algorithm>

#include <iostream>

using namespace std;

state::state(uint16_t len)
{
    this->m = 0;
    this->s_len = len;
    this->s = new uint16_t[len];
    memset(this->s, 0, sizeof(uint16_t)*len);
}

state::state(const state &parent)
{
    this->m = parent.m;
    this->s_len = parent.s_len;
    this->s = new uint16_t[parent.s_len];
    memcpy(this->s, parent.s, sizeof(uint16_t)*parent.s_len);
}

state::~state()
{
    delete [] this->s;
    this->s = NULL;
}


uint32_t state::busy() const
{
    uint32_t ret = 0;
    for(uint_fast16_t k = 0; k < s_len; k++){
        ret += s[k]*(k+1);
    }
    return ret;
}


uint32_t state::apps() const
{
    uint32_t ret = 0;
    for(uint_fast16_t k = 0; k < s_len; k++){
        ret += s[k];
    }
    return ret;
}

void state::add_app(uint16_t serv)
{
    if((this->s_len - this->busy()) >= serv){
        this->s[serv - 1]++;
    }
}

state& state::operator= (const state &right)
{
    if (this == &right)
        return *this;
    
    this->m = right.m;
    if(this->s == NULL){
        this->s = new uint16_t[right.s_len];
        this->s_len = right.s_len;
    } else if(this->s_len != right.s_len){
        delete [] this->s;
        this->s = new uint16_t[right.s_len];
        this->s_len = right.s_len;
    }
    memcpy(this->s, right.s, sizeof(uint16_t)*right.s_len);

    return *this;
}

bool operator< (const state &a, const state &b)
{
    if(a.s_len != b.s_len){
        throw "Length of states must be equal for comparison.";
    }

    if(a.m < b.m){
        return true;
    } else if(a.m == b.m) {
        for(unsigned int k = 0; k < a.s_len; k++){
            if(a.s[k] < b.s[k]){
                return true;
            } else if (a.s[k] > b.s[k]){
                return false;
            }
        }
    }
    return false;
}

bool operator== (const state &a, const state &b)
{
    if(a.s_len != b.s_len){
        throw "Length of states must be equal for comparison.";
    }
    if(a.m == b.m){
        for(unsigned int k = 0; k < a.s_len; k++){
            if(a.s[k] != b.s[k]){
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}

state::operator std::string() const
{
    stringstream ret;
    ret << (unsigned int)m << '|';
    if(this->s_len > 0){
        for(unsigned int k = 0; k < (this->s_len - 1); k++){
            ret << (unsigned int)s[k] << ',';
        }
        ret << (unsigned int)s[this->s_len - 1];
    } else{
        ret << (unsigned int)s[this->s_len];
    }
    return ret.str();
}

void gen_combs(state prev, const server_dist &dist, uint16_t level, uint16_t max_servers, uint16_t last,  vector<state> &ret)
{
    if((prev.busy() == max_servers) || (prev.apps() == level)){
        ret.push_back(prev);
        return;
    } else{
        // Starting with last to remove repetitions:
        for(uint_fast16_t k = last; k < dist.length; k++){
            state tmp(prev);
            tmp.s[dist.serv_count[k]-1]++;
            // Number of servers increased when k increased.
            // If an application of this size does not fit, 
            // then an application of a larger size will also not fit.
            if(tmp.busy() > max_servers){
                // Removing impossible states when distribution is truncated on the right.
                if (dist.max_serv() > dist.serv_count[last]){
                    ret.push_back(prev);
                }
                break;
            }
            gen_combs(tmp, dist, level, max_servers, k, ret);
        }
        //*/
    }
    
}

vector<state> generate_all_states(uint8_t mode_count, uint16_t max_servers, uint16_t level, const server_dist &dist)
{
    vector<state> ret;
    state init(max_servers);
    gen_combs(init, dist, level, max_servers, 0, ret);
    sort(ret.begin(), ret.end());
    unsigned int r = ret.size();
    for(unsigned int j = 1; j < mode_count; j++){
        for(unsigned int k = 0; k < r; k++){
            auto tmp = ret[k];
            tmp.m = j;
            ret.push_back(tmp);
        }   
    }
    return ret;
}