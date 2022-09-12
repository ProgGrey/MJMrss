/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */

#include "dist.h"

#include <cstring>
#include <sstream>

using namespace std;

server_dist::server_dist(uint16_t length)
{
    this->length = length;
    this->serv_count = new uint16_t[length];
    this->prob = new double[length];
    this->mu = new double[length];
}

server_dist::~server_dist()
{
    delete [] this->mu;
    delete [] this->prob;
    delete [] this->serv_count;
    this->mu = NULL;
    this->prob = NULL;
    this->serv_count = NULL;
}

server_dist::server_dist(const server_dist &parent)
{
    this->mu = new double[parent.length];
    this->prob = new double[parent.length];
    this->serv_count = new uint16_t[parent.length];
    this->length = parent.length;
    memcpy(this->mu, parent.mu, sizeof(double)*parent.length);
    memcpy(this->prob, parent.prob, sizeof(double)*parent.length);
    memcpy(this->serv_count, parent.serv_count, sizeof(uint16_t)*parent.length);
}

server_dist& server_dist::operator= (const server_dist &right)
{
    if(&right == this){
        return *this;
    }
    
    if(this->prob == NULL){
        this->mu = new double[right.length];
        this->prob = new double[right.length];
        this->serv_count = new uint16_t[right.length];
    } else if(this->length != right.length){
        delete [] this->mu;
        delete [] this->prob;
        delete [] this->serv_count;
        this->mu = new double[right.length];
        this->prob = new double[right.length];
        this->serv_count = new uint16_t[right.length];
    }
    this->length = right.length;
    memcpy(this->mu, right.mu, sizeof(double)*right.length);
    memcpy(this->prob, right.prob, sizeof(double)*right.length);
    memcpy(this->serv_count, right.serv_count, sizeof(uint16_t)*right.length);

    return *this;
}

server_dist server_dist::cut_prob(uint16_t free_servers)
{
    uint16_t len = 0;
    double possible = 0;
    for(unsigned int k = 0; k < this->length; k++){
        if(serv_count[k] > free_servers){
            len++;
            possible += this->prob[k];
        }
    }
    server_dist cutted(len);
    unsigned int p = 0;
    for(unsigned int k = 0; k < this->length; k++){
        if(this->serv_count[k] > free_servers){
            cutted.serv_count[p] = this->serv_count[k];
            cutted.prob[p] = this->prob[k] / possible;
            cutted.mu[p] = this->mu[k];
            p++;
        }
    }
    return cutted;
}

double server_dist::get_mu(uint16_t serv) const
{
    for(unsigned int k = 0; k < length; k++){
        if(serv == serv_count[k]){
            return this->mu[k];
        }
    }
    return 0;
}

server_dist::operator std::string() const
{
    stringstream ret;
    for(unsigned int k = 0; k < this->length; k++){
        ret << this->serv_count[k] << '\t' << this->prob[k] << '\t' << this->mu[k] << endl;
    }
    return ret.str();
}