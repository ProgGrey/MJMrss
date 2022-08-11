
 # Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 # All rights reserved.
 #
 # This source code is licensed under the BSD-3 clause license.

#' @export build_model

#' @importFrom methods new

#' @title build_model
#' @description Built multiserver job model with random speed scaling
#' @param lambda Input rate.
#' @param N Number of servers.
#' @param mu Job size for each class.
#' @param p_i Discrete distribution of server count required by client.
#' @param f Speed in each modes.
#' @param P_a Mode switching stochastic matrix when arrivals.
#' @param P_d Mode switching stochastic matrix when departures.
#' @return Class that describes model.
build_model = function(lambda, N, mu, p_i, f, P_a, P_d)
{

    if((nrow(P_a) == nrow(P_d)) && (nrow(P_a) == ncol(P_a)) && (nrow(P_a) == ncol(P_d))){
        if(length(mu) != length(p_i)){
            stop("Length of mu and p_i must be equal.")
        }
        if((length(lambda) > 1) || (length(N) > 1)){
            stop("lambda and N must be scalars.")
        }
        if(length(f) != nrow(P_a)){
            stop("Number of velocitites must be equal number of modes.")
        }
        if(N < 1){
            stop("c must be positive integer.")
        }
        if(lambda <= 0){
            stop("Input rate must be greater than zero.")
        }
        return(new(Model, lambda, N, mu, p_i, f, P_a, P_d))
    } else {
        stop("P_a and P_d must be equal size square matrix.")
    }
}