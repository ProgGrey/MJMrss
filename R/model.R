
 # Copyright (c) 2022-2023 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 # All rights reserved.
 #
 # This source code is licensed under the BSD-3 clause license.

#' @export build_model

#' @importFrom methods new

#' @title build_model
#' @description Built multiserver job model with random speed scaling
#' @param lambda Input rate.
#' @param N Number of servers.
#' @param class_desc Description for clients class, i.e. matrix where first row is server count, second row is probability, last row is mu.
#' @param f Speed in each modes.
#' @param P_a Mode switching stochastic matrix when arrivals.
#' @param P_d Mode switching stochastic matrix when departures.
#' @return Class that describes model.
#' @examples
#' library(MJMrss)
#' lambda = 10 # Input rate
#' N = 10 # Number of servers
#' f = c(1, 2.2) # Speed (frequency)
#' P_a = matrix(c(0.1, 0.9,
#'                0.0, 1.0), nrow = 2, byrow = TRUE)
#' P_d = matrix(c(1.0, 0.0,
#'                0.2, 0.8), nrow = 2, byrow = TRUE)
#' classes = matrix(c(1,      2,    3,   4,# j  
#'                    0.5, 0.25, 0.15, 0.1,# p_j
#'                    1,    1.9,  2.5, 3), # mu_j
#'                    nrow = 3, byrow = TRUE)
#' m = build_model(lambda, N, classes, f, P_a, P_d)
#' # mean customers in system and mean queue length:
#' m$mean_clients
#' m$mean_queue
#' # Distribution from level 0 to level 20:
#' dist = m$distribution(20)
build_model = function(lambda, N, class_desc, f, P_a, P_d)
{
    if((nrow(P_a) == nrow(P_d)) && (nrow(P_a) == ncol(P_a)) && (nrow(P_a) == ncol(P_d))){
        if((nrow(class_desc) != 3 ) || (ncol(class_desc) < 1)){
            stop("class_desc must be matrix with 3 rows and at least 1 column.")
        }
        if((length(lambda) > 1) || (length(N) > 1)){
            stop("lambda and N must be scalars.")
        }
        if(length(f) != nrow(P_a)){
            stop("Number of velocitites must be equal number of modes.")
        }
        if(N < 1){
            stop("N must be positive integer.")
        }
        if(lambda <= 0){
            stop("Input rate must be greater than zero.")
        }
        if(min(class_desc[1,]) < 1){
            stop("Server count must be positive integer.")
        }
        if((min(class_desc[2,]) < 0) || (max(class_desc[2,]) > 1)){
            stop("Second row contains probabilities.")
        }
        if(min(class_desc[3,]) < 0){
            stop("mu is a rate and must be greater than zero.")
        }
        return(new(Model, lambda, N, class_desc, f, P_a, P_d))
    } else {
        stop("P_a and P_d must be equal size square matrix.")
    }
}