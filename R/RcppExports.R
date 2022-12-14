# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name Rcpp_ModelTransient
#' @aliases ModelTransient
#' @aliases Rcpp_ModelTransient-class
#' @rdname ModelTransient-class
#' @title ModelTransient Class
#' @description Internal representation for multi-server job model with random speed scaling for transient analysis.
#' @slot h - Returns current step size.
#' @slot mean_clients.
#' \itemize{
#' \item Parameter: max_time - maximum time for which you want to calculate mean customers in system;
#' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
#' \item Returns: mean clients in system.
#' }
#' @slot mean_queue
#' \itemize{
#' \item Parameter: max_time - maximum time for which you want to calculate mean queue length;
#' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
#' \item Returns: mean queue length
#' }
#' @slot distribution 
#' \itemize{
#' \item Parameter: max_time - maximum level for which you want to calculate distribution;
#' \item Parameter: pi_0 - list of vectors contains distribution at zero time.
#' \item Returns: list of lists of vectors containing transient distribution.
#' }
NULL

#' @name Rcpp_Model
#' @aliases Model
#' @aliases Rcpp_Model-class
#' @rdname Model-class
#' @title Model Class
#' @description Internal representation for multi-server job model with random speed scaling.
#' @slot rho Returns rho computated by using Neuts ergodicity criteria.
#' @slot mean_clients Returns mean clients in system.
#' @slot mean_queue Returns mean queue length.
#' @slot sum_pi_from_c_to_inf Returns sum of distribution vectors from level N to infinity level.
#' @slot pi_0_c Returns distribution for first N levels.
#' @slot distribution 
#' \itemize{
#' \item Parameter: max_level - maximum level for which you want to calculate distribution;
#' \item Returns: list of vectors containing stationary distribution.
#' }
#' @slot level_description 
#' \itemize{
#' \item Parameter: level -  level for which you want to get a description;
#' \item Returns: human-readable description of level in form "m|s_1,s_2,..." where "m" is a speed mode and "s_i" is
#' a number of serviced customers of "i" class.
#' }
#' @slot level_busy_servers 
#' \itemize{
#' \item Parameter: level -  level for which you want to get a description;
#' \item Returns: vector of numbers of busy servers for each state at specified level.
#' }
#' @slot level_serviced_clients 
#' \itemize{
#' \item Parameter: level -  level for which you want to get a description;
#' \item Returns: vector of numbers of serviced servers for each state at specified level.
#' }
#' @slot transient_analysis 
#' \itemize{
#' \item Parameter: order - numerical order of Taylor series method
#' \item Parameter: stepsize - If positive, then  step length of the algorithm. If negative - the multiplier of the 
#' maximum step of the algorithm.
#' \item Returns: internal representation for model transient analysis (ModelTransient class).
#' }
NULL

