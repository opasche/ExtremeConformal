
#' Single-sided extreme conformal prediction
#'
#' @description The main function of this package, performing single-sided extreme conformal prediction from extreme quantile regression predictions.
#'   This function computes the conformal correction score-quantile \eqn{\hat{q}^{e}_{\alpha}} (or \eqn{\hat{q}_{\alpha}}) 
#'   to be added to the extreme quantile regression prediction(s) to obtain the desired PIs.
#'   It can optionally perform weighted conformal prediction for nonexchangeable data (e.g. under distribution shift or drift),
#'   by providing calibration-sample weights.
#'
#' @param EQR_pred Vector of extreme quantile predictions for the calibration data, 
#'   ideally at level `coverage_proba` (same length as `y_calibr`).
#' @param y_calibr Vector of observed response values of the calibration data (same length as `EQR_pred`).
#' @param coverage_proba,coverage_alpha Marginal coverage probability (or level alpha) for the conformal prediction interval. 
#'   Only one of `coverage_proba` or `coverage_alpha` must be provided, as `coverage_alpha = 1 - coverage_proba`.
#' @param method Method for constructing the conformal prediction interval. Defaults to `'GPD_safeprofile'`.
#'   See the Details section for more information.
#' @param threshold_lvls Threshold probability level for the GPD-based methods. 
#'   The `'GPD_max'` model instead requires a vector of threshold levels.
#'   Defaults to 0.95 (or to a sequence from 0.8 to 0.99 with the `'GPD_max'` method). 
#' @param alpha_correction Confidence correction method for the CI-based extreme conformal prediction intervals 
#'   used to choose \eqn{\alpha_1} and \eqn{\alpha_2}. Defaults to `'Sidak'`.
#' @param correction_prop Relative (pseudo-)proportion between the values of \eqn{\alpha_1} and \eqn{\alpha_2}, 
#'   for the CI-based extreme conformal prediction intervals. Defaults to 0.5 (equal values of \eqn{\alpha_1} and \eqn{\alpha_2}).
#'   Larger values increase \eqn{\alpha_1} and decrease \eqn{\alpha_2}.
#' @param R Number of bootstrap replicates for the `''GPD_boot'` method.
#' @param min_obs_GPD Minimum number of observations above the threshold for the GPD-based methods.
#' @param profile_init_step_pos,profile_init_step_neg Initial binary-search step size, in the positive and negative directions, 
#'   for the profile-likelihood methods. See [ExtremeCI::GPD_profile_CI()] for more details.
#' @param profile_tol Tolerance for the profile-likelihood search. See [ExtremeCI::GPD_profile_CI()] for more details.
#' @param profile_steps_beyond_conf Number of safety steps beyond the confidence line for the profile-likelihood search. 
#'   See [ExtremeCI::GPD_profile_CI()] for more details.
#' @param profile_max_steps Maximum number of initial search steps for the profile-likelihood CIs.
#'   See [ExtremeCI::GPD_profile_CI()] for more details.
#' @param alpha_profile_naive CI confidence level (i.e., the equivalent of \eqn{\alpha_2}) for the `'GPD_profile_naive'` method.
#' @param obs_weights Optional vector of sample weights (same length as `y_calibr`), to perform weighted conformal prediction 
#'   for nonexchangeable data (e.g. under distribution shift or drift).
#' @param test_weight Weight of the test point for which the weighted conformal prediction is performed. 
#'   Only necessary for the `'Empirical'` method if `obs_weights` is provided (defaults to `max(obs_weights)`).
#' @param return_q_list Boolean indicating whether to return a debug quantile list.
#' @param verbose Verbose level (0 for no messages, 1 for warnings, 2 for warnings duplicated as [cat()] prints).
#' @param .pre_sorted_scores (For development only.) If already computed: Vector of pre-sorted calibration nonconformity scores.
#'
#' @returns A named list containing the following elements.
#' \item{dQ}{The conformal correction \eqn{\hat{q}^{e}_{\alpha}} or \eqn{\hat{q}_{\alpha}} to be added to the extreme quantile prediction(s) to obtain the conformal PI endpoint.}
#' \item{coverage_proba}{The marginal coverage probability for the conformal PI.}
#' \item{coverage_alpha}{The marginal coverage alpha level for the conformal PI.}
#' \item{threshold_lvls}{The threshold level(s) effectively used for the GPD-based methods.}
#' \item{method}{The `method` used to obtain the conformal PI.}
#' \item{dQ_thresholds}{A debugging vector of quantiles at the specified threshold levels. Only for GPD-based methods when `return_q_list==TRUE`.}
#' \item{sigma}{The GPD scale parameter estimate. Currently only for method 'GPD_simple', when `return_q_list==TRUE`.}
#' \item{xi}{The GPD shape parameter estimate. Currently only for method 'GPD_simple', when `return_q_list==TRUE`.}
#' 
#' @export
#'
#' @details The `method` argument specifies the conformalization method used to construct the (extreme) conformal prediction intervals (PIs). 
#' The available method options are:
#' \describe{
#'   \item{'GPD_safeprofile'}{Recommended choice for extreme conformal prediction. 
#'     Tries the method 'GPD_profile' first, and falls back to 'GPD_boot' if the former suffers from numerical instability.}
#'   \item{'GPD_profile'}{GPD-based extreme conformalization using the profile-likelihood CI endpoint of the extreme score quantile.
#'     It captures the asymetric uncertainty of the score-quantile best, yielding the most reliable coverage.
#'     It might sometimes overcover or suffer from numerical convergence issues.}
#'   \item{'GPD_boot'}{GPD-based extreme conformalization using the nonparametric bootstrap percentile CI endpoint of the extreme score quantile.}
#'   \item{'GPD_delta'}{GPD-based extreme conformalization using the Delta method CI endpoint of the extreme score quantile.}
#'   \item{'Empirical'}{The classical (non-extreme) conformalized quantile regression method, relying on the empirical quantile of the nonconformity scores.
#'     Yields infinitely wide PIs if the `coverage_proba` is larger than `1-1/(length(y_calibr)+1`).}
#'   \item{'GPD_profile_naive'}{A naive version of `'GPD_profile'`, without the `alpha_correction` for multiple testing.
#'     Is likely to undercover, only use for comparison.}
#'   \item{'GPD_max'}{A naive GPD-based approach repeating the `'GPD_simple'` method for a range of threshold levels, keeping the most conservative results.
#'     Is likely to undercover, only use for comparison.}
#'   \item{'GPD_simple'}{A naive GPD-based extreme conformalization using a simple extrapolated score quantile estimate instead of a CI endpoint.
#'     Is likely to undercover, only use for comparison.}
#'   \item{'Hill'}{(Not implemented) Extreme conformalization based on the Hill estimator from extreme value analysis.}
#' }
#' See Pasche et al. (2025), referenced below, for the technical details of extreme conformal prediction.
#'
#' @references 
#' Pasche, O. C., Lam, H., and Engelke, S. (2025). "Extreme Conformal Prediction: Reliable Intervals for High-Impact Events." *ArXiv Preprint*. \doi{doi:10.48550/arXiv.2505.08578}.
conformalize_EQR_unilat <- function(EQR_pred, y_calibr, coverage_proba=NULL, coverage_alpha=NULL,
                                    method = c('GPD_safeprofile','GPD_boot','GPD_profile','Empirical',
                                               'GPD_profile_naive','GPD_delta','GPD_max','GPD_simple','Hill'),
                                    threshold_lvls=NULL, alpha_correction=c('Sidak','Bonferroni'), correction_prop=0.5,
                                    R=1004, min_obs_GPD=10, profile_init_step_pos=100, profile_init_step_neg=10,
                                    profile_tol=0.001, profile_steps_beyond_conf=5, profile_max_steps=1e3,
                                    alpha_profile_naive=0.01, obs_weights=NULL, test_weight=NULL,
                                    return_q_list=FALSE, verbose=1, .pre_sorted_scores=NULL) {
  # TODO: option for list of weight vectors and block IDs for seasonal weights, store dQ vector.
  # TODO: Other alpha_corrections?
  
  method <- match.arg(method)
  alpha_correction <- match.arg(alpha_correction)
  if(correction_prop <= 0 || correction_prop >= 1){stop("'correction_prop' should be in (0,1).")}
  if(!is.null(coverage_proba)){
    if(!is.null(coverage_alpha)){
      if(coverage_proba!=(1-coverage_alpha)){stop('coverage_proba should equal 1-coverage_alpha.')}
    }else{
      coverage_alpha <- 1-coverage_proba
    }
  }else if(!is.null(coverage_alpha)){
    coverage_proba <- 1-coverage_alpha
  }else{
    stop("Please provide either 'coverage_proba' or 'coverage_alpha' to 'conformalize_EQR_unilat'.")
  }
  
  if(any(method == c('GPD_safeprofile','GPD_profile','GPD_boot','GPD_delta','GPD_simple','GPD_profile_naive'))){
    if(is.null(threshold_lvls)){
      threshold_lvls <- 0.95
    }
    if(length(threshold_lvls)>1){stop("Provide single 'threshold_lvls' for the selected 'method'.")}
  }
  if(any(method == c('GPD_max'))){
    if(is.null(threshold_lvls)){
      threshold_lvls <- seq(0.8,0.99,0.01)
    }
    threshold_lvls <- unique(threshold_lvls)
  }
  
  if(is.null(.pre_sorted_scores)){
    .pre_sorted <- FALSE
    n <- length(y_calibr)
    if(length(EQR_pred)!=n){stop("'EQR_pred' and 'y_calibr' should have the same length in 'conformalize_*()'.")}
    
    s_cal <- (y_calibr - EQR_pred)
    
  }else{
    .pre_sorted <- TRUE
    n <- length(.pre_sorted_scores)
    
    s_cal <- .pre_sorted_scores
  }
  
  if(any(method == c('GPD_safeprofile','GPD_boot','GPD_profile','GPD_profile_naive','GPD_delta','GPD_max','GPD_simple'))){
    if(n*(1-min(threshold_lvls)) < max(min_obs_GPD,5)){
      stop(paste0("The number of calibration observations is too small for the selected method '", method,"' and 'threshold_lvls', given 'min_obs_GPD'.\n",
                  "For the selected method and 'threshold_lvls', at least ", ceiling(max(min_obs_GPD,5)/(1-min(threshold_lvls))), " observations are required.\n",
                  "Please consider more calibration observations or decrease 'threshold_lvls' (or decrease 'min_obs_GPD' if it is large).\n",
                  "The calibration size n should be > 1000, and typical threshold_lvls values are between 0.9 and 0.99, depending on the nature of the data and on n."))
    }
    
    if(coverage_proba <= min(threshold_lvls)){
      verb_warning(paste0("The 'coverage_proba' is smaller than 'threshold_lvls' (or 'coverage_alpha' > 1-threshold_lvls). The 'method' was changed to 'Empirical' instead of '", method,"'."),
                   verbose=verbose)
      method <- 'Empirical'
    }
    
  }
  
  # Check weights
  if(!is.null(obs_weights)){
    if(length(obs_weights)!=n){stop('obs_weights should be of the same length as Y.')}
    if(any(obs_weights<0)){stop('obs_weights should be non-negative.')}
    if(sum(obs_weights)==0){stop('obs_weights cannot be all zero.')}
    if(is.null(test_weight) & method=='Empirical'){
      test_weight <- max(obs_weights)
      verb_warning(paste0("'test_weight' was set to the max of 'obs_weights': ", roundm(test_weight,4)), verbose=verbose)
    }
  }
  if(!is.null(test_weight)){
    if(test_weight<0){stop('test_weight should be non-negative.')}
    if(test_weight==0){stop('test_weight cannot be zero.')}
  }
  
  #TODO: normalize weights?
  # w <- w / sum(w) * length(w) # sum(w) = length(w)
  
  if(any(method == c('GPD_safeprofile','GPD_profile','GPD_boot','GPD_delta'))){
    if(alpha_correction=='Bonferroni'){
      alpha_profile <- coverage_alpha*(1-correction_prop)
      alpha_profile_qlvl <- coverage_alpha*correction_prop
    }else if(alpha_correction=='Sidak'){
      alpha_profile <- 1-(1-coverage_alpha)^(1-correction_prop)
      alpha_profile_qlvl <- 1-(1-coverage_alpha)^(correction_prop)
    }else{
      stop('alpha_correction not implemented.')
    }
  }
  
  
  if (method == 'Empirical') {
    dQ_res_list <- conformalize_Empirical(s_cal=s_cal, coverage_proba=coverage_proba,
                                          obs_weights=obs_weights, test_weight=test_weight,
                                          verbose=verbose, .pre_sorted=.pre_sorted)
    
  } else if (method == 'GPD_simple') {
    dQ_res_list <- conformalize_GPD_max(s_cal=s_cal, coverage_proba=coverage_proba,
                                        threshold_lvls=threshold_lvls, min_obs_GPD=0, obs_weights=obs_weights,
                                        verbose=verbose, .pre_sorted=.pre_sorted)
    
  } else if (method == 'GPD_max') {
    dQ_res_list <- conformalize_GPD_max(s_cal=s_cal, coverage_proba=coverage_proba,
                                        threshold_lvls=threshold_lvls, min_obs_GPD=min_obs_GPD, obs_weights=obs_weights,
                                        verbose=verbose, .pre_sorted=.pre_sorted)
    threshold_lvls <- dQ_res_list$threshold_lvls
    
  } else if (method=='GPD_safeprofile' | method=='GPD_profile' | method=='GPD_profile_naive'){
    if (method=='GPD_profile_naive'){
      alpha_profile_qlvl <- coverage_alpha
      alpha_profile <- alpha_profile_naive
    }
    dQ_res_list <- conformalize_GPD_profile(s_cal=s_cal, coverage_proba=coverage_proba, threshold_lvls=threshold_lvls,
                                            alpha_profile=alpha_profile, alpha_profile_qlvl=alpha_profile_qlvl,
                                            init_step_pos=profile_init_step_pos, init_step_neg=profile_init_step_neg, tol=profile_tol,
                                            steps_beyond_conf=profile_steps_beyond_conf, max_steps=profile_max_steps,
                                            obs_weights=obs_weights, verbose=verbose, .pre_sorted=.pre_sorted)
    
    if(method=='GPD_safeprofile' & is.infinite(dQ_res_list$dQ)){
      warning(paste0("The profile-based 'dQ' is infinite for the method '", method, "'. The bootstrap-based 'dQ' is used instead."))
      dQ_res_list <- conformalize_GPD_boot(s_cal=s_cal, threshold_lvls=threshold_lvls,
                                           alpha_profile=alpha_profile, alpha_profile_qlvl=alpha_profile_qlvl,
                                           R=R, obs_weights=obs_weights, verbose=verbose, .pre_sorted=.pre_sorted)
    }
    
  } else if (method == 'GPD_boot') {
    dQ_res_list <- conformalize_GPD_boot(s_cal=s_cal, threshold_lvls=threshold_lvls,
                                         alpha_profile=alpha_profile, alpha_profile_qlvl=alpha_profile_qlvl,
                                         R=R, obs_weights=obs_weights, verbose=verbose, .pre_sorted=.pre_sorted)
    
  } else if (method == 'GPD_delta') {
    dQ_res_list <- conformalize_GPD_delta(s_cal=s_cal, threshold_lvls=threshold_lvls,
                                          alpha_profile=alpha_profile, alpha_profile_qlvl=alpha_profile_qlvl,
                                          obs_weights=obs_weights, verbose=verbose, .pre_sorted=.pre_sorted)
    
  } else {
    stop('Conformalization method not implemented.')
  }
  
  dQ <- dQ_res_list$dQ
  dQ_thresholds <- dQ_res_list$dQ_thresholds
  
  EQR_conformalizer <- list(dQ=dQ,
                            coverage_proba=coverage_proba,
                            coverage_alpha=coverage_alpha,
                            threshold_lvls=threshold_lvls,
                            method=method)
  if(return_q_list){
    if(any(method == c('GPD_safeprofile','GPD_boot','GPD_profile','GPD_profile_naive','GPD_delta','GPD_max','GPD_simple'))){
      names(dQ_thresholds) <- paste0(threshold_lvls*100,'%')
    }
    EQR_conformalizer$dQ_thresholds <- dQ_thresholds
    if (method == 'GPD_simple'){
      EQR_conformalizer$sigma <- dQ_res_list$sigma
      EQR_conformalizer$xi <- dQ_res_list$xi
    }
  }
  
  class(EQR_conformalizer) <- c('EQR_conformalizer','unilat_conformalizer','conformalizer')
  
  return(EQR_conformalizer)
}

#' Single-sided conformal prediction interval from conformalizer
#'
#' @param Q_pred Vector of extreme quantile regression predictions for the test data (same length as the number of test points).
#' @param unilat_conformalizer Either a conformalizer object obtained from [conformalize_EQR_unilat()], or, directly, the conformal correction to be added to the extreme quantile regression predictions.
#'   In the latter case, either a single conformal correction value or a vector of conformal corrections of the same length as `Q_pred` is expected.
#' @param return_format Format of the returned prediction interval. Either 'upper_limit' for a numerical vector of upper PI limits (default), 
#'   'interval' for a data frame with columns of lower and upper PI limits, or 'text' for a single string description of the PI (only available for single predictions).
#' @param ymin Lower endpoint for the response distribution (if known). Can be a single value (marginal lower endpoint) or 
#' a vector of the same length as `Q_pred` (conditional lower endpoint). Default is `-Inf`.
#' @param coverage_proba,coverage_alpha (Optional) Marginal coverage probability (or level alpha) for the conformal prediction interval. 
#'   Only one of `coverage_proba` or `coverage_alpha` must be provided, as `coverage_alpha = 1 - coverage_proba`. Only used for certain `return_format` options.
#'
#' @returns Depending on the `return_format` argument, either a numerical vector of upper prediction interval (PI) limits, 
#'   a data frame with lower and upper PI limits as columns (and optionally coverage probability and alpha), 
#'   or a single string description of the PI (for single predictions only).
#' @export
#'
#' @examples conformal_PI_unilat(Q_pred=c(10,12), unilat_conformalizer=0.3, return_format='upper_limit')
conformal_PI_unilat <- function(Q_pred, unilat_conformalizer, return_format=c('upper_limit','interval','text'), 
                                ymin=c(-Inf), coverage_proba=NULL, coverage_alpha=NULL){
  return_format <- match.arg(return_format)
  #TODO: dims Q_pred
  n <- length(Q_pred)
  if(any(class(unilat_conformalizer)=='unilat_conformalizer')){ #todo: class
    dQ <- unilat_conformalizer$dQ
    coverage_proba <- unilat_conformalizer$coverage_proba
    coverage_alpha <- unilat_conformalizer$coverage_alpha
  }else{
    dQ <- unilat_conformalizer
  }
  if(!is.null(coverage_proba)){
    if(!is.null(coverage_alpha)){
      if(coverage_proba!=1-coverage_alpha){stop('coverage_proba should equal 1-coverage_alpha.')}
    }else{
      coverage_alpha <- 1-coverage_proba
    }
  }else if(!is.null(coverage_alpha)){
    coverage_proba <- 1-coverage_alpha
  }
  
  PI_up <- Q_pred + dQ
  
  if(length(ymin)==1){
    ymin <- rep(ymin, n)
  }else{
    if(length(ymin)!=n){stop("'ymin' should be either a single value or a vector of the same length as 'Q_pred'.")}
  }
  # PI_up <- pmax(PI_up, ymin)
  if(any(PI_up < ymin)){
    stop("In 'conformal_PI_unilat': Some upper PI limits are smaller than 'ymin'. 
         There must be an error with the provided 'Q_pred', 'unilat_conformalizer', or 'ymin'.")
  }
  
  if(return_format=='upper_limit'){
    return(PI_up)
  }else if(return_format=='interval'){
    res <- data.frame(PI_down=c(ymin),
                      PI_up=c(PI_up))
    if(!is.null(coverage_proba)){
      res$coverage_proba <- rep(coverage_proba, n)
      res$coverage_alpha <- rep(coverage_alpha, n)
    }
    return(res)
  }else{
    if(n>1){stop("return_format='text' is only available for single predictions.")}
    res <- paste0('Prediction Interval: (', ymin, ', ', PI_up, ')')
    if(!is.null(coverage_proba)){
      res <- paste0(res, ' with ', coverage_proba*100, '% confidence')
    }
    return(res)
  }
}


#' Single-sided block-weighted conformal prediction interval from conformalizer
#'
#' @param Q_pred Vector of extreme quantile regression predictions for the test data (same length as the number of test points).
#' @param block_ids Vector of block IDs for each test point, of the same length as `Q_pred`. 
#'   The block IDs should be integers between 1 and the number of blocks, and correspond to the order of the `dQ_blocks` vector.
#' @param dQ_blocks Vector of block-specific conformal corrections to be added to the extreme quantile regression predictions, of the same length as the number of blocks.
#' @param return_format Format of the returned prediction interval. Either 'upper_limit' for a numerical vector of upper PI limits (default), 
#'   'interval' for a data frame with columns of lower and upper PI limits, or 'text' for a single string description of the PI (only available for single predictions).
#' @param ymin Lower endpoint for the response distribution (if known). Can be a single value (marginal lower endpoint), or 
#'   a vector of the same length as `Q_pred` (conditional lower endpoint), or 
#'   a vector of the same length as `Q_pred` (block-conditional lower endpoint). Default is `-Inf`.
#' @param coverage_proba,coverage_alpha (Optional) Marginal coverage probability (or level alpha) for the conformal prediction interval. 
#'   Only one of `coverage_proba` or `coverage_alpha` must be provided, as `coverage_alpha = 1 - coverage_proba`. Only used for certain `return_format` options.
#'
#' @returns Depending on the `return_format` argument, either a numerical vector of upper prediction interval (PI) limits, 
#'   a data frame with lower and upper PI limits, and the block IDs, as columns (and optionally coverage probability and alpha), 
#'   or (only available for compatibility) a single string description of the PI (for single predictions only).
#' @export
block_weighted_conformal_PI_unilat <- function(Q_pred, block_ids, dQ_blocks, return_format=c('upper_limit','interval','text'),
                                               ymin=c(-Inf), coverage_proba=NULL, coverage_alpha=NULL){
  return_format <- match.arg(return_format)
  
  # Check Q_pred and dQ_blocks are vectors and not matrices/arrays
  if((!is.vector(Q_pred) || is.matrix(Q_pred) || is.array(Q_pred)) && (ncol(Q_pred) > 1)){stop("'Q_pred' should be a vector.")}
  if((!is.vector(dQ_blocks) || is.matrix(dQ_blocks) || is.array(dQ_blocks)) && (ncol(dQ_blocks) > 1 && (nrow(dQ_blocks) > 1))){stop("'dQ_blocks' should be a vector.")}
  if((!is.vector(block_ids) || is.matrix(block_ids) || is.array(block_ids)) && (ncol(block_ids) > 1)){stop("'block_ids' should be a vector.")}
  
  n <- length(Q_pred)
  nb_blocks <- length(dQ_blocks)
  
  if(length(block_ids)!=n){stop("'block_ids' should be of the same length as 'Q_pred'.")}
  if(max(block_ids)>nb_blocks || min(block_ids)<1){stop("'block_ids' entries should be between 1 and the length of 'dQ_blocks'.")}
  
  dQ_obs <- c(dQ_blocks)[c(block_ids)]
  
  PI_up <- c(Q_pred) + dQ_obs
  
  if(length(ymin)==1){
    ymin <- rep(ymin, n)
  }else if (length(ymin)==nb_blocks) {
     ymin <- c(ymin)[c(block_ids)]
  }else if (length(ymin)==n) {
     ymin <- c(ymin)
  }else{
    stop("'ymin' should be either a single value or a vector of the same length as 'Q_pred'.")
  }
  # PI_up <- pmax(PI_up, ymin)
  if(any(PI_up < ymin)){
    stop("In 'block_weighted_conformal_PI_unilat': Some upper PI limits are smaller than 'ymin'. 
         There must be an error with the provided 'Q_pred', 'dQ_blocks', or 'ymin'.")
  }
  
  if(return_format=='upper_limit'){
    return(PI_up)
  }else if(return_format=='interval'){
    res <- data.frame(PI_down=ymin,
                      PI_up=c(PI_up),
                      block_id=c(block_ids))
    if(!is.null(coverage_proba)){
      res$coverage_proba <- rep(coverage_proba, n)
      res$coverage_alpha <- rep(coverage_alpha, n)
    }
    return(res)
  }else{
    if(n>1){stop("return_format='text' is only available for single predictions.")}
    res <- paste0('Prediction Interval: (', ymin, ', ', PI_up, ')')
    if(!is.null(coverage_proba)){
      res <- paste0(res, ' with ', coverage_proba*100, '% confidence')
    }
    return(res)
  }
}



#' Compute extreme quantile from GPD parameters
#'
#' @param p Probability level of the desired extreme quantile.
#' @param p0 Probability level of the (possibly varying) intermediate threshold/quantile.
#' @param t_x0 Value(s) of the (possibly varying) intermediate threshold/quantile.
#' @param sigma Value(s) for the GPD scale parameter.
#' @param xi Value(s) for the GPD shape parameter.
#'
#' @return The quantile value at probability level `p`.
#'
#' @keywords internal
GPD_quantiles <- function(p, p0, t_x0, sigma, xi){
  if(any(xi==0)){
    gpd_qs <- rep(as.double(NA), max(length(t_x0),length(sigma),length(xi)))
    if(length(xi)==1){xi <- rep(xi,length(gpd_qs))}
    gpd_qs[xi!=0] <- ((((1-p)/(1-p0))^{-xi} - 1) * (sigma / xi) + t_x0)[xi!=0]
    gpd_qs[xi==0] <- (log((1-p0)/(1-p)) * sigma + t_x0)[xi==0]
  }else{
    gpd_qs <- (((1-p)/(1-p0))^{-xi} - 1) * (sigma / xi) + t_x0
  }
  return(gpd_qs)
}

#' Weighted empirical quantile estimate for a single probability level
#'
#' @param x .
#' @param prob .
#' @param obs_weights .
#' @param pre_sorted .
#'
#' @returns The weighted empirical quantile estimate.
#'
#' @keywords internal
weighted_quantile <- function(x, prob, obs_weights=NULL, pre_sorted=FALSE){
  if(is.null(obs_weights)){obs_weights <- rep(1,length(x))}
  if(!pre_sorted){
    o <- order(x)
    x <- x[o]
    obs_weights <- obs_weights[o]
  }
  
  inds <- which(cumsum(obs_weights/sum(obs_weights)) >= prob)
  if(length(inds)==0){ # if prob too large or some infinite weights
    return(Inf)
  }else{
    return(x[min(inds)])
  }
}


#' Weighted empirical quantile estimates
#'
#' @param x .
#' @param probs .
#' @param obs_weights .
#' @param pre_sorted .
#'
#' @returns A vector of weighted empirical quantile estimates at the specified probability levels.
#'
#' @keywords internal
weighted_quantiles <- function(x, probs, obs_weights=NULL, pre_sorted=FALSE){
  return(sapply(probs, function(p) weighted_quantile(x=x, prob=p, obs_weights=obs_weights, pre_sorted=pre_sorted)))
}



# INTERNAL CONFORMALIZATION PROCEDURES

#' Classical split conformalization procedure with optional sample weights
#'
#' @param s_cal .
#' @param coverage_proba .
#' @param obs_weights .
#' @param test_weight .
#' @param verbose .
#' @param .pre_sorted .
#'
#' @returns A named list containing the conformal correction `dQ` and a redundant `dQ_thresholds`.
#'
#' @keywords internal
conformalize_Empirical <- function(s_cal, coverage_proba, obs_weights=NULL, test_weight=NULL, verbose=1, .pre_sorted=FALSE) {
  n <- length(s_cal)
  
  if(is.null(obs_weights)){
    
    if(coverage_proba <= 1 - 1/(n+1)){
      if(.pre_sorted){
        dQ <- s_cal[ceiling((n+1)*coverage_proba)]
      }else{
        dQ <- stats::quantile(s_cal, ceiling((n+1)*coverage_proba)/n)
      }
    } else {
      dQ <- Inf
      verb_warning("'coverage_proba' is larger than 1-1/(n+1) with method 'Empirical' in 'conformalize_EQR'. Infitinte intervals are generated.",
                   verbose=verbose)
    }
    
  } else {
    
    dQ <- weighted_quantile(c(s_cal,Inf), coverage_proba, obs_weights=c(obs_weights,test_weight), pre_sorted=.pre_sorted)
    if(is.infinite(dQ)){
      verb_warning("Infitinte intervals are generated with method 'Empirical' in 'conformalize_EQR'. 'coverage_proba' is too large or some weights are infinite.",
                   verbose=verbose)
    }
    
  }
  
  return(list(dQ=dQ, dQ_thresholds=dQ))
}



#' GPD-max conformalization procedure with optional sample weights
#'
#' @param s_cal .
#' @param coverage_proba .
#' @param threshold_lvls .
#' @param min_obs_GPD .
#' @param obs_weights .
#' @param verbose .
#' @param .pre_sorted .
#'
#' @returns A named list containing the conformal correction `dQ`, the vector of `dQ_thresholds`, 
#'   the estimated scale parameters `sigma`, the estimated shape parameters `xi`, and the vector of effectively used thresholds `threshold_lvls`.
#'
#' @keywords internal
conformalize_GPD_max <- function(s_cal, coverage_proba, threshold_lvls=NULL, min_obs_GPD=10, obs_weights=NULL, verbose=1, .pre_sorted=FALSE){
  n <- length(s_cal)
  if(is.null(threshold_lvls)){
    threshold_lvls <-  seq(0.8,0.99,0.01)
  }
  threshold_lvls <- unique(threshold_lvls)
  if(min_obs_GPD>0){
    lvllen1 <- length(threshold_lvls)
    # threshold_lvls <- threshold_lvls[threshold_lvls <= 1-min_obs_GPD/n]
    threshold_lvls[threshold_lvls > 1-min_obs_GPD/n] <- 1-min_obs_GPD/n
    threshold_lvls <- unique(threshold_lvls)
    lvllen2 <- length(threshold_lvls)
    if(lvllen1!=lvllen2){
      verb_warning(paste0("'threshold_lvls' were discarded due to 'min_obs_GPD'. Remaining: ",
                          lvllen1,"/",lvllen2," (",lvllen2/lvllen1*100,"%)"),
                   verbose=verbose)
    }
  }
  # TODO: check if threshold_lvls < coverage_proba for each threshold_lvls, or replace by empirical in loop?
  
  w <- if(is.null(obs_weights)){1}else{obs_weights}
  
  if(is.null(obs_weights)){
    if(.pre_sorted){
      # u <- s_cal[ceiling(n*threshold_lvls)]
      u_seq <- c(stats::quantile(s_cal, threshold_lvls))
    }else{
      # u_seq <- unique(c(stats::quantile(s_cal, threshold_lvls)))
      u_seq <- c(stats::quantile(s_cal, threshold_lvls))
    }
  }else{
    u_seq <- weighted_quantiles(s_cal, threshold_lvls, obs_weights=obs_weights, pre_sorted=.pre_sorted)
  }
  
  dQ_thresholds <- rep(as.double(NA), length(threshold_lvls))
  sigma <- rep(as.double(NA), length(threshold_lvls))
  xi <- rep(as.double(NA), length(threshold_lvls))
  for (i in seq_along(u_seq)) {
    try({ # debug TODO: rm at some point
      u <- u_seq[i]
      if(is.null(obs_weights)){
        # TODO: change to new MLE from ExtremeCI?
        pars <- ismev::gpd.fit(s_cal, u, show = FALSE, maxit=1e6)$mle
        sigma[i] <- pars[1]
        xi[i] <- pars[2]
      }else{
        # TODO: change to new weighted MLE from ExtremeCI?
        GPD_extRemes <- extRemes::fevd(s_cal, threshold=u, type="GP", method='MLE', time.units="years", period.basis="obs",
                                       weights=w)
        pars <- GPD_extRemes$results$par
        sigma[i] <- pars["scale"]
        xi[i] <- pars["shape"]
      }
      dQ_thresholds[i] <- GPD_quantiles(coverage_proba, threshold_lvls[i], u, sigma[i], xi[i])
    })
  }
  dQ <- max(dQ_thresholds, na.rm=TRUE)
  
  return(list(dQ=dQ, dQ_thresholds=dQ_thresholds, sigma=sigma, xi=xi, threshold_lvls=threshold_lvls))
}



#' GPD-boot conformalization procedure with optional sample weights
#'
#' @param s_cal .
#' @param threshold_lvls .
#' @param alpha_profile .
#' @param alpha_profile_qlvl .
#' @param R .
#' @param obs_weights .
#' @param verbose .
#' @param .pre_sorted .
#'
#' @returns A named list containing the conformal correction `dQ`, a redundant `dQ_thresholds`, 
#'   and placeholders for the estimated scale parameters `sigma` and shape parameters `xi` (currently set to NULL).
#'
#' @keywords internal
conformalize_GPD_boot <- function(s_cal, threshold_lvls, alpha_profile, alpha_profile_qlvl, R=1004, obs_weights=NULL, verbose=1, .pre_sorted=FALSE){
  n <- length(s_cal)
  
  w <- if(is.null(obs_weights)){1}else{obs_weights}
  
  if(is.null(obs_weights)){
    if(.pre_sorted){
      u <- s_cal[ceiling(n*threshold_lvls)]
    }else{
      u <- stats::quantile(s_cal, threshold_lvls)
    }
  }else{
    u <- weighted_quantile(s_cal, threshold_lvls, obs_weights=obs_weights, pre_sorted=.pre_sorted)
  }
  
  GPD_extRemes <- extRemes::fevd(s_cal, threshold=u, type="GP", method='MLE', time.units="years", period.basis="obs",
                                 weights=w)
  bci100 <- extRemes::ci.fevd(GPD_extRemes, threshold=u, alpha=alpha_profile, type="return.level",
                              return.period=1./(alpha_profile_qlvl), R=R, method="boot", verbose=FALSE)
  # GPD_extRemes <- extRemes::fevd(s_cal, threshold=u, type="GP", method='MLE', time.units="days", period.basis="year", weights=w)
  # bci100 <- extRemes::ci.fevd(GPD_extRemes, threshold=u, alpha=alpha_profile, type="return.level",
  #                             return.period=1./(alpha_profile_qlvl*365.25), R=R, method="boot", verbose=FALSE)
  dQ <- bci100[[3]]
  return(list(dQ=dQ, dQ_thresholds=dQ, sigma=NULL, xi=NULL)) #TODO: sigma, xi OPTIONAL
}


#' GPD-delta conformalization procedure with optional sample weights
#'
#' @param s_cal .
#' @param threshold_lvls .
#' @param alpha_profile .
#' @param alpha_profile_qlvl .
#' @param obs_weights .
#' @param verbose .
#' @param .pre_sorted .
#'
#' @returns A named list containing the conformal correction `dQ`, a redundant `dQ_thresholds`, 
#'   and placeholders for the estimated scale parameters `sigma` and shape parameters `xi` (currently set to NULL).
#'
#' @keywords internal
conformalize_GPD_delta <- function(s_cal, threshold_lvls, alpha_profile, alpha_profile_qlvl, obs_weights=NULL, verbose=1, .pre_sorted=FALSE){
  n <- length(s_cal)
  
  w <- if(is.null(obs_weights)){1}else{obs_weights}
  
  if(is.null(obs_weights)){
    if(.pre_sorted){
      u <- s_cal[ceiling(n*threshold_lvls)]
    }else{
      u <- stats::quantile(s_cal, threshold_lvls)
    }
  }else{
    u <- weighted_quantile(s_cal, threshold_lvls, obs_weights=obs_weights, pre_sorted=.pre_sorted)
  }
  
  GPD_extRemes <- extRemes::fevd(s_cal, threshold=u, type="GP", method='MLE', time.units="years", period.basis="obs",
                                 weights=w)
  nci100 <- extRemes::ci.fevd(GPD_extRemes, threshold=u, alpha=alpha_profile, type="return.level",
                              return.period=1./(alpha_profile_qlvl), method="normal", verbose=FALSE)
  # GPD_extRemes <- extRemes::fevd(s_cal, threshold=u, type="GP", method='MLE', time.units="days", period.basis="year", weights=w)
  # nci100 <- extRemes::ci.fevd(GPD_extRemes, threshold=u, alpha=alpha_profile, type="return.level",
  #                             return.period=1./(alpha_profile_qlvl*365.25), method="normal", verbose=FALSE)
  dQ <- nci100[[3]]
  return(list(dQ=dQ, dQ_thresholds=dQ, sigma=NULL, xi=NULL)) #TODO: sigma, xi OPTIONAL
}


#' GPD-profile conformalization procedure with optional sample weights
#'
#' @param s_cal .
#' @param coverage_proba .
#' @param threshold_lvls .
#' @param alpha_profile .
#' @param alpha_profile_qlvl .
#' @param init_step_pos .
#' @param init_step_neg .
#' @param tol .
#' @param steps_beyond_conf .
#' @param max_steps .
#' @param obs_weights .
#' @param verbose .
#' @param .pre_sorted .
#'
#' @returns A named list containing the conformal correction `dQ`, a redundant `dQ_thresholds`, 
#'   and placeholders for the estimated scale parameters `sigma` and shape parameters `xi` (currently set to NULL).
#'
#' @keywords internal
conformalize_GPD_profile <- function(s_cal, coverage_proba, threshold_lvls, alpha_profile, alpha_profile_qlvl,
                                     init_step_pos=100, init_step_neg=10, tol=0.001, steps_beyond_conf=5, max_steps=1e3,
                                     obs_weights=NULL, verbose=1, .pre_sorted=FALSE){
  n <- length(s_cal)
  
  if(is.null(obs_weights)){
    if(.pre_sorted){
      u <- s_cal[ceiling(n*threshold_lvls)]
    }else{
      u <- stats::quantile(s_cal, threshold_lvls, names=FALSE)
    }
  }else{
    u <- weighted_quantile(s_cal, threshold_lvls, obs_weights=obs_weights, pre_sorted=.pre_sorted)
  }
  
  pll_rl <- ExtremeCI::GPD_profile_CI(s_cal, threshold=u, threshold_lvl=threshold_lvls, parameter="quantile", subparam_id=0,
                                      alpha=alpha_profile, quantile_lvl=(1-alpha_profile_qlvl),
                                      # X=matrix(df_era5_past$GMST_nonsmooth_ANO), x_rlvl=last_elem(df_era5_all$GMST_nonsmooth_ANO),
                                      scale_cols=NULL, shape_cols=NULL, warmstart_table=NULL,
                                      init_step_pos=init_step_pos, init_step_neg=init_step_neg, tol=tol, steps_beyond_conf=steps_beyond_conf,
                                      initial_MLE_para="classical", max_steps=max_steps,
                                      obs_weights=obs_weights, ill_defined_value=-10^6,
                                      hessian=TRUE, maxit=1e6, method="Nelder-Mead", method_prof="BFGS",
                                      verbose=verbose) #"Nelder-Mead" "Brent" "BFGS"
  # cat("profile 100y return level estimate for 2021:", pll_rl$mle["return_level"], "and CI (", pll_rl$ci, ")\n")
  dQ <- pll_rl$ci[2]-0.02
  return(list(dQ=dQ, dQ_thresholds=dQ, sigma=NULL, xi=NULL)) #TODO: sigma, xi OPTIONAL?
}







