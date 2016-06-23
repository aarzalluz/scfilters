#' Compute a quality metric.
#'
#' Computes a statistical metric to assess differences between the distribution
#' of correlation density in the data and the negative control.
#'
#' The aim of this function is to compare the median correlation value for both the
#' mean control curve (see \code{\link{compute_control}}) and the data correlations
#' (see \code{\link{correlate}}), for a given window, and compute the difference. 
#' 
#' This metric reflects the differences between the two distributions, and assesses
#' the quality of the data in that window. The aim is to determine which windows present
#' a significant similarity to a random distribution, and therefore are highly affected by
#' systematic noise, in order to discard them.
#'
#' \strong{Important note:} the \code{compute_histogram} function was designed to operate
#' on the output of the \code{\link{compute_control}} and \code{\link{correlate}}
#' functions.
#'
#' @param control A data frame containing a subset of random control correlation
#' vectors, computed for a given window.
#' 
#' @param correlations A a list containing correlations of each window against 
#' the top window.
#' 
#' @param window_number An integer indicating the number of the window for which
#' the histogram metric is tobe computed.
#'
#' @return A number of type double, namely the difference between the median 
#' correlation values of the control and the window data.

compute_histogram <- function(control, correlations, window_number){

    # compute absolute correlations of control
    abs_control <- list()
    
    for (i in seq(length(control))){
       abs_control[[i]] <- abs(control[[i]])
    }
    
    # calculate mean absolute correlations of all control windows
    all_abs_control <- do.call(rbind, abs_control) %>% as.data.frame()
    mean_control <- colMeans(all_abs_control)
    
    # calculate absolute correlations for selected window number
    abs_window <- abs(correlations[[window_number]])
    
    # compute median of control and windows
    control_median <- median(mean_control)
    window_median <- median(abs_window)

    # compute difference
    hist_value <- window_median - control_median

    return(hist_value)
}
