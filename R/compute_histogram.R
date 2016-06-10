#' Compute a quality metric.
#'
#' Computes a statistical metric to assess differences between data distribution
#' of correlation density and negative control.
#'
#' The aim of this function is to compute the ratio of correlation density that is non-
#' coincident with the density distribution of the negative control. Both distributions are
#' tested to find the section where data density distribution is greater than the control's,
#' and the differences between both are computed. This metric is divided by the total data
#' density of correlation.
#'
#' Graphically, this is the computation of the area where the correlation density distribution
#' of the data does not overlap with that of the negative control, divided by the total area
#' under the data correlation density curve. The aim is to determine which windows present
#' a significant similarity to a random distribution, and therefore are highly affected by
#' systematic noise, in order to discard them.
#'
#' The metric is computed by the mean control curve, as well as for the maximum and minimum
#' curves, obtaining the error values. The returned data frame containis the following three
#' columns:
#'
#' \itemize{
#'      \item \code{hist_value}: the result of computing the metric for the mean control
#'      curve.
#'
#'      \item \code{error_up}: the value of the metric for the maximum control curve.
#'
#'      \item \code{error_down}: the value of the metric for the minimum control curve.
#'  }
#'
#' \strong{Important note:} the \code{compute_histogram} function was designed to operate
#' on the output of the \code{compute_density} function.
#'
#' @param density_data A data frame containing the density data.
#'
#' @return A data frame containing the value of the metric for each window, and values for
#' error bars to plot them as a histogram.

compute_histogram <- function(density_data){

    # create empty vectors to store calculations
    data_y <- vector()
    mean_y <- vector()
    min_y <- vector()
    max_y <- vector()

    # extract the section of the plot where curves are non-overlapping, only in the right side
    for (i in seq_len(nrow(density_data))){

        if ((density_data$data_y[i] > density_data$control_y[i])){

            data_y <- density_data$data_y[i:nrow(density_data)]
            mean_y <- density_data$control_y[i:nrow(density_data)]
            min_y <- density_data$max_control_y[i:nrow(density_data)]
            max_y <- density_data$min_control_y[i:nrow(density_data)]

            break()
        }
    }

    mean_correlated_area <- sum(data_y - mean_y) # mean non-overlapping area
    min_correlated_area <- sum(data_y - min_y) # min non-overlapping area
    max_correlated_area <- sum(data_y - max_y) # max non-overlapping area
    window_area <- sum(density_data$data_y) # total area under the window curve

    #           4. Calculate ratio:
    hist_value <- mean_correlated_area / window_area
    error_down <- max_correlated_area / window_area
    error_up <- min_correlated_area / window_area

    #           5. Bind results:
    histogram <- data.frame(cbind(hist_value, error_up, error_down))
    return(histogram)
}
