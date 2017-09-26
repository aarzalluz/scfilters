#' Determine until when bin of features should be kept based on the metric table
#'
#' Takes the output of \code{\link{get_mean_median}} and decide until which window to keep
#' based on background level and a threshold.
#'
#' Background level is estimated by averaging correlation coefficient obtained
#' from top window randomisation.
#'
#' Bins of genes are kept until there mean (or median) correlation coefficient fell under
#' a threshold value \code{threshold x background level}.
#'
#' @param metric_table A data frame, usually the output of \code{\link{get_mean_median}}.
#'
#' @param threshold How many time higher than the background should the last bin be? Default to 2.
#'
#' @param metric Which metric to use (i.e. which column from metric_table to work with). Default to \code{mean}.
#'
#' @param random_function_summarisation A function used to agregate the randomised control accross bin. Default to \code{mean}.
#'
#' @return A number, the first bin of features to discard.
#'
#' @seealso \code{\link{get_mean_median}}, \code{\link{plot_metric}}
#'
#' @examples
#' myData <- data_frame(
#'     bin = rep(c(1, 2, 3), each = 3),
#'     window = rep(c("top_window", "shuffled_top_window_1", "shuffled_top_window_2"), 3),
#'     mean = c(0.8, 0.1, 0.11, 0.14, 0.12, 0.09, 0.10, 0.13, 0.08)
#' )
#' determine_bin_cutoff(myData)

determine_bin_cutoff <- function(
    metric_table,
    threshold = 2,
    metric = c("mean", "median", "score"),
    random_function_summarisation = mean
) {
    metric <- match.arg(metric)

    metricsTable <- data_frame(
        bin = unique(metric_table$bin),
        top_window = dplyr::filter(metric_table, window == "top_window") %>%
            dplyr::select_(metric) %>%
            unlist,
        ctrl_window_median = dplyr::filter(metric_table, window != "top_window") %>%
            dplyr::rename_(metric = metric) %>%
            dplyr::select(bin, metric) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(med = median(metric)) %>%
            dplyr::select(med) %>%
            unlist,
        ctrl_window_sd = dplyr::filter(metric_table, window != "top_window") %>%
            dplyr::rename_(metric = metric) %>%
            dplyr::select(bin, metric) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(med = sd(metric)) %>%
            dplyr::select(med) %>%
            unlist,
        diff = top_window - ctrl_window_median
    )
    score_th <- threshold * random_function_summarisation(metricsTable$ctrl_window_median)
    keep_until <- 1
    while(
        keep_until < nrow(metricsTable) &&
        metricsTable$top_window[keep_until] > score_th
    ) {
        keep_until <- keep_until + 1
    }
    return(keep_until)
}

