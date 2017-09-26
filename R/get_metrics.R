#' Extract mean and median correlation coefficient values
#'
#' Takes the output of \code{\link{correlate_windows}} and extract the mean and the median correlation value
#' for each window comparison.
#'
#' @param df A data frame, usually the output of \code{\link{correlate_windows}}.
#'
#' @param absolute_cc Should the function work of absolute value of correlation coefficients?
#' Default to \code{TRUE} to simplify plots and avoid annoying, non-symmetrical, near 0, shifts of distributions.
#'
#' @return A \code{data_frame} with columns \code{bin}, \code{window}, \code{mean} and \code{median}.
#'
#' @examples
#' expMat <- matrix(
#'     c(1, 1, 5,
#'       1, 2, 3,
#'       0, 1, 4,
#'       0, 0, 2),
#'     ncol = 3, byrow = TRUE, dimnames = list(paste("gene", 1:4), paste("cell", 1:3))
#' )
#'
#' calculate_cvs(expMat) %>%
#'     define_top_genes(window_size = 2) %>%
#'     bin_scdata(window_number = 1) %>%
#'     correlate_windows(n_random = 2) %>%
#'     get_mean_median
#'
#' @export
get_mean_median <- function(df, absolute_cc = TRUE) {
    if (absolute_cc) {
        dplyr::group_by(df, bin, window) %>%
            dplyr::summarise(mean = mean(abs(cor_coef)), median = median(abs(cor_coef))) %>%
            dplyr::ungroup()
    } else {
        dplyr::group_by(df, bin, window) %>%
            dplyr::summarise(mean = mean(cor_coef), median = median(cor_coef)) %>%
            dplyr::ungroup()
    }

}