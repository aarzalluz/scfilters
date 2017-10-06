# tiny unneccessary utility function
my_summerise_func <- function(my_mat, my_func) {
    my_func(as.vector(my_mat))
}

# suppress CHECK annoying handling of NSE
utils::globalVariables(c("top_window_size"))
#' Utility plot to choose a top_window size
#'
#' Plot mean autocorrelation value of the features of the top window depending
#' on increasing top window size.
#'
#' @param sc_data A tibble, usually the output of \code{link{calculate_cvs}}.
#'
#' @param from Minimum size of the top window.
#'
#' @param to Maximum size of the top window.
#'
#' @param by Size of the steps to walk form \code{from} to \code{to}. See \code{seq}.
#'
#' @param ... Arguments to be passed to \code{cor}, for example \code{method = "spearman"}
#'
#' @return A \code{ggplot2} plot.
#'
#' @examples
#' plot_top_window_autocor(calculate_cvs(scData_hESC))
#'
#' @export

plot_top_window_autocor <- function(sc_data, from = 10, to = 400, by = 2, ...) {

    cor_mat <- dplyr::arrange(sc_data, dplyr::desc(mean)) %>%
        dplyr::slice(seq_len(to)) %>%
        dplyr::select(-geneName, -mean, -sd, -cv) %>%
        as.matrix %>%
        t %>%
        stats::cor()
    cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA
    cor_mat <- abs(cor_mat)

    plot_table <- tibble::tibble(
        top_window_size = seq(from = from, to = to, by = by),
        mean = vapply(
            top_window_size,
            function(x) {
                my_summerise_func(cor_mat[seq_len(x), seq_len(x)], function(y) {
                    mean(y, na.rm = TRUE)
                })
            },
            0
        )
    )

    top_window_plot <- ggplot(plot_table, aes(x = top_window_size)) +
        geom_line(aes(y = mean)) +
        labs(x = "top window size", y = "mean auto-correlation coefficient")

    return(top_window_plot)
}


