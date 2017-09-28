#' Filter binned expression matrix
#'
#' Takes a binned expression table (the output of \code{\link{bin_scdata}}), a bin number
#' (usually the output of \code{\link{determine_bin_cutoff}}) and returned a filtered expression
#' table or matrix.
#'
#' @param bined_table A \code{tibble}, usually the output of \code{\link{bin_scdata}}.
#'
#' @param bin_cutoff the number of the first bin to be filtered out. Can be the
#' output of \code{\link{determine_bin_cutoff}}).
#'
#' @param as_matrix A boolean. Should the return be a \code{tibble} (\code{FALSE}, the default) or
#' a \code{matrix} (\code{TRUE}).
#'
#' @return A \code{tibble} or a \code{matrix} depending on the value of \code{as_matrix}
#'
#' @seealso \code{\link{bin_scdata}}, \code{\link{determine_bin_cutoff}}
#'
#' @examples
#' myData <- tibble::data_frame(
#'     bin = rep(c(1, 2, 3), each = 3),
#'     mean = 9:1,
#'     sd = runif(9),
#'     cv = runif(9),
#'     cell1 = 8:0 + runif(9),
#'     cell2 = 8:0 + runif(9)
#' )
#' filter_expression_table(myData, bin_cutoff = 2)
#' filter_expression_table(myData, bin_cutoff = 3)
#'
#' @export
filter_expression_table <- function(bined_table, bin_cutoff, as_matrix = FALSE) {

    if (!is.null(bin_cutoff)) {
        bined_table <- dplyr::filter(bined_table, bin < bin_cutoff)
    }

    bined_table <-  dplyr::select(bined_table, -mean, -sd, -cv, -bin)

    if (as_matrix) {
        geneName <- bined_table$geneName
        bined_table <- dplyr::select(bined_table, -geneName) %>%
            as.matrix
        rownames(bined_table) <- geneName
    }

    return(bined_table)
}

