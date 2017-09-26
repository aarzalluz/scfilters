# take a bottom window, a top window, and return a vector of correlation values
.correlate_window <- function(topMatrix, bottomMatrix, ...) {
    correlations <- list()
    for(j in seq_len(nrow(topMatrix))){
        correlations[[j]] <- cor(topMatrix[j,], t(bottomMatrix), ...) %>% as.vector
    }
    return(do.call(c, correlations))
}

#' Calculate correlations against top window.
#'
#' Calculates correlations of each genes in each window against each genes in the top window.
#'
#' This function:
#' \itemize{
#'      \item correlates each gene in each window to each gene in the top window.
#'
#'      \item randomize the top window by shuffling expression value, and correlate each gene in each
#'      window to the randomized top window. This negative control is repeated as many time as specified by
#'      the \code{n_random} parameter.
#' }
#' The input of this function is usually the output of the \code{\link{bin_scdata}} function.
#'
#' @param dataset A data frame containing all the binned genes. Usually the output of \code{\link{bin_scdata}}.
#'
#' @param n_random Number of top window randomization to serve as a negative control. Default to 3.
#'
#' @param ... Additional arguments to be passed to \code{\link{cor}}. Default method is \code{pearson}
#'  which is the fastest.
#'
#' @return A \code{tibble} containing correlation values.
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
#'     correlate_windows
#'
#' @export
correlate_windows <- function(dataset, n_random = 3, ...){

    # extract the top window genes
    top_window <- dplyr::filter(dataset, bin == 1) %>%
        dplyr::select(-geneName, -mean, -sd, -cv, -bin) %>%
        as.matrix
    if(ncol(top_window) <= 1) stop("Needs more than 1 cell!")
    shuffled_top_windows <- lapply(
        seq_len(n_random),
        function(x) t(apply(top_window, 1, sample))
    )

    # iterate bins in the dataset
    corTable <- dplyr::bind_rows(
        lapply(
            unique(dataset$bin),
            function(i) {
                # select the genes in the chosen window using the bin number
                selected_window <- dplyr::filter(dataset, bin == i) %>%
                    dplyr::select(-geneName, -mean, -sd, -cv, -bin) %>%
                    as.matrix

                with_top_window <- tibble::tibble(
                    bin = i,
                    window = "top_window",
                    cor_coef = .correlate_window(top_window, selected_window)
                )

                with_controls <- dplyr::bind_rows(
                    lapply(
                        seq_len(n_random),
                        function(j) {
                            tibble::tibble(
                                bin = i,
                                window = paste0("shuffled_top_window_", j),
                                cor_coef = .correlate_window(shuffled_top_windows[[j]], selected_window)
                            )
                        }
                    )
                )

                return(dplyr::bind_rows(with_top_window, with_controls))
            }

        )
    )

    return(corTable)
}
