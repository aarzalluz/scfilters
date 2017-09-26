#' Filter scRNA-seq expression matrix to keep only highly informative features. Integrated pipeline.
#'
#' This pipeline function takes an expression matrix as an input and filter it to
#' keep only the most highly expressed features, with a threshold estimated through
#' correlations with a reference bin constituted of the very most expressed features
#' (see the vignette for details on the method).
#'
#' The function can optionally produce three plots of \code{print_plots} is \code{TRUE}.
#' It is recommanded to open a graphical device (i.e. through \code{pdf} or \code{png}),
#' to call \code{scFeatureFilter},and then to close the device with \code{dev.off}.
#'
#' @param sc_data A data frame or a matrix, containing expression values for each gene as rows, and
#' expression values for the cells as columns.
#'
#' @param print_plots A boolean. Should the function produce three plots as a side effect?
#' Plots are the output of \code{\link{plot_mean_variance}}, \code{\link{plot_correlations_distributions}}
#' and \code{\link{plot_metric}}.
#'
#' @param max_zeros A number between 0 and 1. Maximum proportion of cells with 0 expression
#' for a feature to be kept.
#'
#' @param threshold A number higher than 1. The higher the more stringent the feature selection
#' will be. See \code{\link{determine_bin_cutoff}}.
#'
#' @param top_window_size Size of the reference bin. See \code{\link{define_top_genes}}
#'
#' @param other_window_size Size of the other bins of feature. See \code{\link{bin_scdata}}
#'
#' @return A \code{matrix} or a \code{tibble}, depending on the type of \code{sc_data},
#' containing only the tope expressed features.
#'
#' @examples
#' sc_feature_filter(sc_data)
#'
#' # with plots
#' \dontrun{
#' pdf("diagnostic.pdf")
#' sc_feature_filter(sc_data, print_plots = TRUE)
#' dev.off()
#' }

sc_feature_filter <- function(
    sc_data,
    print_plots = FALSE,
    max_zeros = 0.75,
    threshold = 2,
    top_window_size = 100,
    other_window_size = 1000
){
    binned_data <- sc_data %>%
        calculate_cvs(max_zeros = max_zeros) %>%
        define_top_genes(window_size = top_window_size) %>%
        bin_scdata(window_size = other_window_size)

    cor_dis <- correlate_windows(binned_data, n_random = 3)

    metrics <- get_mean_median(cor_dis)

    if(is.matrix(sc_data)) {
        is_matrix <- TRUE
    } else {
        is_matrix <- FALSE
    }

    filtered_data <- filter_expression_table(
        binned_data,
        bin_cutoff = determine_bin_cutoff(metrics, threshold = threshold),
        as_matrix = is_matrix
    )

    if(print_plots) {
        print(
            plot_mean_variance(binned_data) + annotation_logticks(sides = "l")
        )
        print(
            plot_correlations_distributions(
                correlations_to_densities(cor_dis),
                metrics = metrics
            )
        )
        print(
            plot_metric(metrics, threshold = threshold)
        )
    }

    return(filtered_data)
}
