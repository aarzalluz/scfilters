#' Compute mean expression levels, standard deviations and coefficients of variation of each feature.
#'
#' Compute mean expression levels, standard deviations and coefficients of variation (CV)
#' of each feature (i.e. gene or transcript) in the supplied data.
#' Filter features with high proportion of 0 expression.
#'
#' Before CV computation, the function
#' removes all rows that have a proportion of zeros above the specified
#' threshold. Genes with many 0s are poorly informative, and would bias the later correlations.
#' Removing them also prevents division by zero when calculating CVs.
#'
#' The data provided must cell/sample names as column names. Feature name can be given either in the
#' first column or as rownames.
#'
#' In the output, mean, standard deviation and CV are incorporated as new columns in the data
#' frame, named \code{mean}, \code{sd} and \code{cv}.
#'
#' @param data A data frame, a matrix or a \code{SingleCellExperiment} object.
#' If data frame or matrix, it should contain expression values for each gene as rows, and
#' expression values for the cells as columns.
#'
#' @param max_zeros A number between 0 and 1 indicating the maximum proportion of zero expression values
#' allowed per row. Features with a higher proportion of 0 will be discarded.
#'
#' @param sce_assay, if \code{data} is an \code{SingleCellExperiment} object,
#' \code{sce_assay} should be one of \code{names(assays(<SingleCellExperiment>))}.
#'
#' @return A data frame, containing the filtered data with additional columns: mean, standard deviation and cv
#' values for each row.
#'
#' @examples
#' expMat <- matrix(
#'     c(1, 1, 1,
#'       1, 2, 3,
#'       0, 1, 2,
#'       0, 0, 2),
#'     ncol = 3, byrow = TRUE, dimnames = list(paste("gene", 1:4), paste("cell", 1:3))
#' )
#' calculate_cvs(expMat)
#' calculate_cvs(expMat, max_zeros = 0.5)
#'
#' @export
calculate_cvs <- function(data, max_zeros = 0.75, sce_assay = NULL){

    if(is.data.frame(data)) {
        data <- .createGeneExpressionMatrixFromDataFrame(data)
    } else if (suppressWarnings(requireNamespace("SingleCellExperiment")) &&
               methods::is(data, "SingleCellExperiment")) {
        if(is.null(sce_assay)) {
            stop("sce_assay needs to be specified for SingleCellExperiment objects.
                 \nShould be one of names(assays(<SingleCellExperiment>))")
        }
        data <- SummarizedExperiment::assay(data, sce_assay)
    }

    # look for negative expression values
    if(any(data < 0)) {
        warning(
            "Genes should not have a negative expression value."
        )
    }

    if(max_zeros < 0 | max_zeros > 1) {
        stop("max_zeros should be between 0 and 1.")
    }

    # remove genes with bigger proportion of zero values than indicated threshold
    data_final <- subset(data, (rowSums(data == 0) / ncol(data)) <= max_zeros)

    # calculate mean, st dev and cv of rows
    mean <- rowMeans(data_final, na.rm = TRUE)
    stdev <- apply(data_final, 1, sd, na.rm = TRUE)
    CV <- stdev/mean

    return(dplyr::bind_cols(
        tibble::tibble(
            geneName = rownames(data_final),
            mean = mean,
            sd = stdev,
            cv = CV
        ),
        tibble::as_tibble(data_final)
    ))
}
