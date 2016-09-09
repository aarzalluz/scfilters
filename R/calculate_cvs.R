#' Calculate coefficient of variation of genes.
#'
#' Calculates the coefficient of variation (CV) of each row in the supplied data table.
#'
#' Before CV computation, the function 
#' removes all rows that have a proportion of zeros above the specified
#' threshold. Genes with many 0s are poorly informative, and would bias the later correlations.
#' Removing them also prevents division by zero when calculating CVs.
#'
#' The data provided must contain gene names as a column in the first position, as output by the
#' \code{read_tsv} function in the \code{dplyr} package, and an cell names as column names.
#' These names will be assigned as the rownames of the output data frame.
#'
#' In the output, mean, standard deviation and CV are incorporated as new columns in the data
#' frame, named \code{mean}, \code{sd} and \code{CV}.
#'
#' @param data A data frame, containing expression values for each gene as rows, and
#' expression values for the cells as columns.
#'
#' @param max_zeros A double indicating the maximum proportion of zero expression values
#' allowed per row. The value of this argument must be \eqn{0 => max_zeros > 1}
#'
#' @return A data frame, containing the filtered data and the mean, standard deviation and cv
#' values for each row.

calculate_cvs <- function(data, max_zeros = 0.75){

    if(is.data.frame(data)) {
        data <- .createGeneExpressionMatrixFromDataFrame(data)
    }

    # look for negative expression values
    if(any(data < 0)) {
        stop(
            "Genes cannot have a negative expression value."
        )
    }
    
    if(max_zeros < 0 | max_zeros >= 1) {
        stop("max_zeros should be between 0 (included) and 1 (excluded).")
    }

    # remove genes with bigger proportion of zero values than indicated threshold
    data_final <- subset(data, (rowSums(data == 0) / ncol(data)) <= max_zeros)

    # calculate mean, st dev and cv of rows
    mean <- rowMeans(data_final, na.rm = T)
    stdev <- apply(data_final, 1, sd, na.rm = T)
    CV <- stdev/mean

    return(as_data_frame(cbind(
        geneName = rownames(data_final),
        mean = mean,
        sd = stdev,
        cv = CV,
        data_final
    )))
}
