#' Calculate coefficient of variation of genes.
#'
#' Calculates the coefficient of variation (CV) of each row in the supplied data table.
#'
#' Before CV computation, the function modifies the raw data in two ways:
#'
#' \itemize{
#'      \item First, it replaces all negative values in the data frame for zeros, as negative
#'      expression values are systematic artifacts caused when genes have very low expression,
#'      and can interfere in CV calculation.
#'
#'      \item Then, it removes all rows that have a proportion of zeros above the specified
#'      threshold. These will be considered highly disrupted by systematic noise. Removing
#'      them prevents division by zero when calculating CVs, which could lead to CV = 0 and to
#'      having missing values in the correlation vectors.
#'
#'      The value of this argument must be \eqn{0 > max_zeros > 1}.
#' }
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
#' allowed per row.
#'
#' @return A data frame, containing the filtered data and the mean, standard deviation and cv
#' values for each row.

calculate_cvs <- function(data, max_zeros){

    # convert dataset provided to data frame
    data <- data.frame(data)

    # use first column as row names and eliminate it to keep only numeric values
    data_labeled <- data.frame(data[, -1], row.names = make.names(data[, 1], unique = T))

    # look for negative expression values and set to zero
    m <- as.matrix(data_labeled)
    m[m < 0] <- 0
    data_labeled <- as.data.frame(m)

    # calculate proportion of zero expression values per gene (row)
    # and remove genes with bigger proportion of zero values than indicated threshold
    data_final <- subset(data_labeled, (rowSums(data_labeled == 0) / ncol(data_labeled)) <= max_zeros)

    # calculate mean and st dev of rows and add them as columns to dataframe
    data_final$mean <- rowMeans(data_final, na.rm = T)
    data_final$stdev <- apply(data_final, 1, sd, na.rm = T)

    # calculate CVs
    data_final$CV <- data_final$stdev/data_final$mean

    # return full dataset with new columns for mean, st dev and CV
    return(data_final)
}
