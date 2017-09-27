#' Define the top highly expressed gene window.
#'
#' Define the group of genes in the dataset that will be considered as reference,
#' the top window, by specifying either a number of genes or an expression threshold.
#'
#' There are three selection methods available:
#'
#' \itemize{
#'
#'      \item \code{window_size}: genes are ranked by mean expression across cells, and the top slice
#'      of the specified size is selected.
#'
#'      \item \code{mean_expression}: the \code{mean} column is checked, and all genes with mean
#'      expression above the threshold indicated are selected.
#'
#'      \item \code{min_expression}: genes where all expression values are above the
#'      expression threshold indicated are selected.
#'
#' }
#'
#' In general, it is advisable to avoid generating top windows larger than 250 genes
#' (100 genes is the recommended value),
#' to prevent excessively long computation time as well as to preserve the quality of the
#' analysis, as the top window should only include a subset of reliable values.
#'
#' @param dataset A data frame, containing genes as rows and cells as columns, and where
#' the mean expression value for each gene has been added as a column. Usually the output of
#' \code{calculate_cvs}.
#'
#' @param window_size Number of genes in the defined top window. Recommended to 100 genes.
#'
#' @param mean_expression A number. Genes with a mean expression across cells higher than the value
#' will be selected. Ignored if \code{window_size} is defined.
#'
#' @param min_expression A number. Genes with a minimum expression across all cells higher than the value
#' will be selected. Ignored if \code{window_size} or \code{mean_expression} is defined.
#'
#' @return A list with two elements, both data frames: the defined top window, and
#' the rest of the genes.
#'
#' @examples
#' library(magrittr)
#' expMat <- matrix(
#'     c(1, 1, 1,
#'       1, 2, 3,
#'       0, 1, 2,
#'       0, 0, 2),
#'     ncol = 3, byrow = TRUE, dimnames = list(paste("gene", 1:4), paste("cell", 1:3))
#' )
#'
#' calculate_cvs(expMat) %>%
#'     define_top_genes(window_size = 2)
#'
#' calculate_cvs(expMat) %>%
#'     define_top_genes(mean_expression = 1.5)
#'
#' @export
#' @importFrom stats cor density median sd
define_top_genes <- function(dataset,
                              window_size = NULL,
                              mean_expression = NULL,
                              min_expression = NULL
                              ){
    divided_data <- list()
    expr_values <- data.frame()
    sorted_values <- data.frame()

    if (is.null(window_size) & is.null(mean_expression) & is.null(min_expression)) stop("Need to provide one of the following parameter: window_size, mean_expression or min_expression")

    if (is.numeric(window_size)){

        # sort data by mean expression - original row names are lost
        sorted_values <- dplyr::arrange(dataset, dplyr::desc(mean))
        # select the top x genes (x=window size selected)
        divided_data$topgenes <- sorted_values[1:window_size, ]
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display mean FPKM value of the last gene in the top window
        message(paste("Mean expression of last top gene:",
                      sorted_values[window_size, ]$mean))
        # store the rest of the genes a the second element of the list
        divided_data$restofgenes <- sorted_values[-(1:nrow(divided_data$topgenes)), ]

    } else if (is.numeric(mean_expression)){

        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset, dataset$mean > mean_expression)
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display number of genes in the window
        message(paste("Number of genes in top window:",
                      nrow(divided_data$topgenes)))
        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset, dataset$mean < mean_expression)

    } else if (is.numeric(min_expression)){

        # extract expression values
        expr_values <- dplyr::select(dataset, -geneName, -mean, -sd, -cv)
        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset,
                                        rowSums(expr_values > min_expression) == ncol(expr_values))
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 1
        # display number of genes in the window
        message(paste("Number of genes in top window:",
                      nrow(divided_data$topgenes)))
        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset,
                                           rowSums(expr_values > min_expression) != ncol(expr_values))

    } else {
        stop("The second parameter should be numeric.")
    }

    return(divided_data)
}

# suppress CHECK annoying handling of NSE
utils::globalVariables(c("geneName", "cv", "bin"))
#' Bin genes by mean expression.
#'
#' Divides the genes that were not included in the top window in windows of the same size with decreasing
#' mean expression levels.
#'
#' Two binning methods are available:
#'
#' \itemize{
#'      \item \code{window_number}: Divides the genes into the number of windows specified.
#'
#'      \item \code{window_size}: Divides the genes into windows of the size specified.
#' }
#'
#' This function adds a bin number column to the data frame.
#'
#' This function is designed to take the list output by the
#' \code{extract_top_window} function as an argument, operating only on the second element
#' of it. Once the genes in it have been binned, both elements of the list are bound
#' together in a data frame and returned. The output contains a new column \code{bin},
#' which indicates the window number assigned to each gene.
#'
#' @param dataset A list, containing the top window generated by \code{extract_top_genes}
#' as the first element, and the rest of undivided genes as the second. Usually the output
#' of \code{define_top_genes}
#'
#' @param window_number An integer, indicating the number of bins to be used.
#'
#' @param window_size An integer, indicating the number of genes to be included
#'  in each window. Ignored if \code{window_size} is defined.
#'
#' @param verbose A boolean. Should the function print a message about window size
#' or the number of windows created?
#'
#' @return A data frame containing the binned genes.
#'
#' @examples
#' library(magrittr)
#' expMat <- matrix(
#'     c(1, 1, 1,
#'       1, 2, 3,
#'       0, 1, 2,
#'       0, 0, 2),
#'     ncol = 3, byrow = TRUE, dimnames = list(paste("gene", 1:4), paste("cell", 1:3))
#' )
#'
#' calculate_cvs(expMat) %>%
#'     define_top_genes(window_size = 1) %>%
#'     bin_scdata(window_number = 2)
#'
#' calculate_cvs(expMat) %>%
#'     define_top_genes(window_size = 1) %>%
#'     bin_scdata(window_size = 1)
#'
#' @export
bin_scdata <- function(dataset, window_number = NULL, window_size = NULL, verbose = TRUE){

    if (is.null(window_number) & is.null(window_size)) stop("Need to provide window_number or window_size.")

    if (is.numeric(window_number)){
        # bin into the selected number of windows
        windows <- dplyr::ntile(dplyr::desc(dataset[[2]]$mean), window_number) + 1
    } else if (is.numeric(window_size)){
        # calculate number of windows of selected window size possible and bin
        windows <- dplyr::ntile(dplyr::desc(dataset[[2]]$mean), trunc(nrow(dataset$restofgenes)/window_size)) + 1
    } else {
        stop("The second parameter should be numeric.")
    }
    dataset$restofgenes <- dplyr::mutate(dataset$restofgenes, bin = windows)

    if (verbose) {
        if (is.numeric(window_number)){
            # print size of the desired number of windows created
            message(paste("Window size:", length(which(dataset[[2]]$bin == 2))))
        } else if (is.numeric(window_size)){
            # print number of windows of desired size created
            message(paste("Number of windows:", max(dataset[[2]]$bin)))
        }
    }

    # bind all the data and correct bin number
    dataset <- dplyr::bind_rows(dataset) %>%
        dplyr::select(geneName, mean, sd, cv, bin, dplyr::everything())
    return(dataset)
}
