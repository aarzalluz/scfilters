#' Select top highly expressed genes.
#'
#' Selects the group of genes in the dataset that will be considered more highly expressed,
#' the top window, using the filtering method chosen by the user.
#'
#' There are three selection methods available:
#'
#' \itemize{
#'
#'      \item \code{"window_size"}: genes are ranked by mean expression, and a subset of the size
#'      indicated in \code{parameter} is selected from the top.
#'
#'      \item \code{"min_expression"}: genes where all expression values are above a minimum
#'      expression threshold indicated in \code{parameter} are selected.
#'
#'      \item \code{"mean_expression"}: the \code{mean} column is checked, and all genes with mean
#'      expression above the threshold indicated in \code{parameter} are selected.
#' }
#'
#' There are no restrictions to the \code{parameter} argument, however, the value should
#' be coherent with the characteristics of the data set provided and the chosen method.
#'
#' In general, it is adviseable to avoid generating top windows much larger than 250 genes,
#' to prevent excessively long computation time as well as to preserve the quality of the
#' analysis, as the top window should only include a subset of reliable values. As a rule,
#' the bigger the top window is, the more likely is that the reliability of the values is
#' compromised, given the characteristics of single cell RNA sequencing data.
#'
#' @param dataset A data frame, containing genes as rows and cells as columns, and where
#' the mean expression value for each gene has been added as a column.
#'
#' @param method A string indicating the method to use when creating the top window. If no
#' method is indicated, \code{"window_size"} will be used.
#'
#' @param parameter An integer. Indicates the numeric parameter to use in the previously
#' chosen method.
#'
#' @return A list with two elements, both data frames: the generated top window, and
#' the rest of the genes.

extract_top_genes <- function(dataset,
                              method = c("window_size", "min_expression", "mean_expression"),
                              parameter){
    divided_data <- list()
    expr_values <- data.frame()
    sorted_values <- data.frame()

    method <- match.arg(method)

    if (method == "window_size"){

        # the parameter will be the window size
        window_size <- parameter

        # save row names
        dataset$rownames <- rownames(dataset)
        # sort data by mean expression - original row names are lost
        sorted_values <- arrange(dataset, desc(mean))
        # assign row names again
        rownames(sorted_values) <- sorted_values$rownames
        # eliminate the row names column
        sorted_values <- select(sorted_values, -rownames)

        # select the top x genes (x=window size selected)
        divided_data$topgenes <- sorted_values[1:window_size, ]
        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 0

        # display mean FPKM value of the last gene in the top window
        message(paste("Generated top window successfully! Mean expression of last top gene:",
                      sorted_values[window_size, ]$mean))

        # store the rest of the genes a the second element of the list
        divided_data$restofgenes <- sorted_values[-(1:nrow(divided_data$topgenes)), ]

    } else if (method == "min_expression"){

        # the parameter will be the min expression threshold
        min_expr_threshold <- parameter

        # extract expression values
        expr_values <- select(dataset, -mean, -stdev, -CV)

        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset,
                                        rowSums(expr_values > min_expr_threshold) == ncol(expr_values))

        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 0

        # display number of genes in the window
        message(paste("Generated top window successfully! Number of genes in top window:",
                      nrow(divided_data$topgenes)))

        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset,
                                           rowSums(expr_values > min_expr_threshold) != ncol(expr_values))

    } else if (method == "mean_expression"){

        # the parameter will be the mean expression threshold
        mean_threshold <- parameter

        # select top genes and store in first position of the list
        divided_data$topgenes <- subset(dataset, dataset$mean > mean_threshold)

        # assign bin 0 to the top window
        divided_data[[1]]$bin <- 0

        # display number of genes in the window
        message(paste("Generated top window successfully! Number of genes in top window:",
                      nrow(divided_data$topgenes)))

        # store rest of genes in second position
        divided_data$restofgenes <- subset(dataset, dataset$mean < mean_threshold)
    }

    return(divided_data)
}

#' Bin genes by mean expression.
#'
#' Divides the genes that were not included in the top window in windows of the same size,
#' by their mean expression.
#'
#' There are two binning methods available:
#'
#' \itemize{
#'      \item \code{"window_number"}: Divides the genes into the number of windows specified in
#'      \code{parameter}, regardless of their size.
#'
#'      \item \code{"window_size"}: Divides the genes into windows of the size specified in
#'      \code{parameter}, regardless of the number of windows generated.
#' }
#'
#' This function uses the \code{ntile} function, in the \code{dplyr} package to assign a bin
#' number to each gene based on the value contained in the \code{mean} column, corresponding
#' to its mean expression. These are then added as a the column \code{bin} using the \code{mutate}
#' function, also in the \code{dplyr} package.
#'
#' \strong{Important note:} This function is designed to take the list output by the
#' \code{extract_top_window} function as an argument, operating only on the second element
#' of it.
#'
#' Once the genes in it have been binned, both elements of the list are bound
#' together in a data frame and returned. The output is similar, but a new column \code{bin}
#' is added, which indicates the window number assigned to each gene.
#'
#' @param dataset A list, containing the top window generated by \code{extract_top_genes}
#' as the first element, and the rest of undivided genes as the second.
#'
#' @param method A string, indicating the method to be used to bin the genes by mean
#' expression.
#'
#' @param parameter An integer. Indicates the numeric parameter to use in the previously
#' chosen method. Values are not restricted, but should be coherent with the method of choice.
#'
#' @return A data frame containing the binned genes.

bin_scdata <- function(dataset,  method = c("window_number", "window_size"), parameter){

    method <- match.arg(method)

    # save row names
    dataset[[2]]$rownames <- rownames(dataset$restofgenes)

    # classify data into windows
    if (method == "window_number"){

        # bin into the selected number of windows
        windows <- ntile(desc(dataset[[2]]$mean), parameter)

    } else if (method == "window_size"){
        # calculate number of windows of selected window size possible and bin
        windows <- ntile(desc(dataset[[2]]$mean), trunc(nrow(dataset$restofgenes)/parameter))
    }

    # add new column indicating the window in which each gene falls
    dataset$restofgenes <- mutate(dataset$restofgenes, bin = windows)
    # set row names again
    rownames(dataset$restofgenes) <- dataset[[2]]$rownames
    # delete row names column
    dataset$restofgenes <- select(dataset$restofgenes, -rownames)

    if (method == "window_number"){

        # print size of the desired number of windows created
        message(paste("Binned data successfully! Size of windows:", nrow(subset(dataset[[2]], bin == 1))))

    } else if (method == "window_size"){

        # print number of windows of desired size created
        message(paste("Binned data successfully! Number of windows:", max(dataset[[2]]$bin)))
    }

    # bind all the data and correct bin number
    dataset <- do.call("rbind", dataset) %>% as.data.frame()
    dataset$bin <- dataset$bin + 1

    return(dataset)
}
