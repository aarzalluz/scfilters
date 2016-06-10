#' Randomize the top window.
#'
#' Randomly sorts the values of each row of the supplied object.
#'
#' \code{randomize_top} generates the result of a single randomization. To create a consistent
#' negative control, this function should be called a sufficient number of times, using the
#' average of all randomizations as such.
#'
#' Randomization takes place one row at a time, sorting the column values and returning the
#' randomized row to a newly created matrix. Since they are no longer relevant, due to the
#' sorting, column names are removed.
#'
#' \strong{Important note:} This function uses the output of \code{\link{extract_top_genes}}
#' as an input, or any list containing the top subset of genes as the first element. If that
#' is the case, it should be noted that it will be subset as \code{$topgenes}.
#'
#' @param dataset A list containing the top window data frame as the first element.
#'
#' @return A matrix containing a sample of randomized values of the top window.

randomize_top <- function(dataset){

    # extract only expression values from the topgenes data frame ("$" subsetting of the dataset)
    top_expression <- data.frame(dataset$topgenes)
    top_expression <- select(top_expression, -mean, -stdev, -CV, -bin)

    # create a new data frame to store the randomized top window
    random_window <- list()

    # iterate the top window row wise, shuffle it and store it in the equivalent row
    # of the empty data frame
    for(i in seq_len(nrow(top_expression))){

        random_window[[i]] <- top_expression[i,][,sample(ncol(top_expression))]
        colnames(random_window[[i]]) <- seq_len(ncol(random_window[[i]]))
    }

    # bind all elements in the list to create a matrix for all the randomized window
    random_window <- do.call(rbind, random_window)

    # return a randomized window
    return(random_window)
}
