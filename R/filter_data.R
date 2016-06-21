#' Quality filtering of data.
#'
#' Filters out of the binned data matrix those windows that have a quality metric value
#' beneath a set threshold.
#'
#' \code{filter_data} compares the quality metric in every window with the threshold value, and
#' finds out which windows are above threshold. Then, it subsets the binned data, keeping
#' only the selected windows.
#'
#' The output data frame contains gene as rows and cells as columns, but not any other
#' information concerning the data, since \code{mean}, \code{sd}, \code{bin} and \code{CV}
#' columns are removed before it is returned.
#'
#' @param binned_data A data frame, containing the binned genes.
#'
#' @param histogram A data frame, containing the quality metric for each window.
#'
#' @param threshold A number of type double indicating the quality metric value used to
#' filter the data.
#'
#' @return A data frame, containing only genes from windows above set minimum quality.

filter_data <- function(binned_data, histogram, threshold){

    # create empty vector to store selected windows
    windows_above <- vector()
    # add window number as a separate column
    histogram$window <- as.numeric(rownames(histogram))

    # iterate histogram values and compare to threshold
    for (i in seq_len(nrow(histogram))){

        if (histogram$hist_value[i] > threshold){
            # store window number of windows that are above threshold
            windows_above <- c(windows_above, histogram$window[i])

        } else {
            break   # stop iterating as soon as a value above threshold is found
        }
    }

    i_list <- vector()

    # store index of the rows that contain genes from the selected windows
    for (i in seq_len(nrow(binned_data))){

        if (binned_data[i,]$bin %in% windows_above){
            i_list <- c(i_list, i)
        }
    }

    # use stored row numbers to retrieve genes and store them in a data frame
    filtered_data <- data.frame(select(binned_data[i_list,], -mean, -CV, -stdev, -bin))

    # tear the first part of the rowname string to leave only ENSEMBL ID
    rownames <- rownames(filtered_data)
    rownames <- gsub("topgenes.", "", rownames)
    rownames <- gsub("restofgenes.", "", rownames)
    rownames(filtered_data) <- rownames

    return(filtered_data)
}
