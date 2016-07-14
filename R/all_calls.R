#' Run the method.
#' 
#' The \code{all_calls} function calls all the functions that form the package
#' pipeline, according to the provided parameters.
#' 
#' Arguments for this function are passed on to the each one of the functions in 
#' the package as indicated in this page, as they are called inside the function 
#' body.
#' 
#' \code{randomizations} is passed as an argument of the \code{\link{randomize_top}} 
#' function. Even though more random windows than specified by default (25) can be 
#' generated, it should be noted that bigger numbers may not make random negative
#' control any more informative.
#' 
#' \code{max_zeros} is passed as an argument of the \code{\link{calculate_cvs}} 
#' function. For details on values and meaning, see this object's documentation.
#'  
#' When run, \code{all_calls} also displays useful messages that help the user
#' track analysis progression, as well as information on the operations performed 
#' (i.e. number or size of bins generated, lowest mean expression value included
#' in top window, etc.). If the retrieved values do not match what is expected, the
#' pipeline can be stopped at any point and restarted from the beginning.
#' 
#' By default, the analysis is run creating the \strong{top window} by size, selecting
#' 100 genes, and the \strong{binning} of the rest of the genes is done creating
#' a specific number of windows (30 by default). For details on these methods and
#' the parameters, see \code{\link{extract_top_genes}} and \code{\link{bin_scdata}}.
#' 
#' \strong{Important note}: This function is used to run the analysis implemented 
#' in this package. Unless highly familiar with the method, it is not advised to 
#' run the different functions individually.
#' 
#' @param data A file name, path or data frame containing the expression value table.
#' 
#' @param imort A logical. Indicates whether raw data should be imported \code{data}
#' from a file or path (\code{TRUE}) or is in a frame already loaded into R 
#' (\code{FALSE}).
#' 
#' @param max_zeros A number of type double. Indicates the maximum proportion
#' of zeros allowed per row in the raw data.
#' 
#' @param randomizations An integer. Number of times the top window will be 
#' randomized to generate a random negative control.
#' 
#' @param top_method A string. Indicates the method for top window selection.
#' 
#' @param filter_parameter An integer. Specifies a parameter for the top window 
#' selection method.
#' 
#' @param bin_method A string. Indicates the method for binning of the genes.
#' 
#' @param bin_parameter An integer. Specifies the parameter for the selected 
#' binning method.
#' 
#' @param cor_method A string. Indicates the correlation method to use.
#' 
#' @return A list containing all the results generated during the analysis. Each
#' element on the list contains the results of a function call, and is named as follows:
#' 
#' \itemize{
#' 
#'      \item \code{raw} contains the raw data table as it was imported.
#'      
#'      \item \code{cv} contains the data frame output by \code{calculate_cvs}, which
#'      includes all the data plus the mean, cv and standard deviation columns.
#'      
#'      \item \code{binned} contains all the data as output by the combined call of
#'      \code{extract_top_genes} and \code{bin_scdata}, which includes the bin number
#'      column in the data frame.
#'      
#'      \item \code{window_correlations} contains the output of \code{correlate}, which
#'      outputs a list containing the vector of correlation coefficients of each window 
#'      to the top window as a separate element.
#'      
#'      \item \code{controls} contains a nested list with the controls -as many as 
#'      randomizations- computed for each window, as output by calling \code{randomize_top} and 
#'      \code{compute_control} once per window.
#'      
#'      \item \code{densities} contains a list where the outputs of \code{compute_density}
#'      for each window are stored as a separate element.
#'      
#'      \item \code{histogram} contains a data frame with all the histogram values.
#'      
#'      \item \code{info} contains a string in which all the necessary information about
#'      the analysis is included. This should be used as the histogram title.
#' }
#' 
#' \strong{Note:} use these names to retrieve information by subsetting using \code{$}.

all_calls <- function(data, import = TRUE, max_zeros = 0.5, randomizations = 25,
                      # top window parameters
                      top_method = "window_size", filter_parameter = 100,
                      # binning paramaters
                      bin_method = "window_number", bin_parameter = 30,
                      # correlation method
                      cor_method = "pearson"){
    
    # start counting computing time
    t0 <- Sys.time()

    # load data
    if (import == TRUE){
        # from a file
        raw_data <- read_tsv(data)
        
    }
    if (import == FALSE){
        # from an existent data frame
        raw_data <- data
    }
    message("Data loaded successfully.")

    # calculate CVs
    if (max_zeros == 0.5){

        message("Did not specify maximum proportion of zero expression values allowed: default (0.5) will be used")
    }

    # calculate CVs and report results of the zero value filtering
    message("Calculating CVs... Performing zero value filtering on raw data...")
    CV_data <- calculate_cvs(raw_data, max_zeros)
    message(paste("Genes remaining after zero value filtering:", nrow(CV_data)))

    # extract top genes using selected method and parameter
    divided_data <- extract_top_genes(CV_data, method = top_method,
                                          parameter = filter_parameter)

    # bin rest of data using selected method and parameter
    all_data <- bin_scdata(divided_data, method = bin_method,
                                            parameter = bin_parameter)

    # calculate correlations
    message("Calculating correlations...")
    cor_data <- correlate(all_data, cor_method)
    message("Calculation finished.")

    # create a random control:
        # generate a specified number of random top windows windows
    message("Generating random negative control...")
    randomization_seq = c(1:randomizations)
    randomized_windows <- list()

    for (i in randomization_seq){
        # use the divided dataset to select the top window
        randomized_windows[[i]] <- randomize_top(divided_data)
    }

    all_controls <- list()

    for (i in unique(all_data$bin)){
        all_controls[[i]] <- compute_control(randomized_windows, all_data, i,
                                             cor_method)
    }
    message("Generated random negative control successfully for all windows.")

    # compute and plot correlations against control
    message("Computing correlations of each window against control...")
    message("Note: simultaneously computing histogram values.")

    dens <- list()          # list to store density computations
    histogram <- list()     # histogram values will be computed simultaneously

    for (i in seq_len(length(all_controls))){

        # generate the density data table for each window
        dens[[i]] <- compute_density(all_controls[[i]], cor_data, i)

        # calculate the histogram values for each window
        histogram[[i]] <- compute_histogram(dens[[i]])
    }

    # bind all histogram values
    all_hist_values <- data.frame(do.call(rbind, histogram))
    message("Histogram values computed successfully!")

    # create a string with all the info about the analysis to use as a histogram title
    info <- paste("Top window:", top_method, "=", filter_parameter,
              ", size:", nrow(divided_data$topgenes),
              "| Bins:", max(all_data$bin), ",",
              nrow(subset(all_data, bin == 2)), "genes/bin",
              "| Cells:", ncol(raw_data))

    # store all results in a list
    results <- list(raw = raw_data, cv = CV_data, binned = all_data, 
                    window_correlations = cor_data, controls = all_controls,
                    densities = dens, histogram = all_hist_values, info = info)
    
    # report computing time
    t <- trunc(Sys.time() - t0)
    message(paste("Total computation time:", t, "min."))

    # return all the generated results
    return(results)
}
