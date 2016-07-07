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
#' \code{histogram_threshold} is passed as an argument of the \code{\link{filter_data}}
#' function.
#' 
#' \code{all_calls} also includes code use to plot graphs for the analysis. These
#' will be saved in R's working directory, unless otherwise specified. Also, if the 
#' indicated directory does not exist in the current working directory, it will be 
#' created inside it. When unspecified, the function will not save the plots, but 
#' only display them.
#' 
#' Plot aesthetics arguments are used as arguments of the \code{\link[ggplot2]{geom_point}} 
#' function, in the \code{ggplot2} package.
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
#' @param file A file or path containing the expression value table.
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
#' @param scatter_plot,bins_plot,correlation_plots,histogram_plot Logicals. 
#' When \code{TRUE}, plots will be generated through the analysis.
#' 
#' @param save Logical. When \code{TRUE}, generated plots will be saved.
#' 
#' @param display Logical. When \code{TRUE}, generated plots will also be displayed.
#' 
#' @param plot_path A string. Supplies a path to save plots to.
#' 
#' @param format A string. Indicates the file format of plot files.
#' 
#' @param scatter_parameters,bins_parameters Numeric vector providing aesthetic
#' parameters for scatter plots: [1] dot shape, [2] dot size, [3] alpha.
#' 
#' @param density Logical. Indicates whether scatter density should be shown
#' in the scatter plot.
#' 
#' @param cor_method A string. Indicates the correlation method to use.
#' 
#' @param histogram_threshold A number of type double. Provides threshold for
#' final data filtering. 
#' 
#' @return A data frame containing expression values for genes remaining after
#' filtering, output by the \code{\link{filter_data}} function.

all_calls <- function(file, max_zeros = 0.5, randomizations = 25,
                      # top window parameters
                      top_method = "window_size", filter_parameter = 100,
                      # binning paramaters
                      bin_method = "window_number", bin_parameter = 30,
                      save = FALSE, plot_path = NULL, format = "png", display = TRUE,
                      # scatter parameters
                      scatter_plot = TRUE, scatter_parameters = c(16, 0.5, 0.5),
                      density = FALSE,
                      # bins plot parameters
                      bins_plot = TRUE, bins_parameters = c(16, 0.5, 0.5),
                      # correlation plots and histogram
                      correlation_plots = TRUE, cor_method = "pearson",
                      histogram_plot = TRUE, histogram_threshold = 0.1){
    t0 <- Sys.time()

    # make sure that plot_path is not null before its value is tested
    if (is.null(plot_path)){

        plot_path <- getwd()

    # if a directory was specified and it doesn't exist, create it
    } else if (!(plot_path %in% list.files())){

            dir.create(plot_path)
    }

    # load data from specified path
    raw_data <- read_tsv(file)
    message("Data loaded successfully.")

    # calculate CVs
    if (max_zeros == 0.5){

        message("Did not specify maximum proportion of zero expression values allowed: default (0.5) will be used")
    }

    # calculate CVs and report results of the zero value filtering
    message("Calculating CVs... Performing zero value filtering on raw data...")
    CV <- calculate_cvs(raw_data, max_zeros)
    message(paste("Genes remaining after zero value filtering:", nrow(CV)))

    # scatter plot
    if (scatter_plot == TRUE){

        pl <- ggplot(CV, aes(x = CV, y = mean + 1)) +
            # point shape and characteristics
            geom_point(shape = scatter_parameters[1], size = scatter_parameters[2],
                       alpha = scatter_parameters[3]) +
            # set log10 scale in the y axis
            scale_y_log10() +
            ylab("Mean expression + 1") + xlab("Coefficient of variation")

        if(density == TRUE){
            pl <- pl + geom_density_2d()
        }

        # displaying condition
        if (display == TRUE){

            plot(pl)
            message("Generated scatter plot.")
        }

        # saving condition
        if (save == TRUE){

            ggsave(paste("scatter_plot.", format, sep=""), path = plot_path)
            message("Saved scatter plot in selected directory as 'scatter_plot'.")
        }
    }

    # extract top genes using selected method and parameter
    divided_data <- extract_top_genes(CV, method = top_method,
                                          parameter = filter_parameter)

    # bin rest of data using selected method and parameter
    all_data <- bin_scdata(divided_data, method = bin_method,
                                            parameter = bin_parameter)

    # bins plot
    if (bins_plot == TRUE){

        # plot data assigning colors by window
        pl <- ggplot(all_data, aes(x = CV, y = mean + 1, colour = factor(all_data$bin+1))) +
            geom_point(shape = bins_parameters[1], alpha = bins_parameters[2], size = bins_parameters[3]) +
            scale_y_log10() +
            ylab("Mean expression") + xlab("Coefficient of variation")

        # displaying condition
        if (display == TRUE){
            plot(pl)
            message("Generated plot of binned data.")
        }

        # saving condition
        if (save == TRUE){

            ggsave(paste("bins_plot.", format, sep = ""), path = plot_path)
            message("Save binned data plot as 'bins_plot' in selected directory.")
        }
    }

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

    histogram <- list()     # histogram values will be computed simultaneously

    for (i in seq_len(length(all_controls))){

        # generate the density data table for each window
        dens <- compute_density(all_controls[[i]], cor_data, i)

        # calculate the histogram values for each window
        histogram[[i]] <- compute_histogram(dens)

        # plotting condition
        if (correlation_plots == TRUE){

            # create label for each window
            window_label <- paste("Window", i)

            # plot the data
            pl <- ggplot() +
                geom_ribbon(aes(x = dens$control_x, ymin = dens$min_control_y,
                                ymax = dens$max_control_y), color = "grey", alpha = 0.2) +
                geom_line(data = dens, aes(x = control_x, y = control_y, color = "Control")) +
                geom_line(data = dens, aes(x = data_x, y = data_y, color = window_label)) +
                ylab("Density") + xlab("Absolute value of correlation") +
                theme_bw() + scale_color_manual(values = c("black","red"),
                                                labels = c("Control", window_label))

             # displaying condition
            if (display == TRUE){
                plot(pl)
            }

            # saving condition
            if (save == TRUE){

                ggsave(paste("Random negative control - WIndow", i, ".", format, sep = ""),
                    path = plot_path)
                message(paste("Correlation plot", i, "saved in selected directory."))
            }

            if (display == TRUE){
                message("See correlation plots displayed.")
            }
        }
    }

    # bind all histogram values
    all_hist_values <- data.frame(do.call(rbind, histogram))
    message("Histogram values computed successfully!")

    # plot the histogram
    if (histogram_plot == TRUE){
        
        # generate plot title
        hist_title <- paste("Top window:", top_method, "=", filter_parameter,
              ", size:", nrow(divided_data$topgenes),
              "| Bins:", max(all_data$bin), ",",
              nrow(subset(all_data, bin == 2)), "genes/bin",
              "| Cells:", ncol(raw_data))
        
        pl <- ggplot(all_hist_values, aes(x = factor(as.numeric(rownames(all_hist_values))),
                                          y = hist_value, fill = "")) +
            geom_bar(position = position_dodge(), stat="identity") + 
            ylab("Proportion of correlated genes") + xlab("Window number") +
            ggtitle(hist_title) + ylim(0, 0.8) +
            geom_errorbar(aes(ymax = error_up, ymin = error_down), width = 0.2)

        # display condition
        if (display == TRUE){
            plot(pl)
            message("Plotted histogram.")
        }

        # saving condition
        if (save == TRUE){

            ggsave(paste("histogram.", format, sep=""), path = plot_path)
            message("Saved histogram as 'histogram' in selected directory.")
        }
    }

    # filter the data using the histogram values
    message("Filtering data...")
    quality_filtered_data <- filter_data(all_data, all_hist_values, histogram_threshold)
    message("Data was filtered according to set threshold. Analysis finished.")

    t <- trunc(Sys.time() - t0)
    message(paste("Total computation time:", t, "min."))

    # return the filtered data
    return(quality_filtered_data)
}
