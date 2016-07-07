#' Calculate density of correlation.
#'
#' Computes the density of correlation of the random negative control and the data correlations.
#'
#' The density data is returned only for the selected window, therefore, \code{compute_density}
#' needs to be called once per window.
#'
#' For both the control and the data, \code{compute_density} converts correlations to absolute
#' values before computing their density.
#'
#' Density is computed using the \code{density} function, from the \code{stats} package. However,
#' this function returns a density object, therefore, x and y values of density are subset
#' using \code{$}.
#'
#' Since there is more than one control correlation vector for a single window, given the
#' necessity to randomize the top window several times, the mean of density y values will
#' be returned instead. Standard error of the mean is also calculated and added and substracted
#' from the mean to compute the maximum and minimum estimated density.
#'
#' It should be noted that all control densities are assumed to provide approximately similar
#' x values, and therefore, only those for the first control density are returned.
#'
#' All these, together with x and y density values for the data correlations, are bound in
#' a data frame and returned.
#'
#' Column names of the returned data frame and their meaning:
#'
#' \itemize{
#'
#'      \item \code{control_x}: x values of the density object, taken from density of one
#'      the random window's correlations.
#'
#'      \item \code{control_y}: mean y values of density for all random correlations.
#'
#'      \item \code{min_control_y}: \code{control_y} minus standard error of mean.
#'      \item \code{max_control_y}: \code{control_y} plus standard error of mean
#'
#'      \item \code{data_y}: y values of the density object returned after computing the density
#'      of the window correlations.
#'
#'      \item \code{data_x}: x values of the density object returned after computing the density
#'      of the window correlations.
#' }
#'
#'
#' @param control A list containing random control correlation vectors for a selected window.
#'
#' @param correlations A list containing correlations of each window against the top window.
#'
#' @param window_number An integer indicating the number of the window for which the density
#' is to be computed.
#'
#' @return A data frame containing all density data from the control and the window correlations
#' as named columns.

compute_density <- function(control, correlations, window_number){

    # calculate densities for the control correlations
    control_dens <- list()

    for (i in seq_len(length(control))){
        control_dens[[i]] <- data.frame(density(abs(control[[i]]), cut = 0, from = 0, to = 1)$x,
                                        density(abs(control[[i]]), cut=0, from = 0, to = 1)$y)
        colnames(control_dens[[i]]) <- c("x", "y")
    }

    # extract y column of each element in the control_dens list (lapply)
    dens_y <- do.call(rbind, lapply(control_dens, function(my_dens) my_dens$y))

    # calculate standard deviation of the data
    sd_y <- apply(dens_y, 2, sd)
    sd_y <- sd_y/sqrt(nrow(dens_y))

    # calculate the column wise mean of all the density curves
    mean_control_dens <- colMeans(dens_y)

    # calculate min and max control curves
    max_control_y <- mean_control_dens + sd_y
    min_control_y <- mean_control_dens - sd_y

    # compute density for actual data
    data_density <- density(abs(correlations[[window_number]]), from = 0, to = 1, cut = 0)

    # create a label for the data line using the given window number
    window_label <- paste("Window", window_number)

    # bind all in one data frame and return
    density_data <- data.frame(control_x = control_dens[[1]]$x, control_y = mean_control_dens,
                               min_control_y = min_control_y, max_control_y = max_control_y,
                               data_y = data_density$y, data_x = data_density$x)

    return(density_data)
}
