# suppress CHECK annoying handling of NSE
utils::globalVariables(c("cor_coef", "ctrl_window_sd"))
#' Produce a mean expression x coefficient of variation scatter plot.
#'
#' Use the output of \code{\link{calculate_cvs}} or \code{\link{bin_scdata}} and plot a feature
#' mean expression x coefficient of variation scatter plot. Mean expression is represented as
#' \code{log10(mean + 1)}. Each dot represents a feature.
#' Means and coefficient of variations were obtained across single cells.
#' Optionally, colours each dot according to the defined bins of features.
#' Optionally, adds a density2d geom.
#'
#' @param df A \code{tibble}, usually the output of \code{\link{calculate_cvs}} or \code{\link{bin_scdata}}.
#'
#' @param density A boolean. Should a \code{density2d} geom be added to the plot?
#'
#' @param colourByBin A boolean. Should feature be coloured by bin? Need a \code{bin} column in \code{df}
#'  (i.e. the output of \code{\link{bin_scdata}}).
#'
#' @param density_color Colour of the density2d curves.
#'
#' @param ... Further arguments are passed to \code{geom_point} such as \code{size}.
#'
#' @return A ggplot2 plot.
#'
#' @seealso \code{\link{calculate_cvs}}, \code{\link{bin_scdata}}
#'
#' @examples
#' library(magrittr)
#' scData_hESC %>%
#'    calculate_cvs %>%
#'    plot_mean_variance(colourByBin = FALSE)
#'
#' scData_hESC %>%
#'    calculate_cvs %>%
#'    define_top_genes(window_size = 100) %>%
#'    bin_scdata(window_size = 1000) %>%
#'    plot_mean_variance
#'
#' @export
#' @import ggplot2
#' @importFrom magrittr "%>%"
plot_mean_variance <- function(df, density = TRUE, colourByBin = TRUE, density_color = "blue", ...){
    if(colourByBin) {
        pl <- ggplot(df, aes(x = cv, y = mean + 1, colour = factor(bin))) +
            labs(x ="Coefficient of variation", y = "Mean expression + 1" , color = "Bin")

    } else {
        pl <- ggplot(df, aes(x = cv, y = mean + 1)) +
            labs(x ="Coefficient of variation", y = "Mean expression + 1")
    }

    pl <- pl +
        geom_point(...) +
        scale_y_log10()
    # add scatter density to plot
    if (density == TRUE) {
        pl <- pl + geom_density_2d(color = density_color)
    }

    return(pl)
}


#' Produce a density plot of correlation values for each window of feature
#'
#' Feature by feature correlation values between every windows and the reference to window
#' of features are visualized as density lines, one facet per comparison.
#' Two density lines are drown in each facets:
#' \itemize{
#'      \item A thin colored line, the correlations between the bin and the reference top bin of features
#'      \item A thicker blue line with grey error area, the correlations between the bin and the \strong{randomized}
#'      top bin of features. The lines are not shown if \code{n_random = 0} in \code{\link{correlate_windows}}.
#' }
#'
#' @param df A \code{tibble}, usually the output of \code{\link{correlations_to_densities}}.
#'
#' @param metrics Optional. The output of \code{\link{get_mean_median}}. Dashed line will represent
#' mean or median of the correlation coefficient distributions.
#'
#' @param vlines A string, either "mean" or "median". Should the dashed line represent the mean or the median
#' of the correlation coefficient distributions? Ignored if \code{metrics} is \code{NULL}.
#'
#' @param facet_ncol The number of columns to arrange the plots.
#'
#' @return A ggplot2 plot.
#'
#' @seealso \code{\link{correlations_to_densities}}, \code{\link{get_mean_median}}
#'
#' @examples
#' library(magrittr)
#' myData <- scData_hESC %>%
#' calculate_cvs %>%
#'     define_top_genes(window_size = 100) %>%
#'     bin_scdata(window_size = 1000)
#'
#' corDistrib <- correlate_windows(myData, n_random = 3)
#'
#' corDens <- correlations_to_densities(corDistrib)
#'
#' plot_correlations_distributions(corDens)
#'
#' metrics <- get_mean_median(corDistrib)
#'
#' plot_correlations_distributions(corDens, metrics = metrics)
#'
#' @export
plot_correlations_distributions <- function(df, metrics = NULL, vlines = c("mean", "median" ), facet_ncol = 4) {

    vlines <- match.arg(vlines)

    pl <- ggplot() +
        geom_smooth(data = dplyr::filter(df, window != "top_window"),
                    aes(x = cor_coef, y = density), span = 0.8) +
        geom_line(  data = dplyr::filter(df, window == "top_window"),
                    aes(x = cor_coef, y = density, color = factor(bin))) +
        facet_wrap(~bin, labeller = "label_both", ncol = facet_ncol) +
        labs(x = "correlation coefficient") +
        guides(color = FALSE)
    if (!is.null(metrics)) {

        metrics <- dplyr::select_(metrics, "bin", "window", vlines) %>%
            dplyr::rename(metric = !!vlines) %>%
            dplyr::mutate(window = sub("_[0-9]+$", "", window)) %>%
            dplyr::group_by(bin, window) %>%
            dplyr::summarise(metric = median(abs(metric)))

        pl <- pl +
            geom_vline(
                data = dplyr::filter(metrics, window == "top_window"),
                aes(xintercept = metric, color = factor(bin)),
                linetype = "dashed"
            ) +
            geom_vline(
                data = dplyr::filter(metrics, window != "top_window"),
                aes(xintercept = metric),
                linetype = "dashed"
            )
    }

    return(pl)
}


#' Produce a bar chart of mean (or median) correlation coefficient per bin of feature.
#'
#' Use the output of \code{\link{get_mean_median}} and produce a bar chart of mean
#' (or median) correlation coefficient per bin of features. Correlations against the
#' randomised top window are shown as dot-and-whiskers, and are used to estimate a
#' background level.
#'
#' @param metric_table A \code{tibble}, usually the output of \code{\link{get_mean_median}}.
#'
#' @param selected_metric Which column in \code{metricsTable} to use? Default to \code{mean}.
#'
#' @param show_ctrl A boolean. Should a dashed line indicate the estimated background level?
#'
#' @param control_color The colour of the background dashed line (default to blue).
#'
#' @param show_threshold A boolean. Should a dashed line indicate the estimated threshold level?
#'
#' @param threshold How many times the background level should be multiplies do determine a threshold?
#' Default to 2. The higher the more stringent.
#'
#' @param threshold_color The colour of the threshold dashed line (default to blue).
#'
#' @param line_size Thickness of the dashed lines.
#'
#' @param annotate_lines A boolean. Should the dashed lines be annotated?
#'
#' @return A ggplot2 plot.
#'
#' @seealso \code{\link{get_mean_median}}
#'
#' @examples
#' library(magrittr)
#' scData_hESC %>%
#'     calculate_cvs %>%
#'     define_top_genes(window_size = 100) %>%
#'     bin_scdata(window_size = 1000) %>%
#'     correlate_windows(n_random = 3) %>%
#'     get_mean_median %>%
#'     plot_metric
#'
#' @export
plot_metric <- function(
    metric_table,
    selected_metric = c("mean", "median", "score"),
    show_ctrl = TRUE,
    control_color = "blue",
    show_threshold = TRUE,
    threshold = 2,
    threshold_color = "red",
    line_size = 1,
    annotate_lines = TRUE

) {
    selected_metric <- match.arg(selected_metric)
    eq_selected_metric <- enquo(selected_metric)

    metricsTable <- tibble::tibble(
        bin = unique(metric_table$bin),
        top_window = dplyr::filter(metric_table, window == "top_window") %>%
            dplyr::select(!!eq_selected_metric) %>%
            unlist,
        ctrl_window_median = dplyr::filter(metric_table, window != "top_window") %>%
            dplyr::rename(metric = !!selected_metric) %>%
            dplyr::select(bin, metric) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(med = median(metric)) %>%
            dplyr::select(med) %>%
            unlist,
        ctrl_window_sd = dplyr::filter(metric_table, window != "top_window") %>%
            dplyr::rename(metric = !!selected_metric) %>%
            dplyr::select(bin, metric) %>%
            dplyr::group_by(bin) %>%
            dplyr::summarise(med = sd(metric)) %>%
            dplyr::select(med) %>%
            unlist,
        diff = top_window - ctrl_window_median
    )

    p <- ggplot(metricsTable, aes(x = bin, fill = factor(bin))) +
        geom_bar(aes(y = top_window), stat = "identity")

    if(show_ctrl) {
        p <- p +
            geom_hline(yintercept = mean(metricsTable$ctrl_window_median), linetype = "dashed",
                       color = control_color, size = line_size)
        if(annotate_lines) {
            p <- p +
                annotate("text", label = "Background", x = 0.6, y = mean(metricsTable$ctrl_window_median) + 0.01,
                         color = control_color, hjust = 0, vjust = 0)
        }
    }
    if(show_threshold) {
        p <- p +
            geom_hline(yintercept = threshold * mean(metricsTable$ctrl_window_median), linetype = "dashed",
                       color = threshold_color, size = line_size)
        if(annotate_lines) {
            p <- p +
                annotate("text", label = paste0(threshold, " x Background"),
                         x = 0.6, y = threshold * mean(metricsTable$ctrl_window_median) + 0.01,
                         color = threshold_color, hjust = 0, vjust = 0)
        }
    }

    p <- p +
        geom_errorbar(
            aes(ymin = ctrl_window_median - ctrl_window_sd, ymax = ctrl_window_median + ctrl_window_sd),
            width = 0.5, color = "gray30"
        ) +
        geom_point(aes(y = ctrl_window_median), color = "gray30") +
        guides(fill = FALSE) +
        labs(y = paste0("Correlation coefficient (", selected_metric, ")"))

    return(p)
}