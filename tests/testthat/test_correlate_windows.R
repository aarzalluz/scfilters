context("correlate_windows")

topM <- matrix(c(1, 2, 3, 1, 1, 2), byrow = TRUE, ncol = 3)
botM <- matrix(c(1, 2, 3, 1, 3, 2, 1, 3, 2), byrow = TRUE, ncol = 3)

result_pearson <- c(1, 0.5, 0.5, 0.86602540378, 0, 0)
result_kendall <- c(1, 1/3, 1/3, 0.81649658092, 0, 0)

test_that(".correlate_window", {
    expect_equal(
        .correlate_window(topM, botM, method = "spearman"),
        result_pearson
    )
    expect_equal(
        .correlate_window(topM, botM, method = "kendall"),
        result_kendall
    )
})

input <- suppressMessages(
    data.frame(gene = letters[1:5],
               x1 = c(4, 4, 1, 1, 1),
               x2 = c(8, 4, 2, 3, 3),
               x3 = c(12, 8, 3, 2, 2),
               stringsAsFactors = TRUE) %>%
        calculate_cvs %>%
        define_top_genes(window_size = 2) %>%
        bin_scdata(window_number = 2)
)

test_that("correlate_window", {
    expect_equal(
        dim(correlate_windows(input, n_random = 2)),
        c(30, 3)
    )
    expect_equal(
        unique(correlate_windows(input, n_random = 2)$bin),
        c(1, 2, 3)
    )
    expect_equal(
        unique(correlate_windows(input, n_random = 2)$window),
        c("top_window", "shuffled_top_window_1", "shuffled_top_window_2")
    )
    expect_is(
        correlate_windows(input, n_random = 2)$cor_coef,
        "numeric"
    )
    expect_equal(
        dplyr::filter(correlate_windows(input, n_random = 2), bin == 1 & window == "top_window")$cor_coef,
        c(1, 0.86602540378, 0.86602540378, 1)
    )
    expect_equal(
        dplyr::filter(correlate_windows(input, n_random = 2), bin == 2 & window == "top_window")$cor_coef,
        c(1, 0.5, 0.86602540378, 0)
    )
    expect_equal(
        dplyr::filter(correlate_windows(input, n_random = 2), bin == 3 & window == "top_window")$cor_coef,
        c(0.5, 0)
    )
})
