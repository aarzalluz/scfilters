context("bin_scdata")

preinput <- suppressMessages(
    data.frame(gene = letters[1:6],
               x1 = c(1, 2, 3, 4, 5, 6),
               x2 = c(2, 1, 4, 3, 6, 5),
               stringsAsFactors = TRUE) %>%
        calculate_cvs
)
input <- suppressMessages(define_top_genes(preinput, window_size = 2))
result <- dplyr::arrange(preinput, desc(mean)) %>%
    dplyr::mutate(bin = rep(seq(1.0, 3.0, by = 1.0), each = 2)) %>%
    dplyr::select(geneName, mean, sd, cv, bin, dplyr::everything())

test_that("Giving the expected result", {
    expect_equal(
        suppressMessages(bin_scdata(input, window_number = 2)),
        result
    )
    expect_equal(
        suppressMessages(bin_scdata(input, window_size = 2)),
        result
    )
})

