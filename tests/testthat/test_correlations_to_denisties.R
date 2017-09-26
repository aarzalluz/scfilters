library(scfilters)
context("correlation to density")

myData1 <- data_frame(
    bin = rep(c(1, 2), each = 4),
    window = "top_window",
    cor_coef = c(1, 1, 1, 1, 0, 0, 1, 1)
)

myData2 <- data_frame(
    bin = rep(1, each = 3),
    window = "top_window",
    cor_coef = c(-0.5, 0.5, 1)
)

test_that("correlations_to_densities", {
    # n is working
    expect_equal(
        correlations_to_densities(myData1, n = 3)$cor_coef,
        c(0, 0.5, 1, 0, 0.5, 1)
    )
    expect_equal(
        correlations_to_densities(myData1, n = 2)$cor_coef,
        c(0, 1, 0, 1)
    )
    # absolute_cc is working
    expect_equal(
        correlations_to_densities(myData2, n = 2)$cor_coef,
        c(0.0, 1.0)
    )
    expect_equal(
        correlations_to_densities(myData2, n = 2, absolute_cc = FALSE)$cor_coef,
        c(-1.0, 1.0)
    )

    # weak test of density columns
    expect_is(
        correlations_to_densities(myData1)$density,
        "numeric"
    )
    expect_is(
        correlations_to_densities(myData2)$density,
        "numeric"
    )

})
