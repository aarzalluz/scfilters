context("filter expression table")

myData <- data_frame(
    bin = rep(c(1, 2, 3), each = 3),
    mean = 9:1,
    sd = runif(9),
    cv = runif(9),
    cell1 = 8:0 + runif(9),
    cell2 = 8:0 + runif(9)
)

test_that("filter_expression_table", {
    expect_equal(
        filter_expression_table(myData, bin_cutoff = 2),
        myData[1:3, c("cell1", "cell2")]
    )
    expect_equal(
        filter_expression_table(myData, bin_cutoff = 3),
        myData[1:6, c("cell1", "cell2")]
    )
})
