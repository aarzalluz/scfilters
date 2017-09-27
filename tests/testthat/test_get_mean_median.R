context("Extract metrics")

myData1 <- data_frame(
    bin = rep(c(1, 1, 2), each = 3),
    window = rep(c("top_window", "shuffled_top_window_1", "top_window"), each = 3),
    cor_coef = c(0, 0.5, 1, 0, 0.6, 1, 0, 0.7, 1)
)

expectedTbl <- data_frame(
    bin = c(1, 1, 2),
    window = c("shuffled_top_window_1", "top_window", "top_window"),
    mean = c(1.6/3, 0.5, 1.7/3),
    median = c(0.6, 0.5, 0.7)
)

test_that("get_mean_median", {
    expect_equal(
        get_mean_median(myData1),
        expectedTbl
    )
})