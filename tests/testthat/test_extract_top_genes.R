library(scFeatureFilter)
context("define_top_genes")

myDF1 <- data.frame(gene = letters[1:6],
                    x1 = c(1, 2, 3, 4, 5, 6),
                    x2 = c(2, 1, 4, 3, 6, 5),
                    stringsAsFactors = TRUE)
input1 <- suppressMessages(calculate_cvs(myDF1))

result1 <- list(
    topgenes = input1[c(5, 6), ],
    restofgenes = input1[c(3, 4, 1, 2), ]
)
result1$topgenes$bin <- 1

result2 <- list(
    topgenes = input1[c(5, 6, 3), ],
    restofgenes = input1[c(4, 1, 2), ]
)
result2$topgenes$bin <- 1

test_that("Giving the expected result", {
    expect_equal(
        suppressMessages(define_top_genes(input1, window_size = 2)),
        result1
    )
    expect_equal(
        suppressMessages(define_top_genes(input1, mean_expression = 4)),
        result1
    )
    expect_equal(
        suppressMessages(define_top_genes(input1, min_expression = 4)),
        result1
    )
})

test_that("Working with ties", {
    expect_equal(
        suppressMessages(define_top_genes(input1, window_size = 3)),
        result2
    )
})