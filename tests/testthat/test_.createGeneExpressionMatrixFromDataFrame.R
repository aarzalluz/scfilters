library(scFeatureFilter)
context(".createGeneExpressionMatrixFromDataFrame")

test_that("working with different types:", {
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data.frame(a = letters[1:5], b = 1:5)
        ),
        matrix(1:5, dimnames = list(letters[1:5], "b"))
    )
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data.frame(a = letters[1:5], b = 1:5, stringsAsFactors = FALSE)
        ),
        matrix(1:5, dimnames = list(letters[1:5], "b"))
    )
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data_frame(a = letters[1:5], b = 1:5)
        ),
        matrix(1:5, dimnames = list(letters[1:5], "b"))
    )
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data.frame(b = 1:5)
        ),
        matrix(1:5, dimnames = list(NULL, "b"))
    )
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data.frame(b = 1:5, stringsAsFactors = FALSE)
        ),
        matrix(1:5, dimnames = list(NULL, "b"))
    )
    expect_equal(
        .createGeneExpressionMatrixFromDataFrame(
            data_frame(b = 1:5)
        ),
        matrix(1:5, dimnames = list(NULL, "b"))
    )
})

