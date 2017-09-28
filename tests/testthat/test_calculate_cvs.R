context("calculate_cvs")

myDF1 <- data.frame(gene = letters[1:6], x1 = 0:5, x2 = c(0, 5:1), stringsAsFactors = TRUE )
myDF2 <- data.frame(gene = letters[1:6], x1 = 0:5, x2 = c(0, 5:1), stringsAsFactors = FALSE)
myMatrix <- matrix(c(0:5, 0, 5:1), ncol = 2, dimnames = list(letters[1:6], c("x1", "x2")))
myTibble <- tibble::tibble(gene = letters[1:6], x1 = 0:5, x2 = c(0, 5:1))

result <- tibble::tibble(
    geneName = letters[2:6],
    mean = rep(3, 5),
    sd = c(sd(c(1, 5)), sd(c(2, 4)), 0, sd(c(2, 4)), sd(c(1, 5))),
    cv = sd / 3,
    x1 = seq(1.0, 5.0, by = 1.0),
    x2 = seq(5.0, 1.0, by = -1.0)
)

test_that("giving the right result", {
    expect_equal(
        calculate_cvs(myDF1),
        result
    )
    expect_equal(
        calculate_cvs(myDF2),
        result
    )
    expect_equal(
        calculate_cvs(myMatrix),
        result
    )
    expect_equal(
        calculate_cvs(myTibble),
        result
    )
})


