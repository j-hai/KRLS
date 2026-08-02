## First differences for binary predictors must land on the binary predictor's
## own column when cat_columns has one-hot expanded a factor.
##
## fdskrls() locates binary columns in object$X (the raw input) but writes into
## derivatives / avgderivatives / binaryindicator, which are in the processed
## kernel design. Once a factor is expanded the two spaces differ in width, so
## raw positions used as processed positions land on the wrong column: with the
## factor listed before the binary treatment, the treatment's first difference
## was stored under the first factor dummy.

make_cat_data <- function(n = 200, seed = 77L) {
  set.seed(seed)
  g <- factor(sample(c("a", "b", "c"), n, TRUE))
  d <- rbinom(n, 1, 0.5)
  y <- 2.1 * d + c(a = 0, b = 1, c = -1)[as.character(g)] + rnorm(n, sd = 0.3)
  list(g = g, d = d, y = y)
}

test_that("the first difference is attributed to the binary column, not a factor dummy", {
  dat <- make_cat_data()
  fit <- krls(X = data.frame(g = dat$g, d = dat$d), y = dat$y,
              cat_columns = "g", print.level = 0)

  flagged <- colnames(fit$avgderivatives)[as.logical(fit$binaryindicator)]
  expect_identical(flagged, "d")

  # and the stored value is the treatment's first difference, near the truth
  expect_equal(unname(fit$avgderivatives[1, "d"]), 2.1, tolerance = 0.15)
})

test_that("results do not depend on the column order of X", {
  dat <- make_cat_data()
  a <- krls(X = data.frame(g = dat$g, d = dat$d), y = dat$y,
            cat_columns = "g", print.level = 0)
  b <- krls(X = data.frame(d = dat$d, g = dat$g), y = dat$y,
            cat_columns = "g", print.level = 0)

  cols <- c("d", "ga", "gb", "gc")
  expect_equal(a$avgderivatives[1, cols], b$avgderivatives[1, cols],
               tolerance = 1e-8)
  expect_identical(as.logical(a$binaryindicator), as.logical(b$binaryindicator))
})

test_that("attribution is correct for an unnamed matrix with numeric cat_columns", {
  # X supplied as a matrix has no column names, so the raw/processed spaces have
  # to be reconciled by origin index rather than by name.
  set.seed(303)
  n <- 200
  gcode <- sample(1:3, n, TRUE)
  d <- rbinom(n, 1, 0.5)
  y <- 2.1 * d + c(0, 1, -1)[gcode] + rnorm(n, sd = 0.3)

  M_gd <- cbind(gcode, d); dimnames(M_gd) <- NULL
  M_dg <- cbind(d, gcode); dimnames(M_dg) <- NULL
  a <- krls(X = M_gd, y = y, cat_columns = 1L, print.level = 0)
  b <- krls(X = M_dg, y = y, cat_columns = 2L, print.level = 0)

  # in both fits the continuous block comes first, so the binary column is
  # processed column 1 regardless of where it sat in the raw input
  expect_identical(which(as.logical(a$binaryindicator)), 1L)
  expect_identical(which(as.logical(b$binaryindicator)), 1L)
  expect_equal(unname(a$avgderivatives[1, 1]), 2.1, tolerance = 0.15)
  expect_equal(unname(a$avgderivatives[1, ]), unname(b$avgderivatives[1, ]),
               tolerance = 1e-8)
})

test_that("a two-level factor in cat_columns is skipped rather than misattributed", {
  # Such a column expands to two one-hot columns and so has no single processed
  # counterpart. It must not claim a first difference on an arbitrary column.
  set.seed(404)
  n <- 150
  sex <- factor(sample(c("m", "f"), n, TRUE))
  x1 <- rnorm(n)
  y <- 1.2 * (sex == "m") + 0.5 * x1 + rnorm(n, sd = 0.3)

  fit <- krls(X = data.frame(sex = sex, x1 = x1), y = y,
              cat_columns = "sex", print.level = 0)
  expect_false(any(as.logical(fit$binaryindicator)))
})

test_that("several factors plus a binary treatment stay order-invariant", {
  set.seed(505)
  n <- 150
  g1 <- factor(sample(c("a", "b", "c"), n, TRUE))
  g2 <- factor(sample(c("p", "q", "r"), n, TRUE))
  d <- rbinom(n, 1, 0.5)
  y <- 2.1 * d + c(a = 0, b = 1, c = -1)[as.character(g1)] +
    c(p = 0, q = 0.5, r = -0.5)[as.character(g2)] + rnorm(n, sd = 0.3)

  a <- krls(X = data.frame(g1 = g1, g2 = g2, d = d), y = y,
            cat_columns = c("g1", "g2"), print.level = 0)
  b <- krls(X = data.frame(d = d, g1 = g1, g2 = g2), y = y,
            cat_columns = c("g1", "g2"), print.level = 0)

  expect_identical(colnames(a$avgderivatives)[as.logical(a$binaryindicator)], "d")
  expect_identical(colnames(b$avgderivatives)[as.logical(b$binaryindicator)], "d")
  cols <- c("d", "g1a", "g1b", "g1c", "g2p", "g2q", "g2r")
  expect_equal(a$avgderivatives[1, cols], b$avgderivatives[1, cols],
               tolerance = 1e-8)
})

test_that("a fit with no categorical columns is unaffected", {
  set.seed(11)
  n <- 120
  x1 <- rnorm(n)
  d <- rbinom(n, 1, 0.5)
  y <- 1.5 * d + 0.5 * x1 + rnorm(n, sd = 0.3)
  # cat_columns = integer(0) acknowledges the binary d column, as elsewhere
  fit <- krls(X = cbind(x1 = x1, d = d), y = y, cat_columns = integer(0),
              print.level = 0)

  flagged <- colnames(fit$avgderivatives)[as.logical(fit$binaryindicator)]
  expect_identical(flagged, "d")
})
