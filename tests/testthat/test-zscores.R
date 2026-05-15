test_that("in-sample z-scores have mean near 0 and sd near 1", {
  set.seed(123)
  x <- runif(300, 10, 20)
  y <- 8 + 2 * sin(x) + rnorm(300, 0, x / 11)
  fit <- lmsqreg.fit(y, x, maxit = 25)

  expect_true(fit[[1]]$converged)

  z <- fit[[1]]$finalz
  expect_lt(abs(mean(z)), 0.1)
  expect_lt(abs(sd(z) - 1), 0.15)
})

test_that("out-of-sample zscores() has mean near 0 and sd near 1", {
  set.seed(123)
  x <- runif(300, 10, 20)
  y <- 8 + 2 * sin(x) + rnorm(300, 0, x / 11)
  fit <- lmsqreg.fit(y, x, maxit = 25)

  set.seed(456)
  xnew <- runif(150, 10, 20)
  ynew <- 8 + 2 * sin(xnew) + rnorm(150, 0, xnew / 11)
  z <- zscores(ynew, xnew, fit)

  expect_lt(abs(mean(z)), 0.2)
  expect_lt(abs(sd(z) - 1), 0.25)
})

test_that("non-convergence warning mentions maxit", {
  set.seed(123)
  x <- runif(300, 10, 20)
  y <- 8 + 2 * sin(x) + rnorm(300, 0, x / 11)
  expect_warning(
    lmsqreg.fit(y, x, maxit = 3),
    regexp = "maxit"
  )
})
