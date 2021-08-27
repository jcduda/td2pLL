

test_that("Data must be a data.frame", {
  expect_error(TDR(data = "a"))
})

test_that("No time-dependency: Always plot result if strict_stop = FALSE", {

  test_data <- data.frame(time = rep(c(1, 3, 7), each = 12),
                          dose = rep(seq(0, 1, length = 4), each = 3))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 3, gamma = 0, c0 = 0.2, delta = 0.5
  )

  n_test <- 100
  set.seed(1905)
  for(i in 1:n_test){

    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)

    res <- suppressMessages(TDR(curr_data))

    expect_error(plot(res$fit), NA)

  }

})

test_that("With time-dependency: Always plot result if strict_stop = FALSE", {

  test_data <- data.frame(time = rep(c(1, 3, 7), each = 16),
                          dose = rep(seq(0, 1, length = 4), each = 4))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 2, gamma = 1, c0 = 0.2, delta = 0.3
  )

  n_test <- 100
  set.seed(1905)
  for(i in 1:n_test){

    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)

    res <- suppressMessages(TDR(curr_data))

    expect_error(plot(res$fit, xaxis_scale = "linear"), NA)

  }

})


