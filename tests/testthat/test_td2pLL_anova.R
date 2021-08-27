

test_that("Data must be a data.frame", {
  expect_error(td2pLL_anova(data = "a"))

  expect_error(td2pLL_anova(data = list(time = 1:10,
                                        dose = 1:10,
                                        resp = 1:10)))
})

test_that("Do not reject too often if there is no time-effect", {

  test_data <- data.frame(time = rep(c(1, 3, 7), each = 12),
                          dose = rep(seq(0, 1, length = 4), each = 3))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 3, gamma = 0, c0 = 0.2, delta = 0.5
  )

  reject_no_effect <- 0
  n_test <- 100
  set.seed(1905)
  for(i in 1:n_test){

    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)

    res <-suppressMessages(td2pLL_anova(curr_data))

    if(res$signif == TRUE)
      reject_no_effect <- reject_no_effect + 1

  }

  expect_lt(reject_no_effect / n_test, 0.1)

})

test_that("Detect time-effect if present", {

  test_data <- data.frame(time = rep(c(1, 3, 7), each = 16),
                          dose = rep(seq(0, 1, length = 4), each = 4))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 2, gamma = 1, c0 = 0.2, delta = 0.3
  )

  # plot_td2pLL(td2pLL_coefs = c(h = 2, gamma = 1, c0 = 0.2, delta = 0.3),
  #            add_ext_data = test_data, xaxis_scale = "linear")

  reject_no_effect <- 0
  n_test <- 100
  set.seed(1905)
  for(i in 1:n_test){

    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)

    plot_td2pLL(td2pLL_coefs = c(h = 2, gamma = 1, c0 = 0.2, delta = 0.3),
                add_ext_data = curr_data, xaxis_scale = "linear")

    res <-suppressMessages(td2pLL_anova(curr_data))

    if(res$signif == TRUE)
      reject_no_effect <- reject_no_effect + 1

  }

  expect_gt(reject_no_effect / n_test, 0.9)

})

