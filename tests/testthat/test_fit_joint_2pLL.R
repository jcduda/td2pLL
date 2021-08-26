


test_that("Data must be a data.frame", {
  expect_error(fit_joint_2pLL(data = "a"))
})

test_that("At least two doses", {
  test_data = cytotox[cytotox$compound == "ASP" & cytotox$dose == 0, c("expo", "dose", "resp")]
  expect_error(fit_joint_2pLL(data = test_data), "Data must contain at least two different doses.")
})

test_that("Convergence if sparse, noisy and no time-effect", {
  test_data <- data.frame(time = rep(c(1, 3, 7), each = 12),
                          dose = rep(seq(0, 1, length = 4), each = 3))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 3, gamma = 0, c0 = 0.2, delta = 0.5
  )

  non_conv_count <- 0

  set.seed(1905)
  for(i in 1:100){
    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)
    res <- tryCatch({
      suppressWarnings(fit_joint_2pLL(curr_data))
    },
    error = function(cond) return(0))

    if(!is.list(res)) non_conv_count <- non_conv_count + 1
  }

  expect_equal(non_conv_count, 0)
})


