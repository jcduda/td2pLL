

test_that("Data must be a data.frame", {
  expect_error(fit_td2pLL(data = "a"),
               'data must be a numeric data.frame with colnames "time", "dose" and "resp"')
})

test_that("At least two doses and times", {
  test_data = cytotox[cytotox$compound == "ASP" & cytotox$expo == 2, c("expo", "dose", "resp")]
  colnames(test_data)[1] <- "time"
  expect_error(fit_td2pLL(data = test_data), "Data must contain at least 3 different times.")
  test_data = cytotox[cytotox$compound == "ASP" & cytotox$dose == 0, c("expo", "dose", "resp")]
  colnames(test_data)[1] <- "time"
  expect_error(fit_td2pLL(data = test_data), "Data must contain at least two different doses.")
})

test_that("Correct column names of data frame", {
  test_data = cytotox[cytotox$compound == "ASP",]
  expect_error(fit_td2pLL(test_data),
               'data must be a numeric data.frame with colnames "time", "dose" and "resp"')
})



test_that("Convergence for small linear dose grid and no time-effect", {
  test_data <- data.frame(time = rep(c(1, 3, 7), each = 33),
                          dose = rep(seq(0, 1, 0.1), each = 3))
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
      suppressWarnings(fit_td2pLL(curr_data))
    },
    error = function(cond) return(0))

    if(!is.list(res)) non_conv_count <- non_conv_count + 1
  }

  expect_equal(non_conv_count, 0)

})


test_that("Convergence for big log dose grid and no time-effect", {
  test_data <- data.frame(time = rep(c(1, 3, 7), each = 33),
                          dose = rep(c(0, sqrt(10)^(-3:6)), each = 3))
  test_data$resp <- td2pLL(
    time = test_data$time,
    dose = test_data$dose,
    h = 3, gamma = 0, c0 = 0.2, delta = 0.5
  )

  # take a look, if you want:
  # plot_td2pLL(td2pLL_coefs = c(h = 3, gamma = 0, c0 = 0.2, delta = 0.5),
  #             add_ext_data = test_data)

  non_conv_count <- 0

  set.seed(1905)
  for(i in 1:100){
    curr_data <- test_data
    curr_data$resp <- pmax(curr_data$resp + rnorm(nrow(curr_data), sd = 10), 0)
    res <- tryCatch({
      suppressWarnings(fit_td2pLL(curr_data))
    },
    error = function(cond) return(0))

    if(!is.list(res)) non_conv_count <- non_conv_count + 1
  }

  expect_equal(non_conv_count, 0)

})


