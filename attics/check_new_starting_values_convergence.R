
# Quick check: Does the fitting work for each compound?

all_compounds <- cytotox$compound %>% unique

for(compound in all_compounds){
  print(compound)
  data_subset <- cytotox[cytotox$compound == compound, c("expo", "dose", "resp")]
  colnames(data_subset)[1] <- "time"
  fit <- fit_td2pLL(data = data_subset)
  print(plot(fit, add_data = data_subset, zaxis_title = compound))
}

# First problem: MePA (20th compound)

for(compound in all_compounds[21:length(all_compounds)]){
  print(compound)
  data_subset <- cytotox[cytotox$compound == compound, c("expo", "dose", "resp")]
  colnames(data_subset)[1] <- "time"
  fit <- fit_td2pLL(data = data_subset)
  print(plot(fit, add_data = data_subset, zaxis_title = compound))
}



# apparently, only problem with MePA

data <- cytotox[cytotox$compound == "MePA", c("expo", "dose", "resp")]
colnames(data)[1] <- "time"
coefs <- td2pLL:::get_starting_values(data = data) %>% unlist()

plot_td2pLL(td2pLL_coefs = coefs, add_ext_data = data,
            dose_lim = c(0.001, max(data$dose)), time_lim = range(data$time))


# does it work now?


compound <. "MePA"
data_subset <- cytotox[cytotox$compound == compound, c("expo", "dose", "resp")]
colnames(data_subset)[1] <- "time"
fit <- fit_td2pLL(data = data_subset)
print(plot(fit, add_data = data_subset, zaxis_title = compound, xaxis_scale = "linear"))

