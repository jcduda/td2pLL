library(rds)
library(dplyr)
data <- readRDS("C:/Users/duda/Projekte/td2pLL/data/cytotox.rds")

names(data) <- tolower(names(data))
colnames(data)[6] <- "sample_id"

cytotox <- data

save(cytotox, file ="C:/Users/duda/Projekte/td2pLL/data/cytotox.RData")

cytotox$compound %>% unique %>% length
# 30 compounds

cytotox %>% select(compound, donor) %>% distinct %>%
  group_by(compound) %>% tally() %>% select(n) %>% table()
# always 3 donors per compound

cytotox %>% select(compound, donor, expo) %>% distinct %>%
  group_by(compound, donor) %>% tally %>% ungroup %>% select(n) %>% table()
# for each donor-compund combination: 3 different exposure times

cytotox %>% select(compound, donor, expo, dose) %>% distinct() %>%
  group_by(compound, donor, expo) %>% tally() %>% ungroup() %>% select(n) %>% table()
# per compound, donor and exposure time: mostly 6 doses (240 cases), rarely
# 7 doses (27 cases) or 8 doses (3 cases)

cytotox %>% group_by(compound) %>% tally() %>% ungroup %>% select(n) %>% table()
# total number of observations per compound hence varies between 214 and 276

cytotox %>% group_by(compound, donor, expo, dose) %>% tally() %>%
  ungroup() %>% select(n) %>% table()
# for at each dose-time level (given a compound and donor), there are
# mostly 3 measurements (1536 cases) and sometimes 8 (114 cases).
# Rarely, there are just 2 (1 case) or 3 (2 cases) measurements.
# The 8 measurements are alwways for dose = 0,
# whereas sometimes, for dose = 0, there are also just 4 measurements:

cytotox %>% group_by(compound, donor, expo, dose) %>% tally() %>%
  ungroup() %>% filter(dose == 0) %>% select(n) %>% table()

# sample_id: enumerates the replicates within
# compound, donor, expo and dose

# control_mean: within compound, donor and expo: a control_mean is calculated
# as mean of the raw responses at dose=0.

# resp = raw_resp / control_mean

# the control_mean is the arithmetic mean of raw_resp measurements for
# a compound, donor and time combination.
# Sometimes, there are two sets of such control measurements, leading to
# two control_mean's for such a setting.
# by_control indiciates, which control_mean is then used of the two for the
# normalization.

# by_control: factor with levels "0", "1" and "2"
#     - "0": that means that for a compound, donor expo combination,
#               there was just one set of control-measurements that are used
#               for all doses for normalization
#     - "1": There are two sets of control measurements fot a compound,
#               donor and expo combination and for this observation,
#               the first one (control_1 in original excel files) is used
#     - "2": There are two sets of control measurements for a compound,
#               donor and expo combination and for this observation,
#               the second  one (control_2 in original excel files) is used


# The already pre-processed data  will go thourgh
# a secondary, optional pre-processing step, the refitting.
# The new responses are stored in resp_refit.

#
# The refitting algorithm does the follwing:
#
# For a given compound, exposure time and donor, a dose-response curve (4pLL)
# is fitted.
# The corresponding data is divided through the resulting left (upper) asymptote
# (left_asymp)
# and again multiplied by 100.
# This way, the data is expected to better follow the assumption of having a
# left asymptote at 100 [percent].

# left_asymp: Left asymptote calculated in refitting.
# resp_refit: Response values after refitting.

