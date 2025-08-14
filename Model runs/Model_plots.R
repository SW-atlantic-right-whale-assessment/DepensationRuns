library(StateSpaceSIR)
library(EnvStats)
source("R/plot_density.R", echo=TRUE)

################################################################################
# Base model ----
################################################################################
file_name <- "Model runs/Base/Base"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base[[1]],  file_name = file_name)
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))
plot_abs_abundance(sir_base[[1]],  file_name = file_name)
plot_density(SIR = sir_base,  file_name = file_name, posteriors_lwd = c(3,1), posteriors_lty = rep(1, 2), posteriors_col = c(1,1))
plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_base[[1]],  file_name = file_name)


################################################################################
# Sensitivity 1 - lognormal prior on Rmax  ----
################################################################################
file_name <- "Model runs/Sensitivity_1/Sensitivity_1"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_1[[1]],  file_name = file_name)
plot_trajectory(sensitivity_1[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_1[[1]], sensitivity_1[[2]], sir_base[[1]], sensitivity_1[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_1[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_1[[1]],  file_name = file_name)


################################################################################
# Sensitivity 2 - smaller CV on lognromal prior on Rmax  ----
################################################################################
file_name <- "Model runs/Sensitivity_2/Sensitivity_2"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_2[[1]],  file_name = file_name)
plot_trajectory(sensitivity_2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_2[[1]], sensitivity_2[[2]], sir_base[[1]], sensitivity_2[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_2[[1]],  file_name = file_name)

################################################################################
# sensitivity_3 ----
################################################################################
# - Prior on rmax of lognormal(-2.65, 0.5)T(0.2, 0.11)
# - Using rlnormTrunc fron "EnvStats" package
file_name <- "Model runs/sensitivity_3/sensitivity_3"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_3[[1]],  file_name = file_name)
plot_trajectory(sensitivity_3[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_3[[1]], sensitivity_3[[2]], sir_base[[1]], sensitivity_3[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_3[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_3[[1]],  file_name = file_name)



################################################################################
# Sensitivity_4 ----
################################################################################
file_name <- "Model runs/Sensitivity_4/Sensitivity_4"
# Upper bound of process error is 100 * lower bound
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_4[[1]],  file_name = file_name)
plot_trajectory(sensitivity_4[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_4[[1]], sensitivity_4[[2]], sir_base[[1]], sensitivity_4[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_4[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_4[[1]],  file_name = file_name)



################################################################################
# Sensitivity_5 ----
################################################################################
file_name <- "Model runs/sensitivity_5/sensitivity_5"
# Upper bound of process error variance is 5 times lower bound
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_5[[1]],  file_name = file_name)
plot_trajectory(sensitivity_5[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_5[[1]], sensitivity_5[[2]], sir_base[[1]], sensitivity_5[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_5[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_5[[1]],  file_name = file_name)


################################################################################
# sensitivity_6 ----
################################################################################
file_name <- "Model runs/sensitivity_6/sensitivity_6"
# Nrecent is 2004
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_6[[1]],  file_name = file_name)
plot_abs_abundance(sensitivity_6[[1]],  file_name = file_name)
plot_trajectory(sensitivity_6[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_6[[1]], sensitivity_6[[2]], sir_base[[1]], sensitivity_6[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1), target = FALSE) 
plot_ioa(sensitivity_6[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_6[[1]],  file_name = file_name)


################################################################################
# sensitivity_7 ----
################################################################################
file_name <- "Model runs/sensitivity_7/sensitivity_7"
# No struck and loss rates
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_7[[1]],  file_name = file_name)
plot_trajectory(sensitivity_7[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_7[[1]], sensitivity_7[[2]], sir_base[[1]], sensitivity_7[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_7[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_7[[1]],  file_name = file_name)



################################################################################
# sensitivity_8 ----
################################################################################
file_name <- "Model runs/sensitivity_8/sensitivity_8"
# - Catch time series is only low
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_8[[1]],  file_name = file_name)
plot_trajectory(sensitivity_8[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_8[[1]], sensitivity_8[[2]], sir_base[[1]], sensitivity_8[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_8[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_8[[1]],  file_name = file_name)



################################################################################
# sensitivity_9 ----
################################################################################
file_name <- "Model runs/sensitivity_9/sensitivity_9"
# -- Catch time series is high
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_9[[1]],  file_name = file_name)
plot_trajectory(sensitivity_9[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_9[[1]], sensitivity_9[[2]], sir_base[[1]], sensitivity_9[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_9[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_9[[1]],  file_name = file_name)


################################################################################
# sensitivity_10 min haplotypes = 0 ----
################################################################################
file_name <- "Model runs/sensitivity_10/sensitivity_10"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_10[[1]],  file_name = file_name)
plot_trajectory(sensitivity_10[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_10[[1]], sensitivity_10[[2]], sir_base[[1]], sensitivity_10[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_10[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_10[[1]],  file_name = file_name)


################################################################################
# sensitivity_11 min haplotypes = 25 ----
################################################################################
file_name <- "Model runs/sensitivity_11/sensitivity_11"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_11[[1]],  file_name = file_name)
plot_trajectory(sensitivity_11[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_11[[1]], sensitivity_11[[2]], sir_base[[1]], sensitivity_11[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_11[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_11[[1]],  file_name = file_name)


################################################################################
# Sensitivity_12 min haplotypes = 37 ----
################################################################################
file_name <- "Model runs/sensitivity_12/sensitivity_12"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_12[[1]],  file_name = file_name)
plot_trajectory(sensitivity_12[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_12[[1]], sensitivity_12[[2]], sir_base[[1]], sensitivity_12[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_12[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_12[[1]],  file_name = file_name)



################################################################################
# Sensitivity_13 additional cv ----
################################################################################
file_name <- "Model runs/sensitivity_13/sensitivity_13"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_13[[1]],  file_name = file_name)
plot_trajectory(sensitivity_13[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_13[[1]], sensitivity_13[[2]], sir_base[[1]], sensitivity_13[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_13[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_13[[1]],  file_name = file_name)



################################################################################
# Sensitivity 14 - power equation q analytical q ----
################################################################################
file_name <- "Model runs/sensitivity_14/sensitivity_14"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_14[[1]],  file_name = file_name)
plot_trajectory(sensitivity_14[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_14[[1]], sensitivity_14[[2]], sir_base[[1]], sensitivity_14[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_14[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_14[[1]],  file_name = file_name)



################################################################################
# Sensitivity 15 - power equation q explicit q prior ----
################################################################################
file_name <- "Model runs/sensitivity_15/sensitivity_15"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_15[[1]],  file_name = file_name)
plot_trajectory(sensitivity_15[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_15[[1]], sensitivity_15[[2]], sir_base[[1]], sensitivity_15[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sensitivity_15[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_15[[1]],  file_name = file_name)

