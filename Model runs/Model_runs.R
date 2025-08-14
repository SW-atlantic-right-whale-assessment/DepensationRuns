library(StateSpaceSIR)
library(EnvStats)
library(plyr)



################################################################################
# Read in data ----
################################################################################
# -- Catch
sw_right_data<-read.delim("Data/datosModeloBallenasmiles2020Miles1648to2019.csv", sep=";",header=FALSE)   
names(sw_right_data)<- c("Year","CatchMin","CatchMax","Nt")
sw_right_data <- rbind(sw_right_data, data.frame(Year = 2022:2023, CatchMin = 0, CatchMax = 0, Nt=0))


# Four periods of SLRs
# - Period 1: 1648-1770: SLR = 1
# - Period 2: 1771-1850: SLR ~ N(1.6, 0.04)
# - Period 3: 1851-1973: SLR ~ N(1.09, 0.04)
# - Period 4: 1974-Present: SLR = 1
catch_list <- list(sw_right_data[which(sw_right_data$Year < 1771),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1771 & sw_right_data$Year <= 1850),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1851 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])


# -- Absolute abundance
Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 4029, CV.obs = NA) # FIXME: not used as of 4/24/21
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725).


# -- Relative abundance
# - Index 1: Accumulated number of whales
sw_right_rel_abundance<-read.csv("Data/Accumulated_n_whales_1999_to_2023.csv") 

Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_rel_abundance)), 
                                    Year = sw_right_rel_abundance$Year, 
                                    IA.obs = sw_right_rel_abundance$A_xy_mu_sim)
var_covar <- sw_right_rel_abundance[,paste0("X",1:20)]
colnames(var_covar) <- 1:20
Rel.Abundance.SWRight = cbind(Rel.Abundance.SWRight, var_covar)

for(i in 1:15){
  dir.create(paste0("Model runs/Sensitivity_",i))
}

################################################################################
# Base model ----
################################################################################
file_name <- "Model runs/Base/Base"

sir_base <- list()
for(i in 1:2){
  sir_base[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_base, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base[[1]],  file_name = file_name)
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_base[[1]]),  file_name = file_name,   priors = list(sir_base[[2]]), inc_reference = FALSE)
# plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_base[[1]],  file_name = file_name)


################################################################################
# Sensitivity 1 - lognormal prior on Rmax ----
################################################################################
file_name <- "Model runs/Sensitivity_1/Sensitivity_1"

sensitivity_1 <- list()
for(i in 1:2){
  sensitivity_1[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(rlnorm, -2.67, 0.5),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_1[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_1[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_1, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_1[[1]],  file_name = file_name)
plot_trajectory(sensitivity_1[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_1[[1]]),  file_name = file_name,   priors = list(sensitivity_1[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_1[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_1[[1]],  file_name = file_name)


################################################################################
# Sensitivity 2 - smaller CV on lognromal prior on Rmax  ----
################################################################################
file_name <- "Model runs/Sensitivity_2/Sensitivity_2"

sensitivity_2 <- list()
for(i in 1:2){
  sensitivity_2[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(rlnorm, -2.67, 0.3),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_2[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_2[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_2, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_2[[1]],  file_name = file_name)
plot_trajectory(sensitivity_2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_2[[1]]),  file_name = file_name,   priors = list(sensitivity_2[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_2[[1]],  file_name = file_name)

################################################################################
# Sensitivity 3 ----
################################################################################
# - Prior on rmax of lognormal(-2.65, 0.5)T(0.2, 0.11)
# - Using rlnormTrunc fron "EnvStats" package

file_name <- "Model runs/sensitivity_3/sensitivity_3"

sensitivity_3 <- list()
for(i in 1:2){
  sensitivity_3[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(rlnormTrunc, -2.67, 0.5, 0.02, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_3[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_3[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_3, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_3[[1]],  file_name = file_name)
plot_trajectory(sensitivity_3[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_3[[1]]),  file_name = file_name,   priors = list(sensitivity_3[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_3[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_3[[1]],  file_name = file_name)



################################################################################
# Sensitivity 4 ----
################################################################################
file_name <- "Model runs/Sensitivity_4/Sensitivity_4"
# Upper bound of process error is 100 * lower bound

sensitivity_4 <- list()
for(i in 1:2){
  sensitivity_4[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 100),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 0.5e-3, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_4[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_4[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_4, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_4[[1]],  file_name = file_name)
plot_trajectory(sensitivity_4[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_4[[1]]),  file_name = file_name,   priors = list(sensitivity_4[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_4[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_4[[1]],  file_name = file_name)



################################################################################
# Sensitivity 5 ----
################################################################################
file_name <- "Model runs/sensitivity_5/sensitivity_5"
# Upper bound of process error variance is 5 times lower bound

sensitivity_5 <- list()
for(i in 1:2){
  sensitivity_5[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 2),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_5[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_5[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_5, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_5[[1]],  file_name = file_name)
plot_trajectory(sensitivity_5[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_5[[1]]),  file_name = file_name,   priors = list(sensitivity_5[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_5[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_5[[1]],  file_name = file_name)


################################################################################
# Sensitivity 6 ----
################################################################################
file_name <- "Model runs/sensitivity_6/sensitivity_6"
# Nrecent is 2004

sensitivity_6 <- list()
for(i in 1:2){
  sensitivity_6[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2004,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_6[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_6[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_6, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_6[[1]],  file_name = file_name)
plot_trajectory(sensitivity_6[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_6[[1]]),  file_name = file_name,   priors = list(sensitivity_6[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_6[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_6[[1]],  file_name = file_name)


################################################################################
# Sensitivity 7 ----
################################################################################
file_name <- "Model runs/sensitivity_7/sensitivity_7"
# No struck and loss rates

sensitivity_7 <- list()
for(i in 1:2){
  sensitivity_7[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(1), 
      make_prior(1),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_7[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_7[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_7, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_7[[1]],  file_name = file_name)
plot_trajectory(sensitivity_7[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_7[[1]]),  file_name = file_name,   priors = list(sensitivity_7[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_7[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_7[[1]],  file_name = file_name)



################################################################################
# Sensitivity 8 ----
################################################################################
file_name <- "Model runs/sensitivity_8/sensitivity_8"
# - Catch time series is only low

sensitivity_8 <- list()
for(i in 1:2){
  sensitivity_8[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             catch_sample = make_prior(0) # Set to 0 used low catch only
                             ),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_8[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_8[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_8, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_8[[1]],  file_name = file_name)
plot_trajectory(sensitivity_8[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_8[[1]]),  file_name = file_name,   priors = list(sensitivity_8[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_8[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_8[[1]],  file_name = file_name)



################################################################################
# Sensitivity 9 ----
################################################################################
file_name <- "Model runs/sensitivity_9/sensitivity_9"
# -- Catch time series is high

sensitivity_9 <- list()
for(i in 1:2){
  sensitivity_9[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             catch_sample = make_prior(1) # Set to 1 used high catch only
                             ),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_9[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_9[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_9, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_9[[1]],  file_name = file_name)
plot_trajectory(sensitivity_9[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_9[[1]]),  file_name = file_name,   priors = list(sensitivity_9[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_9[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_9[[1]],  file_name = file_name)


################################################################################
# Sensitivity 10 min haplotypes = 0 ----
################################################################################
file_name <- "Model runs/sensitivity_10/sensitivity_10"

sensitivity_10 <- list()
for(i in 1:2){
  sensitivity_10[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 0,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_10[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_10[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_10, file = paste0(file_name, ".Rdata"))



load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_10[[1]],  file_name = file_name)
plot_trajectory(sensitivity_10[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_10[[1]]),  file_name = file_name,   priors = list(sensitivity_10[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_10[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_10[[1]],  file_name = file_name)


################################################################################
# Sensitivity 11 min haplotypes = 25 ----
################################################################################
file_name <- "Model runs/sensitivity_11/sensitivity_11"

sensitivity_11 <- list()
for(i in 1:2){
  sensitivity_11[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 25,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_11[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_11[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_11, file = paste0(file_name, ".Rdata"))



load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_11[[1]],  file_name = file_name)
plot_trajectory(sensitivity_11[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_11[[1]]),  file_name = file_name,   priors = list(sensitivity_11[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_11[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_11[[1]],  file_name = file_name)


################################################################################
# Sensitivity12 min haplotypes = 37 ----
################################################################################
file_name <- "Model runs/sensitivity_12/sensitivity_12"

sensitivity_12 <- list()
for(i in 1:2){
  sensitivity_12[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 37,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_12[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_12[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_12, file = paste0(file_name, ".Rdata"))



load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_12[[1]],  file_name = file_name)
plot_trajectory(sensitivity_12[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_12[[1]]),  file_name = file_name,   priors = list(sensitivity_12[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_12[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_12[[1]],  file_name = file_name)



################################################################################
# Sensitivity 13 additional cv ----
################################################################################
file_name <- "Model runs/sensitivity_13/sensitivity_13"

sensitivity_13 <- list()
for(i in 1:2){
  sensitivity_13[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max = make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             add_VAR_IA = make_prior(rlnorm, log(0.2), 0.5),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 37,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_13[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_13[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_13, file = paste0(file_name, ".Rdata"))



load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_13[[1]],  file_name = file_name)
plot_trajectory(sensitivity_13[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_13[[1]]),  file_name = file_name,   priors = list(sensitivity_13[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_13[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_13[[1]],  file_name = file_name)



################################################################################
# Sensitivity 14 - power equation q analytical ----
################################################################################
file_name <- "Model runs/sensitivity_14/sensitivity_14"

sensitivity_14 <- list()
for(i in 1:2){
  sensitivity_14[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             q_IA1 = make_prior(use = FALSE),
                             q_IA2 = make_prior(rnorm, 0, 0.1),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_14[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_14[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_14, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_14[[1]],  file_name = file_name)
plot_trajectory(sensitivity_14[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_14[[1]]),  file_name = file_name,   priors = list(sensitivity_14[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_14[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_14[[1]],  file_name = file_name)


################################################################################
# Sensitivity 15 - power equation explicit q prior ----
################################################################################
file_name <- "Model runs/sensitivity_15/sensitivity_15"

sensitivity_15 <- list()
for(i in 1:2){
  sensitivity_15[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             q_IA1 = make_prior(rlnorm, log(0.5), 0.2),
                             q_IA2 = make_prior(rnorm, 0, 0.1),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2023, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sensitivity_15[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sensitivity_15[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sensitivity_15, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sensitivity_15[[1]],  file_name = file_name)
plot_trajectory(sensitivity_15[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sensitivity_15[[1]]),  file_name = file_name,   priors = list(sensitivity_15[[2]]), inc_reference = FALSE)
plot_ioa(sensitivity_15[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sensitivity_15[[1]],  file_name = file_name)
