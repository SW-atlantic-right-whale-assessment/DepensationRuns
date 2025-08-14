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

# - Index 2: Cooke et al 2001 mark-recapture females
# - Not used 10/05/2023
# cooke_ioa <-read.csv("Data/Cooke_et_al_2001.csv") 
# Rel.Abundance.SWRight.Cooke <- data.frame(Index = rep(2, nrow(cooke_ioa)), 
#                                           Year = cooke_ioa$Year, 
#                                           IA.obs = cooke_ioa$Population)
# Rel.Abundance.SWRight.Cooke = cbind(Rel.Abundance.SWRight.Cooke, diag(sqrt(log(cooke_ioa$PopSE/cooke_ioa$Population^2 + 1)))) #CV to lognmoral sigma
# 
# # - Combine
# Rel.Abundance.SWRight <- rbind.fill(Rel.Abundance.SWRight, Rel.Abundance.SWRight.Cooke)

# -- Set up directories
for(i in 1:8){
  dir.create(paste0("Model runs/Depensation_",i))
}

################################################################################
# Base model ----
################################################################################
file_name <- "Model runs/Base2/Base2"

sir_base2 <- list()
for(i in 1:2){
  sir_base2[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 0,
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_base2[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base2[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_base2, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base2[[1]],  file_name = file_name)
plot_trajectory(sir_base2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_base2[[1]]),  file_name = file_name,   priors = list(sir_base2[[2]]), inc_reference = FALSE)
plot_ioa(sir_base2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_base2[[1]],  file_name = file_name)



################################################################################
# Depensation model 1 - Hilborn et al 2014 ----
################################################################################
file_name <- "Model runs/Depensation_1/Depensation_1"

sir_depensation1 <- list()
for(i in 1:2){
  sir_depensation1[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 1,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             P50 = make_prior(runif, 0, .2)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation1[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation1[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation1, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation1[[1]],  file_name = file_name)
plot_trajectory(sir_depensation1[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation1[[1]]),  file_name = file_name,   priors = list(sir_depensation1[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation1[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation1[[1]],  file_name = file_name)



################################################################################
# Depensation model 2 - Logistic ----
################################################################################
file_name <- "Model runs/Depensation_2/Depensation_2"

sir_depensation2 <- list()
for(i in 1:2){
  sir_depensation2[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 2,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             P50 = make_prior(runif, 0, .2)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation2[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation2[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation2, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation2[[1]],  file_name = file_name)
plot_trajectory(sir_depensation2[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation2[[1]]),  file_name = file_name,   priors = list(sir_depensation2[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation2[[1]],  file_name = file_name)



################################################################################
# Depensation model 3 - Lin & Li 2002 ----
################################################################################
file_name <- "Model runs/Depensation_3/Depensation_3"

sir_depensation3 <- list()
for(i in 1:2){
  sir_depensation3[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 3,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8),
                             P50 = make_prior(runif, 0, .2)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation3[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation3[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation3, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation3[[1]],  file_name = file_name)
plot_trajectory(sir_depensation3[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation3[[1]]),  file_name = file_name,   priors = list(sir_depensation3[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation3[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation3[[1]],  file_name = file_name)



################################################################################
# Depensation model 4 - Haider et al 2017 ----
################################################################################
file_name <- "Model runs/Depensation_4/Depensation_4"

sir_depensation4 <- list()
for(i in 1:2){
  sir_depensation4[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 4,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8),
                             P50 = make_prior(runif, 0.001, .2)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation4[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation4[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation4, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation4[[1]],  file_name = file_name)
plot_trajectory(sir_depensation4[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation4[[1]]),  file_name = file_name,   priors = list(sir_depensation4[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation4[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation4[[1]],  file_name = file_name)



################################################################################
# Depensation model 5 - Hilborn et al 2014 w beta prior ----
################################################################################
file_name <- "Model runs/Depensation_5/Depensation_5"

sir_depensation5 <- list()
for(i in 1:2){
  sir_depensation5[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 1,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             P50 = make_prior(rbeta, 1, 10)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation5[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation5[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation5, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation5[[1]],  file_name = file_name)
plot_trajectory(sir_depensation5[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation5[[1]]),  file_name = file_name,   priors = list(sir_depensation5[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation5[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation5[[1]],  file_name = file_name)



################################################################################
# Depensation model 6 - Logistic w/ beta prior ----
################################################################################
file_name <- "Model runs/Depensation_6/Depensation_6"

sir_depensation6 <- list()
for(i in 1:2){
  sir_depensation6[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 2,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8),
                             P50 = make_prior(rbeta, 1, 10)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation6[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation6[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation6, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation6[[1]],  file_name = file_name)
plot_trajectory(sir_depensation6[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation6[[1]]),  file_name = file_name,   priors = list(sir_depensation6[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation6[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation6[[1]],  file_name = file_name)



################################################################################
# Depensation model 7 - Lin & Li 2002 w/ beta prior ----
################################################################################
file_name <- "Model runs/Depensation_7/Depensation_7"

sir_depensation7 <- list()
for(i in 1:2){
  sir_depensation7[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 3,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8),
                             P50 = make_prior(rbeta, 1, 10)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation7[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation7[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation7, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation7[[1]],  file_name = file_name)
plot_trajectory(sir_depensation7[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation7[[1]]),  file_name = file_name,   priors = list(sir_depensation7[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation7[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation7[[1]],  file_name = file_name)



################################################################################
# Depensation model 8 - Haider et al 2017 w/ beta prior ----
################################################################################
file_name <- "Model runs/Depensation_8/Depensation_8"

sir_depensation8 <- list()
for(i in 1:2){
  sir_depensation8[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    allee_model = 4,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.6, 0.8),
                             P50 = make_prior(rbeta, 1, 10)), # curve(dbeta(x, 1, 10), from = 0 ,to = 1)
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
    control = sir_control(threshold = 1e-05, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_depensation8[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_depensation8[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_depensation8, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation8[[1]],  file_name = file_name)
plot_trajectory(sir_depensation8[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_depensation8[[1]]),  file_name = file_name,   priors = list(sir_depensation8[[2]]), inc_reference = FALSE)
#plot_ioa(sir_depensation8[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation8[[1]],  file_name = file_name)



