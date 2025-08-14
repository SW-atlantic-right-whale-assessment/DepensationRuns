library(StateSpaceSIR)
library(EnvStats)
library(dplyr)



################################################################################
# Read in data
################################################################################
# -- Catch
nz_right_data<-read.csv("Data/Jackson 2016/Catch_inputs.csv")   # "Year","CatchMin","CatchMax","Nt"
catch.names=c("Year","CoastNZlow","CoastNZhigh","CoastEA","AmericanBaylow","AmericanBayhigh","AmericanoffNZ","AmericanoffNZSE","AmericanoffEA","AmericanoffEASE","French","SovietEA","SovietNZ")
colnames(nz_right_data)<-catch.names

AmericanBaymean=(nz_right_data$AmericanBaylow+nz_right_data$AmericanBayhigh)/2 # average of low & high American bay whaling catches, used to convert French total catches into 'bay' and 'offshore' components

France.offEA=(French.catch*offshore.AmericanEA3)/(AmericanBaymean+offshore.AmericanEA3+offshore.AmericanNZ3) # this breaks down the French catch series into component types of whaling, based on relative proportion estimated for the American campaign
France.offNZ=(French.catch-France.offEA)*(offshore.AmericanNZ3/(offshore.AmericanNZ3+AmericanBaymean)) # this derives the French catch off NZ and breaks down into relative bay and offshore components, as above.	
France.bayNZ=(French.catch-France.offEA)*(AmericanBaymean/(offshore.AmericanNZ3+AmericanBaymean))
# French catch is already smeared

# American offshore needs to be drawn from normal dist given SE and then smeared

# --- Low 
nzlow_catch_list <- list(
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$CoastNZlow + nz_right_data$AmericanBaylow, CatchMax = CatchMin), # Coastal NZ
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$Soviet.NZ, CatchMax = nz_right_data$Soviet.NZ) # Soviet catches (no SLR)
)


swpaclow_catch_list <- list(
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$CoastNZlow + nz_right_data$AmericanBaylow + nz_right_data$CoastEA, CatchMax = CatchMin), # Coastal SW Pacific
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$Soviet.NZ + nz_right_data$SovietEA, CatchMax = CatchMin) # Soviet catches SW Pacific (no SLR)
)


# --- High
nz_high_catch_list <- list(
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$CoastNZhigh + nz_right_data$AmericanBayhigh, CatchMax = CatchMin), # Coastal NZ
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$Soviet.NZ, CatchMax = nz_right_data$Soviet.NZ) # Soviet catches (no SLR)
)


swpac_high_catch_list <- list(
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$CoastNZhigh + nz_right_data$AmericanBayhigh + nz_right_data$CoastEA, CatchMax = CatchMin), # Coastal SW Pacific
  data.frame(Year = nz_right_data$Year, CatchMin = nz_right_data$Soviet.NZ + nz_right_data$SovietEA, CatchMax = CatchMin) # Soviet catches SW Pacific (no SLR)
)

# -- Absolute abundance
Abs.Abundance.2009 <- data.frame(Year = 2009, N.obs = 2148, CV.obs = 0.2) 

# -- Relative abundance
Rel.Abundance.NZRight <- data.frame(Index = rep(1, 4), 
                                    Year = c(1995, 1998, 2006, 2009), 
                                    IA.obs = c(533, 619, 910, 1074)) 
Rel.Abundance.NZRight <- cbind(Rel.Abundance.NZRight, diag(0.2,4))

for(i in 1:15){
  dir.create(paste0("NZ model runs/Sensitivity_",i))
}

################################################################################
# Base model
################################################################################
file_name <- "NZ model runs/Base/Base"

sir_base <- list()
for(i in 1:2){
  sir_base[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 1000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.12),
                             N_obs = make_prior(runif, 100, 20000),
                             var_N = make_prior(0),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(0.6)),
    catch_multipliers = make_multiplier_list(
      make_prior(rnorm, 1.27, 0.05), 
      make_prior(rnorm, 1.45, 0.054),
      make_prior(1)
      ),
    target.Yr = 2009,
    num.haplotypes = 12,
    output.Yrs = c(2021, 2030),
    abs.abundance = Abs.Abundance.2009,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.NZRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-3, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, FALSE, TRUE))
}
resample_summary_reference <- summary_sir(sir_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
save(sir_base, file = paste0(file_name, ".Rdata"))


load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base[[1]],  file_name = file_name)
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))
plot_density(SIR = list(sir_base[[1]]),  file_name = file_name,   priors = list(sir_base[[2]]), inc_reference = FALSE)
plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_base[[1]],  file_name = file_name)

