library(StateSpaceSIR)

# Load all the models
file_names <- c("Base/Base",
                "Sensitivity_1/Sensitivity_1",
                "Sensitivity_2/Sensitivity_2",
                "Sensitivity_3/Sensitivity_3",
                "Sensitivity_4/Sensitivity_4",
                "Sensitivity_5/Sensitivity_5",
                "Sensitivity_6/Sensitivity_6",
                "Sensitivity_7/Sensitivity_7",
                "Sensitivity_8/Sensitivity_8",
                "Sensitivity_9/Sensitivity_9",
                "Sensitivity_10/Sensitivity_10",
                "Sensitivity_11/Sensitivity_11",
                "Sensitivity_12/Sensitivity_12",
                "Sensitivity_13/Sensitivity_13",
                "sensitivity_14/sensitivity_14",
                "Model_average/Model_average")

for(i in 1:length(file_names)){
  load(file = paste0("Model runs/",file_names[i], ".Rdata"))
}

# Densities and trajectories with reference
sir_list <- list(sir_base,
                 sensitivity_1, 
                 sensitivity_2, 
                 sensitivity_3,
                 sensitivity_4,
                 sensitivity_5,
                 sensitivity_6,
                 sensitivity_7,
                 sensitivity_8,
                 sensitivity_9,
                 sensitivity_10,
                 sensitivity_11,
                 sensitivity_12,
                 sensitivity_13,
                 sensitivity_14)

for(i in 1:length(sir_list)){
  plot_abs_abundance(sir_list[[i]][[1]],  file_name = paste0("Model runs/",file_names[i]))
}

for(i in 2:length(sir_list)){
  sir_list_tmp <- list(sir_list[[i]][[1]], sir_base[[1]], sir_list[[i]][[1]])
  plot_density(SIR = sir_list_tmp,  file_name = paste0("Model runs/",file_names[i]),   priors = list(sir_list[[i]][[2]]), inc_reference = TRUE, target = ifelse(i == 7, FALSE, TRUE))
  plot_trajectory( SIR = sir_list[[i]][[1]], Reference = sir_list[[1]][[1]],  file_name = paste0("Model runs/",file_names[i]))
}


#############################################################
#### Model averaging
#############################################################
# Get bayes factor for models with comparable likelihoods
bayes_f <- bayes_factor(SIR = list(sir_base[[1]],
                                   sensitivity_1[[1]], 
                                   sensitivity_2[[1]], 
                                   sensitivity_3[[1]],
                                   #sensitivity_4[[1]],
                                   #sensitivity_5[[1]],
                                   sensitivity_6[[1]],
                                   sensitivity_7[[1]],
                                   #sensitivity_8[[1]],
                                   #sensitivity_9[[1]],
                                   sensitivity_10[[1]],
                                   sensitivity_11[[1]],
                                   sensitivity_12[[1]],
                                   sensitivity_13[[1]],
                                   sensitivity_14[[1]]))


# Create a new model based on bayes factors
model_average <- weight_model(SIR = list(sir_base[[1]],
                                         sensitivity_1[[1]], 
                                         sensitivity_2[[1]], 
                                         sensitivity_3[[1]],
                                         #sensitivity_4[[1]],
                                         #sensitivity_5[[1]],
                                         sensitivity_6[[1]],
                                         sensitivity_7[[1]],
                                         #sensitivity_8[[1]],
                                         #sensitivity_9[[1]],
                                         sensitivity_10[[1]],
                                         sensitivity_11[[1]],
                                         sensitivity_12[[1]],
                                         sensitivity_13[[1]],
                                         sensitivity_14[[1]]), 
                              bayes_factor = bayes_f)
model_average$inputs$output.Years <- c(1999, 2021, 2023, 2030)

# For plotting make a vector of bayes factors, set NA for models that cant be compared (different likelihood)
bayes_vec <- c(bayes_f[1:4], NA,NA, bayes_f[5:6],  NA, NA, bayes_f[7:11], NA)
model_names <-  c( "B", paste0("S-", 1:14), "MA")
table2 <- data.frame(Model = model_names, BayesFactor = round(bayes_vec,4))
write.csv(table2, file = paste0(paste0("Model runs/",file_names[16],"_bayes_factors.csv")))

# Compare Aposteriors of all
compare_posteriors(
  reference_sir = TRUE, 
  SIR = list(sir_base[[1]],
             sensitivity_1[[1]], 
             sensitivity_2[[1]], 
             sensitivity_3[[1]],
             sensitivity_4[[1]],
             sensitivity_5[[1]],
             sensitivity_6[[1]],
             sensitivity_7[[1]],
             sensitivity_8[[1]],
             sensitivity_9[[1]],
             sensitivity_10[[1]],
             sensitivity_11[[1]],
             sensitivity_12[[1]],
             sensitivity_13[[1]],
             sensitivity_14[[1]],
             model_average), 
  model_names = model_names, 
  bayes_factor = round(bayes_vec,2),
  file_name = paste0("Model runs/",file_names[16]),
  years = c(2021, 2030))

# Plot and get parameter values from Model Average
file_name <-paste0("Model runs/",file_names[16])
trajectory_summary_reference <- summary_sir(model_average$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
plot_trajectory(model_average, Reference = sir_base[[1]],  file_name = file_name)
plot_abs_abundance(model_average,  file_name = file_name)
sir_list_ma <- list(model_average, sir_base[[2]], sir_base[[1]], model_average)


source("R/plot_density_model_average.R", echo=TRUE)
plot_density_ma(SIR = sir_list_ma, posteriors_lwd = c(3, 1,3,3), posteriors_lty = c(1,1,1,1), posteriors_col = c(1, "grey45", "grey45", 1) ,  file_name = file_name)
plot_ioa(model_average,  file_name = file_name, ioa_names = NULL)
summary_table(model_average,  file_name = file_name)
save(model_average, file = paste0(file_name, ".Rdata"))


## Compare indices of relative abundance
sw_right_rel_abundance_2023 <- read.csv("Data/Accumulated_n_whales_1999_to_2023.csv") 
sw_right_rel_abundance_2019 <- read.csv("Data/Accumulated_n_whales_1999_to_2019.csv") 

# - Against one another
plot(x = sw_right_rel_abundance_2019$A_xy_mu_sim, y = sw_right_rel_abundance_2023$A_xy_mu_sim[1:17],
     xlab = "1999-2019 Index", ylab = "1999-2023 Index", pch = 16)
abline(1,1, col = "grey", lwd = 1.5)

# - On time-series
sw_right_rel_abundance_2019$Upper95 <- qlnorm(0.975, mean = log(sw_right_rel_abundance_2019$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2019$ln_A_xy_var_sim))
sw_right_rel_abundance_2019$Lower95 <- qlnorm(0.025, mean = log(sw_right_rel_abundance_2019$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2019$ln_A_xy_var_sim))

sw_right_rel_abundance_2023$Upper95 <- qlnorm(0.975, mean = log(sw_right_rel_abundance_2023$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2023$ln_A_xy_var_sim))
sw_right_rel_abundance_2023$Lower95 <- qlnorm(0.025, mean = log(sw_right_rel_abundance_2023$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2023$ln_A_xy_var_sim))

xlim <- range(c(sw_right_rel_abundance_2019$Year, sw_right_rel_abundance_2023$Year))
ylim <- range(c(sw_right_rel_abundance_2019$Upper95, sw_right_rel_abundance_2023$Upper95, sw_right_rel_abundance_2019$Lower95, sw_right_rel_abundance_2023$Lower95))

# - Plot
plot(NA, NA, ylim = ylim, xlim = xlim, xlab = "Year", ylab = "Index")

# Relative abundance
# - 2019
points( x = sw_right_rel_abundance_2019$Year-0.15,
        y = sw_right_rel_abundance_2019$A_xy_mu_sim,
        col = 1, pch = 16, cex = 2)
arrows( x0 = sw_right_rel_abundance_2019$Year-0.15,
        y0 = sw_right_rel_abundance_2019$Lower95,
        x1 = sw_right_rel_abundance_2019$Year-0.15,
        y1 = sw_right_rel_abundance_2019$Upper95,
        length=0.05, angle=90, code=3, lwd = 3, col = 1)

points( x = sw_right_rel_abundance_2023$Year+0.15,
        y = sw_right_rel_abundance_2023$A_xy_mu_sim,
        col = "grey", pch = 16, cex = 2)
arrows( x0 = sw_right_rel_abundance_2023$Year+0.15,
        y0 = sw_right_rel_abundance_2019$Lower95,
        x1 = sw_right_rel_abundance_2023$Year+0.15,
        y1 = sw_right_rel_abundance_2023$Upper95,
        length=0.05, angle=90, code=3, lwd = 3, col = "grey")

legend("topleft", c("1999-2019 Index", "1999-2023 Index"), bty = "n", col = c(1, "grey"), pch = 16, pt.cex = 2)

