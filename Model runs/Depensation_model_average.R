# Updated from Romero et al 2022 with models that include Allee effects
library(StateSpaceSIR)
source("R/plot_ioa_medians.R")
source("R/plot_surplus_prod_function.R")
source("R/plot_density_model_average.R")
source("R/plot_density_depensation.R")

# Load all the models
file_names <- c("Base2/Base2",
                "Depensation_1/Depensation_1",
                "Depensation_2/Depensation_2",
                "Depensation_3/Depensation_3",
                "Depensation_4/Depensation_4",
                "Depensation_5/Depensation_5",
                "Depensation_6/Depensation_6",
                "Depensation_7/Depensation_7",
                "Depensation_8/Depensation_8")


for(i in 1:length(file_names)){
  load(file = paste0("Model runs/",file_names[i], ".Rdata"))
}





# Plot surplus production function
plot_suplus_prod(SIRlist = list(sir_base2[[1]],
                                sir_depensation1[[1]],
                                sir_depensation5[[1]],
                                sir_depensation2[[1]],
                                sir_depensation6[[1]]),
                 coolors = c("#FB8B24", "#D90368", "#820263", "#291720", "#04A777"),
                 model_names = c("Base", "Dep-1", "Dep-1b", "Dep-2", "Dep-2b"),
                 file_name = "Model runs/Model_average_2/Model_average_2_set1")


# Plot surplus production function
plot_suplus_prod(SIRlist = list(sir_base2[[1]],
                                sir_depensation3[[1]],
                                sir_depensation7[[1]],
                                sir_depensation4[[1]],
                                sir_depensation8[[1]]),
                 coolors = c("#FB8B24", "#D90368", "#820263", "#291720", "#04A777"),
                 model_names = c("Base", "Dep-3", "Dep-3b", "Dep-4", "Dep-4b"),
                 file_name = "Model runs/Model_average_2/Model_average_2_set2")



# Densities and trajectories with reference
sir_list <- list(sir_base2,
                 sir_depensation1,
                 sir_depensation2,
                 sir_depensation3,
                 sir_depensation4,
                 sir_depensation5,
                 sir_depensation6,
                 sir_depensation7,
                 sir_depensation8)

for(i in 1:length(sir_list)){
  plot_abs_abundance(sir_list[[i]][[1]],  file_name = paste0("Model runs/",file_names[i]))
}

for(i in 2:length(sir_list)){
  sir_list_tmp <- list(sir_list[[i]][[1]], sir_base2[[1]], sir_list[[i]][[1]])
  plot_density_depensation(SIR = sir_list_tmp,  file_name = paste0("Model runs/",file_names[i]),   priors = list(sir_list[[i]][[2]]), inc_reference = TRUE, target = ifelse(i == 7, FALSE, TRUE))
  plot_trajectory( SIR = sir_list[[i]][[1]], Reference = sir_list[[1]][[1]],  file_name = paste0("Model runs/",file_names[i]))
}


#############################################################
#### Model averaging
#############################################################
# Get bayes factor for models with comparable likelihoods
bayes_f <- bayes_factor(SIR = list(sir_base2[[1]],
                                   sir_depensation1[[1]],
                                   # sir_depensation5[[1]],
                                   sir_depensation2[[1]],
                                   # sir_depensation6[[1]],
                                   sir_depensation3[[1]],
                                   # sir_depensation7[[1]],
                                   sir_depensation4[[1]]
                                   # sir_depensation8[[1]]
))

waic <- waic(SIR = list(sir_base2[[1]],
                        sir_depensation1[[1]],
                        # sir_depensation5[[1]],
                        sir_depensation2[[1]],
                        # sir_depensation6[[1]],
                        sir_depensation3[[1]],
                        # sir_depensation7[[1]],
                        sir_depensation4[[1]]
                        # sir_depensation8[[1]]
                        ))


# Create a new model based on bayes factors
model_average <- weight_model(SIR = list(sir_base2[[1]],
                                         sir_depensation1[[1]],
                                         #sir_depensation5[[1]],
                                         sir_depensation2[[1]],
                                         #sir_depensation6[[1]],
                                         sir_depensation3[[1]],
                                         #sir_depensation7[[1]],
                                         sir_depensation4[[1]]
                                         #sir_depensation8[[1]]
), 
bayes_factor = bayes_f)


# For plotting make a vector of bayes factors, set NA for models that cant be compared (different likelihood)
bayes_vec <- c(bayes_f, NA)
waic <- c(waic, NA)
model_names <-  c("Base", "Dep-1", "Dep-1b", "Dep-2", "Dep-2b", "Dep-3", "Dep-3b", "Dep-4", "Dep-4b", "MA")
model_names_short <-  c("Base", "Dep-1", "Dep-2", "Dep-3", "Dep-4", "MA")

table2 <- data.frame(Model = model_names_short, BayesFactor = round(bayes_vec,4), WAIC = round(waic, 4))
write.csv(table2, file = paste0(paste0("Model runs/",file_names[10],"_bayes_factors_and_waic.csv")))


# Compare posteriors of all
compare_posteriors(
  reference_sir = TRUE, 
  SIR = list(sir_base2[[1]],
             sir_depensation1[[1]],
             sir_depensation5[[1]],
             sir_depensation2[[1]],
             sir_depensation6[[1]],
             sir_depensation3[[1]],
             sir_depensation7[[1]],
             sir_depensation4[[1]],
             sir_depensation8[[1]],
             model_average), 
  model_names = model_names, 
  bayes_factor = round(c(bayes_vec[1], bayes_vec[2], NA, bayes_vec[3], NA, bayes_vec[4], NA, bayes_vec[5], NA, NA),2),
  file_name = paste0("Model runs/",file_names[10]),
  years = c(2021, 2030))


# Plot IOA medians
plot_ioa_medians(SIR = list(sir_base2[[1]],
                            sir_depensation1[[1]],
                            sir_depensation5[[1]],
                            sir_depensation2[[1]],
                            sir_depensation6[[1]],
                            sir_depensation3[[1]],
                            sir_depensation7[[1]],
                            sir_depensation4[[1]],
                            sir_depensation8[[1]],
                            model_average),
                 model_names = model_names,
                 file_name = paste0("Model runs/",file_names[10]))


# Plot and get parameter values from Model Average
file_name <-paste0("Model runs/",file_names[10])
trajectory_summary_reference <- summary_sir(model_average$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)
plot_trajectory(model_average, Reference = sir_base2[[1]],  file_name = file_name)
plot_abs_abundance(model_average,  file_name = file_name)
sir_list_ma <- list(model_average, sir_base2[[2]], sir_base2[[1]], model_average)

plot_density_ma(SIR = sir_list_ma, posteriors_lwd = c(3, 1,3,3), posteriors_lty = c(1,1,1,1), posteriors_col = c(1, "grey45", "grey45", 1) ,  file_name = file_name)
plot_ioa(model_average,  file_name = file_name, ioa_names = NULL)
StateSpaceSIR::summary_table(model_average,  file_name = file_name)
save(model_average, file = paste0(file_name, ".Rdata"))
