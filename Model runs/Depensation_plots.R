library(StateSpaceSIR)
library(EnvStats)
source("R/plot_surplus_prod_function.R")
source("R/plot_density_depensation.R")
source("R/plot_density.R")


################################################################################
# Base model
################################################################################
file_name <- "Model runs/Base2/Base2"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_base2[[1]],  file_name = file_name)
plot_trajectory(sir_base2[[2]],  file_name = paste0(file_name, "prior"))
# plot_abs_abundance(sir_base2[[1]],  file_name = file_name)
plot_density(SIR = sir_base2,  file_name = file_name, posteriors_lwd = c(3,1), posteriors_lty = rep(1, 2), posteriors_col = c(1,1))
plot_ioa(sir_base2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_base2[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_base2[[1]]), file_name = file_name, model_names = "Base", coolors = c("#104F55", "#2F0A28", "#3185FC"))


################################################################################
# Depensation model 1 - Hilborn et al 2014
################################################################################
file_name <- "Model runs/Depensation_1/Depensation_1"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation1[[1]],  file_name = file_name)
plot_trajectory(sir_depensation1[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation1[[1]], sir_depensation1[[2]], sir_base2[[1]], sir_depensation1[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation1[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation1[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation1[[1]], sir_depensation1[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-1", "Dep-1 prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))

################################################################################
# Depensation model 2 - Logistic
################################################################################
file_name <- "Model runs/Depensation_2/Depensation_2"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation2[[1]],  file_name = file_name)
plot_trajectory(sir_depensation2[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation2[[1]], sir_depensation2[[2]], sir_base2[[1]], sir_depensation2[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation2[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation2[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation2[[1]], sir_depensation2[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-2", "Dep-2 prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))



################################################################################
# Depensation model 3 - Lin & Li 2002
################################################################################
file_name <- "Model runs/Depensation_3/Depensation_3"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation3[[1]],  file_name = file_name)
plot_trajectory(sir_depensation3[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation3[[1]], sir_depensation3[[2]], sir_base2[[1]], sir_depensation3[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation3[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation3[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation3[[1]], sir_depensation3[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-3", "Dep-3 prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))



################################################################################
# Depensation model 4 - Haider et al 2017
################################################################################
file_name <- "Model runs/Depensation_4/Depensation_4"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation4[[1]],  file_name = file_name)
plot_trajectory(sir_depensation4[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation4[[1]], sir_depensation4[[2]], sir_base2[[1]], sir_depensation4[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation4[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation4[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation4[[1]], sir_depensation4[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-4", "Dep-4 prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))


################################################################################
# Depensation model 5 - Hilborn et al 2014 w beta prior
################################################################################
file_name <- "Model runs/Depensation_5/Depensation_5"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation5[[1]],  file_name = file_name)
# plot_abs_abundance(sir_depensation5[[1]],  file_name = file_name)
plot_trajectory(sir_depensation5[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation5[[1]], sir_depensation5[[2]], sir_base2[[1]], sir_depensation5[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation5[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation5[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation5[[1]], sir_depensation5[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-1b", "Dep-1b prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))


################################################################################
# Depensation model 6 - Logistic w/ beta prior
################################################################################
file_name <- "Model runs/Depensation_6/Depensation_6"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation6[[1]],  file_name = file_name)
plot_trajectory(sir_depensation6[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation6[[1]], sir_depensation6[[2]], sir_base2[[1]], sir_depensation6[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation6[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation6[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation6[[1]], sir_depensation6[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-2b", "Dep-2b prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))



################################################################################
# Depensation model 7 - Lin & Li 2002 w/ beta prior
################################################################################
file_name <- "Model runs/Depensation_7/Depensation_7"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation7[[1]],  file_name = file_name)
plot_trajectory(sir_depensation7[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation7[[1]], sir_depensation7[[2]], sir_base2[[1]], sir_depensation7[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation7[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation7[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation7[[1]], sir_depensation7[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-3b", "Dep-3b prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))



################################################################################
# Depensation model 8 - Haider et al 2017 w/ beta prior
################################################################################
file_name <- "Model runs/Depensation_8/Depensation_8"
load(file = paste0(file_name, ".Rdata"))
plot_trajectory(sir_depensation8[[1]],  file_name = file_name)
plot_trajectory(sir_depensation8[[2]],  file_name = paste0(file_name, "prior"))
plot_density_depensation(SIR = list(sir_depensation8[[1]], sir_depensation8[[2]], sir_base2[[1]], sir_depensation8[[1]]),  file_name = file_name,  posteriors_lwd = c(3,1,3,3), posteriors_lty = c(1,1,2,1), posteriors_col = c(1,1,"grey45",1))
plot_ioa(sir_depensation8[[1]],  file_name = file_name, ioa_names = NULL )
summary_table(sir_depensation8[[1]],  file_name = file_name)
plot_suplus_prod(list(sir_depensation8[[1]], sir_depensation8[[2]], sir_base2[[1]]), file_name = file_name, model_names = c("Dep-4", "Dep-4b prior", "Base"), coolors = c("#104F55", "#2F0A28", "#3185FC"))

