#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Compare indices of relative abundance ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
sw_right_rel_abundance_2024 <- read.csv("Data/Accumulated_n_whales_1999_to_2024.csv") 
sw_right_rel_abundance_2019 <- read.csv("Data/Accumulated_n_whales_1999_to_2019.csv") 

# - Against one another
plot(x = sw_right_rel_abundance_2019$A_xy_mu_sim, y = sw_right_rel_abundance_2024$A_xy_mu_sim[1:17],
     xlab = "1999-2019 Index", ylab = "1999-2024 Index", pch = 16)
abline(1,1, col = "grey", lwd = 1.5)

# - On time-series
sw_right_rel_abundance_2019$Upper95 <- qlnorm(0.975, mean = log(sw_right_rel_abundance_2019$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2019$ln_A_xy_var_sim))
sw_right_rel_abundance_2019$Lower95 <- qlnorm(0.025, mean = log(sw_right_rel_abundance_2019$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2019$ln_A_xy_var_sim))

sw_right_rel_abundance_2024$Upper95 <- qlnorm(0.975, mean = log(sw_right_rel_abundance_2024$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2024$ln_A_xy_var_sim))
sw_right_rel_abundance_2024$Lower95 <- qlnorm(0.025, mean = log(sw_right_rel_abundance_2024$A_xy_mu_sim), sd = sqrt(sw_right_rel_abundance_2024$ln_A_xy_var_sim))

xlim <- range(c(sw_right_rel_abundance_2019$Year, sw_right_rel_abundance_2024$Year))
ylim <- range(c(sw_right_rel_abundance_2019$Upper95, sw_right_rel_abundance_2024$Upper95, sw_right_rel_abundance_2019$Lower95, sw_right_rel_abundance_2024$Lower95))

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

points( x = sw_right_rel_abundance_2024$Year+0.15,
        y = sw_right_rel_abundance_2024$A_xy_mu_sim,
        col = "grey", pch = 16, cex = 2)
arrows( x0 = sw_right_rel_abundance_2024$Year+0.15,
        y0 = sw_right_rel_abundance_2019$Lower95,
        x1 = sw_right_rel_abundance_2024$Year+0.15,
        y1 = sw_right_rel_abundance_2024$Upper95,
        length=0.05, angle=90, code=3, lwd = 3, col = "grey")

legend("topleft", c("1999-2019 Index", "1999-2024 Index"), bty = "n", col = c(1, "grey"), pch = 16, pt.cex = 2)

