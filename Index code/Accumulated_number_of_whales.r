#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load libraries ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(StateSpaceSIR)
require(MASS)
library(matrixStats)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Read in and clean data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
balle<- read.table( "Data/SRW_aerial_survey_1999-2024.txt", header = TRUE)


balle<- balle [c(-11,-12),] # 11 y 12 (2004) vuelos con el aerocomander// esta es la seleccion de datos que hice para la los paràmetros para el nùmero de ballenas que dan la vuelta anualmente por PV
balle <- subset (balle, Year<2025) # Subset years

# * Variable Respuesta ----
balle$RTA <- (balle$T) # Total of observed whales
balle$Year <- as.factor(balle$Year)
nyrs <- length(unique(balle$Year))


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# First stage - regression model ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#	Regression model selected (up to  2024)
regresion.balle.nb.jul.cuad <- glm.nb(RTA ~ Year + Juliano + I(Juliano^2), data = balle, link = log)
summary (regresion.balle.nb.jul.cuad)
nb.Jul.cuad <- cbind(Estimate = coef(regresion.balle.nb.jul.cuad))
vcov.nb.Jul.cuad <- vcov(regresion.balle.nb.jul.cuad) # Variance-covariance matrix

nb.Jul.cuad # Parameter estimates
vcov.nb.Jul.cuad # NOTE: dispersion parameter is not included


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Second Stage - accumulated number of whales ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# -- Run for years of interest and using nb regression model
A_xy <- data.frame(Year = sort(unique(balle$Year)), B = c(0, nb.Jul.cuad[2:nyrs])) # Years which we want to calculate the accumulated number of whales from years which we have data and the associated regression parameters. NOTE: 1990 is the intercept 

# Get A_xy using MASS and R
A_xy$A_xy <- NA
for(i in 1:nrow(A_xy)){
  A_xy$A_xy[i] <- accum_fun(a = nb.Jul.cuad[1],      # Intercept
                            b = A_xy$B[i],           # Year parameter
                            c = nb.Jul.cuad[nyrs+1], # Julian day parameter
                            d = nb.Jul.cuad[nyrs+2], # Julian day^2 parameter
                            mu = 60,                 # mu from manuscript
                            sigma = 8.66,            # sigma from manuscript
                            x = 320,
                            ReportList = FALSE)      # Calculate until day 320
  
}
A_xy


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Variance of second stage via numerical simulation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
ndraw = 100000
set.seed(123)
Param_draws <- mvrnorm(n = ndraw, mu = as.numeric(nb.Jul.cuad), Sigma = vcov.nb.Jul.cuad) # Simulate 100,000 parameter sets from MLE and variance-covariance matrix, assuming asymptotic normality
A_xy_mat <- matrix(NA, nrow = nrow(A_xy), ncol = ndraw) # Matrix to save A_xy estimates from each parameter draw

# * Get A_xy using MASS and R ----
for(draw in 1:ndraw){ # Loop through draws
  for(i in 1:nrow(A_xy)){ # Loop through years
    A_xy_mat[i,draw] <- accum_fun(a = Param_draws[draw, 1], # Intercept
                              b = c(0,Param_draws[draw, 2:nyrs])[i], # Year parameter (1999 is 0 because intercept)
                              c = Param_draws[draw,nyrs+1], # Julian day parameter
                              d = Param_draws[draw,nyrs+2], # Julian day^2 parameter
                              mu = 60,                      # mu from manuscript
                              sigma = 8.66,                 # sigma from manuscript
                              x = 320,
                              ReportList = FALSE) # Calculate until day 320
    
  }
}
A_xy$ln_A_xy_mu_sim <- rowMeans(log(A_xy_mat)) # Expected value of A_xy via numerical simulation.
A_xy$A_xy_mu_sim <- rowMeans((A_xy_mat))       # Expected value of A_xy via numerical simulation.
A_xy$ln_A_xy_var_sim <- rowVars(log(A_xy_mat)) # Variance of A_xy via numerical simulation.


# * Build variance covariance ----
A_xy_vcov <- diag(A_xy$ln_A_xy_var_sim)
for(i in 1:nrow(A_xy_vcov)){
  for(j in 1:nrow(A_xy_vcov)){
    if(i!=j){
      A_xy_vcov[i,j] <- cov(log(A_xy_mat[i,]), log(A_xy_mat[j,]))
    }
  }
}
A_xy_vcov
A_xy <- cbind(A_xy, A_xy_vcov)

write.csv(A_xy, file = "Data/Accumulated_n_whales_1999_to_2024.csv")
