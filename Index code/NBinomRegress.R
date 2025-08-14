################################################################################
# Load libraries
################################################################################
rm(list=ls())
require(MASS)
library(glmmTMB)
library(TMB)
library(StateSpaceSIR)

################################################################################
# Read in and clean data
################################################################################
balle<- read.csv ( "Data/ballena3.csv" , sep=";" , dec=".")


balle<- balle [c(-11,-12),] #11 y 12 (2004) vuelos con el aerocomander// esta es la seleccion de datos que hice para la los paràmetros para el nùmero de ballenas que dan la vuelta anualmente por PV
balle <- subset (balle, Year<2020) # Subset years

####Variable Respuesta ############
balle$RTA <- (balle$T) # Total of observed whales
balle$Year <- as.factor(balle$Year)


################################################################################
# First stage - regression model
################################################################################
#	Regresion model selected (up to  2019) using Mass package
regresion.balle.nb.jul.cuad <- glm.nb(RTA ~ Year + Juliano, data = balle, link = log)
regresion.balle.glmmTMB.jul.cuad <- glmmTMB(RTA ~ Year + Juliano, data = balle, family = nbinom2(link = "log")) # Doesnt converge
summary (regresion.balle.nb.jul.cuad)
nb.Jul.cuad <- cbind(Estimate = coef(regresion.balle.nb.jul.cuad))
nb.Jul.cuad # Estimates

# Run in TMB to get SD of stage two
setwd("R code")
compile("NBinomRegress.cpp")
dyn.load(dynlib("NBinomRegress"))
setwd("../")

# Set up predictions matrix'
form <- formula(~ Year + Juliano) # Whats the model look like
Xhat <- merge(data.frame(Juliano = 1:365), data.frame(Year = unique(balle$Year)), all = TRUE)
XhatMat <- model.matrix(form, data = Xhat) # Make into matrix

# Set up data, params and estimate
X = model.matrix(form, data = balle) # Make into matrix
data <- list(yobs = balle$RTA, # Response
             Xd=matrix(1, nrow = nrow(balle), ncol = 1), # Dispersion matrix
             X = X, # Design matrix
             Xhat = XhatMat, # Prediction matrix for calculating A_xy
             Nyr = length(unique(balle$Year)), # Number of years
             P_t = dnorm(1:365, 60, 8.66, FALSE), # Probability whale remains in area
             Iday = 100, # Day beginning Ingress
             Eday = 320 # Day ending ingress/egress
             )
parameters <- list(beta=rep(0, ncol(X)), betad=1)
obj <- MakeADFun(data, parameters, DLL="NBinomRegress", hessian = TRUE)
opt <- do.call("optim", obj)


# Compare parameter estimates
Param_est <- data.frame(Name = rownames(nb.Jul.cuad), Mass = round(as.numeric(nb.Jul.cuad), 6), TMB = round(opt$par[1:18], 6))
Param_est$Dif = Param_est$Mass - Param_est$TMB
Param_est # Great! Very close



################################################################################
# Second Stage - accumulated number of whales
################################################################################
# -- Run for years of interest and using nb regression model
A_xy <- data.frame(Year = sort(unique(balle$Year)), B = c(0, nb.Jul.cuad[2:17])) # Years which we want to calculate the accumulated number of whales from years which we have data and the associated regression parameters. NOTE: 1990 is the intercept 

# Get A_xy using MASS and R
A_xy$A_xy <- NA
for(i in 1:nrow(A_xy)){
  A_xy$A_xy[i] <- accum_fun(a = nb.Jul.cuad[1], # Intercept
                            b = A_xy$B[i], # Year parameter
                            c = nb.Jul.cuad[18], # Julian day parameter
                            d = 0, # Julian day^2 parameter
                            mu = 60, # mu from manuscript
                            sigma = 8.66, # sigma from manuscript
                            x = 320,
                            ReportList = FALSE) # Calculate until day 320
  
}
A_xy


# Get A_xy using TMB
report <- obj$report(obj$env$last.par.best)
rep <- sdreport(obj)
Xhat$A_xy <- report$A_xyLong
A_xyTMB <- Xhat[which(Xhat$Juliano == 320),]
A_xyTMB
