

#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated posterior densities of parameters from  SIR  model.
#'
#' @param SIR A fit SIR model or list of SIR models. Plots in the order provided.
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param lower Vector of lower bounds for x-axis
#' @param upper Vector of upper bounds for x-axis
#' @param probs Lower and upper quantiles to use for plot limits if lower and upper are not specified.
#' @param posteriors_lwd Line width for models
#' @param posteriors_lty Line type for models
#' @param posteriors_col Line color for models
#'
#' @return Returns and saves a figure with the posterior densities of parameters.
#' @export
plot_density <- function(SIR, posteriors_lwd = rep(3, length(SIR)), posteriors_lty = rep(1, length(SIR)), posteriors_col = rep(1, length(SIR)),  file_name = NULL, lower = NULL, upper = NULL, probs = c(0.025, 0.975), target = TRUE){
  
  # Make into list
  if(class(SIR) == "SIR"){
    SIR <- list(SIR)
  }
  
  # Vars of interest
  num.IA <- sort(unique(c( sapply(SIR, function(x) x$inputs$rel.abundance$Index))))
  if(target){
    years <- sort(unique(unlist(c( 
      sapply(SIR, function(x) x$inputs$target.Yr),
      sapply(SIR, function(x) x$inputs$output.Years)))))
  }
  if(!target){
    years <- sort(unique(unlist(c( 
      #sapply(SIR, function(x) x$inputs$target.Yr),
      sapply(SIR, function(x) x$inputs$output.Years)))))
  }
  vars <- c("r_max", "K", "Pmsy", "var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years))#, paste0("q_IA1", num.IA), paste0("q_IA2", num.IA), "add_VAR_IA", paste0("catch_multiplier_",2:3), "catch_parameter")
  vars_latex <- c("$r_{max}$", "$K$", "$P_{MSY}$", "$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "$P_{min}$", paste0("$P_{", years,"}$"))#, paste0("$q_{flt", num.IA, "}$"), paste0("$\beta_{q_{flt", num.IA,"}}$"), "$tau_q$", paste0("$SLR_",1:2,"$"), "$pi$")
  
  
  vars <- c("r_max", "K", "Pmsy", "var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", num.IA), paste0("q_IA2", num.IA), "add_VAR_IA", paste0("catch_multiplier_",2:3), "catch_parameter")
  vars_latex <- c("$r_{max}$", "$K$", "$P_{MSY}$", "$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "$P_{min}$", paste0("$P_{", years,"}$"), paste0("$q_{flt", num.IA, "}$"), paste0("$\beta_{q_{flt", num.IA,"}}$"), "$tau_q$", paste0("$SLR_",1:2,"$"), "$pi$")
  
  
  # Only select vars that have multiple unique parameters
  x <- SIR[[1]]$resamples_output[,vars]
  unique_vars <- sapply(x, function(x) length(unique(x)))
  
  for(i in 1:length(SIR)){
    SIR[[i]]$resamples_output$var_N <- sqrt(SIR[[i]]$resamples_output$var_N)
  }
  
  vars <- vars[which(unique_vars > 1)]
  vars_latex <- vars_latex[which(unique_vars > 1)]
  
  # -- Get posterior of q
  # Get posterior of q
  q_posteriors <- list() # Each layer is an index
  
  for(k in 1:length(SIR)){
    q_posteriors[[k]] <- list()
    
    # -- Determining the number of Indices of Abundance available
    rel.abundance <- SIR[[k]]$inputs$rel.abundance
    indices <- unique(rel.abundance$Index)
    IA.yrs <- rel.abundance$Year
    N_hat <- SIR[[k]]$resamples_trajectories[, paste0("N_", IA.yrs)] # Estimates of N within IOA years
    
    # -- Q2 for exponent
    q1_cols <- grep("q_IA1", colnames(SIR[[k]]$resamples_output)) # Columns of resample Q estimates
    q1_est <- SIR[[k]]$resamples_output[, q1_cols]
    q1_est <- as.matrix(q1_est, ncol = length(q1_cols))
    
    q2_cols <- grep("q_IA2", colnames(SIR[[k]]$resamples_output)) # Columns of resample Q estimates
    q2_est <- SIR[[k]]$resamples_output[, q2_cols]
    q2_est <- as.matrix(q2_est, ncol = length(q2_cols))
    
    # -- Make var-covar into wide and tall with cov = 0 for different indices
    rel.var.covar.tall <-  subset(rel.abundance, select = -c(Index,Year,IA.obs,IndYear))
    rel.var.covar.wide <- rel.var.covar.tall[which(rel.abundance$Index == 1),]
    rel.var.covar.wide <- rel.var.covar.wide[1:nrow(rel.var.covar.wide),1:nrow(rel.var.covar.wide)]
    
    rel.hess.wide <- solve(rel.var.covar.wide[1:nrow(rel.var.covar.wide), 1: nrow(rel.var.covar.wide)])
    
    if(num.IA>1){
      for(i in 2:length(unique(rel.abundance$Index))){
        var.cov.tmp <- as.matrix(rel.var.covar.tall[which(rel.abundance$Index == i),])
        var.cov.tmp <- var.cov.tmp[1:nrow(var.cov.tmp), 1:nrow(var.cov.tmp)]
        colnames(var.cov.tmp) <- NULL
        rownames(var.cov.tmp) <- NULL
        rel.var.covar.wide <- Matrix::bdiag(as.matrix(rel.var.covar.wide), var.cov.tmp)
        rel.hess.tall <- plyr::rbind.fill.matrix(rel.hess.tall, solve(var.cov.tmp))
      }
    }
    rel.var.covar.wide <- as.matrix(rel.var.covar.wide)
    
    # -- Loop through posterior draws
    for(j in 1:nrow(SIR[[k]]$resamples_trajectories)){
      # -- Sample q
      q_posteriors_tmp <- exp(MASS::mvrnorm(
        n = 5,
        mu = as.numeric(log(rel.abundance$IA.obs/N_hat[j,] ^ (q2_est[j,rel.abundance$Index] + 1)) - diag(rel.var.covar.wide)/2),
        Sigma = rel.var.covar.wide))
      
      # q_est <- exp(sum(rel.hess.wide %*% as.numeric(log(rel.abundance$IA.obs/N_hat[j,] ^ (q2_est[j,rel.abundance$Index] + 1))))/(sum(rel.hess.wide))) # q_i
      
      # -- Assign to list
      for(i in indices){
        if(j == 1){
          q_posteriors[[k]][[i]] <- c(q_posteriors_tmp[,which(rel.abundance$Index == i)])
        } else {
          q_posteriors[[k]][[i]] <- c(q_posteriors[[k]][[i]], c(q_posteriors_tmp[,which(rel.abundance$Index == i)]))
        }
      }
    }
  }
  
  
  # Plot
  for(j in 1:(1 + as.numeric(!is.null(file_name)) * 2)){
    
    # PNG
    if(j == 2){
      filename <- paste0(file_name, "_posterior_density", ".png")
      png( file = filename , width=10, height = 110 / 25.4, family = "serif", units = "in", res = 300)
    }
    
    # PDF
    if(j == 3){
      filename <- paste0(file_name, "_posterior_density", ".pdf")
      pdf( file = filename , width=10, height = 110 / 25.4, family = "serif")
    }
    
    par(mfrow = c(4,ceiling(length(vars)/4) + 2))
    par( mar=c(3, 0.05 , 0.5 , 0.55) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    
    plot.new()
    
    
    # Loop through vars
    for(i in 1:length(vars)){
      
      # Extract posterio densities
      posterior_dens <- list()
      for(k in 1:length(SIR)){
        posterior_dens[[k]] <- density(as.numeric(as.character(SIR[[k]]$resamples_output[,vars[i]])))
        if(vars[i] %in% paste0("q_IA1", num.IA)){
          posterior_dens[[k]] <- density(q_posteriors[[k]][[as.numeric(unlist(strsplit(paste0("q_IA1", num.IA), "q_IA1"))[2])]])
        }
      }
      
      # Get x range
      if(is.null(lower[i])){
        xlow <- quantile(sapply(posterior_dens, "[", "x")$x, probs= probs[1])
      } else if(is.na(lower[i])){
        xlow <- quantile(sapply(posterior_dens, "[", "x")$x, probs= probs[1])
      } else{
        xlow <- lower[i]
      }
      
      if(is.null(upper[i])){
        xup <- quantile(sapply(posterior_dens, "[", "x")$x, probs= probs[2])
      }
      else if(is.na(upper[i])){
        xup <- quantile(sapply(posterior_dens, "[", "x")$x, probs= probs[2])
      } else{
        xup <- upper[i]
      }
      
      
      # Plot them
      plot(NA,
           xlim = c(xlow, xup),
           ylim = c(0, range(sapply(posterior_dens, "[", "y"))[2]),
           ylab = NA, xlab = latex2exp::TeX(vars_latex[i]), yaxt = "n")
      mapply(lines, posterior_dens, lwd = posteriors_lwd, lty = posteriors_lty, col = posteriors_col[1:length(posterior_dens)])
      
      if(i %in% c(ceiling(length(vars)/4) * 1:3))  {
        plot.new()
        plot.new()
      }
      
      
      if(i %in% c(1, ceiling(length(vars)/4) * 1:4 + 1) ) {
        mtext(side = 2, "Density", line = 1, cex= 0.75)
      }
      
      
    }
    
    if(j > 1){ dev.off()}
  }
}
