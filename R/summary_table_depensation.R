summary_table <- function( SIR, file_name = NULL){
  
  num.IA <- max(SIR$inputs$rel.abundance$Index)
  
  # Vars of interest
  years <- sort(c( SIR$inputs$target.Yr, SIR$inputs$output.Years))
  
  if(SIR$inputs$allee_model == 0){
    vars <- c("r_max", "K", "z", "Pmsy", "var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")
    vars_latex <- c("$r_{max}$", "$K$", "$z$", "$Pmsy$","$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years), paste0("$q_{flt", 1:num.IA, "}$"), paste0("$\beta_{q_{flt", 1:num.IA,"}}$"), "$sigma_q$")
  } else{
    
    vars <- c("r_max", "K", "z", "Pmsy", "P50","var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")
    vars_latex <- c("$r_{max}$", "$K$", "$z$", "$Pmsy$", "$P_{50}$","$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years), paste0("$q_{flt", 1:num.IA, "}$"), paste0("$\beta_{q_{flt", 1:num.IA,"}}$"), "$sigma_q$")
  }
  
  
  pop_vars <- c("K", "Nmin", paste0("N", years))
  depletion_vars <- c("Max_Dep", paste0("status", years),"P50","sigma", paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")
  
  results <- data.frame(matrix(NA, nrow = length(vars), ncol = 8))
  colnames(results) <- c("Parameter","Mean", "Median", "2.5% CI", "25% CI", "75% CI", "97.5% CI", "Unique")
  
  x <- SIR$resamples_output[,vars]
  x$var_N <- sqrt(x$var_N)
  
  # Get posterior of q
  q_posteriors <- list() # Each layer is an index
  
  # -- Determining the number of Indices of Abundance available
  rel.abundance <- SIR$inputs$rel.abundance
  indices <- unique(rel.abundance$Index)
  IA.yrs <- rel.abundance$Year
  N_hat <- SIR$resamples_trajectories[, paste0("N_", IA.yrs)] # Estimates of N within IOA years
  
  # -- Q2 for exponent
  q1_cols <- grep("q_IA1", colnames(SIR$resamples_output)) # Columns of resample Q estimates
  q1_est <- SIR$resamples_output[, q1_cols]
  q1_est <- as.matrix(q1_est, ncol = length(q1_cols))
  
  q2_cols <- grep("q_IA2", colnames(SIR$resamples_output)) # Columns of resample Q estimates
  q2_est <- SIR$resamples_output[, q2_cols]
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
  for(j in 1:nrow(SIR$resamples_trajectories)){
    # -- Sample q
    q_posteriors_tmp <- exp(MASS::mvrnorm(
      n = 5,
      mu = as.numeric(log(rel.abundance$IA.obs/(N_hat[j,] ^ (1+q2_est[j,rel.abundance$Index]))) - diag(rel.var.covar.wide)/2),
      Sigma = rel.var.covar.wide))
    
    # q_est <- exp(sum(rel.hess.wide %*% as.numeric(log(rel.abundance$IA.obs/N_hat[j,] ^ (q2_est[j,rel.abundance$Index] + 1))))/(sum(rel.hess.wide))) # q_i
    
    # -- Assign to list
    for(i in indices){
      if(j == 1){
        q_posteriors[[i]] <- c(q_posteriors_tmp[,which(rel.abundance$Index == i)])
      } else {
        q_posteriors[[i]] <- c(q_posteriors[[i]], c(q_posteriors_tmp[,which(rel.abundance$Index == i)]))
      }
    }
  }
  
  
  # Get summary statistics
  results[,1] <- vars_latex
  results[,2] <- sapply(x, mean)
  results[,3:7] <- t(sapply(x, quantile, probs= c(0.5, 0.025, 0.25, 0.75, 0.975)))
  results[,8] <- sapply(x, function(x) length(unique(x)))
  
  # Update q for posteriors
  posterior_q_results <- data.frame(matrix(NA, nrow = num.IA, ncol = 8))
  posterior_q_results[,1] <- paste0("$p(q)_{flt", 1:num.IA, "}$")
  posterior_q_results[,2] <- round(sapply(q_posteriors, mean),3)
  posterior_q_results[,3:7] <- round(t(sapply(q_posteriors, quantile, probs= c(0.5, 0.025, 0.25, 0.75, 0.975))), 3)
  colnames(posterior_q_results) <- c("Parameter","Mean", "Median", "2.5% CI", "25% CI", "75% CI", "97.5% CI", "Unique")
  
  # Format things
  results[c(1,3:6),2:7] <- round(results[c(1,3:6),2:7], 3)
  results[which(vars %in% depletion_vars),2:7] <- round(results[which(vars %in% depletion_vars),2:7], 3)
  results[which(vars %in% pop_vars),2:7] <- format(round(results[which(vars %in% pop_vars),2:7], 0),big.mark=",",scientific=FALSE)
  results[,8] <- format(round(results[,8], 0),big.mark=",",scientific=FALSE)
  results <- rbind(results, posterior_q_results)
  
  if(!is.null(file_name)){
    write.csv(results, file = paste0(file_name, "_summary_table.csv"))
  }
  return(results)
}