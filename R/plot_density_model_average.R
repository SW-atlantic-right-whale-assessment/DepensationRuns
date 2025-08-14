

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
plot_density_ma <- function(SIR, posteriors_lwd = rep(3, length(SIR)), posteriors_lty = rep(1, length(SIR)), posteriors_col = rep(1, length(SIR)),  file_name = NULL, lower = NULL, upper = NULL, probs = c(0.025, 0.975) ){
  
  # Make into list
  if(class(SIR) == "SIR"){
    SIR <- list(SIR)
  }
  
  # Vars of interest
  num.IA <- sort(unique(c( sapply(SIR, function(x) x$inputs$rel.abundance$Index))))
  years <- sort(unique(unlist(c( 
    sapply(SIR, function(x) x$inputs$target.Yr),
    sapply(SIR, function(x) x$inputs$output.Years)))))
  vars_resamples <- c("r_max", "K", "Pmsy", "var_N", "Nmin", "Max_Dep")#, paste0("q_IA1", num.IA), paste0("q_IA2", num.IA), "add_VAR_IA", paste0("catch_multiplier_",2:3), "catch_parameter")
  vars_resamples_latex <- c("$r_{max}$", "$K$", "$P_{MSY}$", "$sigma$", "$N_{min}$", "$P_{min}$")#, paste0("$q_{flt", num.IA, "}$"), paste0("$\beta_{q_{flt", num.IA,"}}$"), "$tau_q$", paste0("$SLR_",1:2,"$"), "$pi$")
  
  traj_vars <- c(paste0("N_", years), paste0("status", years))
  traj_vars_latex <- c( paste0("$N_{", years, "}$"),  paste0("$P_{", years,"}$"))
  # vars <- c("r_max", "K", "Pmsy", "var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", num.IA), paste0("q_IA2", num.IA), "add_VAR_IA", paste0("catch_multiplier_",2:3), "catch_parameter")
  # vars_latex <- c("$r_{max}$", "$K$", "$P_{MSY}$", "$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "$P_{min}$", paste0("$P_{", years,"}$"), paste0("$q_{flt", num.IA, "}$"), paste0("$\beta_{q_{flt", num.IA,"}}$"), "$tau_q$", paste0("$SLR_",1:2,"$"), "$pi$")
  
  
  # Only select vars that have multiple unique parameters
  x <- SIR[[1]]$resamples_output[,vars_resamples]
  unique_vars <- sapply(x, function(x) length(unique(x)))
  
  for(i in 1:length(SIR)){
    SIR[[i]]$resamples_output$var_N <- sqrt(SIR[[i]]$resamples_output$var_N)
  }
  
  vars_resamples <- vars_resamples[which(unique_vars > 1)]
  vars_resamples_latex <- vars_resamples_latex[which(unique_vars > 1)]
  
  # - Combine
  vars <- c(vars_resamples, traj_vars)
  vars_latex <- c(vars_resamples_latex, traj_vars_latex)
  
  # -- Get posterior of q
  # Get posterior of q
  q_posteriors <- list() # Each layer is an index
  


  
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
    
    par(mfrow = c(2,ceiling(length(vars)/2) + 2))
    par( mar=c(3, 0.05 , 0.5 , 0.55) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    
    plot.new()
    
    
    # Loop through vars
    for(i in 1:length(vars)){
      
      # Extract posterio densities
      posterior_dens <- list()
      for(k in 1:length(SIR)){
        if(vars[i] %in% paste0("q_IA1", num.IA)){
          posterior_dens[[k]] <- density(q_posteriors[[k]][[as.numeric(unlist(strsplit(paste0("q_IA1", num.IA), "q_IA1"))[2])]])
        } else if(vars[i] %in% paste0("N_", years)){
          posterior_dens[[k]] <- density(as.numeric(as.character(SIR[[k]]$resamples_trajectories[,vars[i]])))
        } else if(vars[i] %in% paste0("status", years)){
          posterior_dens[[k]] <- density(as.numeric(as.character(SIR[[k]]$resamples_trajectories[,vars[i-length(years)]]))/
                                           as.numeric(as.character(SIR[[k]]$resamples_output$K))
                                           )
        } else {
          posterior_dens[[k]] <- density(as.numeric(as.character(SIR[[k]]$resamples_output[,vars[i]])))
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
      
      if(i %in% c(ceiling(length(vars)/2)))  {
        plot.new()
        plot.new()
      }
      
      
      if(i %in% c(1, ceiling(length(vars)/2) + 1) ) {
        mtext(side = 2, "Density", line = 1, cex= 0.75)
      }
      
      
    }
    
    if(j > 1){ dev.off()}
  }
}
