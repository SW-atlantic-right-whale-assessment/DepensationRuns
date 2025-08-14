#' OUTPUT FUNCTION
#'
#' Function that provides a plot of the estimated indices of abundance a SIR  model including: median, 95%
#' credible interval, 90% credible interval, catch, and observed indices of abundance abundance.
#'
#' @param SIR A fit SIR model
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, does not save.
#' @param line_col 
#' @param model_names 
#' @param ioa_names names of indices of abundance used.
#'
#' @return Returns and saves a figure with the IOA trajectories.
#'
#' @export
plot_ioa_medians <- function(SIR, file_name = NULL, ioa_names = NULL, line_col = NULL, model_names = NULL){
  
  rel.abundance <- SIR[[1]]$inputs$rel.abundance
  row_names <- c("mean", "median",
                 "2.5%PI", "97.5%PI",
                 "5%PI", "95%PI",
                 "min", "max", "n")
  
  if(is.null(rel.abundance)){
    stop("SIR model did not include an IOA")
  }
  
  if(is.null(line_col)){
      line_col <- rev(oce::oce.colorsViridis(length(SIR)+1))[-1]
  }
  
  ## Determining the number of Indices of Abundance available
  num.IA <- max(rel.abundance$Index)
  
  ## Make var-covar into wide and tall with cov = 0 for different indices
  rel.var.covar.tall <-  subset(rel.abundance, select = -c(Index,Year,IA.obs))
  rel.var.covar.wide <- rel.var.covar.tall[which(rel.abundance$Index == 1),]
  rel.var.covar.wide <- rel.var.covar.wide[1:nrow(rel.var.covar.wide),1:nrow(rel.var.covar.wide)]
  
  rel.hess.tall <- solve(rel.var.covar.wide[1:nrow(rel.var.covar.wide), 1: nrow(rel.var.covar.wide)])
  
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
  
  rel.abundance$Upper95 <- qlnorm(0.975, mean = log(rel.abundance$IA.obs), sd = sqrt(diag(rel.var.covar.wide)))
  rel.abundance$Lower95 <- qlnorm(0.025, mean = log(rel.abundance$IA.obs), sd = sqrt(diag(rel.var.covar.wide)))
  
  # Predict IOA
  ymax <- c()           # Maximum predicted IOA
  IA.yr.range <- list() # Year range for each IOA
  IA_summary <- list()
  
  q1_cols <- list()
  q1_est <- list()
  q2_cols <- list()
  q2_est <- list()
  indices <- unique(rel.abundance$Index)
  IA.yrs <- rel.abundance$Year
  
  for(k in 1:length(SIR)){
    q1_cols[[k]] <- grep("q_IA1", colnames(SIR[[k]]$resamples_output)) # Columns of resample Q estimates
    q1_est[[k]] <- SIR[[k]]$resamples_output[, q1_cols[[k]]]
    q1_est[[k]] <- as.matrix(q1_est[[k]], ncol = length(q1_cols[[k]]))
    
    q2_cols[[k]] <- grep("q_IA2", colnames(SIR[[k]]$resamples_output)) # Columns of resample Q estimates
    q2_est[[k]] <- SIR[[k]]$resamples_output[, q2_cols[[k]]]
    q2_est[[k]] <- as.matrix(q2_est[[k]], ncol = length(q2_cols[[k]]))
    
    IA_summary[[k]] <- list()
    
    # Predict and calculate summary
    for(i in 1:length(q1_cols[[k]])){ # Loop across indices
      
      # Get IOA specifications
      rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]
      rel.var.covar.wide.sub <- as.matrix(rel.var.covar.tall[which(rel.abundance$Index == i),])
      IA.yrs <- rel.abundance.sub$Year
      IA.yr.range[[i]] <- c((min(IA.yrs)):(max(IA.yrs))) # Range +- 1 of IOA years
      
      # Get N and Q2
      N_hat <- SIR[[k]]$resamples_trajectories[, paste0("N_", IA.yr.range[[i]])] # Estimates of N within IOA years
      q2_tmp <- matrix(q2_est[[k]][,i], ncol = 1)
      
      # Predict
      IA_pread <- q1_est[[k]][,i] * N_hat ^ (1 + q2_tmp) # Use analytic q
      
      # Summarize
      IA_summary[[k]][[i]] <-  matrix(nrow = length(row_names), ncol = dim(IA_pread)[2])
      IA_summary[[k]][[i]][1, ] <- sapply(IA_pread, mean)
      IA_summary[[k]][[i]][2:6, ] <- sapply(IA_pread, quantile, probs= c(0.5, 0.025, 0.975, 0.25, 0.75))
      IA_summary[[k]][[i]][7, ] <- sapply(IA_pread, min)
      IA_summary[[k]][[i]][8, ] <- sapply(IA_pread, max)
      IA_summary[[k]][[i]][9, ] <- sapply(IA_pread, length)
      
      IA_summary[[k]][[i]] <- data.frame(IA_summary[[k]][[i]])
      names(IA_summary[[k]][[i]]) <- paste0("IA", i, "_", IA.yr.range[[i]])
      row.names(IA_summary[[k]][[i]]) <- row_names
      
      # Plot limits
      ymax[i] <- max(unlist(c(IA_summary[[k]][[i]][2:6,], rel.abundance.sub$Lower95, rel.abundance.sub$Upper95, ymax[i]))) # Max of posterior
    }
  }
  
  ymin <- 0
  
  
  for(j in 1:(1 + !is.null(file_name))){
    if(j == 2){
      filename <- paste0(file_name, "_IOA_median_fits", ".png")
      png( file = filename , width=169 / 25.4, height = 100 / 25.4, family = "serif", units = "in", res = 300)
    }
    
    par(mfrow = c(1,length(IA_summary[[1]])), mar=c(3, 3 , 0.5 , 0.3) , oma=c(0 , 0 , 0 , 0), tcl = -0.35, mgp = c(1.75, 0.5, 0))
    
    # Loop through IOA indices
    for(i in 1:length(IA_summary[[1]])){
      rel.abundance.sub <- rel.abundance[which(rel.abundance$Index == i),]
      
      # Plot configuration
      plot(y = NA, x = NA,
           ylim = c(ymin, ymax[i]),
           xlim = c(min(IA.yr.range[[i]]), max(IA.yr.range[[i]])),
           xlab = "Year", ylab = "Relative abundance")
      
      # Relative abundance
      points( x = rel.abundance.sub$Year,
              y = rel.abundance.sub$IA.obs,
              col = 1, pch = 16, cex = 2)
      arrows( x0 = rel.abundance.sub$Year,
              y0 = rel.abundance.sub$Lower95,
              x1 = rel.abundance.sub$Year,
              y1 = rel.abundance.sub$Upper95,
              length=0.05, angle=90, code=3, lwd = 3, col = 1)
      
      # Median
      for(k in 1:length(IA_summary)){
        lines( x = IA.yr.range[[i]], y = IA_summary[[k]][[i]][2, ], col = line_col[k], lwd = 3) # Median
      }
      
      if(!is.null(model_names) & i == 1 ){
        legend(x = 2001, y = 3400, legend = model_names, lwd = 3, col = line_col, bty = "n", cex = .7)
      }
      
      if(!is.null(ioa_names)){
        legend("topright", legend = ioa_names[i] ,bty = "n")
      }
      
    }
    if(j == 2){ dev.off()}
  }
}
