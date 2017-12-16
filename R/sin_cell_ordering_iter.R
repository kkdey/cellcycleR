#' @title Fitting sinusoidal model
#'
#' @param cycle_data a N x G matrix, where N is number of cells, G number of genes
#' @param cell_times a N x 1 vector of cell times (not ordered)
#' @param fix.phase TRUE if cell_times is an ordered vector, i.e., the
#'                  relative position of cells in time is known, and FALSE otherwise.
#'
#' @description Compute sinusoidal model estimates.
#'
#' @author Kushal K Dey, Chiaowen Joyce Hsiao
#'
#' @export
#'
#' @example
# G <- 100;
# num_cells <- 256;
# amp_genes <- rep(10, G);
# phi_genes <- rep(c(2,4,6,8), each=G/4);
# sigma_genes <- rchisq(G, 0.01);
# cell_times_sim <- seq(0,2*pi, length.out=num_cells);
# cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);
# celltime_levels <- 200
# celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels - 1))
# cell_times_iter <- sample(celltimes_choice, dim(cycle_data)[1], replace=TRUE)
# fit <- sin_cell_ordering_iter(cycle_data, celltime_levels,
#                        cell_times_iter,
#                        parallel=FALSE,
#                        fix.phase = FALSE, phase_in=NULL)

sin_cell_ordering_iter <- function(cycle_data,
                                   celltime_levels,
                                   cell_times_iter, freq=1,
                                   fix.phase=FALSE, phase_in=NULL,
                                   parallel = FALSE, n_cores=5)
{
  if(fix.phase==TRUE & is.null(phase_in))
    stop("fix.phase=TRUE and phase not provided")
  if(fix.phase==FALSE & !is.null(phase_in))
    stop("fix.phase=FALSE and phase provided")
  if(length(unique(cell_times_iter))==1)
    stop("All the points have converged at same point on cycle");

  if(parallel==TRUE & is.null(n_cores))
    stop("parallel=TRUE and number of cores not provided")

#  if(n_cores==NULL) {n_cores <- parallel::detectCores()}

  G <- dim(cycle_data)[2]
  numcells <- dim(cycle_data)[1]
  sigma <- array(0,G)
  amp <- array(0,G)
  phi <- array(0,G)

  # Fit linear models for each gene $g$ given the cell times [ linear model depends on fix.phase]

  if(!fix.phase){

    fit_lm_notfixed_phase <- function(g) {
      fit <- lm(cycle_data[,g]  ~ sin(freq*cell_times_iter) + cos(freq*cell_times_iter) -1);
      out_sigma <- sd(fit$residuals);
      beta1 <- fit$coefficients[1];
      beta2 <- fit$coefficients[2];
      if(beta1==0 & beta2==0){
        stop(paste0("You have a gene with all 0 counts at gene",g));
      }
      out_amp <- sqrt(beta1^2 + beta2^2);
      out_phi <- atan3(as.numeric(beta2), as.numeric(beta1));
      ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma)
      return(ll)
    }

    if (parallel) {
      fit_lm_all <- parallel::mclapply(1:G, function(g) fit_lm_notfixed_phase(g), mc.cores=n_cores)
    }
    if (!parallel) {
      fit_lm_all <- lapply(1:G, function(g) fit_lm_notfixed_phase(g))
    }

    amp <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_amp))))
    phi <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_phi))))
    sigma <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_sigma))))
  }

  if(fix.phase){
    phi <- phase_in

    fit_lm_fixed_phase <- function(g) {
      fit <- lm(cycle_data[,g]  ~ sin(freq*cell_times_iter+phi[g]) -1);
      out_sigma <- sd(fit$residuals);
      out_amp <- abs(fit$coefficients[1]);
      out_phi <- phi;
      ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma)
      return(ll)
    }

    if (parallel) {
      fit_lm_all <- parallel::mclapply(1:G, function(g) fit_lm_fixed_phase(g), mc.cores=n_cores)
    }
    if (!parallel) {
      fit_lm_all <- lapply(1:G, function(g) fit_lm_fixed_phase(g))
    }

    amp <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_amp))))
    phi <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_phi))))
    sigma <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_sigma))))
  }


  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  num_celltime_class <- length(cell_times_class);

  sin_class_times <- sin(freq*cell_times_class);
  cos_class_times <- cos(freq*cell_times_class);
  sin_phi_genes <- sin(phi);
  cos_phi_genes <- cos(phi);
  sinu_signal <- cbind(sin_class_times, cos_class_times) %*% rbind(amp*cos_phi_genes, amp*sin_phi_genes);
  options(digits=12)
  signal_intensity_per_class <- matrix(0, numcells, num_celltime_class)

  signal_intensity_per_class <- do.call(rbind, lapply(1:numcells, function(cell)
  {
    res_error <- sweep(sinu_signal,2,cycle_data[cell,]);
    res_error_adjusted <- -(res_error^2);
    res_error_adjusted <- sweep(res_error_adjusted, 2, 2*sigma^2, '/');
    out <- rowSums(sweep(res_error_adjusted,2,log(sigma)) - 0.5*log(2*pi));
    return(out)
  }))


  signal_intensity_class_exp <- do.call(rbind,lapply(1:dim(signal_intensity_per_class)[1], function(x)
  {
    out <- exp(signal_intensity_per_class[x,]- max(signal_intensity_per_class[x,]));
    return(out)
  }));

  cell_times <- cell_times_class[unlist(lapply(1:dim(signal_intensity_class_exp)[1], function(x)
  {
    temp <- signal_intensity_class_exp[x,];
    if(length(unique(signal_intensity_class_exp[x,]))==1)
      out <- sample(1:dim(signal_intensity_class_exp)[2],1)
    else
      out <- which(rmultinom(1,1,signal_intensity_class_exp[x,])==1);
    return(out)
  }))];

  out <- list("cell_times_iter"=cell_times, "amp_iter"=amp, "phi_iter"=phi, "sigma_iter"=sigma, "signal_intensity_iter"=signal_intensity_per_class);
  return(out)
}
