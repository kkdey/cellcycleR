bump_cell_ordering_iter <- function(cycle_data, celltime_levels, cell_times_iter)
{
  if(length(unique(cell_times_iter))==1)
    stop("All the points have converged at same point on cycle");

  # cycle_data: a N \times G matrix, where N is number of cells, G number of genes
  # cell_times_iter:  the vector of cell times taken as input (a N \times 1)

  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  sigma <- array(0,G);
  amp <- array(0,G);
  phi <- array(0,G);

  # Fit linear models for each gene $g$ given the cell times [ linear model depends on fix.phase]

    lmfit_list <- parallel::mclapply(1:G, function(g)
    {
      temp1 <- scales::rescale(cycle_data[,g], to=c(0,1))
      fit1 <- lm(temp1  ~ sin(0.5*cell_times_iter) + cos(0.5*cell_times_iter) -1);
      temp2 <- scales::rescale(cycle_data[,g], to = c(-1, 0))
      fit2 <- lm(temp2  ~ sin(0.5*cell_times_iter) + cos(0.5*cell_times_iter) -1);

      if(sum(fit1$residuals^2) < sum(fit2$residuals^2)){
        fit <- fit1
        scale <- 0
      }else{
        fit <- fit2
        scale <- 1
      }

      out_sigma <- sd(fit$residuals);
      beta1 <- fit$coefficients[1];
      beta2 <- fit$coefficients[2];
      if(beta1==0 & beta2==0){
        stop(paste0("You have a gene with all 0 counts at gene",g));
      }
      out_amp <- sqrt(beta1^2 + beta2^2);
      out_phi <- atan3(as.numeric(beta2), as.numeric(beta1));
      ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma,
                 "out_scale" = scale)
      return(ll)
    }, mc.cores=parallel::detectCores())

    amp <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_amp))));
    phi <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_phi))));
    sigma <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_sigma))));
    scale <- as.numeric(unlist(lapply(1:length(lmfit_list), function(n) return(lmfit_list[[n]]$out_scale))));

    cycle_data_scaled <- matrix(0, dim(cycle_data)[1], dim(cycle_data)[2])
    for(l in 1:dim(cycle_data_scaled)[2]){
      if(scale[l] == 0){
        cycle_data_scaled[,l] <- scales::rescale(cycle_data[,l], to=c(0,1))
      }else{
        cycle_data_scaled[,l] <- scales::rescale(cycle_data[,l], to=c(-1,0))
      }
    }

  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  num_celltime_class <- length(cell_times_class);

  sin_class_times <- sin(0.5*cell_times_class);
  cos_class_times <- cos(0.5*cell_times_class);
  sin_phi_genes <- sin(phi);
  cos_phi_genes <- cos(phi);
  sinu_signal <- cbind(sin_class_times, cos_class_times) %*% rbind(amp*cos_phi_genes, amp*sin_phi_genes);
  options(digits=12)
  signal_intensity_per_class <- matrix(0, numcells, num_celltime_class)

  signal_intensity_per_class <- do.call(rbind,parallel::mclapply(1:numcells, function(cell)
  {
    res_error <- sweep(sinu_signal,2,cycle_data_scaled[cell,]);
    res_error_adjusted <- -(res_error^2);
    res_error_adjusted <- sweep(res_error_adjusted, 2, 2*sigma^2, '/');
    out <- rowSums(sweep(res_error_adjusted,2,log(sigma)) - 0.5*log(2*pi));
    return(out)
  }, mc.cores=parallel::detectCores()));


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

  out <- list("cell_times_iter"=cell_times, "amp_iter"=amp,
              "phi_iter"=phi, "sigma_iter"=sigma, "scaled_data" = cycle_data_scaled,
              "signal_intensity_iter"=signal_intensity_per_class);
  return(out)
}
