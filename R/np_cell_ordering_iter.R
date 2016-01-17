nonparam_cell_ordering_iter <- function(cycle_data, celltime_levels, cell_times_iter, method=c("LOESS", "B-spline", "Wavelet"))
{
  # cycle_data: a N \times G matrix, where N is number of cells, G number of genes
  # cell_times_iter:  the vector of cell times taken as input (a N \times 1)

  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));

  np_signal <- do.call(cbind, parallel::mclapply(1:G, function(g)
                                  {

                                      if(method=="LOESS"){
                                              ordered_vec <- as.numeric(tapply(cycle_data[order(cell_times_iter),g], factor(sort(cell_times_iter)), mean));
                                              ordered_vec_out <- approx(unique(sort(cell_times_iter)), ordered_vec, xout = cell_times_class, ties = "ordered")$y
                                              fit <- loess(ordered_vec_out ~ cell_times_class)$fitted
                                      }
                                      if(method=="B-spline"){
                                              ordered_vec <- as.numeric(tapply(cycle_data[order(cell_times_iter),g], factor(sort(cell_times_iter)), mean));
                                              ordered_vec_out <- approx(unique(sort(cell_times_iter)), ordered_vec, xout = cell_times_class, ties = "ordered")$y
                                              fit <- smooth.spline(cell_times_class, ordered_vec_out)$y
                                      }
                                      if(method=="Wavelet"){
                                              ordered_vec <- as.numeric(tapply(cycle_data[order(cell_times_iter),g], factor(sort(cell_times_iter)), mean));
                                              ordered_vec_out <- approx(unique(sort(cell_times_iter)), ordered_vec, xout = cell_times_class, ties = "ordered")$y
                                              fit <-  wr(threshold(wd(ordered_vec_out), type="soft"));
                                      }
                                      return(fit)
  }, mc.cores=parallel::detectCores()))

  options(digits=12)
  signal_intensity_per_class <- matrix(0, numcells, celltime_levels)

  signal_intensity_per_class <- do.call(rbind,parallel::mclapply(1:numcells, function(cell)
  {
    res_error <- sweep(np_signal,2,cycle_data[cell,]);
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

  out <- list("cell_times_iter"=cell_times, "signal_intensity_iter"=signal_intensity_per_class, "fitted_signal"=np_signal);
  return(out)
}





