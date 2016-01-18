#' @title Nonparametric cell ordering into different phases on the cell cycle
#'
#' @description The function runs a Gibbs sampler for a number of iterations to obtain the estimates
#'             of cell times on cell cycle (modulo rotation) as well as the gene features (smoothing fit)
#'
#' @param cycle_data: a N x G matrix, where N is number of cells, G number of genes
#' @param celltime_levels:  The number of phase classes chosen (The deafult is 100). It splits up
#'                          circle or the range (0, 2pi) into celltime_levels uniform splits.
#' @param num_iter The number of iterations to run for the Gibbs sampler scheme
#' @param save_path The file path to save the RDA file containing the information on gene characteristics
#'                   and the cell ordering times. Default is NULL in which case, it will not save output
#' @return Returns a list containing the following items
#'  \item{cell_times}{estimated cell times from Gibbs sampler}
#'  \item{sigma}{The estimated non-signal or noise variation of the genes}
#'  \item{loglik}{The model log likelihood of the fit}
#'  \item{signal_intensity}{The intensity of each cell in each of the cell times classes}
#'
#' @author  Kushal K Dey
#' @export


np_cell_ordering_class <- function(cycle_data, celltime_levels, num_iter, method=c("LOESS", "B-spline", "Wavelet"), save_path=NULL)
{
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];

  celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  cell_times_init <- sample(celltimes_choice, numcells, replace=TRUE);

  cell_times_iter <- cell_times_init;

  for(iter in 1:num_iter)
  {
    fun <- np_cell_ordering_iter(cycle_data, celltime_levels, cell_times_iter, method=c("LOESS", "B-spline", "Wavelet"));
    cell_times_iter <- fun$cell_times_iter;
    fitted_signal <- fun$fitted_signal;
    signal_intensity_iter <- fun$signal_intensity_iter;
    sigma_iter <- fun$sigma_iter;
    loglik_iter <- np_loglik_cellcycle(cycle_data, cell_times_iter, fitted_signal, sigma_iter);
    cat("The loglikelihood after iter", iter, "is:", loglik_iter,"\n")
  }

  out <- list("cell_times"=cell_times_iter,
              "fitted_signal"=fitted_signal,
              "sigma"=sigma_iter,
              "loglik"=loglik_iter,
              "signal_intensity"=signal_intensity_iter)

  if(!is.null(save_path)){
    save(out,file=save_path);
  }
  return(out)
}

