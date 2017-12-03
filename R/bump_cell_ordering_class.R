#' @title Cell ordering into different phases on the cell cycle using bump function features (genes)
#'
#' @description The function runs a Gibbs sampler for a number of iterations to obtain the estimates
#'             of cell times (modulo rotation) when the genes vary as per a bump functin.
#'             The function also returns gene features (amplitude, phase and noise variation)
#'
#' @param cycle_data: a N x G matrix, where N is number of cells, G number of genes
#' @param celltime_levels:  The number of phase classes chosen (The deafult is 100). It splits up
#'                          circle or the range (0, 2pi) into celltime_levels uniform splits.
#' @param num_iter The number of iterations to run for the Gibbs sampler scheme
#' @param save_path The file path to save the RDA file containing the information on gene characteristics
#'                   and the cell ordering times. Default is NULL in which case, it will not save output
#' @param start the starting time slots of the cells.
#' @param verbose if TRUE, prints the loglikelihood at each step/iteration. If FALSE,just prints final loglikelihood.
#'
#' @return Returns a list containing the following items
#'  \item{cell_times}{estimated cell times from Gibbs sampler}
#'  \item{amp}{The estimated amplitudes of the genes}
#'  \item{phi}{The estimated phase angles of the genes}
#'  \item{sigma}{The estimated non sinusoid signal variation of the genes}
#'  \item{loglik}{The model log likelihood of the fit}
#'  \item{signal_intensity}{The intensity of each cell in each of the cell times classes}
#'
#' @author  Kushal K Dey
#' @export


bump_cell_ordering_class <- function(cycle_data, celltime_levels=100, num_iter=100,
                                    save_path=NULL,
                                    start = NULL,
                                    verbose = FALSE)
{

  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];

  celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  if(is.null(start)){
    cell_times_init <- sample(celltimes_choice, numcells, replace=TRUE);
  }else{
    cell_times_init <- start
  }

  cell_times_iter <- cell_times_init;

  for(iter in 1:num_iter)
  {
    fun <- bump_cell_ordering_iter(cycle_data, celltime_levels, cell_times_iter);
    cell_times_iter <- fun$cell_times_iter;
    amp_iter <- fun$amp_iter;
    phi_iter <- fun$phi_iter;
    sigma_iter <- fun$sigma_iter;
    signal_intensity_iter <- fun$signal_intensity_iter;
    cycle_scaled_data <- fun$scaled_data
    loglik_iter <- bump_loglik_cellcycle(cycle_scaled_data,
                                        cell_times_iter, amp_iter, phi_iter,
                                        sigma_iter);
    if (verbose == TRUE) {
      cat("The loglikelihood after iter", iter, "is:", loglik_iter,"\n")
    }
  }

  if (verbose == FALSE) {
    cat("Final loglikelihood, iter", num_iter, ":", loglik_iter,"\n")
  }
  out <- list("cell_times"=cell_times_iter,
              "amp"=amp_iter,
              "phi"=phi_iter,
              "sigma"=sigma_iter,
              "loglik"=loglik_iter,
              "signal_intensity"=signal_intensity_iter)

  if(!is.null(save_path)){
    save(out,file=save_path);
  }
  return(out)
}
