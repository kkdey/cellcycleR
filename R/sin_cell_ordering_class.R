#' @title Cell ordering into different phases on the cell cycle using sinusoidal features (genes)
#'
#' @description This function fits a sinsoidal model to gene expression data and 
#'              computes cell-order estimates on cell-cycle (modulo rotation),
#'              as well as the gene-specific cell-cycle features, including amplitude,
#'              phase and noise variation.
#'
#' @param cycle_data a N x G matrix, where N is number of cells, G number of genes
#' @param celltime_levels  The number of phase classes chosen (The deafult is 100). It splits up
#'                          circle or the range (0, 2pi) into celltime_levels uniform splits.
#' @param num_iter Number of iterations. If not specified, number of iteration
#'                 is set to be maxiter (the maximum number of iteratoins).
#' @param save_path The file path to save the RDA file containing the information on gene characteristics
#'                   and the cell ordering times. Default is NULL in which case, it will not save output
#' @param fix.phase if TRUE, the phase will be fixed in inference for the genes, default is FALSE
#' @param phase_in if fix.phase is TRUE, then phase_in is G x 1 vector of user input gene phases.
#'        Default is NULL as is the case if fix.phase=FALSE.
#' @param maxiter The maximum number of iterations. Default = 500. 
#'        
#' @param freq The frequency of the sinusoidal genes. The default is 1.
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
#' @author  Kushal K Dey, Joyce Hsiao
#' @export

sin_cell_ordering_class <- function(cycle_data, celltime_levels, 
                                    num_iter = NULL, save_path=NULL,
                                    fix.phase=FALSE, phase_in=NULL, 
                                    verbose = FALSE, tol = .01, maxiter = 500)
{
  G <- dim(cycle_data)[2]
  numcells <- dim(cycle_data)[1]

  if (is.null(num_iter)) num_iter <- maxiter
  
  # intialize cell-phase
  # randomly assign cells to (0, 2*pi) with 
  # intervals fixed at 2*pi/(celltime_levels - 1)
  celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels - 1))
  cell_times_previous <- sample(celltimes_choice, numcells, replace=TRUE)

  # initialize converge check
  loglik_previous <- .Machine$double.xmax
  eps <- tol+1
  iter <- 0
  
  while (TRUE) {
    estimates <- sin_cell_ordering_iter(cycle_data, 
                                        celltime_levels, 
                                        cell_times_previous, 
                                        fix.phase, phase_in)
    cell_times_iter <- estimates$cell_times_iter
    amp_iter <- estimates$amp_iter
    phi_iter <- estimates$phi_iter
    sigma_iter <- estimates$sigma_iter
    signal_intensity_iter <- estimates$signal_intensity_iter
    loglik_iter <- sin_loglik_cellcycle(cycle_data, 
                                        cell_times_iter, 
                                        amp_iter, 
                                        phi_iter, 
                                        sigma_iter)
    if (verbose) message("log-likelihood:", loglik_iter)

    eps <- abs(loglik_iter - loglik_previous)/abs(loglik_previous)
    # loop out if converged
    if (!(eps > tol & iter < maxiter)) break
    iter <- iter + 1
    loglik_previous <- loglik_iter
    
    if (verbose) {
        message("Iteration: ", loglik_iter, " eps: ", eps) }
  
    for(iter in 1:num_iter) {
      fun <- sin_cell_ordering_iter(cycle_data, celltime_levels, cell_times_iter,
                                  fix.phase, phase_in, freq);
      cell_times_iter <- fun$cell_times_iter;
      amp_iter <- fun$amp_iter;
      phi_iter <- fun$phi_iter;
      sigma_iter <- fun$sigma_iter;
      signal_intensity_iter <- fun$signal_intensity_iter;
      loglik_iter <- sin_loglik_cellcycle(cycle_data, cell_times_iter, amp_iter, phi_iter, sigma_iter, freq);
      if (verbose == TRUE) {
        cat("The loglikelihood after iter", iter, "is:", loglik_iter,"\n")
      }
    
      cell_times_previous <- cell_times_iter
    }
    }
    output <- list(cell_times = cell_times_iter,
                   amp = amp_iter,
                   phi = phi_iter,
                   sigma = sigma_iter,
                   loglik = loglik_iter,
                   signal_intensity = signal_intensity_iter)

  if( !is.null(save_path) ) { save(output, file=save_path) }

  return(output)
}
