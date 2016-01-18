#' @title Cell ordering into different phases on the cell cycle
#'
#' @description The function runs a Gibbs sampler for a number of iterations to obtain the estimates
#'             of cell times on cell cycle (modulo rotation) as well as the gene features (amplitude,
#'              phase and noise variation)
#'
#' @param cycle_data: a N x G matrix, where N is number of cells, G number of genes
#' @param celltime_levels:  The number of phase classes chosen (The deafult is 100). It splits up
#'                          circle or the range (0, 2pi) into celltime_levels uniform splits.
#' @param num_iter The number of iterations to run for the Gibbs sampler scheme
#' @param save_path The file path to save the RDA file containing the information on gene characteristics
#'                   and the cell ordering times. Default is NULL in which case, it will not save output
#' @param fix.phase if TRUE, the phase will be fixed in inference for the genes, default is FALSE
#' @param phase_in if fix.phase is TRUE, then phase_in is G x 1 vector of user input gene phases.
#'        Default is NULL as is the case if fix.phase=FALSE.
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
#' @examples
#'
#' G <- 500;
#' num_cells <- 400;
#' amp_genes <- rep(10, G);
#' phi_genes <- runif(G, 0, 2*pi)
#' sigma_genes <- rchisq(G, 4);
#' cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);
#' cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);
#' celltime_levels <- 100;
#' out <- sin_cell_ordering_class(cycle_data, celltime_levels = 100, num_iter=100)

sin_cell_ordering_class <- function(cycle_data, celltime_levels, num_iter, save_path=NULL,
                                  fix.phase=FALSE, phase_in=NULL)
{
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];

  celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  cell_times_init <- sample(celltimes_choice, numcells, replace=TRUE);

  cell_times_iter <- cell_times_init;

  for(iter in 1:num_iter)
  {
    fun <- sin_cell_ordering_iter(cycle_data, celltime_levels, cell_times_iter, fix.phase, phase_in);
    cell_times_iter <- fun$cell_times_iter;
    amp_iter <- fun$amp_iter;
    phi_iter <- fun$phi_iter;
    sigma_iter <- fun$sigma_iter;
    signal_intensity_iter <- fun$signal_intensity_iter;
    loglik_iter <- loglik_cell_cycle(cycle_data, cell_times_iter, amp_iter, phi_iter, sigma_iter);
    cat("The loglikelihood after iter", iter, "is:", loglik_iter,"\n")
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
