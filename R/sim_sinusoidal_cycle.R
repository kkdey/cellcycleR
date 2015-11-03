#' @title Simulate sinusoidal gene expression along cell phases on cell cycle
#'
#'@description The function simulates sinusoidal gene patterns over a cell cycle with the cell phases
#'              provided by the user.
#' @param num_genes The number of sinusoidal genes to simulate
#' @param amp_genes The amplitude vector of the sinusoidal genes (of length equal to num_genes)
#' @param phi_genes The phase angle vector of the sinusoidal genes (of length equal to num_genes)
#' @param sigma_genes The noise variation vector of the sinusoidal genes (of length equal to num_genes)
#' @param cell_times The phases of the cellson the cell cycle (a vector of values between 0 to 2 pi degrees)
#'        The length of cell_times corresponds to number of single cells considered for simulation (say N).
#'
#'
#' @return It returns a matrix of size N x num_genes, of sinusoidal gene expression patterns
#'         across the cells.
#' @author  Kushal K Dey
#' @export
#' @examples
#'  G <- 500;
#'  num_cells <- 400;
#'  amp_genes <- rep(10, G);
#'  phi_genes <- runif(G, 0, 2*pi)
#'  sigma_genes <- rchisq(G, 4);
#'  cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);
#'  cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);
#'




sim_sinusoidal_cycle <- function(num_genes, amp_genes, phi_genes, sigma_genes, cell_times)
{
  signal <- matrix(0, length(cell_times), G);
  for(s in 1:length(cell_times))
  {
    signal[s,] <- mapply(rnorm, 1, amp_genes * sin(cell_times[s] + phi_genes), sigma_genes);
  }
  return(signal)
}
