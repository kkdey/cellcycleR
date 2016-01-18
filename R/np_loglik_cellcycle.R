#' @title Loglikelihood of of sinusoidal genes in cell cycle
#'
#' @param cycle_data: a N x G matrix, where N is number of cells, G number of genes
#' @param cell_times : a N x 1 vector of cell times
#' @param fit: a celltime_levels(C) x G matrix of model fit from nonparametric cellcycleR
#' @param sigma: the G x 1 vector of gene variation
#'
#' @description Computes the loglikelihood of all the cells in the cycle under nonparametric smoothing
#'
#'  @author  Kushal K Dey
#'
#'  @export
#'

np_loglik_cellcycle <- function(cycle_data, cell_times, fit, sigma)
{
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  celltime_levels <- dim(fit)[1];

  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));

  sum <- 0;

  for(s in 1:numcells)
  {
    ind <- which(cell_times[s]==cell_times_class);
    sum <- sum + sum(mapply(dnorm, cycle_data[s,],fit[ind,], sigma, log=TRUE));
  }

  return(sum)
}
