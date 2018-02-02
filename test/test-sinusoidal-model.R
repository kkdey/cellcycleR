# Description
#   Test all functions involved in fitting sin_cell_ordering_class.
# Goal is to check if all run without errors.
#
# Date: 2017-12-16

devtools::install_github("jhsiao999/cellcycleR", ref="0.1.3")

library(cellcycleR)

G <- 500
num_cells <- 400
amp_genes <- rep(10, G)
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4)
cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE)

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim)



# fitting sin_cell_ordering_iter
celltime_levels <- 100
celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels - 1))
cell_times_previous <- sample(celltimes_choice, num_cells, replace=TRUE)
fit_iter <- sin_cell_ordering_iter(cycle_data=cycle_data,
                                   celltime_levels=celltime_levels,
                                   cell_times_iter=cell_times_previous,
                                   freq=1,
                                   fix.phase=FALSE, phase_in=NULL,
                                   n_cores=3)

# fitting sin_cell_ordering_class
# no parallel computing
fit_lapply <- sin_cell_ordering_class(cycle_data=cycle_data,
                               celltime_levels=celltime_levels,
                               num_shuffle=1,
                               fix.phase=FALSE,
                               verbose = TRUE, freq = 1,
                               tol = 1e-6, maxiter = 500,
                               n_cores=3)
str(fit)

# parallele computing
fit_mclapply <- sin_cell_ordering_class(cycle_data=cycle_data,
                               celltime_levels=celltime_levels,
                               num_iter=2,
                               verbose = FALSE, freq = 1,
                               tol = .01, maxiter = 500,
                               parallel = TRUE, n_cores=2)
str(fit_mclapply)
