# cellcycleR

cellcycleR is a R package for cell cycle phasing and analysis

To install the package, follow the commands below 

```
install.packages("devtools")
library(devtools)
install_github('kkdey/cellcycleR')
```

To load the package

```
library(cellcycleR)
```
For the main workhorse function, the cell ordering based on classes, check

```
?cell_ordering_class
```

Example fit on a simulated cell cycle data. The simulation model

```
G <- 500;
num_cells <- 400;
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);
cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);
```

The cell ordering into time classes on the cell cycle is given by 

```
out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100)
```

The full order of the cells can then be obtained by 

```
cell_ordering_full(out$signal_intensity, G);
```


