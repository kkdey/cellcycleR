
## testing the code for nonparametric cell cycle re-ordering

library(cellcycleR)
library(wavethresh)

G <- 100;
num_cells <- 256;
amp_genes <- rep(10, G);
phi_genes <- rep(c(2,5), each=G/2);
sigma_genes <- rchisq(G, 4);
cell_times_sim <- sort(sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE));
cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);

celltime_levels <- 128;

system.time(out <- np_cell_ordering_class(cycle_data, celltime_levels = 100, num_iter=100, method="B-spline"))


## Post processing

plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_sim)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_sim),radial.pos=sort(cell_times_sim),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)

