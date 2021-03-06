
## Testing the package cellcycleR

library(devtools)
install_github('kkdey/cellcycleR')
library(cellcycleR)

G <- 500;
num_cells <- 400;
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);

celltime_levels <- 100;

system.time(out <- sin_cell_ordering_class(cycle_data, celltime_levels = 100, num_iter=100))


## Post processing

plot(amp_genes, out$amp, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, out$phi, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_sim)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_sim),radial.pos=sort(cell_times_sim),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)

