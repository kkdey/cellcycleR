
library(cellcycleR)
library(wavethresh)

G <- 100;
num_cells <- 256;
phi_genes <- seq(1, pi, length.out=G);
sigma_genes <- rchisq(G, 0.01);
cell_times_sim <- sort(sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE));
cycle_data <- matrix(0,num_cells,G)
# base <- c(1:(num_cells/4), (num_cells/4):1, 1:(num_cells/4), (num_cells/4):1);
base <- c(1:(num_cells/4), (num_cells/4):1, rep(0,num_cells/4), rep(0,num_cells/4));
for(g in 1:G){
cycle_data[,g] <- rnorm(num_cells,shift(base,2*g, dir="right"),sigma_genes[g]);
}

plot(cycle_data[,1], type="l")
plot(cycle_data[,60], type="l")
plot(cycle_data[,30], type="l")


celltime_levels <- 256;

sample_reorder <- sample(1:num_cells,num_cells, replace=FALSE);
cell_times_reorder <- cell_times_sim[sample_reorder];
cycle_data_reorder <- cycle_data[sample_reorder,];

plot(cycle_data_reorder[,60], type="l")

###  B-spline smoothing

system.time(out <- np_cell_ordering_class(cycle_data_reorder, celltime_levels = 256, num_iter=500, method="B-spline"))


###### Post processing

plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_reorder)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_reorder),radial.pos=sort(cell_times_reorder),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_reorder)), lwd=2)

plot(cycle_data_reorder[order(out$cell_times),1], type="l")
plot(cycle_data[,1],type="l")
plot(cycle_data_reorder[order(out$cell_times),7], type="l")
plot(cycle_data[,7],type="l")
plot(cycle_data_reorder[order(out$cell_times),60], type="l")
plot(cycle_data[,60],type="l")
plot(cycle_data_reorder[order(out$cell_times),80], type="l")
plot(cycle_data[,80],type="l")



### LOWESS smoothing

system.time(out <- np_cell_ordering_class(cycle_data_reorder, celltime_levels = 128, num_iter=100, method="LOESS"))


###### Post processing

plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_reorder)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_reorder),radial.pos=sort(cell_times_reorder),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_reorder)), lwd=2)

plot(cycle_data_reorder[order(out$cell_times),1], type="l")
plot(cycle_data[,1],type="l")

### Wavelet smoothing

system.time(out <- np_cell_ordering_class(cycle_data_reorder, celltime_levels = 256, num_iter=500, method="Wavelet"))


###### Post processing

plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_reorder)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_reorder),radial.pos=sort(cell_times_reorder),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_reorder)), lwd=2)

plot(cycle_data_reorder[order(out$cell_times),1], type="l")
plot(cycle_data[,1],type="l")
plot(cycle_data_reorder[order(out$cell_times),7], type="l")
plot(cycle_data[,7],type="l")
plot(cycle_data_reorder[order(out$cell_times),60], type="l")
plot(cycle_data[,60],type="l")
plot(cycle_data_reorder[order(out$cell_times),80], type="l")
plot(cycle_data[,80],type="l")

### Sinusoidal smoothing

system.time(out <- sin_cell_ordering_class(cycle_data_reorder, celltime_levels = 256, num_iter=500))

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_reorder)],
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_reorder),radial.pos=sort(cell_times_reorder),
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_reorder)), lwd=2)

plot(cycle_data_reorder[order(out$cell_times),1], type="l")
plot(cycle_data[,1],type="l")
plot(cycle_data_reorder[order(out$cell_times),7], type="l")
plot(cycle_data[,7],type="l")
plot(cycle_data_reorder[order(out$cell_times),60], type="l")
plot(cycle_data[,60],type="l")
plot(cycle_data_reorder[order(out$cell_times),80], type="l")
plot(cycle_data[,80],type="l")




