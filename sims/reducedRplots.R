# Imports
library("ggplot2")
require(gridExtra)

# Relevant parameters
p <- 0.95
nReals <- 2500


# ====================================================================================================
# Simulation 2: Moving circles closer together
# ====================================================================================================
# Plot: Number of observations vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim2_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim2/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim2_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim2_data <- sim2_data[c("n","Distance",toString(p))]

# Sort the unique n
n <- sort(unique(sim2_data$n))
reduced_n <-c(100,300,500)

# Sort the unique distances
d <- sort(unique(sim2_data$Distance))
reduced_d <- c(0,20,40)

# Binomial confidence line
bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = n, y = rep(p,length(n)))
uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Reduce to just for some n
tmp <- sim2_data[(sim2_data$Distance==reduced_d[1]),]

# sort by distance
tmp <- tmp[order(tmp$n),]

xmin <- 0
xmax <- 500

ymin <- 1-2*(1-p)
ymax <- 1
  
# Loop through and add the remaining n
for (d in reduced_d[2:length(reduced_d)]){
  
  # Reduce to just for some n
  tmp2 <- sim2_data[(sim2_data$Distance==d),]
  
  # sort by distance
  tmp2 <- tmp2[order(tmp2$n),]
  
  # sort by distance
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$Distance <- as.factor(tmp$Distance)

# Save line parameters 
tmp$truep <- p
tmp$lowp <- p-bin_conf
tmp$uppp <- p+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
sim2_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
  geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(40,500) + ylim(0.925,0.99) + 
  scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
  geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  labs(title = 'Simulation 1: Varying Distance, Circle Signal', x = 'Number of Observations', y = 'Observed Coverage') +
  theme(plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"), plot.title = element_text(size=12), plot.subtitle = element_text(size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# ====================================================================================================
# Simulation 8: Moving squares closer together
# ====================================================================================================
# Plot: Number of observations vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim8_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim8/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim8_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim8_data <- sim8_data[c("n","Distance",toString(p))]

# Sort the unique n
n <- sort(unique(sim8_data$n))
reduced_n <-c(100,300,500)

# Sort the unique distances
d <- sort(unique(sim8_data$Distance))
reduced_d <- c(0,20,40)

# Binomial confidence line
bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = n, y = rep(p,length(n)))
uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Reduce to just for some n
tmp <- sim8_data[(sim8_data$Distance==reduced_d[1]),]

# sort by distance
tmp <- tmp[order(tmp$n),]

xmin <- 0
xmax <- 500

ymin <- 1-2*(1-p)
ymax <- 1
  
# Loop through and add the remaining n
for (d in reduced_d[2:length(reduced_d)]){
  
  # Reduce to just for some n
  tmp2 <- sim8_data[(sim8_data$Distance==d),]
  
  # sort by distance
  tmp2 <- tmp2[order(tmp2$n),]
  
  # sort by distance
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$Distance <- as.factor(tmp$Distance)

# Save line parameters 
tmp$truep <- p
tmp$lowp <- p-bin_conf
tmp$uppp <- p+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
sim8_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
  geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(40,500) + ylim(0.925,0.99) + 
  scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
  geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  labs(title = 'Simulation 4: Squares Moving', x = 'Number of Observations', y = 'Observed Coverage') +
  theme(plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"), plot.title = element_text(size=12), plot.subtitle = element_text(size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# ====================================================================================================
# Save Sim 2 plot
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim1.png', width = 700, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim2_n_vs_cov, ncol=1, nrow=1)
dev.off()

# ====================================================================================================
# Save Sim 8 plot
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim4.png', width = 700, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim8_n_vs_cov, ncol=1, nrow=1)
dev.off()

# ====================================================================================================
# Combine Sim 2 and 8 plots
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim1and4.png', width = 1800, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim2_n_vs_cov, sim8_n_vs_cov, ncol=2, nrow=1)
dev.off()

# ====================================================================================================
# Simulation 18: Varying corrleation, square signal
# ====================================================================================================
# Plot: Correlation vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim18_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim18/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim18_data) <- c('cfgId','n', 'Correlation','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim18_data <- sim18_data[c("n","Correlation",toString(p))]

# Sort the unique n
n <- sort(unique(sim18_data$n))
reduced_n <-c(100,300,500)

# Sort the unique Correlations
corr <- sort(unique(sim18_data$Correlation))
reduced_corr <- c(-1,-0.5,0, 0.5, 1)

# Binomial confidence line
bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = corr, y = rep(p,length(corr)))
uppline <- data.frame( x = corr, y = rep(p+bin_conf,length(corr)))
lowline <- data.frame( x = corr, y = rep(p-bin_conf,length(corr)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Need to round the correlations due to casting issues
sim18_data['Correlation'] <- round(sim18_data['Correlation'],digits=2)

# Reduce to just for some n
tmp <- sim18_data[(sim18_data$n==reduced_n[1]),]

# sort by Correlation
tmp <- tmp[order(tmp$Correlation),]

xmin <- 0
xmax <- 50

ymin <- 1-2*(1-p)
ymax <- 1
  
# Loop through and add the remaining n
for (n in reduced_n[2:length(reduced_n)]){
  
  # Reduce to just for some n
  tmp2 <- sim18_data[(sim18_data$n==n),]
  
  # sort by Correlation
  tmp2 <- tmp2[order(tmp2$Correlation),]
  
  # sort by Correlation
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$n <- as.factor(tmp$n)

# Save line parameters 
tmp$truep <- p
tmp$lowp <- p-bin_conf
tmp$uppp <- p+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
sim18_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=`0.95`, group=`n`, color=`n`)) + 
  geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(-1,1) + ylim(0.925,0.99) + 
  scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
  geom_line(aes(x=`Correlation`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  labs(title = 'Simulation 2: Varying Correlation, Square Signal', x = 'Correlation Between Noise Fields', y = 'Observed Coverage') + 
  theme(plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"), plot.title = element_text(size=12), plot.subtitle = element_text(size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))

# ====================================================================================================
# Simulation 20: Varying Signal Gradient, Ramp signal
# ====================================================================================================
# Plot: Gradient vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim20_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim20/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim20_data) <- c('cfgId','n', 'Gradient','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim20_data <- sim20_data[c("n","Gradient",toString(p))]

# Sort the unique n
n <- sort(unique(sim20_data$n))
reduced_n <-c(100,300,500)

# Sort the unique Gradients
grad <- sort(unique(sim20_data$Gradient))
reduced_grad <- c(2,4,6)/4

# Binomial confidence line
bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = grad, y = rep(p,length(grad)))
uppline <- data.frame( x = grad, y = rep(p+bin_conf,length(grad)))
lowline <- data.frame( x = grad, y = rep(p-bin_conf,length(grad)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Need to round the Gradients due to casting issues
sim20_data['Gradient'] <- round(sim20_data['Gradient'],digits=3)

# Reduce to just for some n
tmp <- sim20_data[(sim20_data$n==reduced_n[1]),]

# sort by Gradient
tmp <- tmp[order(tmp$Gradient),]

xmin <- 0
xmax <- 50

ymin <- 1-2*(1-p)
ymax <- 1
  
# Loop through and add the remaining n
for (n in reduced_n[2:length(reduced_n)]){
  
  # Reduce to just for some n
  tmp2 <- sim20_data[(sim20_data$n==n),]
  
  # sort by Gradient
  tmp2 <- tmp2[order(tmp2$Gradient),]
  
  # sort by Gradient
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$n <- as.factor(tmp$n)

# Save line parameters 
tmp$truep <- p
tmp$lowp <- p-bin_conf
tmp$uppp <- p+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
sim20_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=`0.95`, group=`n`, color=`n`)) + 
  geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(0.5,3.5) + ylim(0.925,0.99) + 
  scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
  geom_line(aes(x=`Gradient`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal', x = 'Gradient (per 50 voxels)', y = 'Observed Coverage') +
  theme(plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm"), plot.title = element_text(size=12), plot.subtitle = element_text(size=12), axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))


# ====================================================================================================
# Save Sim 18 plot
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim2.png', width = 700, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim18_corr_vs_cov, ncol=1, nrow=1)
dev.off()

# ====================================================================================================
# Save Sim 20 plot
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim3.png', width = 700, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim20_grad_vs_cov, ncol=1, nrow=1)
dev.off()

# ====================================================================================================
# Combine Sim 18 and 20 plots
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim2and3.png', width = 1800, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim18_corr_vs_cov, sim20_grad_vs_cov, ncol=2, nrow=1)
dev.off()


# ====================================================================================================
# Combine three plots
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim_all_col.png', width = 700, height = 1100,
    units = "px", pointsize = 10, bg = "white", res = 120)
par(
  cex.axis = 2,
  cex.lab  = 2
)
grid.arrange(sim2_n_vs_cov, sim18_corr_vs_cov, sim20_grad_vs_cov, ncol=1, nrow=3)
dev.off()


# ====================================================================================================
# Combine three plots
# ====================================================================================================

png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim_all_row.png', width = 1950, height = 400,
    units = "px", pointsize = 10, bg = "white", res = 120)
grid.arrange(sim2_n_vs_cov, sim18_corr_vs_cov, sim20_grad_vs_cov, ncol=3, nrow=1)
dev.off()