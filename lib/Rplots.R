# Imports
library("ggplot2")


# ====================================================================================================
# Simulations 1: Moving circles closer together (High SNR)
# ====================================================================================================
# Plot: Distance vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim1_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/trueBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim1_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim1_data <- sim1_data[c("n","Distance","0.95")]

# Sort the unique n
n <- sort(unique(sim1_data$n))
reduced_n <-c(100,300,500)

# Sort the unique distances
d <- sort(unique(sim1_data$Distance))
reduced_d <- c(0,20,40)

# Relevant parameters
alpha <- 0.95
p <- 0.95
nReals <- 2500

# Binomial confidence line
bin_conf <- qnorm(alpha)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = d, y = rep(0.95,length(d)))
uppline <- data.frame( x = d, y = rep(0.95+bin_conf,length(d)))
lowline <- data.frame( x = d, y = rep(0.95-bin_conf,length(d)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Reduce to just for some n
tmp <- sim1_data[(sim1_data$n==reduced_n[1]),]

# sort by distance
tmp <- tmp[order(tmp$Distance),]

xmin <- 0
xmax <- 50

ymin <- 0.9
ymax <- 1
  
# Loop through and add the remaining n
for (n in reduced_n[2:length(reduced_n)]){
  
  # Reduce to just for some n
  tmp2 <- sim1_data[(sim1_data$n==n),]
  
  # sort by distance
  tmp2 <- tmp2[order(tmp2$Distance),]
  
  # sort by distance
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$n <- as.factor(tmp$n)

# Save line parameters 
tmp$truep <- 0.95
tmp$lowp <- 0.95-bin_conf
tmp$uppp <- 0.95+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
current_plot <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
  geom_ribbon(aes(ymin=0.95-bin_conf, ymax=0.95+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(0,50) + ylim(0.9,1) + 
  scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
  geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
  geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
  labs(title = 'Simulaton 1: Circles Moving (High SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')

# Show plot
current_plot


# ====================================================================================================
# Simulations 1: Moving circles closer together (High SNR)
# ====================================================================================================
# Plot: Number of observations vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim1_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/trueBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim1_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim1_data <- sim1_data[c("n","Distance","0.95")]

# Sort the unique n
n <- sort(unique(sim1_data$n))
reduced_n <-c(100,300,500)

# Sort the unique distances
d <- sort(unique(sim1_data$Distance))
reduced_d <- c(0,20,40)

# Relevant parameters
alpha <- 0.95
p <- 0.95
nReals <- 2500

# Binomial confidence line
bin_conf <- qnorm(alpha)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = n, y = rep(0.95,length(n)))
uppline <- data.frame( x = n, y = rep(0.95+bin_conf,length(n)))
lowline <- data.frame( x = n, y = rep(0.95-bin_conf,length(n)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Reduce to just for some n
tmp <- sim1_data[(sim1_data$Distance==reduced_d[1]),]

# sort by distance
tmp <- tmp[order(tmp$n),]

xmin <- 0
xmax <- 500

ymin <- 0.9
ymax <- 1
  
# Loop through and add the remaining n
for (d in reduced_d[2:length(reduced_n)]){
  
  # Reduce to just for some n
  tmp2 <- sim1_data[(sim1_data$Distance==d),]
  
  # sort by distance
  tmp2 <- tmp2[order(tmp2$n),]
  
  # sort by distance
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$Distance <- as.factor(tmp$Distance)

# Save line parameters 
tmp$truep <- 0.95
tmp$lowp <- 0.95-bin_conf
tmp$uppp <- 0.95+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
current_plot <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
  geom_ribbon(aes(ymin=0.95-bin_conf, ymax=0.95+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(40,500) + ylim(0.9,1) + 
  scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
  geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
  geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
  labs(title = 'Simulaton 1: Circles Moving (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

# Show plot
current_plot

# ====================================================================================================
# Simulations 2: Moving circles closer together (Low SNR)
# ====================================================================================================
# Plot: Distance vs coverage
# ----------------------------------------------------------------------------------------------------

# Read in data
sim1_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim2/FullResults/trueBdry_intrp.csv',sep=',', header=FALSE)

# Name data
names(sim1_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

# Reduce data
sim1_data <- sim1_data[c("n","Distance","0.95")]

# Sort the unique n
n <- sort(unique(sim1_data$n))
reduced_n <-c(100,300,500)

# Sort the unique distances
d <- sort(unique(sim1_data$Distance))
reduced_d <- c(0,10,20,30,40,50)

# Relevant parameters
alpha <- 0.95
p <- 0.95
nReals <- 2500

# Binomial confidence line
bin_conf <- qnorm(alpha)*sqrt(p*(1-p)/nReals)

# Lines for band
midline <- data.frame( x = d, y = rep(0.95,length(d)))
uppline <- data.frame( x = d, y = rep(0.95+bin_conf,length(d)))
lowline <- data.frame( x = d, y = rep(0.95-bin_conf,length(d)))

# -------------------------------------------------------
# Reformat data
# -------------------------------------------------------

# Reduce to just for some n
tmp <- sim1_data[(sim1_data$n==reduced_n[1]),]

# sort by distance
tmp <- tmp[order(tmp$Distance),]

xmin <- 0
xmax <- 50

ymin <- 0.9
ymax <- 1
  
# Loop through and add the remaining n
for (n in reduced_n[2:length(reduced_n)]){
  
  # Reduce to just for some n
  tmp2 <- sim1_data[(sim1_data$n==n),]
  
  # sort by distance
  tmp2 <- tmp2[order(tmp2$Distance),]
  
  # sort by distance
  tmp <- rbind(tmp,tmp2)
  
}

# Save n
tmp$n <- as.factor(tmp$n)

# Save line parameters 
tmp$truep <- 0.95
tmp$lowp <- 0.95-bin_conf
tmp$uppp <- 0.95+bin_conf

# -------------------------------------------------------
# Make plot
# -------------------------------------------------------

# Create plot
current_plot <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
  geom_ribbon(aes(ymin=0.95-bin_conf, ymax=0.95+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
  xlim(0,50) + ylim(0.9,1) + 
  scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
  geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
  geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
  geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
  labs(title = 'Simulaton 2: Circles Moving (Low SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')

# Show plot
current_plot

  
  