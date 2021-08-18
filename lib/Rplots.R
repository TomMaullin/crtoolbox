# Imports
library("ggplot2")
require(gridExtra)

# Simulations we're looking at
simNumbers = '19-20'

# Relevant parameters
p <- 0.95
nReals <- 2500

if (simNumbers=='1-2'){

  # ====================================================================================================
  # Simulation 1: Moving circles closer together (High SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------
  # Read in data
  sim1_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim1_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim1_data <- sim1_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim1_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim1_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim1_data[(sim1_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
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
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim1_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 2: Circles Moving (High SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 1: Moving circles closer together (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim1_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim1_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim1_data <- sim1_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim1_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim1_data$Distance))
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
  tmp <- sim1_data[(sim1_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
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
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim1_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 2: Circles Moving (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 2: Moving circles closer together (Low SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
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
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim2_data[(sim2_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim2_data[(sim2_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim2_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 1: Circles Moving (Low SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 2: Moving circles closer together (Low SNR)
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
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 1: Circles Moving (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')


  # ====================================================================================================
  # Combine Sim 1 and 2 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim1and2.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim2_n_vs_cov, sim2_d_vs_cov, sim1_n_vs_cov, sim1_d_vs_cov, ncol=2, nrow=2)
  dev.off()

}

if (simNumbers=='7-8'){

  # ====================================================================================================
  # Simulation 7: Moving squares closer together (High SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim7_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim7/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim7_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim7_data <- sim7_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim7_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim7_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim7_data[(sim7_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim7_data[(sim7_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim7_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 4: Squares Moving (High SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 7: Moving squares closer together (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim7_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim7/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim7_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim7_data <- sim7_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim7_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim7_data$Distance))
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
  tmp <- sim7_data[(sim7_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
    # Reduce to just for some n
    tmp2 <- sim7_data[(sim7_data$Distance==d),]
    
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
  sim7_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 4: Squares Moving (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 8: Moving squares closer together (Low SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
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
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim8_data[(sim8_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim8_data[(sim8_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim8_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 3: Squares Moving (Low SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 8: Moving squares closer together (Low SNR)
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
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 3: Squares Moving (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 7 and 8 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim3and4.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim8_n_vs_cov, sim8_d_vs_cov, sim7_n_vs_cov, sim7_d_vs_cov, ncol=2, nrow=2)
  dev.off()

}

if (simNumbers=='11-12'){
  # ====================================================================================================
  # Simulation 11: Moving squares closer together (High SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim11_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim11/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim11_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim11_data <- sim11_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim11_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim11_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim11_data[(sim11_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim11_data[(sim11_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim11_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 6: Squares Moving, Heterogeneous Noise (High SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 11: Moving squares closer together (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim11_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim11/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim11_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim11_data <- sim11_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim11_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim11_data$Distance))
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
  tmp <- sim11_data[(sim11_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
    # Reduce to just for some n
    tmp2 <- sim11_data[(sim11_data$Distance==d),]
    
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
  sim11_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 6: Squares Moving, Heterogeneous Noise (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 12: Moving squares closer together (Low SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim12_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim12/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim12_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim12_data <- sim12_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim12_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim12_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim12_data[(sim12_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim12_data[(sim12_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim12_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,50) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 5: Squares Moving, Heterogeneous Noise (Low SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 12: Moving squares closer together (Low SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim12_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim12/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim12_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim12_data <- sim12_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim12_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim12_data$Distance))
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
  tmp <- sim12_data[(sim12_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
    # Reduce to just for some n
    tmp2 <- sim12_data[(sim12_data$Distance==d),]
    
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
  sim12_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 5: Squares Moving, Heterogeneous Noise (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 11 and 12 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim5and6.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim12_n_vs_cov, sim12_d_vs_cov, sim11_n_vs_cov, sim11_d_vs_cov, ncol=2, nrow=2)
  dev.off()

}


if (simNumbers=='15-16'){
  # ====================================================================================================
  # Simulation 15: Moving squares, one smaller (High SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim15_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim15/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim15_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim15_data <- sim15_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim15_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim15_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim15_data[(sim15_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim15_data[(sim15_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim15_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,46) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 8: Squares Moving, One Smaller (High SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 15: Moving squares, one smaller (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim15_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim15/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim15_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim15_data <- sim15_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim15_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim15_data$Distance))
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
  tmp <- sim15_data[(sim15_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
    # Reduce to just for some n
    tmp2 <- sim15_data[(sim15_data$Distance==d),]
    
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
  sim15_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 8: Squares Moving, One Smaller (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 16: Moving squares, one smaller (Low SNR)
  # ====================================================================================================
  # Plot: Distance vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim16_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim16/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim16_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim16_data <- sim16_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim16_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim16_data$Distance))
  reduced_d <- c(0,20,40)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = d, y = rep(p,length(d)))
  uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
  lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Reduce to just for some n
  tmp <- sim16_data[(sim16_data$n==reduced_n[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$Distance),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim16_data[(sim16_data$n==n),]
    
    # sort by distance
    tmp2 <- tmp2[order(tmp2$Distance),]
    
    # sort by distance
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
  sim16_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(0,46) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Distance`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Distance`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 7: Squares Moving, One Smaller (Low SNR)', subtitle = 'Distance vs Observed Coverage', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 16: Moving squares, one smaller (Low SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim16_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim16/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim16_data) <- c('cfgId','n', 'Distance','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim16_data <- sim16_data[c("n","Distance",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim16_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique distances
  d <- sort(unique(sim16_data$Distance))
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
  tmp <- sim16_data[(sim16_data$Distance==reduced_d[1]),]

  # sort by distance
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (d in reduced_d[2:length(reduced_d)]){
    
    # Reduce to just for some n
    tmp2 <- sim16_data[(sim16_data$Distance==d),]
    
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
  sim16_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Distance`, color=`Distance`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 7: Squares Moving, One Smaller (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 15 and 16 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim7and8.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim16_n_vs_cov, sim16_d_vs_cov, sim15_n_vs_cov, sim15_d_vs_cov, ncol=2, nrow=2)
  dev.off()

}



if (simNumbers=='17-18'){
  # ====================================================================================================
  # Simulation 17: Varying corrleation, square signal (High SNR)
  # ====================================================================================================
  # Plot: Correlation vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim17_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim17/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim17_data) <- c('cfgId','n', 'Correlation','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim17_data <- sim17_data[c("n","Correlation",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim17_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Correlations
  corr <- sort(unique(sim17_data$Correlation))
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
  sim17_data['Correlation'] <- round(sim17_data['Correlation'],digits=2)

  # Reduce to just for some n
  tmp <- sim17_data[(sim17_data$n==reduced_n[1]),]

  # sort by Correlation
  tmp <- tmp[order(tmp$Correlation),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim17_data[(sim17_data$n==n),]
    
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
  sim17_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(-1,1) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Correlation`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Correlation`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Correlation`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 10: Varying Correlation, Square Signal (High SNR)', subtitle = 'Correlation vs Observed Coverage', x = 'Correlation Between Noise Fields', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 17: Varying corrleation, square signal (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim17_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim17/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim17_data) <- c('cfgId','n', 'Correlation','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim17_data <- sim17_data[c("n","Correlation",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim17_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Correlations
  corr <- sort(unique(sim17_data$Correlation))
  reduced_corr <- c(-1,-0.5,0, 0.5, 1)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the correlations due to casting issues
  sim17_data['Correlation'] <- round(sim17_data['Correlation'],digits=2)

  # Reduce to just for some n
  tmp <- sim17_data[(sim17_data$Correlation==reduced_corr[1]),]

  # sort by Correlation
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (corr in reduced_corr[2:length(reduced_corr)]){
    
    # Reduce to just for some n
    tmp2 <- sim17_data[(sim17_data$Correlation==corr),]
    
    # sort by Correlation
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Correlation
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Correlation <- as.factor(tmp$Correlation)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim17_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Correlation`, color=`Correlation`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'green', '1' = 'slategray'), name = 'Correlation') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 10: Varying Correlation, Square Signal (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 18: Varying corrleation, square signal (Low SNR)
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
    xlim(-1,1) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Correlation`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Correlation`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Correlation`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 9: Varying Correlation, Square Signal (Low SNR)', subtitle = 'Correlation vs Observed Coverage', x = 'Correlation Between Noise Fields', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 18: Varying corrleation, square signal (Low SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
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
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the correlations due to casting issues
  sim18_data['Correlation'] <- round(sim18_data['Correlation'],digits=2)

  # Reduce to just for some n
  tmp <- sim18_data[(sim18_data$Correlation==reduced_corr[1]),]

  # sort by Correlation
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (corr in reduced_corr[2:length(reduced_corr)]){
    
    # Reduce to just for some n
    tmp2 <- sim18_data[(sim18_data$Correlation==corr),]
    
    # sort by Correlation
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Correlation
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Correlation <- as.factor(tmp$Correlation)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim18_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Correlation`, color=`Correlation`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'green', '1' = 'slategray'), name = 'Correlation') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 9: Varying Correlation, Square Signal (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 17 and 18 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim9and10.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim18_n_vs_cov, sim18_corr_vs_cov, sim17_n_vs_cov, sim17_corr_vs_cov, ncol=2, nrow=2)
  dev.off()

}

if (simNumbers=='19-20'){
  # ====================================================================================================
  # Simulation 19: Varying Signal Gradient, Ramp signal (High SNR)
  # ====================================================================================================
  # Plot: Gradient vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim19_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim19/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim19_data) <- c('cfgId','n', 'Gradient','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim19_data <- sim19_data[c("n","Gradient",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim19_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Gradients
  grad <- sort(unique(sim19_data$Gradient))
  reduced_grad <- c(2,4,6)

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
  sim19_data['Gradient'] <- round(sim19_data['Gradient'],digits=3)

  # Reduce to just for some n
  tmp <- sim19_data[(sim19_data$n==reduced_n[1]),]

  # sort by Gradient
  tmp <- tmp[order(tmp$Gradient),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim19_data[(sim19_data$n==n),]
    
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
  sim19_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(2,6) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Gradient`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Gradient`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Gradient`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 12: Varying Signal Gradient, Ramp Signal (High SNR)', subtitle = 'Gradient vs Observed Coverage', x = 'Gradient (per 50 voxels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 19: Varying Signal Gradient, Ramp signal (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim19_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim19/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim19_data) <- c('cfgId','n', 'Gradient','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim19_data <- sim19_data[c("n","Gradient",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim19_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Gradients
  grad <- sort(unique(sim19_data$Gradient))
  reduced_grad <- c(2,4,6)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Gradients due to casting issues
  sim19_data['Gradient'] <- round(sim19_data['Gradient'],digits=3)

  # Reduce to just for some n
  tmp <- sim19_data[(sim19_data$Gradient==reduced_grad[1]),]

  # sort by Gradient
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (grad in reduced_grad[2:length(reduced_grad)]){
    
    # Reduce to just for some n
    tmp2 <- sim19_data[(sim19_data$Gradient==grad),]
    
    # sort by Gradient
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Gradient
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Gradient <- as.factor(tmp$Gradient)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim19_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Gradient`, color=`Gradient`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('2' = 'salmon','4' = 'darkorchid','6' = 'slategray'), name = 'Gradient') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 12: Varying Signal Gradient, Ramp Signal (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
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
    xlim(0.5,1.5) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Gradient`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Gradient`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Gradient`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 11: Varying Signal Gradient, Ramp Signal (Low SNR)', subtitle = 'Gradient vs Observed Coverage', x = 'Gradient (per 50 voxels)', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
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
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Gradients due to casting issues
  sim20_data['Gradient'] <- round(sim20_data['Gradient'],digits=3)

  # Reduce to just for some n
  tmp <- sim20_data[(sim20_data$Gradient==reduced_grad[1]),]

  # sort by Gradient
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (grad in reduced_grad[2:length(reduced_grad)]){
    
    # Reduce to just for some n
    tmp2 <- sim20_data[(sim20_data$Gradient==grad),]
    
    # sort by Gradient
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Gradient
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Gradient <- as.factor(tmp$Gradient)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim20_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Gradient`, color=`Gradient`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('0.5' = 'salmon','1' = 'darkorchid','1.5' = 'slategray'), name = 'Gradient') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 11: Varying Signal Gradient, Ramp Signal (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 19 and 20 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim11and12.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim20_n_vs_cov, sim20_grad_vs_cov, sim19_n_vs_cov, sim19_grad_vs_cov, ncol=2, nrow=2)
  dev.off()

}

if (simNumbers=='21-22'){
  # ====================================================================================================
  # Simulation 21: Varying noise magnitude, square signal (High SNR)
  # ====================================================================================================
  # Plot: Magnitude vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim21_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim21/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim21_data) <- c('cfgId','n', 'Magnitude','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim21_data <- sim21_data[c("n","Magnitude",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim21_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Magnitudes
  mag <- sort(unique(sim21_data$Magnitude))
  reduced_mag <- c(1,2,3)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = mag, y = rep(p,length(mag)))
  uppline <- data.frame( x = mag, y = rep(p+bin_conf,length(mag)))
  lowline <- data.frame( x = mag, y = rep(p-bin_conf,length(mag)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Magnitudes due to casting issues
  sim21_data['Magnitude'] <- round(sim21_data['Magnitude'],digits=2)

  # Reduce to just for some n
  tmp <- sim21_data[(sim21_data$n==reduced_n[1]),]

  # sort by Magnitude
  tmp <- tmp[order(tmp$Magnitude),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim21_data[(sim21_data$n==n),]
    
    # sort by Magnitude
    tmp2 <- tmp2[order(tmp2$Magnitude),]
    
    # sort by Magnitude
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
  sim21_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(1,3) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Magnitude`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Magnitude`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Magnitude`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 14: Varying Noise Magnitude, Square Signal (High SNR)', subtitle = 'Magnitude vs Observed Coverage', x = 'Magnitude of Second Noise Field', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 21: Varying noise magnitude, square signal (High SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim21_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim21/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim21_data) <- c('cfgId','n', 'Magnitude','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim21_data <- sim21_data[c("n","Magnitude",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim21_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Magnitudes
  mag <- sort(unique(sim21_data$Magnitude))
  reduced_mag <- c(1,2,3)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Magnitudes due to casting issues
  sim21_data['Magnitude'] <- round(sim21_data['Magnitude'],digits=2)

  # Reduce to just for some n
  tmp <- sim21_data[(sim21_data$Magnitude==reduced_mag[1]),]

  # sort by Magnitude
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (mag in reduced_mag[2:length(reduced_mag)]){
    
    # Reduce to just for some n
    tmp2 <- sim21_data[(sim21_data$Magnitude==mag),]
    
    # sort by Magnitude
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Magnitude
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Magnitude <- as.factor(tmp$Magnitude)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim21_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Magnitude`, color=`Magnitude`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Magnitude') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 14: Varying Noise Magnitude, Square Signal (High SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')

  # ====================================================================================================
  # Simulation 22: Varying noise magnitude, square signal (Low SNR)
  # ====================================================================================================
  # Plot: Magnitude vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim22_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim22/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim22_data) <- c('cfgId','n', 'Magnitude','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim22_data <- sim22_data[c("n","Magnitude",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim22_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Magnitudes
  mag <- sort(unique(sim22_data$Magnitude))
  reduced_mag <- c(1,2,3)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = mag, y = rep(p,length(mag)))
  uppline <- data.frame( x = mag, y = rep(p+bin_conf,length(mag)))
  lowline <- data.frame( x = mag, y = rep(p-bin_conf,length(mag)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Magnitudes due to casting issues
  sim22_data['Magnitude'] <- round(sim22_data['Magnitude'],digits=2)

  # Reduce to just for some n
  tmp <- sim22_data[(sim22_data$n==reduced_n[1]),]

  # sort by Magnitude
  tmp <- tmp[order(tmp$Magnitude),]

  xmin <- 0
  xmax <- 50

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (n in reduced_n[2:length(reduced_n)]){
    
    # Reduce to just for some n
    tmp2 <- sim22_data[(sim22_data$n==n),]
    
    # sort by Magnitude
    tmp2 <- tmp2[order(tmp2$Magnitude),]
    
    # sort by Magnitude
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
  sim22_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=`0.95`, group=`n`, color=`n`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(1,3) + ylim(0.9,1) + 
    scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
    geom_line(aes(x=`Magnitude`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`Magnitude`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`Magnitude`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 13: Varying Noise Magnitude, Square Signal (Low SNR)', subtitle = 'Magnitude vs Observed Coverage', x = 'Magnitude of Second Noise Field', y = 'Observed Coverage')


  # ====================================================================================================
  # Simulation 22: Varying noise magnitude, square signal (Low SNR)
  # ====================================================================================================
  # Plot: Number of observations vs coverage
  # ----------------------------------------------------------------------------------------------------

  # Read in data
  sim22_data <- read.csv(file = '/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim22/FullResults/estBdry_intrp.csv',sep=',', header=FALSE)

  # Name data
  names(sim22_data) <- c('cfgId','n', 'Magnitude','0.00','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00')

  # Reduce data
  sim22_data <- sim22_data[c("n","Magnitude",toString(p))]

  # Sort the unique n
  n <- sort(unique(sim22_data$n))
  reduced_n <-c(100,300,500)

  # Sort the unique Magnitudes
  mag <- sort(unique(sim22_data$Magnitude))
  reduced_mag <- c(1,2,3)

  # Binomial confidence line
  bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)

  # Lines for band
  midline <- data.frame( x = n, y = rep(p,length(n)))
  uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
  lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))

  # -------------------------------------------------------
  # Reformat data
  # -------------------------------------------------------

  # Need to round the Magnitudes due to casting issues
  sim22_data['Magnitude'] <- round(sim22_data['Magnitude'],digits=2)

  # Reduce to just for some n
  tmp <- sim22_data[(sim22_data$Magnitude==reduced_mag[1]),]

  # sort by Magnitude
  tmp <- tmp[order(tmp$n),]

  xmin <- 0
  xmax <- 500

  ymin <- 1-2*(1-p)
  ymax <- 1
    
  # Loop through and add the remaining n
  for (mag in reduced_mag[2:length(reduced_mag)]){
    
    # Reduce to just for some n
    tmp2 <- sim22_data[(sim22_data$Magnitude==mag),]
    
    # sort by Magnitude
    tmp2 <- tmp2[order(tmp2$n),]
    
    # sort by Magnitude
    tmp <- rbind(tmp,tmp2)
    
  }

  # Save n
  tmp$Magnitude <- as.factor(tmp$Magnitude)

  # Save line parameters 
  tmp$truep <- p
  tmp$lowp <- p-bin_conf
  tmp$uppp <- p+bin_conf

  # -------------------------------------------------------
  # Make plot
  # -------------------------------------------------------

  # Create plot
  sim22_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`0.95`, group=`Magnitude`, color=`Magnitude`)) + 
    geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
    xlim(40,500) + ylim(0.9,1) + 
    scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Magnitude') +
    geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
    geom_line(aes(x=`n`,y=`lowp`),linetype="dashed",color="gray65",size=0.1) +
    geom_line(aes(x=`n`,y=`uppp`),linetype="dashed",color="gray65",size=0.1) +
    labs(title = 'Simulation 13: Varying Noise Magnitude, Square Signal (Low SNR)', subtitle = 'Number of Observations vs Observed Coverage', x = 'Number of Observations', y = 'Observed Coverage')



  # ====================================================================================================
  # Combine Sim 21 and 22 plots
  # ====================================================================================================
  
  png(filename = '/home/tommaullin/Documents/ConfRes/FinalSims/sim13and14.png', width = 1800, height = 1000,
      units = "px", pointsize = 12, bg = "white", res = 120)
  grid.arrange(sim22_n_vs_cov, sim22_corr_vs_cov, sim21_n_vs_cov, sim21_corr_vs_cov, ncol=2, nrow=2)
  dev.off()

}