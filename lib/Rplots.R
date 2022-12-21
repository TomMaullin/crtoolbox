  # Imports
  library("ggplot2")
  require(gridExtra)
  library(dplyr)
  
  # All simulations
  simNoArray = c('1-2','7-8','11-12','15-16','17-18','19-20','21-22','23-24','25-26','27-28','29-30','31-32','33-34','M1-M2','M3-M4')
  
  # All p values we're interested in
  pArray <- c(0.80,0.90,0.95)
  
  # Relevant parameters
  nReals <- 2500
  
  # Boundary type
  bdryType <- 'True'
  
  if (bdryType=='Est') {
    bdryStr = 'Estimated Boundary'
  } else if (bdryType=='True') {
    bdryStr = 'True Boundary'
  }
  
  simDirs <- c('2smp/Sim1','2smp/Sim2','2smp/Sim7','2smp/Sim8','2smp/Sim11','2smp/Sim12','2smp/Sim15','2smp/Sim16','2smp/Sim17','2smp/Sim18','2smp/Sim19','2smp/Sim20','2smp/Sim21','2smp/Sim22','2smp/Sim23','2smp/Sim24','2smp/Sim25','2smp/Sim26','2smp/Sim27','2smp/Sim28','2smp/Sim29','2smp/Sim30','2smp/Sim31','2smp/Sim32','Msmp/Sim1','Msmp/Sim2','Msmp/Sim3','Msmp/Sim4')
  
  # Relevant parameters
  nReals <- 2500
  
  # Boundary type
  bdryTypes <- c('True','Est')
  
  # Current max
  y80 <- 0
  y90 <- 0
  y95 <- 0
  
  for(bdry in bdryTypes){
        
      for(simDir in simDirs){
        
        # Read in data
        sim_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/',simDir,'/FullResults/', tolower(bdry), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim_data) <- c('cfgId','n', 'Varying','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim_data <- sim_data[c("n","Varying",'0.8','0.9','0.95')]
        sim_data <- subset(sim_data,n>60)
        
        if ((simDir == '2smp/Sim15') | (simDir == '2smp/Sim16')){
          
          sim_data <- subset(sim_data,Varying<46)
          
        }
        
        if ((simDir == '2smp/Sim20')){
          
          sim_data <- subset(sim_data,Varying>0.9)
          
        }
        
        if ((simDir == '2smp/Sim23') | (simDir == '2smp/Sim24')){
          
          sim_data <- subset(sim_data,Varying>2)
          
        }
        
        # max and min
        y80 <- max(max(sim_data['0.8'])-0.8,0.8-min(sim_data['0.8']),y80)
        y90 <- max(max(sim_data['0.9'])-0.9,0.9-min(sim_data['0.9']),y90)
        y95 <- max(max(sim_data['0.95'])-0.95,0.95-min(sim_data['0.95']),y95)
        
      }
  }
  
  y80 <- ceiling(100*y80)/100
  y90 <- ceiling(100*y90)/100
  y95 <- ceiling(100*y95)/100
  
  
  for(p in pArray){
  
    if (p=='0.8'){
  
      ylimit <- y80
  
    } else if(p=='0.9'){
  
      ylimit <- y90
  
    } else{
  
      ylimit <- y95
  
    }
  
    for(simNumbers in simNoArray){
  
      if (simNumbers=='1-2'){
  
        # ====================================================================================================
        # Simulation 1: Moving circles closer together (High SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","Distance",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim1_data$Distance))
        reduced_d <- c(0,20,40)
  
        # Binomial confidence line
        bin_conf <- qnorm(1-0.5*(1-p))*sqrt(p*(1-p)/nReals)
  
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
        sim1_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (High SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 1: Moving circles closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim1/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim1_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 2: Moving circles closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim2_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (Low SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 2: Moving circles closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim2_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Combine Sim 1 and 2 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim1.png', sep = ''), width = 1800, height = 1000,
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
        sim7_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim7/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim7_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim7_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (High SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 7: Moving squares closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim7_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim7/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim7_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim7_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 8: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim8_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim8/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim8_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim8_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (Low SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 8: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim8_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim8/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim8_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim8_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 7 and 8 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim4.png', sep = ''), width = 1800, height = 1000,
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
        sim11_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim11/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim11_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim11_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (High SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 11: Moving squares closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim11_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim11/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim11_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim11_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 12: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim12_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim12/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim12_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim12_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (Low SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 12: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim12_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim12/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim12_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim12_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 11 and 12 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim5.png', sep = ''), width = 1800, height = 1000,
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
        sim15_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim15/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim15_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim15_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,46) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (High SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 15: Moving squares, one smaller (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim15_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim15/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim15_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim15_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 16: Moving squares, one smaller (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim16_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim16/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim16_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim16_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,46) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Distance`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (Low SNR)', subtitle = paste('Distance vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Distance Between Square Centers (Number of Pixels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 16: Moving squares, one smaller (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim16_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim16/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim16_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim16_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Distance`, color=`Distance`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 15 and 16 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim6.png', sep = ''), width = 1800, height = 1000,
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
        sim17_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim17/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim17_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim17_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(-1,1) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Correlation`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (High SNR)', subtitle = paste('Correlation vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Correlation Between Noise Fields', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 17: Varying corrleation, square signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim17_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim17/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim17_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim17_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Correlation`, color=`Correlation`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'dodgerblue3', '1' = 'slategray'), name = 'Correlation') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 18: Varying corrleation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Correlation vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim18_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim18/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim18_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim18_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(-1,1) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Correlation`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (Low SNR)', subtitle = paste('Correlation vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Correlation Between Noise Fields', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 18: Varying corrleation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim18_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim18/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim18_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim18_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Correlation`, color=`Correlation`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'dodgerblue3', '1' = 'slategray'), name = 'Correlation') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 17 and 18 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim2.png', sep = ''), width = 1800, height = 1000,
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
        sim19_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim19/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim19_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim19_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(2,14) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Gradient`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (High SNR)', subtitle = paste('Gradient vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Gradient (per 50 voxels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 19: Varying Signal Gradient, Ramp signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim19_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim19/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim19_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim19_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Gradient`, color=`Gradient`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'salmon','4' = 'darkorchid','6' = 'slategray'), name = 'Gradient') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
        # ====================================================================================================
        # Plot: Gradient vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim20_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim20/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim20_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim20_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0.5,3.5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Gradient`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (Low SNR)', subtitle = paste('Gradient vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Gradient (per 50 voxels)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim20_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim20/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim20_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim20_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Gradient`, color=`Gradient`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0.5' = 'salmon','1' = 'darkorchid','1.5' = 'slategray'), name = 'Gradient') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 19 and 20 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim3.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim20_n_vs_cov, sim20_grad_vs_cov, sim19_n_vs_cov, sim19_grad_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }
  
      if (simNumbers=='21-22'){
        # ====================================================================================================
        # Simulation 21: Varying noise Standard Deviation, square signal (High SNR)
        # ====================================================================================================
        # Plot: Magnitude vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim21_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim21/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim21_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim21_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(1,3) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Magnitude`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (High SNR)', subtitle = paste('Standard Deviation vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Standard Deviation of Second Noise Field', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 21: Varying noise Standard Deviation, square signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim21_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim21/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim21_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim21_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Magnitude`, color=`Magnitude`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Std.') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 22: Varying noise Standard Deviation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Magnitude vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim22_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim22/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim22_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim22_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(1,3) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Magnitude`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (Low SNR)', subtitle = paste('Standard Deviation vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Standard Deviation of Second Noise Field', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 22: Varying noise Standard Deviation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim22_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim22/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim22_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
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
        sim22_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Magnitude`, color=`Magnitude`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Std.') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 21 and 22 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim7.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim22_n_vs_cov, sim22_corr_vs_cov, sim21_n_vs_cov, sim21_corr_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }
  
      if (simNumbers=='M1-M2'){
        # ====================================================================================================
        # Simulation 1: Varying Number of Study Conditions (High SNR)
        # ====================================================================================================
        # Plot: Number of Study Conditions vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim1/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim1_data$M))
        reduced_M <- c(2,3,4,5)
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = M, y = rep(p,length(M)))
        uppline <- data.frame( x = M, y = rep(p+bin_conf,length(M)))
        lowline <- data.frame( x = M, y = rep(p-bin_conf,length(M)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Need to round the Ms due to casting issues
        sim1_data['M'] <- round(sim1_data['M'],digits=2)
  
        # Reduce to just for some n
        tmp <- sim1_data[(sim1_data$n==reduced_n[1]),]
  
        # sort by M
        tmp <- tmp[order(tmp$M),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- 1-2*(1-p)
        ymax <- 1
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim1_data[(sim1_data$n==n),]
          
          # sort by M
          tmp2 <- tmp2[order(tmp2$M),]
          
          # sort by M
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
        sim1_m_vs_cov <- ggplot(tmp, aes(x=`M`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(2,5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`M`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (High SNR)', subtitle = paste('Number of Study Conditions vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Study Conditions (M)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 1: Varying Number of Study Conditions (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim1/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim1_data$M))
        reduced_M <- c(2,3,4,5)
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = n, y = rep(p,length(n)))
        uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
        lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Need to round the Ms due to casting issues
        sim1_data['M'] <- round(sim1_data['M'],digits=2)
  
        # Reduce to just for some n
        tmp <- sim1_data[(sim1_data$M==reduced_M[1]),]
  
        # sort by M
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- 1-2*(1-p)
        ymax <- 1
          
        # Loop through and add the remaining n
        for (M in reduced_M[2:length(reduced_M)]){
          
          # Reduce to just for some n
          tmp2 <- sim1_data[(sim1_data$M==M),]
          
          # sort by M
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by M
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$M <- as.factor(tmp$M)
  
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim1_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`M`, color=`M`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'indianred1', '3' = 'orange', '4' = 'darkorchid', '5' = 'slategray'), name = 'M') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
        # ====================================================================================================
        # Simulation 2: Varying Number of Study Conditions (Low SNR)
        # ====================================================================================================
        # Plot: Number of Study Conditions vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim2_data$M))
        reduced_M <- c(2,3,4,5)
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = M, y = rep(p,length(M)))
        uppline <- data.frame( x = M, y = rep(p+bin_conf,length(M)))
        lowline <- data.frame( x = M, y = rep(p-bin_conf,length(M)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Need to round the Ms due to casting issues
        sim2_data['M'] <- round(sim2_data['M'],digits=2)
  
        # Reduce to just for some n
        tmp <- sim2_data[(sim2_data$n==reduced_n[1]),]
  
        # sort by M
        tmp <- tmp[order(tmp$M),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- 1-2*(1-p)
        ymax <- 1
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim2_data[(sim2_data$n==n),]
          
          # sort by M
          tmp2 <- tmp2[order(tmp2$M),]
          
          # sort by M
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
        sim2_m_vs_cov <- ggplot(tmp, aes(x=`M`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(2,5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`M`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (Low SNR)', subtitle = paste('Number of Study Conditions vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Study Conditions (M)', y = 'Observed Coverage')
  
  
        # ====================================================================================================
        # Simulation 2: Varying Number of Study Conditions (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim2_data$M))
        reduced_M <- c(2,3,4,5)
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = n, y = rep(p,length(n)))
        uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
        lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Need to round the Ms due to casting issues
        sim2_data['M'] <- round(sim2_data['M'],digits=2)
  
        # Reduce to just for some n
        tmp <- sim2_data[(sim2_data$M==reduced_M[1]),]
  
        # sort by M
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- 1-2*(1-p)
        ymax <- 1
          
        # Loop through and add the remaining n
        for (M in reduced_M[2:length(reduced_M)]){
          
          # Reduce to just for some n
          tmp2 <- sim2_data[(sim2_data$M==M),]
          
          # sort by M
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by M
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$M <- as.factor(tmp$M)
  
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim2_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`M`, color=`M`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'indianred1', '3' = 'orange', '4' = 'darkorchid', '5' = 'slategray'), name = 'M') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim M1 and M2 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim8.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim2_n_vs_cov, sim2_m_vs_cov, sim1_n_vs_cov, sim1_m_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }
  
      if (simNumbers=='M3-M4'){
  
        # ====================================================================================================
        # Simulation 3/4: Nested Squares (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim3_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim3/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim3_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim3_data <- sim3_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim3_data$n))
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = n, y = rep(p,length(n)))
        uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
        lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Save copy
        tmp <- sim3_data
  
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
        tmp$Sim <- 'High SNR'
  
        # Read in data
        sim4_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim4/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
  
        # Name data
        names(sim4_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
  
        # Reduce data
        sim4_data <- sim4_data[c("n","M",toString(p))]
  
        # Sort the unique n
        n <- sort(unique(sim4_data$n))
  
        # Binomial confidence line
        bin_conf <- qnorm(p)*sqrt(p*(1-p)/nReals)
  
        # Lines for band
        midline <- data.frame( x = n, y = rep(p,length(n)))
        uppline <- data.frame( x = n, y = rep(p+bin_conf,length(n)))
        lowline <- data.frame( x = n, y = rep(p-bin_conf,length(n)))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Save copy
        tmp2 <- sim4_data
  
        # Save line parameters 
        tmp2$truep <- p
        tmp2$lowp <- p-bin_conf
        tmp2$uppp <- p+bin_conf
        tmp2$Sim <- 'Low SNR'
  
        # Combine data for each simulation
        finaldata <- rbind(tmp,tmp2)
  
        # Change to factor
        finaldata$Sim <- as.factor(finaldata$Sim)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim34_n_vs_cov <- ggplot(finaldata, aes(x=`n`,y=.data[[toString(p)]], group=`Sim`, color=`Sim`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('Low SNR' = 'salmon','High SNR' = 'darkorchid'), name = 'Simulation') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 9: Nested Squares', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
  
  
  
        # ====================================================================================================
        # Combine Sim 3 and 4 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim9.png', sep = ''), width = 900, height = 500,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim34_n_vs_cov, ncol=1, nrow=1)
        dev.off()
  
      }
      
      if (simNumbers=='23-24'){
        
        # ====================================================================================================
        # Simulation 23: Varying smoothness of signal
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim23_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim23_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim23_data <- sim23_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim23_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesses
        d <- sort(unique(sim23_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
        # Binomial confidence line
        bin_conf <- qnorm(1-0.5*(1-p))*sqrt(p*(1-p)/nReals)
        
        # Lines for band
        midline <- data.frame( x = d, y = rep(p,length(d)))
        uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
        lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))
        
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
        
        # Reduce to just for some n
        tmp <- sim23_data[(sim23_data$n==reduced_n[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
        
        xmin <- 0
        xmax <- 7.8
        
        ymin <- 1-3*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim23_data[(sim23_data$n==n),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$Smoothness),]
          
          # sort by smoothness
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
        sim23_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Smoothness`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (High SNR)', subtitle = paste('Smoothness vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Smoothness of Signals (FWHM)', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Simulation 23: Varying smoothness of signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim23_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim23_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim23_data <- sim23_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim23_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesses
        d <- sort(unique(sim23_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim23_data[(sim23_data$Smoothness==reduced_d[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
        
        xmin <- 0
        xmax <- 500
        
        ymin <- 1-3*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim23_data[(sim23_data$Smoothness==d),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by smoothness
          tmp <- rbind(tmp,tmp2)
          
        }
        
        # Save n
        tmp$Smoothness <- as.factor(tmp$Smoothness)
        
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
        
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
        
        # Create plot
        sim23_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Smoothness`, color=`Smoothness`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values =  c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
        
        # ====================================================================================================
        # Simulation 24:  Varying Signal Smoothness (Low SNR)
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim24_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim24_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim24_data <- sim24_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim24_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesss
        d <- sort(unique(sim24_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim24_data[(sim24_data$n==reduced_n[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
        
        xmin <- 0
        xmax <- 7.8
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim24_data[(sim24_data$n==n),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$Smoothness),]
          
          # sort by smoothness
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
        sim24_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Smoothness`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (Low SNR)', subtitle = paste('Smoothness vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Smoothness of Signals (FWHM)', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Simulation 24: Varying Signal Smoothness(Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim24_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim24_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim24_data <- sim24_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim24_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesss
        d <- sort(unique(sim24_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim24_data[(sim24_data$Smoothness==reduced_d[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
        
        xmin <- 0
        xmax <- 500
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim24_data[(sim24_data$Smoothness==d),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by smoothness
          tmp <- rbind(tmp,tmp2)
          
        }
        
        # Save n
        tmp$Smoothness <- as.factor(tmp$Smoothness)
        
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
        
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
        
        # Create plot
        sim24_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Smoothness`, color=`Smoothness`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Combine Sim 23 and 24 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim10.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim24_n_vs_cov, sim24_d_vs_cov, sim23_n_vs_cov, sim23_d_vs_cov, ncol=2, nrow=2)
        dev.off()
        
      }

      if (simNumbers=='25-26'){
        
        # ====================================================================================================
        # Simulation 25: Varying smoothness of noise
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim25_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim25/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim25_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim25_data <- sim25_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim25_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesses
        d <- sort(unique(sim25_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
        # Binomial confidence line
        bin_conf <- qnorm(1-0.5*(1-p))*sqrt(p*(1-p)/nReals)
        
        # Lines for band
        midline <- data.frame( x = d, y = rep(p,length(d)))
        uppline <- data.frame( x = d, y = rep(p+bin_conf,length(d)))
        lowline <- data.frame( x = d, y = rep(p-bin_conf,length(d)))
        
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
        
        # Reduce to just for some n
        tmp <- sim25_data[(sim25_data$n==reduced_n[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
        
        xmin <- 0
        xmax <- 7.8
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim25_data[(sim25_data$n==n),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$Smoothness),]
          
          # sort by smoothness
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
        sim25_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Smoothness`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (High SNR)', subtitle = paste('Smoothness vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Smoothness of Noise (FWHM)', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Simulation 25: Varying smoothness of noise (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim25_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim25/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim25_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim25_data <- sim25_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim25_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesses
        d <- sort(unique(sim25_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim25_data[(sim25_data$Smoothness==reduced_d[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
        
        xmin <- 0
        xmax <- 500
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim25_data[(sim25_data$Smoothness==d),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by smoothness
          tmp <- rbind(tmp,tmp2)
          
        }
        
        # Save n
        tmp$Smoothness <- as.factor(tmp$Smoothness)
        
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
        
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
        
        # Create plot
        sim25_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Smoothness`, color=`Smoothness`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values =  c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (High SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
        
        # ====================================================================================================
        # Simulation 26:  Varying Noise Smoothness (Low SNR)
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim26_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim26_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim26_data <- sim26_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim26_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesss
        d <- sort(unique(sim26_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim26_data[(sim26_data$n==reduced_n[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
        
        xmin <- 0
        xmax <- 7.8
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim26_data[(sim26_data$n==n),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$Smoothness),]
          
          # sort by smoothness
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
        sim26_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=.data[[toString(p)]], group=`n`, color=`n`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          geom_line(aes(x=`Smoothness`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (Low SNR)', subtitle = paste('Smoothness vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Smoothness of Noise (FWHM)', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Simulation 26: Varying Noise Smoothness(Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim26_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Name data
        names(sim26_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0')
        
        # Reduce data
        sim26_data <- sim26_data[c("n","Smoothness",toString(p))]
        
        # Sort the unique n
        n <- sort(unique(sim26_data$n))
        reduced_n <-c(100,300,500)
        
        # Sort the unique smoothnesss
        d <- sort(unique(sim26_data$Smoothness))
        reduced_d <- c(0,2,4,6)
        
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
        tmp <- sim26_data[(sim26_data$Smoothness==reduced_d[1]),]
        
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
        
        xmin <- 0
        xmax <- 500
        
        ymin <- 1-2*(1-p)
        ymax <- 1
        
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim26_data[(sim26_data$Smoothness==d),]
          
          # sort by smoothness
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by smoothness
          tmp <- rbind(tmp,tmp2)
          
        }
        
        # Save n
        tmp$Smoothness <- as.factor(tmp$Smoothness)
        
        # Save line parameters 
        tmp$truep <- p
        tmp$lowp <- p-bin_conf
        tmp$uppp <- p+bin_conf
        
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
        
        # Create plot
        sim26_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=.data[[toString(p)]], group=`Smoothness`, color=`Smoothness`)) + 
          geom_ribbon(aes(ymin=p-bin_conf, ymax=p+bin_conf), alpha=0.05, fill='turquoise4', colour = NA) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          geom_line(aes(x=`n`,y=`truep`),linetype="dashed",color="black",size=0.1) +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (Low SNR)', subtitle = paste('Number of Observations vs Observed Coverage, p = ', toString(p), ', ', bdryStr, sep=''), x = 'Number of Observations', y = 'Observed Coverage')
        
        
        # ====================================================================================================
        # Combine Sim 25 and 26 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/', toString(100*p), '/', bdryType ,'/sim11.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 12, bg = "white", res = 120)
        grid.arrange(sim26_n_vs_cov, sim26_d_vs_cov, sim25_n_vs_cov, sim25_d_vs_cov, ncol=2, nrow=2)
        dev.off()
        
      }
  
    }
  
  }
