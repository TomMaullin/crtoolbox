  # Imports
  library("ggplot2")
  require(gridExtra)
  library(dplyr)
  
  # All simulations
  simNoArray = c('1-2','7-8','11-12','15-16','17-18','19-20','21-22','M1-M2','M3-M4')
  
  # All p values we're interested in
  pArray <- c(0.80,0.90,0.95)
  
  # Relevant parameters
  nReals <- 2500
  
  # Boundary type
  bdryType <- 'Est'
  
  if (bdryType=='Est') {
    bdryStr = 'Estimated Boundary'
  } else if (bdryType=='True') {
    bdryStr = 'True Boundary'
  }
  
  for(p in pArray){
  
    for(simNumbers in simNoArray){
  
      if (simNumbers=='1-2'){
  
        # ====================================================================================================
        # Simulation 1: Moving circles closer together (High SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim1/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim1_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim1/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim1_data <- cbind(sim1_data, sim1_data_times)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim1_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim1_data[(sim1_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim1_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (High SNR)', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 1: Moving circles closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim1/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim1_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim1/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim1_data <- cbind(sim1_data, sim1_data_times)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim1_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim1_data[(sim1_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim1_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 2: Moving circles closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim2/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim2_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim2/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim2_data <- cbind(sim2_data, sim2_data_times)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim2_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim2_data[(sim2_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim2_data$times)-10)
        ymax <- max(sim2_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim2_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (Low SNR)', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 2: Moving circles closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim2/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim2_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim2/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim2_data <- cbind(sim2_data, sim2_data_times)
  
        # Name data
        names(sim2_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim2_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim2_data[(sim2_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim2_data$times)-10)
        ymax <- max(sim2_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim2_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 1: Varying Distance, Circle Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Combine Sim 1 and 2 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim1.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim7_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim7/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim7_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim7/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim7_data <- cbind(sim7_data, sim7_data_times)
  
        # Name data
        names(sim7_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim7_data <- sim7_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim7_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim7_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim7_data[(sim7_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim7_data$times)-10)
        ymax <- max(sim7_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim7_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (High SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 7: Moving squares closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim7_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim7/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim7_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim7/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim7_data <- cbind(sim7_data, sim7_data_times)
  
        # Name data
        names(sim7_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim7_data <- sim7_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim7_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim7_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim7_data[(sim7_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim7_data$times)-10)
        ymax <- max(sim7_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim7_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 8: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim8_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim8/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim8_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim8/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim8_data <- cbind(sim8_data, sim8_data_times)
  
        # Name data
        names(sim8_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim8_data <- sim8_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim8_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim8_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim8_data[(sim8_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim8_data$times)-10)
        ymax <- max(sim8_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim8_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (Low SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 8: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim8_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim8/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim8_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim8/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim8_data <- cbind(sim8_data, sim8_data_times)
  
        # Name data
        names(sim8_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim8_data <- sim8_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim8_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim8_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim8_data[(sim8_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim8_data$times)-10)
        ymax <- max(sim8_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim8_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 4: Varying Distance, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 7 and 8 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim4.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim11_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim11/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim11_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim11/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim11_data <- cbind(sim11_data, sim11_data_times)
  
        # Name data
        names(sim11_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim11_data <- sim11_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim11_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim11_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim11_data[(sim11_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim11_data$times)-10)
        ymax <- max(sim11_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim11_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (High SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 11: Moving squares closer together (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim11_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim11/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim11_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim11/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim11_data <- cbind(sim11_data, sim11_data_times)
  
        # Name data
        names(sim11_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim11_data <- sim11_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim11_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim11_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim11_data[(sim11_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim11_data$times)-10)
        ymax <- max(sim11_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim11_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 12: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim12_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim12/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim12_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim12/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim12_data <- cbind(sim12_data, sim12_data_times)
  
        # Name data
        names(sim12_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim12_data <- sim12_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim12_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim12_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim12_data[(sim12_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim12_data$times)-10)
        ymax <- max(sim12_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim12_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (Low SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 12: Moving squares closer together (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim12_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim12/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim12_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim12/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim12_data <- cbind(sim12_data, sim12_data_times)
  
        # Name data
        names(sim12_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim12_data <- sim12_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim12_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim12_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim12_data[(sim12_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim12_data$times)-10)
        ymax <- max(sim12_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim12_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 5: Varying Distance, Square Signal, Heterogeneous Noise (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 11 and 12 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim5.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim15_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim15/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim15_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim15/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim15_data <- cbind(sim15_data, sim15_data_times)
  
        # Name data
        names(sim15_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim15_data <- sim15_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim15_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim15_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim15_data[(sim15_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim15_data$times)-10)
        ymax <- max(sim15_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim15_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,46) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (High SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 15: Moving squares, one smaller (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim15_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim15/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim15_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim15/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim15_data <- cbind(sim15_data, sim15_data_times)
  
        # Name data
        names(sim15_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim15_data <- sim15_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim15_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim15_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim15_data[(sim15_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim15_data$times)-10)
        ymax <- max(sim15_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim15_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 16: Moving squares, one smaller (Low SNR)
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim16_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim16/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim16_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim16/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim16_data <- cbind(sim16_data, sim16_data_times)
  
        # Name data
        names(sim16_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim16_data <- sim16_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim16_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim16_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim16_data[(sim16_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim16_data$times)-10)
        ymax <- max(sim16_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim16_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,46) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (Low SNR)', x = 'Distance Between Square Centers (Number of Pixels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 16: Moving squares, one smaller (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim16_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim16/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim16_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim16/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim16_data <- cbind(sim16_data, sim16_data_times)
  
        # Name data
        names(sim16_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim16_data <- sim16_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim16_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim16_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim16_data[(sim16_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim16_data$times)-10)
        ymax <- max(sim16_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim16_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 6: Varying Distance, Square Signal, One Smaller (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 15 and 16 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim6.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim17_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim17/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim17_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim17/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim17_data <- cbind(sim17_data, sim17_data_times)
  
        # Name data
        names(sim17_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim17_data <- sim17_data[c("n","Correlation","times")]
  
        # Sort the unique n
        n <- sort(unique(sim17_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Correlations
        corr <- sort(unique(sim17_data$Correlation))
        reduced_corr <- c(-1,-0.5,0, 0.5, 1)
  
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
  
        ymin <- max(0,min(sim17_data$times)-10)
        ymax <- max(sim17_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim17_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(-1,1) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (High SNR)', x = 'Correlation Between Noise Fields', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 17: Varying corrleation, square signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim17_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim17/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim17_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim17/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim17_data <- cbind(sim17_data, sim17_data_times)
  
        # Name data
        names(sim17_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim17_data <- sim17_data[c("n","Correlation","times")]
  
        # Sort the unique n
        n <- sort(unique(sim17_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Correlations
        corr <- sort(unique(sim17_data$Correlation))
        reduced_corr <- c(-1,-0.5,0, 0.5, 1)
  
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
  
        ymin <- max(0,min(sim17_data$times)-10)
        ymax <- max(sim17_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim17_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Correlation`, color=`Correlation`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'dodgerblue3', '1' = 'slategray'), name = 'Correlation') +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 18: Varying corrleation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Correlation vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim18_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim18/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim18_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim18/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim18_data <- cbind(sim18_data, sim18_data_times)
  
        # Name data
        names(sim18_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim18_data <- sim18_data[c("n","Correlation","times")]
  
        # Sort the unique n
        n <- sort(unique(sim18_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Correlations
        corr <- sort(unique(sim18_data$Correlation))
        reduced_corr <- c(-1,-0.5,0, 0.5, 1)
  
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
  
        ymin <- max(0,min(sim18_data$times)-10)
        ymax <- max(sim18_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim18_corr_vs_cov <- ggplot(tmp, aes(x=`Correlation`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(-1,1) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (Low SNR)', x = 'Correlation Between Noise Fields', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 18: Varying corrleation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim18_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim18/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim18_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim18/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim18_data <- cbind(sim18_data, sim18_data_times)
  
        # Name data
        names(sim18_data) <- c('cfgId','n', 'Correlation','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim18_data <- sim18_data[c("n","Correlation","times")]
  
        # Sort the unique n
        n <- sort(unique(sim18_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Correlations
        corr <- sort(unique(sim18_data$Correlation))
        reduced_corr <- c(-1,-0.5,0, 0.5, 1)
  
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
  
        ymin <- max(0,min(sim18_data$times)-10)
        ymax <- max(sim18_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim18_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Correlation`, color=`Correlation`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('-1' = 'indianred1', '-0.5' = 'orange', '0' = 'darkorchid', '0.5' = 'dodgerblue3', '1' = 'slategray'), name = 'Correlation') +
          labs(title = 'Simulation 2: Varying Correlation, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 17 and 18 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim2.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim19_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim19/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim19_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim19/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim19_data <- cbind(sim19_data, sim19_data_times)
  
        # Name data
        names(sim19_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim19_data <- sim19_data[c("n","Gradient","times")]
  
        # Sort the unique n
        n <- sort(unique(sim19_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Gradients
        grad <- sort(unique(sim19_data$Gradient))
        reduced_grad <- c(2,4,6)
  
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
  
        ymin <- max(0,min(sim19_data$times)-10)
        ymax <- max(sim19_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim19_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(2,14) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (High SNR)', x = 'Gradient (per 50 voxels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 19: Varying Signal Gradient, Ramp signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim19_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim19/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim19_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim19/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim19_data <- cbind(sim19_data, sim19_data_times)
  
        # Name data
        names(sim19_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim19_data <- sim19_data[c("n","Gradient","times")]
  
        # Sort the unique n
        n <- sort(unique(sim19_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Gradients
        grad <- sort(unique(sim19_data$Gradient))
        reduced_grad <- c(2,4,6)
  
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
  
        ymin <- max(0,min(sim19_data$times)-10)
        ymax <- max(sim19_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim19_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Gradient`, color=`Gradient`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'salmon','4' = 'darkorchid','6' = 'slategray'), name = 'Gradient') +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
        # ====================================================================================================
        # Plot: Gradient vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim20_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim20/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim20_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim20/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim20_data <- cbind(sim20_data, sim20_data_times)
  
        # Name data
        names(sim20_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim20_data <- sim20_data[c("n","Gradient","times")]
  
        # Sort the unique n
        n <- sort(unique(sim20_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Gradients
        grad <- sort(unique(sim20_data$Gradient))
        reduced_grad <- c(2,4,6)/4
  
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
  
        ymin <- max(0,min(sim20_data$times)-10)
        ymax <- max(sim20_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim20_grad_vs_cov <- ggplot(tmp, aes(x=`Gradient`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0.5,3.5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (Low SNR)', x = 'Gradient (per 50 voxels)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 20: Varying Signal Gradient, Ramp signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim20_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim20/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim20_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim20/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim20_data <- cbind(sim20_data, sim20_data_times)
  
        # Name data
        names(sim20_data) <- c('cfgId','n', 'Gradient','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim20_data <- sim20_data[c("n","Gradient","times")]
  
        # Sort the unique n
        n <- sort(unique(sim20_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Gradients
        grad <- sort(unique(sim20_data$Gradient))
        reduced_grad <- c(2,4,6)/4
  
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
  
        ymin <- max(0,min(sim20_data$times)-10)
        ymax <- max(sim20_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim20_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Gradient`, color=`Gradient`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0.5' = 'salmon','1' = 'darkorchid','1.5' = 'slategray'), name = 'Gradient') +
          labs(title = 'Simulation 3: Varying Signal Gradient, Ramp Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 19 and 20 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim3.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        sim21_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim21/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim21_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim21/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim21_data <- cbind(sim21_data, sim21_data_times)
  
        # Name data
        names(sim21_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim21_data <- sim21_data[c("n","Magnitude","times")]
  
        # Sort the unique n
        n <- sort(unique(sim21_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Magnitudes
        mag <- sort(unique(sim21_data$Magnitude))
        reduced_mag <- c(1,2,3)
  
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
  
        ymin <- max(0,min(sim21_data$times)-10)
        ymax <- max(sim21_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim21_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(1,3) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (High SNR)', x = 'Standard Deviation of Second Noise Field', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 21: Varying noise Standard Deviation, square signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim21_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim21/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim21_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim21/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim21_data <- cbind(sim21_data, sim21_data_times)
  
        # Name data
        names(sim21_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim21_data <- sim21_data[c("n","Magnitude","times")]
  
        # Sort the unique n
        n <- sort(unique(sim21_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Magnitudes
        mag <- sort(unique(sim21_data$Magnitude))
        reduced_mag <- c(1,2,3)
  
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
  
        ymin <- max(0,min(sim21_data$times)-10)
        ymax <- max(sim21_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim21_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Magnitude`, color=`Magnitude`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Std.') +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 22: Varying noise Standard Deviation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Magnitude vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim22_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim22/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim22_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim22/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim22_data <- cbind(sim22_data, sim22_data_times)
  
        # Name data
        names(sim22_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim22_data <- sim22_data[c("n","Magnitude","times")]
  
        # Sort the unique n
        n <- sort(unique(sim22_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Magnitudes
        mag <- sort(unique(sim22_data$Magnitude))
        reduced_mag <- c(1,2,3)
  
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
  
        ymin <- max(0,min(sim22_data$times)-10)
        ymax <- max(sim22_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim22_corr_vs_cov <- ggplot(tmp, aes(x=`Magnitude`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(1,3) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (Low SNR)', x = 'Standard Deviation of Second Noise Field', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 22: Varying noise Standard Deviation, square signal (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim22_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim22/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim22_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/Sim22/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim22_data <- cbind(sim22_data, sim22_data_times)
  
        # Name data
        names(sim22_data) <- c('cfgId','n', 'Magnitude','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim22_data <- sim22_data[c("n","Magnitude","times")]
  
        # Sort the unique n
        n <- sort(unique(sim22_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Magnitudes
        mag <- sort(unique(sim22_data$Magnitude))
        reduced_mag <- c(1,2,3)
  
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
  
        ymin <- max(0,min(sim22_data$times)-10)
        ymax <- max(sim22_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim22_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Magnitude`, color=`Magnitude`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('1' = 'salmon','2' = 'darkorchid','3' = 'slategray'), name = 'Std.') +
          labs(title = 'Simulation 7: Varying Noise Standard Deviation, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 21 and 22 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim7.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        
        # Halve the times as we want to average over est and true and times include both
        sim1_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim1/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim1_data <- cbind(sim1_data, sim1_data_times)
  
        # Name data
        names(sim1_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim1_data$M))
        reduced_M <- c(2,3,4,5)
  
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
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim1_m_vs_cov <- ggplot(tmp, aes(x=`M`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(2,5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (High SNR)', x = 'Number of Study Conditions (M)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 1: Varying Number of Study Conditions (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim1_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim1/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim1_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim1/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim1_data <- cbind(sim1_data, sim1_data_times)
        
        # Name data
        names(sim1_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim1_data <- sim1_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim1_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim1_data$M))
        reduced_M <- c(2,3,4,5)
  
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
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim1_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`M`, color=`M`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'indianred1', '3' = 'orange', '4' = 'darkorchid', '5' = 'slategray'), name = 'M') +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 2: Varying Number of Study Conditions (Low SNR)
        # ====================================================================================================
        # Plot: Number of Study Conditions vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim2_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim2_data <- cbind(sim2_data, sim2_data_times)
        
        # Name data
        names(sim2_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim2_data$M))
        reduced_M <- c(2,3,4,5)
  
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
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim2_m_vs_cov <- ggplot(tmp, aes(x=`M`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(2,5) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (Low SNR)', x = 'Number of Study Conditions (M)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 2: Varying Number of Study Conditions (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim2_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim2_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim2/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim2_data <- cbind(sim2_data, sim2_data_times)
        
        # Name data
        names(sim2_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim2_data <- sim2_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim2_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique Ms
        M <- sort(unique(sim2_data$M))
        reduced_M <- c(2,3,4,5)
  
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
  
        ymin <- max(0,min(sim1_data$times)-10)
        ymax <- max(sim1_data$times)*1.2
          
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
  
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim2_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`M`, color=`M`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('2' = 'indianred1', '3' = 'orange', '4' = 'darkorchid', '5' = 'slategray'), name = 'M') +
          labs(title = 'Simulation 8: Varying Number of Study Conditions (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim M1 and M2 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim8.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        
        # Halve the times as we want to average over est and true and times include both
        sim3_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim3/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim3_data <- cbind(sim3_data, sim3_data_times)
        
        # Name data
        names(sim3_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim3_data <- sim3_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim3_data$n))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Save copy
        tmp <- sim3_data
  
        tmp$Sim <- 'High SNR'
  
        # Read in data
        sim4_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim4/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim4_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Msmp/Sim4/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2
        
        sim4_data <- cbind(sim4_data, sim4_data_times)
        
        # Name data
        names(sim4_data) <- c('cfgId','n', 'M','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim4_data <- sim4_data[c("n","M","times")]
  
        # Sort the unique n
        n <- sort(unique(sim4_data$n))
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Save copy
        tmp2 <- sim4_data
  
        # Save line parameters 
        tmp2$Sim <- 'Low SNR'
  
        # Combine data for each simulation
        finaldata <- rbind(tmp,tmp2)
  
        # Change to factor
        finaldata$Sim <- as.factor(finaldata$Sim)
        
        
        ymin <- max(0,min(sim3_data$times)-10)
        ymax <- max(sim3_data$times)*1.2
        
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim34_n_vs_cov <- ggplot(finaldata, aes(x=`n`,y=`times`, group=`Sim`, color=`Sim`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('Low SNR' = 'salmon','High SNR' = 'darkorchid'), name = 'Simulation') +
          labs(title = 'Simulation 9: Nested Squares', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
  
        # ====================================================================================================
        # Combine Sim 3 and 4 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim9.png', sep = ''), width = 900, height = 500,
            units = "px", pointsize = 14, bg = "white", res = 120)
        grid.arrange(sim34_n_vs_cov, ncol=1, nrow=1)
        dev.off()
  
      }
  
    }
  
  }
