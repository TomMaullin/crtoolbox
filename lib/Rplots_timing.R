  # Imports
  library("ggplot2")
  require(gridExtra)
  library(dplyr)
  
  # All simulations
  simNoArray = c('1-2','7-8','11-12','15-16','17-18','19-20','21-22','23-24','25-26','27-28','29-30','31-32','33-34','35-36','M1-M2','M3-M4')
  
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



      if (simNumbers=='23-24'){
  
        # ====================================================================================================
        # Simulation 23: Varying smoothness of signal
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim23_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)       
        
        # Halve the times as we want to average over est and true and times include both
        sim23_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim23_data <- cbind(sim23_data, sim23_data_times)
  
        # Name data
        names(sim23_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','times')
     
        # Reduce data
        sim23_data <- sim23_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim23_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim23_data$Smoothness))
        reduced_d <- c(0,2,4,6)
    
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim23_data[(sim23_data$n==reduced_n[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
  
        xmin <- 0
        xmax <- 7.8
  
        ymin <- max(0,min(sim23_data$times)-10)
        ymax <- max(sim23_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim23_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (High SNR)', x = 'Smoothness of Signals (FWHM)', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 23: Varying smoothness of signal (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim23_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim23_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim23/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim23_data <- cbind(sim23_data, sim23_data_times)
  
        # Name data
        names(sim23_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim23_data <- sim23_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim23_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim23_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim23_data[(sim23_data$Smoothness==reduced_d[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim23_data$times)-10)
        ymax <- max(sim23_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim23_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Smoothness`, color=`Smoothness`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 24:  Varying Signal Smoothness (Low SNR)
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim24_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim24_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim24_data <- cbind(sim24_data, sim24_data_times)
  
        # Name data
        names(sim24_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim24_data <- sim24_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim24_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim24_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim24_data[(sim24_data$n==reduced_n[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
  
        xmin <- 0
        xmax <- 7.8
  
        ymin <- max(0,min(sim24_data$times)-10)
        ymax <- max(sim24_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim24_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (Low SNR)', x = 'Smoothness of Signals (FWHM)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 24: Varying Signal Smoothness(Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim24_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim24_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim24/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim24_data <- cbind(sim24_data, sim24_data_times)
  
        # Name data
        names(sim24_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim24_data <- sim24_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim24_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim24_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim24_data[(sim24_data$Smoothness==reduced_d[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim24_data$times)-10)
        ymax <- max(sim24_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim24_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Smoothness`, color=`Smoothness`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          labs(title = 'Simulation 10: Varying Signal Smoothness, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Combine Sim 23 and 24 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim10.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
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
        
        # Halve the times as we want to average over est and true and times include both
        sim25_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim25/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim25_data <- cbind(sim25_data, sim25_data_times)
  
        # Name data
        names(sim25_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','times')
     
        # Reduce data
        sim25_data <- sim25_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim25_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim25_data$Smoothness))
        reduced_d <- c(0,2,4,6)
    
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim25_data[(sim25_data$n==reduced_n[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
  
        xmin <- 0
        xmax <- 7.8
  
        ymin <- max(0,min(sim25_data$times)-10)
        ymax <- max(sim25_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim25_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (High SNR)', x = 'Smoothness of Noise (FWHM)', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 25: Varying smoothness of noise (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim25_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim25/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim25_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim25/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim25_data <- cbind(sim25_data, sim25_data_times)
  
        # Name data
        names(sim25_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim25_data <- sim25_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim25_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim25_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim25_data[(sim25_data$Smoothness==reduced_d[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim25_data$times)-10)
        ymax <- max(sim25_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim25_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Smoothness`, color=`Smoothness`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 26:  Varying Noise Smoothness (Low SNR)
        # ====================================================================================================
        # Plot: Smoothness vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim26_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim26_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim26_data <- cbind(sim26_data, sim26_data_times)
  
        # Name data
        names(sim26_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim26_data <- sim26_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim26_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim26_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim26_data[(sim26_data$n==reduced_n[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$Smoothness),]
  
        xmin <- 0
        xmax <- 7.8
  
        ymin <- max(0,min(sim26_data$times)-10)
        ymax <- max(sim26_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim26_d_vs_cov <- ggplot(tmp, aes(x=`Smoothness`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,7.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (Low SNR)', x = 'Smoothness of Noise (FWHM)', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Simulation 26: Varying Noise Smoothness(Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim26_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim26_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim26/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim26_data <- cbind(sim26_data, sim26_data_times)
  
        # Name data
        names(sim26_data) <- c('cfgId','n', 'Smoothness','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim26_data <- sim26_data[c("n","Smoothness","times")]
  
        # Sort the unique n
        n <- sort(unique(sim26_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique smoothnesses
        d <- sort(unique(sim26_data$Smoothness))
        reduced_d <- c(0,2,4,6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim26_data[(sim26_data$Smoothness==reduced_d[1]),]
  
        # sort by smoothness
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim26_data$times)-10)
        ymax <- max(sim26_data$times)*1.2
          
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
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim26_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Smoothness`, color=`Smoothness`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '2' = 'orange', '4' = 'darkorchid', '6' = 'dodgerblue3'), name = 'Smoothness') +
          labs(title = 'Simulation 11: Varying Noise Smoothness, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Combine Sim 25 and 26 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim11.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
        grid.arrange(sim26_n_vs_cov, sim26_d_vs_cov, sim25_n_vs_cov, sim25_d_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }



      if (simNumbers=='27-28'){
        
        # ====================================================================================================
        # Simulation 27: Varying Threshold
        # ====================================================================================================
        # Plot: Threshold vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim27_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim27/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)       
        
        # Halve the times as we want to average over est and true and times include both
        sim27_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim27/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim27_data <- cbind(sim27_data, sim27_data_times)
  
        # Name data
        names(sim27_data) <- c('cfgId','n', 'Threshold','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','times')
     
        # Reduce data
        sim27_data <- sim27_data[c("n","Threshold","times")]
  
        # Sort the unique n
        n <- sort(unique(sim27_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique thresholds
        d <- sort(unique(sim27_data$Threshold))
        reduced_d <- c(0,0.8,1.6,2.4)
    
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim27_data[(sim27_data$n==reduced_n[1]),]
  
        # sort by threshold
        tmp <- tmp[order(tmp$Threshold),]
  
        xmin <- 0
        xmax <- 3.2
  
        ymin <- max(0,min(sim27_data$times)-10)
        ymax <- max(sim27_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim27_data[(sim27_data$n==n),]
          
          # sort by threshold
          tmp2 <- tmp2[order(tmp2$Threshold),]
          
          # sort by threshold
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$n <- as.factor(tmp$n)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim27_d_vs_cov <- ggplot(tmp, aes(x=`Threshold`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,3.2) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 12: Varying Threshold, Square Signal (High SNR)', x = 'Threshold (c)', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 27: Varying Threshold (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim27_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim27/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim27_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim27/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim27_data <- cbind(sim27_data, sim27_data_times)
  
        # Name data
        names(sim27_data) <- c('cfgId','n', 'Threshold','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim27_data <- sim27_data[c("n","Threshold","times")]
  
        # Sort the unique n
        n <- sort(unique(sim27_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique thresholds
        d <- sort(unique(sim27_data$Threshold))
        reduced_d <- c(0,0.8,1.6,2.4)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim27_data[(abs(sim27_data$Threshold-reduced_d[1]) < 1e-6),]
  
        # sort by threshold
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim27_data$times)-10)
        ymax <- max(sim27_data$times)*1.2
          
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n # sim27_data$Threshold-0.6 < 1e-6
          tmp2 <- sim27_data[(abs(sim27_data$Threshold-d) < 1e-6),]
          
          # sort by threshold
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by threshold
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$Threshold <- as.factor(tmp$Threshold)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim27_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Threshold`, color=`Threshold`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '0.8' = 'orange', '1.6' = 'darkorchid', '2.4' = 'dodgerblue3'), name = 'Threshold') +
          labs(title = 'Simulation 12: Varying Threshold, Square Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 28:  Varying Threshold (Low SNR)
        # ====================================================================================================
        # Plot: Threshold vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim28_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim28/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim28_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim28/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim28_data <- cbind(sim28_data, sim28_data_times)
  
        # Name data
        names(sim28_data) <- c('cfgId','n', 'Threshold','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim28_data <- sim28_data[c("n","Threshold","times")]
  
        # Sort the unique n
        n <- sort(unique(sim28_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique thresholds
        d <- sort(unique(sim28_data$Threshold))
        reduced_d <- c(0,0.2,0.4,0.6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim28_data[(sim28_data$n==reduced_n[1]),]
  
        # sort by threshold
        tmp <- tmp[order(tmp$Threshold),]
  
        xmin <- 0
        xmax <- 0.8
  
        ymin <- max(0,min(sim28_data$times)-10)
        ymax <- max(sim28_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim28_data[(sim28_data$n==n),]
          
          # sort by threshold
          tmp2 <- tmp2[order(tmp2$Threshold),]
          
          # sort by threshold
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$n <- as.factor(tmp$n)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim28_d_vs_cov <- ggplot(tmp, aes(x=`Threshold`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,0.8) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 12: Varying Threshold, Square Signal (Low SNR)', x = 'Threshold (c)', y = 'Time (Seconds)')
  
        
        # ====================================================================================================
        # Simulation 28: Varying Threshold(Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
        
        # Read in data
        sim28_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim28/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim28_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim28/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim28_data <- cbind(sim28_data, sim28_data_times)
  
        # Name data
        names(sim28_data) <- c('cfgId','n', 'Threshold','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim28_data <- sim28_data[c("n","Threshold","times")]
  
        # Sort the unique n
        n <- sort(unique(sim28_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique thresholds
        d <- sort(unique(sim28_data$Threshold))
        reduced_d <- c(0,0.2,0.4,0.6)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim28_data[(abs(sim28_data$Threshold-reduced_d[1]) < 1e-6),]

        # sort by threshold
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim28_data$times)-10)
        ymax <- max(sim28_data$times)*1.2
          
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim28_data[(abs(sim28_data$Threshold-d) < 1e-6),]
          
          # sort by threshold
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by threshold
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$Threshold <- as.factor(tmp$Threshold)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim28_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Threshold`, color=`Threshold`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'indianred1', '0.2' = 'orange', '0.4' = 'darkorchid', '0.6' = 'dodgerblue3'), name = 'Threshold') +
          labs(title = 'Simulation 12: Varying Threshold, Square Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
  
        # ====================================================================================================
        # Combine Sim 27 and 28 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim12.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
        grid.arrange(sim28_n_vs_cov, sim28_d_vs_cov, sim27_n_vs_cov, sim27_d_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }

      if (simNumbers=='33-34'){

        # ====================================================================================================
        # Simulation 33: Moving circles closer together (High SNR), Niave Intersections
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim33_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim33/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)       
        
        # Halve the times as we want to average over est and true and times include both
        sim33_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim33/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim33_data <- cbind(sim33_data, sim33_data_times)
  
        # Name data
        names(sim33_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','times')
     
        # Reduce data
        sim33_data <- sim33_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim33_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim33_data$Distance))
        reduced_d <- c(0,20,40)
    
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim33_data[(sim33_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim33_data$times)-10)
        ymax <- max(sim33_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim33_data[(sim33_data$n==n),]
          
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
        sim33_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 13: Varying Distance, Circle Signal, Niave Intersections (High SNR)', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 33: Moving circles closer together (High SNR), Niave intersections
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim33_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim33/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim33_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim33/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim33_data <- cbind(sim33_data, sim33_data_times)
  
        # Name data
        names(sim33_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim33_data <- sim33_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim33_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim33_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim33_data[(sim33_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim33_data$times)-10)
        ymax <- max(sim33_data$times)*1.2
          
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim33_data[(sim33_data$Distance==d),]
          
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
        sim33_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 13: Varying Distance, Circle Signal, Niave Intersections (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 34: Moving circles closer together (Low SNR), Niave Intersections
        # ====================================================================================================
        # Plot: Distance vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim34_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim34/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim34_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim34/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim34_data <- cbind(sim34_data, sim34_data_times)
  
        # Name data
        names(sim34_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim34_data <- sim34_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim34_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim34_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim34_data[(sim34_data$n==reduced_n[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$Distance),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim34_data$times)-10)
        ymax <- max(sim34_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim34_data[(sim34_data$n==n),]
          
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
        sim34_d_vs_cov <- ggplot(tmp, aes(x=`Distance`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0,50) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 13: Varying Distance, Circle Signal, Niave Intersections (Low SNR)', x = 'Distance Between Circle Centers (Number of Pixels)', y = 'Time (Seconds)')
  
        
        # ====================================================================================================
        # Simulation 34: Moving circles closer together (Low SNR), Niave Intersections
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim34_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim34/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim34_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim34/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim34_data <- cbind(sim34_data, sim34_data_times)
  
        # Name data
        names(sim34_data) <- c('cfgId','n', 'Distance','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim34_data <- sim34_data[c("n","Distance","times")]
  
        # Sort the unique n
        n <- sort(unique(sim34_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique distances
        d <- sort(unique(sim34_data$Distance))
        reduced_d <- c(0,20,40)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim34_data[(sim34_data$Distance==reduced_d[1]),]
  
        # sort by distance
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim34_data$times)-10)
        ymax <- max(sim34_data$times)*1.2
          
        # Loop through and add the remaining n
        for (d in reduced_d[2:length(reduced_d)]){
          
          # Reduce to just for some n
          tmp2 <- sim34_data[(sim34_data$Distance==d),]
          
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
        sim34_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`Distance`, color=`Distance`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0' = 'salmon','20' = 'darkorchid','40' = 'slategray'), name = 'Distance') +
          labs(title = 'Simulation 13: Varying Distance, Circle Signal, Niave Intersections (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Combine Sim 33 and 34 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim13.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
        grid.arrange(sim34_n_vs_cov, sim34_d_vs_cov, sim33_n_vs_cov, sim33_d_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }


    if (simNumbers=='35-36'){

        # ====================================================================================================
        # Simulation 35: Varying Noise Mixture (High SNR)
        # ====================================================================================================
        # Plot: Noise Mixture vs coverage
        # ----------------------------------------------------------------------------------------------------
        # Read in data
        sim35_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim35/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)       
        
        # Halve the times as we want to average over est and true and times include both
        sim35_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim35/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim35_data <- cbind(sim35_data, sim35_data_times)
  
        # Name data
        names(sim35_data) <- c('cfgId','n', 'NoiseMix','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0','times')
     
        # Reduce data
        sim35_data <- sim35_data[c("n","NoiseMix","times")]
  
        # Sort the unique n
        n <- sort(unique(sim35_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique noise mixtures
        nm <- sort(unique(sim35_data$NoiseMix))
        reduced_nm <- c(0.8,2,3.2,4)
    
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim35_data[(sim35_data$n==reduced_n[1]),]
  
        # sort by noise mixture
        tmp <- tmp[order(tmp$NoiseMix),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim35_data$times)-10)
        ymax <- max(sim35_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim35_data[(sim35_data$n==n),]
          
          # sort by noise mixture
          tmp2 <- tmp2[order(tmp2$NoiseMix),]
          
          # sort by noise mixture
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$n <- as.factor(tmp$n)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim35_nm_vs_cov <- ggplot(tmp, aes(x=`NoiseMix`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0.8,4) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 14: Varying Noise Mixture, Circle Signal (High SNR)', x = 'Second Noise Component', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 35: Varying Noise Mixture (High SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim35_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim35/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim35_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim35/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim35_data <- cbind(sim35_data, sim35_data_times)
  
        # Name data
        names(sim35_data) <- c('cfgId','n', 'NoiseMix','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim35_data <- sim35_data[c("n","NoiseMix","times")]
  
        # Sort the unique n
        n <- sort(unique(sim35_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique noise mixtures
        nm <- sort(unique(sim35_data$NoiseMix))
        reduced_nm <- c(0.8,2,3.2,4)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim35_data[(sim35_data$NoiseMix==reduced_nm[1]),]
  
        # sort by noise mixture
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim35_data$times)-10)
        ymax <- max(sim35_data$times)*1.2
          
        # Loop through and add the remaining n
        for (nm in reduced_nm[2:length(reduced_nm)]){
          
          # Reduce to just for some n
          tmp2 <- sim35_data[(sim35_data$NoiseMix==nm),]
          
          # sort by noise mixture
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by noise mixture
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$NoiseMix <- as.factor(tmp$NoiseMix)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim35_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`NoiseMix`, color=`NoiseMix`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0.8' = 'salmon','2' = 'darkorchid','3.2' = 'slategray', '4' = 'indianred1'), name = 'NoiseMix') +
          labs(title = 'Simulation 14: Varying Noise Mixture, Circle Signal (High SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Simulation 36: Varying Noise Mixture (Low SNR)
        # ====================================================================================================
        # Plot: Noise Mixture vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim36_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim36/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim36_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim36/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim36_data <- cbind(sim36_data, sim36_data_times)
  
        # Name data
        names(sim36_data) <- c('cfgId','n', 'NoiseMix','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim36_data <- sim36_data[c("n","NoiseMix","times")]
  
        # Sort the unique n
        n <- sort(unique(sim36_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique noise mixtures
        nm <- sort(unique(sim36_data$NoiseMix))
        reduced_nm <- c(0.8,2,3.2,4)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim36_data[(sim36_data$n==reduced_n[1]),]
  
        # sort by noise mixture
        tmp <- tmp[order(tmp$NoiseMix),]
  
        xmin <- 0
        xmax <- 50
  
        ymin <- max(0,min(sim36_data$times)-10)
        ymax <- max(sim36_data$times)*1.2
          
        # Loop through and add the remaining n
        for (n in reduced_n[2:length(reduced_n)]){
          
          # Reduce to just for some n
          tmp2 <- sim36_data[(sim36_data$n==n),]
          
          # sort by noise mixture
          tmp2 <- tmp2[order(tmp2$NoiseMix),]
          
          # sort by noise mixture
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$n <- as.factor(tmp$n)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim36_nm_vs_cov <- ggplot(tmp, aes(x=`NoiseMix`,y=`times`, group=`n`, color=`n`)) + geom_line() + 
          xlim(0.8,4) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('100' = 'salmon','300' = 'darkorchid','500' = 'slategray'), name = 'n') +
          labs(title = 'Simulation 14: Varying Noise Mixture, Circle Signal (Low SNR)', x = 'Second Noise Component', y = 'Time (Seconds)')
  
        
        # ====================================================================================================
        # Simulation 36: Varying Noise Mixture (Low SNR)
        # ====================================================================================================
        # Plot: Number of observations vs coverage
        # ----------------------------------------------------------------------------------------------------
  
        # Read in data
        sim36_data <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim36/FullResults/', tolower(bdryType), 'Bdry_intrp.csv', sep='') ,sep=',', header=FALSE)
        
        # Halve the times as we want to average over est and true and times include both
        sim36_data_times <- read.csv(file = paste('/home/tommaullin/Documents/ConfRes/FinalSims/2smp/Sim36/FullResults/', 'times.csv', sep='') ,sep=',', header=FALSE)/2

        sim36_data <- cbind(sim36_data, sim36_data_times)
  
        # Name data
        names(sim36_data) <- c('cfgId','n', 'NoiseMix','0.0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.9','0.95','1.0', 'times')
  
        # Reduce data
        sim36_data <- sim36_data[c("n","NoiseMix","times")]
  
        # Sort the unique n
        n <- sort(unique(sim36_data$n))
        reduced_n <-c(100,300,500)
  
        # Sort the unique noise mixtures
        nm <- sort(unique(sim36_data$NoiseMix))
        reduced_nm <- c(0.8,2,3.2,4)
  
        # -------------------------------------------------------
        # Reformat data
        # -------------------------------------------------------
  
        # Reduce to just for some n
        tmp <- sim36_data[(sim36_data$NoiseMix==reduced_nm[1]),]
  
        # sort by noise mixture
        tmp <- tmp[order(tmp$n),]
  
        xmin <- 0
        xmax <- 500
  
        ymin <- max(0,min(sim36_data$times)-10)
        ymax <- max(sim36_data$times)*1.2
          
        # Loop through and add the remaining n
        for (nm in reduced_nm[2:length(reduced_nm)]){
          
          # Reduce to just for some n
          tmp2 <- sim36_data[(sim36_data$NoiseMix==nm),]
          
          # sort by noise mixture
          tmp2 <- tmp2[order(tmp2$n),]
          
          # sort by noise mixture
          tmp <- rbind(tmp,tmp2)
          
        }
  
        # Save n
        tmp$NoiseMix <- as.factor(tmp$NoiseMix)
  
        # -------------------------------------------------------
        # Make plot
        # -------------------------------------------------------
  
        # Create plot
        sim36_n_vs_cov <- ggplot(tmp, aes(x=`n`,y=`times`, group=`NoiseMix`, color=`NoiseMix`)) + geom_line() + 
          xlim(40,500) + ylim(ymin,ymax) + 
          scale_color_manual(values = c('0.8' = 'salmon','2' = 'darkorchid','3.2' = 'slategray', '4' = 'indianred1'), name = 'NoiseMix') +
          labs(title = 'Simulation 14: Varying Noise Mixture, Circle Signal (Low SNR)', x = 'Number of Observations', y = 'Time (Seconds)')
  
        # ====================================================================================================
        # Combine Sim 35 and 36 plots
        # ====================================================================================================
        
        png(filename = paste('/home/tommaullin/Documents/ConfRes/FinalSims/Times/sim14.png', sep = ''), width = 1800, height = 1000,
            units = "px", pointsize = 14, bg = "white", res = 120)
        grid.arrange(sim36_n_vs_cov, sim36_nm_vs_cov, sim35_n_vs_cov, sim35_nm_vs_cov, ncol=2, nrow=2)
        dev.off()
  
      }
  
    }
  
  }
