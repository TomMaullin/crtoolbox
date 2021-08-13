# Load libraries
library(ggplot2)
library(fields)

# Handy function, taken from https://stackoverflow.com/questions/68772335/color-palette-for-overlapping-binary-images-in-r/68773493#68773493

mask_to_df <- function(mask, alpha, fillcol)
{ 
  mask <- t(apply(mask * alpha, 2, rev)) 
  data.frame(values = as.vector(mask),
             fillcol = fillcol,
             X = as.vector(row(mask)),
             Y = as.vector(col(mask)))
}

# Set directory with mask files (outdir in python script)
maskdir <- '/home/tommaullin/Documents/ConfRes'

# Set simulation number
simNo <- 17

if (simNo %in% c(1,3,11)){
  mlist = c(2)
} else if (simNo %in% c(15)){
  mlist = c(2:5)
} else if (simNo %in% c(17)){
  mlist = c(3)
}

# Set simulation directory
simdir <- paste(maskdir,paste('sim',toString(simNo),sep=''),sep=.Platform$file.sep)

# Loop through range of m values
for(m in mlist){

  # Get figure directory
  figdir <- paste(simdir,paste(toString(m),'sample',sep=''),sep=.Platform$file.sep)

  # Read in masks
  mask1 <- read.csv(file = paste(figdir,'mask1.csv',sep=.Platform$file.sep),sep=',', header=FALSE)
  mask2 <- read.csv(file = paste(figdir,'mask2.csv',sep=.Platform$file.sep),sep=',', header=FALSE)
  maskall <- read.csv(file = paste(figdir,'maskall.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

  # Read in boundaries
  bdry1 <- read.csv(file = paste(figdir,'bdry1.csv',sep=.Platform$file.sep),sep=',', header=FALSE)
  bdry2 <- read.csv(file = paste(figdir,'bdry2.csv',sep=.Platform$file.sep),sep=',', header=FALSE)
  bdryall <- read.csv(file = paste(figdir,'bdryall.csv',sep=.Platform$file.sep),sep=',', header=FALSE)


  # Add a box
  box <- matrix(data=0,nrow=100,ncol=100)
  box[1,]=1
  box[100,]=1
  box[,1]=1
  box[,100]=1

  # Create plot with two excursion sets
  img <- ggplot(mask_to_df(mask1, 0.4, "goldenrod1"), aes(x = X, y = Y, alpha = values, fill = fillcol)) +
          geom_raster() +
          geom_raster(data = mask_to_df(bdry1, 0.7, "goldenrod1")) +
          geom_raster(data = mask_to_df(mask2, 0.4, "coral2")) +
          geom_raster(data = mask_to_df(bdry2, 0.7, "coral2")) +
          scale_fill_identity() +
          scale_alpha_identity() +
          theme_light() +
          coord_equal() + 
          theme_void()

  # Add the third set if we have it
  if(m >= 3){

    # Read in mask
    mask3 <- read.csv(file = paste(figdir,'mask3.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Read in boundary
    bdry3 <- read.csv(file = paste(figdir,'bdry3.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Add to image
    img <- img +
           geom_raster(data = mask_to_df(mask3, 0.4, "seagreen3")) +
           geom_raster(data = mask_to_df(bdry3, 0.7, "seagreen3"))

  } 

  # Add the fourth set if we have it
  if(m >= 4){

    # Read in mask
    mask4 <- read.csv(file = paste(figdir,'mask4.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Read in boundary
    bdry4 <- read.csv(file = paste(figdir,'bdry4.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Add to image
    img <- img +
           geom_raster(data = mask_to_df(mask4, 0.4, "dodgerblue")) +
           geom_raster(data = mask_to_df(bdry4, 0.7, "dodgerblue"))


  } 

  # Add the fifth set if we have it
  if(m >= 5){

    # Read in mask
    mask5 <- read.csv(file = paste(figdir,'mask5.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Read in boundary
    bdry5 <- read.csv(file = paste(figdir,'bdry5.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

    # Add to image
    img <- img +
           geom_raster(data = mask_to_df(mask5, 0.4, "darkorange1")) +
           geom_raster(data = mask_to_df(bdry5, 0.7, "darkorange1"))

  } 

  #orchid3
  #orchid3

  # Add the mask for Fc and the boundary box. Then scale and theme settings.
  img <- img +        
         geom_raster(data = mask_to_df(maskall, 0.8, "plum3")) +
         geom_raster(data = mask_to_df(bdryall, 1.0, "plum4")) +
         geom_raster(data = mask_to_df(box, 0.8, "gray39")) 

  png(paste(figdir,'DemoFig',sep=.Platform$file.sep))
  print(img)
  dev.off()

}

# Make figure for heterogeneous noise
hetnoise <-read.csv(file=paste(maskdir,'heterogen_noise','data.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

# Save to png
png(paste(maskdir,'heterogen_noise','DemoFig',sep=.Platform$file.sep))
image.plot(t(hetnoise), xaxt='n', yaxt='n', ann=FALSE)
dev.off()

# Make figure for homogeneous noise
homnoise <-read.csv(file=paste(maskdir,'homogen_noise','data.csv',sep=.Platform$file.sep),sep=',', header=FALSE)

# Save to png
png(paste(maskdir,'homogen_noise','DemoFig',sep=.Platform$file.sep))
image.plot(t(homnoise), xaxt='n', yaxt='n', ann=FALSE)
dev.off()