require(gstat)
require(geoR)
library(readr)


#dat <- read_csv("E:/NMT MS/Fall 21/MATH-586-01-Spatial Variability & Geostats/Project/Assays_Pen_Bx_4m_comp_Before2010/Assays_Pen_Bx_4m_comp_Before2010.csv")

# dat$z_cat <- ifelse(dat$mid_z>0 & dat$mid_z<1300,"till 1300",
#                     ifelse(dat$mid_z>1300 & dat$mid_z<1500,"1300 to 1500","Above 1500"))
# dat$z_cat <- factor(dat$z_cat , levels=c("till 1300", "1300 to 1500", "Above 1500"))
# 
# 
# boundaries <- boxplot(AG_PPM ~z_cat,
#                       data=dat,
#                       main= paste("Concentration of ","AG_PPM", " by depth", sep = ""),
#                       xlab="Depth",
#                       ylab="Concentration",
#                       col="orange",
#                       border="brown",outline=FALSE
# )
# 
# 
# # Add sample size on top
# nbGroup <- nlevels(dat$z_cat)
# text( 
#   x=c(1:nbGroup), 
#   y=boundaries$stats[nrow(boundaries$stats),] + 2, 
#   paste("n = ",table(dat$z_cat),sep="")  
# )
# 
# 
# 

org_dat <- read_csv("E:/NMT MS/Fall 21/MATH-586-01-Spatial Variability & Geostats/Project/Assays_Pen_Bx_4m_comp_Before2010/Assays_Pen_Bx_4m_comp_Before2010.csv")

org_data <- data.frame(org_data)
org_data <- org_data[order(org_data$mid_z),]


################## Correlation and scatter plot between depth and AG_PPM

par(mfrow=c(1,1))
plot(org_data$mid_z,org_data$AG_PPM,main="Scatter plot between AG_PPM and depth",xlab= "Depth",ylab="AG_PPM")

yout=org_data$AG_PPM>1000   # some logical rule for identifying an outlier
plot(org_data$mid_z[!yout],org_data$AG_PPM[!yout])

cor(org_data$mid_z,org_data$AG_PPM,use="complete.obs")

par(mfrow=c(1,2))
hist(org_dat$AG_PPM,main = "Original AG_PPM",xlab = "AG_PPM")
box()
hist(log(org_dat$AG_PPM),main="Log transformed AG_PPM",xlab = "Log(AG_PPM)")
box()


####################################################


#define number of data frames to split into
n <- 10
#split data frame into n equal-sized data frames
data_split <- split(org_data, factor(sort(rank(row.names(org_data))%%n)))





plot3d_layers <- function(splt){
  require(rgl)    #   require(rgl,lib="c:/temp")
  
  # Data preperation
  dat <- data_prep(splt)
  
  x = dat[,1]
  y = dat[,2]    
  V = dat[,3]
  
  plot3d(x,y,V)
}

plot3d_layers(splt=1)

## Data preperation function
data_prep <- function(splt){
  
  sub_dat <- data.frame(data_split[splt])
  colnames(sub_dat) <- colnames(org_data)
  dat <- sub_dat
  dat <- dat[,c(10,11,4)]
  dat <- na.omit(dat)
  coordinates(dat) <- c("mid_x", "mid_y")
  zd <- zerodist(dat)
  dat2 <- dat[-zd[,2], ]
  dat3 <- subset(dat, !(1:nrow(dat) %in% zd[,2]))
  dat <- as.data.frame(dat3)
  
  dat$AG_PPM <- log(dat$AG_PPM)
  
  dat
}

### Variogram function
variogram_AG_PPM <- function(splt){
  
  
  # Data preperation
  dat <- data_prep(splt)
  hp = as.geodata(dat, coords.col = c(1,2), data.col = 3)   # data.names = c(east,north, Welev)) 
  
  
  
  vario2a = variog(hp, max.dist=500, direction = 0)
  plot(vario2a, type="o",main= paste("Variogram for data split ", splt))
  
  vario2b = variog(hp, max.dist=500, direction = pi/2)
  lines(vario2b, lty=2)
  
  
}

#### Kriging function
kriging_mining_AG_PPM <- function(splt){
  
  require(geoR)
  library(ggplot2)
  library(sf)
  library(grid)
  library(gridExtra)
  
  
  # Data preperation
  dat <- data_prep(splt)
  hp = as.geodata(dat, coords.col = c(1,2), data.col = 3)   # data.names = c(east,north, Welev)) 
  
  
  x = dat[,1]
  y = dat[,2]    
  V = dat[,3]
  
  
  
  # obtain kriging maps
  meanx = mean(x)
  meany = mean(y)
  x = x - meanx;  y = y - meany  
  
  
  # obtain kriging maps
  xc = seq(min(dat$mid_x),max(dat$mid_x),l=100) 
  yc = seq(min(dat$mid_y), max(dat$mid_y),l=100)
  nx = length(xc)
  ny = length(yc)
  
  gridxy = expand.grid(x = xc, y = yc)
  
  
  kc = krige.conv(hp, loc = gridxy,
                  krige = krige.control(cov.pars = c(2.5,400), nugget = 0.5, type.krige ="OK",cov.model = "sph") 
  )
  
  
  
  ok_kmean = kc$predict 
  ok_ksd = sqrt(kc$krige.var)
  
  # create the plot, the geom_contour may not be needed, but I find it helpful
  
  p1 <- ggplot(gridxy) + 
    aes(x = x, y = y, z = ok_kmean, fill = ok_kmean) + 
    geom_tile() + 
    geom_contour(color = "white", alpha = .2) + 
    scale_fill_gradientn(colours = sf.colors(3))+ 
    theme(panel.grid.major = element_line(colour = "white"))
  
  p2 <- ggplot(gridxy) + 
    aes(x = x, y = y, z = ok_ksd, fill = ok_ksd) + 
    geom_tile() + 
    geom_contour(color = "white", alpha = .2) + 
    scale_fill_gradientn(colours = sf.colors(3))+ 
    theme(panel.grid.major = element_line(colour = "white"))
  
  ####### Simple kriging ##########
  
  sk = krige.conv(hp, loc = gridxy,
                  krige = krige.control(cov.pars = c(2.5,400), nugget = .5, type.krige ="SK",cov.model = "sph",beta = mean(dat$AG_PPM)) 
  )
  
  
  
  sk_kmean = sk$predict 
  sk_ksd = sqrt(sk$krige.var)
  
  # create the plot, the geom_contour may not be needed, but I find it helpful
  
  #par(mfrow=c(2,1))
  p3 <- ggplot(gridxy) + 
    aes(x = x, y = y, z = sk_kmean, fill = sk_kmean) + 
    geom_tile() + 
    geom_contour(color = "white", alpha = .2) + 
    scale_fill_gradientn(colours = sf.colors(3))+ 
    theme(panel.grid.major = element_line(colour = "white"))
  
  p4 <- ggplot(gridxy) + 
    aes(x = x, y = y, z = sk_ksd, fill = sk_ksd) + 
    geom_tile() + 
    geom_contour(color = "white", alpha = .2) + 
    scale_fill_gradientn(colours = sf.colors(3))+ 
    theme(panel.grid.major = element_line(colour = "white"))

  # Plotting together:  
  grid.arrange(p1, p2,p3,p4, ncol=2,top = textGrob(paste("Spatial kriging data split ", splt),gp=gpar(fontsize=20,font=3)))
  
  
}




# Testing the function:
par(mfrow=c(2,2))
variogram_AG_PPM(splt=1)
variogram_AG_PPM(splt=2)
variogram_AG_PPM(splt=3)
variogram_AG_PPM(splt=4)
variogram_AG_PPM(splt=5)
variogram_AG_PPM(splt=6)
variogram_AG_PPM(splt=7)
variogram_AG_PPM(splt=8)
variogram_AG_PPM(splt=9)
variogram_AG_PPM(splt=10)


#### Kriging Plots ########

for(i in 1:1){
  library("jpeg")
  mypath <- file.path("E:/NMT MS/Fall 21/MATH-586-01-Spatial Variability & Geostats/Project/Output_file/Plots/Kriging/",paste("kriging_data_",i, ".jpg", sep = ""))
  
  jpeg(file=mypath)
  #mytitle = paste("my title is", names[i])
  kriging_mining_AG_PPM(splt=i)
  dev.off()
}



#vario2b = variog(hp, max.dist=200, direction = pi/2)
#lines(vario2b, lty=2)


