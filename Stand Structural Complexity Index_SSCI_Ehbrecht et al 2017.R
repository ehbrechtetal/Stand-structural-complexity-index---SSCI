##---------------------------------------------------------------------------------------------------------------------
## Algorithm to compute the Stand Structural Complexity Index - SSCI after Ehbrecht et al. (2017)
##---------------------------------------------------------------------------------------------------------------------
## The current version is tailored to point clouds that were scanned with a FARO Focus 120 or Faro M70 scanner and exported 
## with the scanner specific software FARO Scene (v.7.1). For data acquisition, the scanner should be placed on a tripod 
## in ~1.3 m above ground and set to scan a field of view of 300° degrees vertically and 360° horizontally with an angular 
## step width of 0.03515625°. Each scan should then be imported to the hardware-specific software FARO SCENE 
## (Faro Technologies Inc., Lake Mary, USA, v.7.1.1.81) and standard filter algorithms should be applied to each scan file 
## to erase stray and erroneous points from the point cloud. Point clouds should be exported as 5-columned .xyz files, 
## incl. the vertical and horizontal laser beam direction stored in the first two columns, respectively. Column 3,4 and 5 
## should be x, y, and z coordinates. The current version of the algorithm works with only every 4th vertical and horizontal 
## beam direction (16th of original scan resolution).
##---------------------------------------------------------------------------------------------------------------------

#install.package(data.table)
#install.package(sphereplot)
#install.package(pracma)
#install.package(plot3D)
#install.package(rgl)
#install.package(sp)
#install.package(spatialEco)
#install.package(animation)

start_time <- Sys.time()

library(data.table)
library(sphereplot)
library(pracma)
library(plot3D)
library(rgl)
library(sp)
library(spatialEco)
library(animation)

directory<-"..." # specify directory of .xyz files here
files<-list.files(directory,full.names = T,pattern="\\.xyz$")
j<-files[1]

#for(j in files){
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## import point cloud
  pointcloud <- fread(j, sep=" ")
  pointcloud<-data.frame(pointcloud)
  pointcloud<-pointcloud[,1:5]
  names(pointcloud) <-c("beam_direction_vertical","beam_direction_horizontal","x","y","z")
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## remove points below scanner height
  
  pointcloud<-pointcloud[pointcloud[,"z"]>0,] 
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## convert x,y,z coordinates to spherical coordinates
  pointcloud<-data.frame(pointcloud,cart2sph(cbind(pointcloud[,"x"],pointcloud[,"y"],pointcloud[,"z"])))
  
  ## convert to degrees
  pointcloud$theta_grad<-pointcloud$theta*180/pi
  pointcloud$phi_grad<-pointcloud$phi*180/pi
  
  ## determine angle based on laser beam direction, which is more precise than actual coordinates
  pointcloud$angle_vertical<-pointcloud$beam_direction_vertical*0.03515625
  pointcloud$angle_horizontal<-pointcloud$beam_direction_horizontal*0.03515625
  
  summary(pointcloud)
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## Constructing cross-sectional polygons
  ## In total, there are 10236 horizontal beam directions (360°/0.03515625°)
  ## The current version works with only every 4th horizontal and vertical laser beam (as exported from FARO Scene)

  horizontal_1<-seq(0,5116,4)
  horizontal_2<-seq(5120,10236,4)
  opposite<-data.frame(horizontal_1,horizontal_2)

  
  pointcloud<-pointcloud[pointcloud$beam_direction_horizontal<=10236,] 
  pointcloud$angle_vertical_mod <- ifelse(pointcloud$angle_horizontal>=180,180-pointcloud$angle_vertical,pointcloud$angle_vertical)
  pointcloud$angle_horizontal_mod <- 180-pointcloud$angle_horizontal 
  pointcloud<-pointcloud[order(pointcloud$beam_direction_horizontal,pointcloud$angle_vertical_mod),]

  stopifnot(length(unique(pointcloud$angle_horizontal))/2 == 1280)
  
  out<-data.frame(sph2cart(cbind(ifelse(pointcloud$angle_vertical_mod>90,180,-0)/180*pi, pointcloud$angle_vertical/180*pi,pointcloud$r)))
 
  pointcloud$x2D<-out$x
  pointcloud$z2D<-out$z
  pointcloud$y2D<-out$y
  
  rm(out)
  gc()
  

  ##---------------------------------------------------------------------------------------------------------------------
  ## Computing perimeter, area and fractal dimension for each cross-sectional polygon
  
  xylim<-range(c((pointcloud$x),(pointcloud$y),(pointcloud$x2D),(pointcloud$y2D)))
  zlim<-range(c(((pointcloud$z))))
  #i=1
  result<-lapply(seq(1,nrow(opposite)),function(i){
    
    test <- pointcloud[pointcloud$beam_direction_horizontal %in% c(opposite$horizontal_2[i],opposite$horizontal_1[i]),]
    test$XX<-test$x2D
    test$YY<-test$z2D
    coordinates(test)<-c("XX","YY")
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    
 
    area<-round(sapply(slot(sps, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area")),1)
    peri<-round(spatialEco::polyPerimeter(sps),1)
    frac<-2*log(0.25*peri)/log(area)
    
   
    out<-apply(test@data,2,mean)
    out<-data.frame(rbind(out))
    out$area<-area
    out$peri<-peri
    out$frac<-frac
    rownames(out)<-NULL

    
    
    return(list(out,sps,test))
    
    
  }) 
  

  ##---------------------------------------------------------------------------------------------------------------------
  ## Computing the Effective number of layers - ENL (Ehbrecht et al. 2016)

  
  pointcloud_20<-pointcloud[pointcloud$x>-20&pointcloud$x<20 & pointcloud$y>-20 & pointcloud$y<20,]
  pointcloud_20 <- subset(pointcloud_20, pointcloud_20$z >= 0)
  ##-------------------------------------------
  ## Point cloud voxelization
  pointcloud_20$x<-round(pointcloud_20$x*5)/5
  pointcloud_20$y<-round(pointcloud_20$y*5)/5
  pointcloud_20$z<-round(pointcloud_20$z*5)/5+0.5
  head(pointcloud_20[,c("x","y","z")])
  pointcloud_20<- unique(pointcloud_20[,c("x","y","z")])
  #pointcloud_20<-pointcloud_20[order(pointcloud_20$x,pointcloud_20$y,pointcloud_20$z,decreasing = F),]
  ##-------------------------------------------
  ## Stratifying point cloud intro 1m thick vertical layers
  pointcloud_20$z_layer<-round(pointcloud_20$z)
  #pointcloud_20<-pointcloud_20[order(pointcloud_20$z_layer),]

  tabelle<-table(pointcloud_20$z_layer)
  sum_voxel<-sum(tabelle)
  dat$top.height<-max(pointcloud_20$z_layer) # Maximum canopy height
  dat$height<-length(tabelle) # Number of layers
  dat$enl<-1/sum((tabelle/sum_voxel)^2) # Inverse Simpson Index --> Effective number of layers
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## Stand structural complexity Index - SSCI
  dat$ssci <- dat$frac^log(dat$enl)
  dput(names(dat))

  write.csv(dat[,c("area","peri", "frac", "name", "file", "top.height", "height", 
                     "enl", "ssci")],file= gsub(".xyz$",".csv",j),row.names = F)
  

#} 

end_time <- Sys.time()
end_time - start_time
warnings()
