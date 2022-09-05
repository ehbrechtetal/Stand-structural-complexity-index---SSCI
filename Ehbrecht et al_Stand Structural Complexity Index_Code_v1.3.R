##---------------------------------------------------------------------------------------------------------------------
## Algorithm to compute the Stand Structural Complexity Index - SSCI after Ehbrecht et al. (2017) - version 1.3
##---------------------------------------------------------------------------------------------------------------------
## The current version is tailored to point clouds that were scanned with a FARO Focus 120 or Faro M70 scanner and exported 
## with the scanner specific software FARO Scene (v.7.1). For data acquisition, the scanner should be placed on a tripod 
## in ~1.3 m above ground and set to scan a field of view of 300 degrees vertically and 360 horizontally with an angular 
## step width of 0.03515625. Each scan should then be imported to the hardware-specific software FARO SCENE 
## (Faro Technologies Inc., Lake Mary, USA, v.7.1.1.81) and standard filter algorithms should be applied to each scan file 
## to erase stray and erroneous points from the point cloud. Point clouds should be exported as 5-columned .xyz files, 
## incl. the vertical and horizontal laser beam direction stored in the first two columns, respectively. Column 3,4 and 5 
## should be x, y, and z coordinates. The current version of the algorithm works with only every 4th vertical and horizontal 
## beam direction (16th of original scan resolution).
##
## Update to version 1.3 now includes the computation of canopy openness for an opening angle of 60° above the scanner. 
## Canopy openness is computed by simulating a hemispherical photo and calculated as the percentage of sky pixels 
## in the simulated hemispherical photograph
##---------------------------------------------------------------------------------------------------------------------

start_time <- Sys.time()

library(data.table)
library(sphereplot)
library(pracma)
library(plot3D)
library(rgl)
library(sp)
library(spatialEco)
library(animation)
library(rgeos) 
library(raster)


##---------------------------------------------------------------------------------------------------------------------
# Specify directory where .xyz files are located
##---------------------------------------------------------------------------------------------------------------------

directory<-"..." 

files<-list.files(directory,full.names = T,pattern="\\.xyz$")
j<-files[1]


for(j in files){
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## import point cloud
  ##---------------------------------------------------------------------------------------------------------------------
  pointcloud <- fread(j, sep=" ")
  pointcloud<-data.frame(pointcloud)
  pointcloud<-pointcloud[,1:5]
  names(pointcloud) <-c("beam_direction_vertical","beam_direction_horizontal","x","y","z")
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## convert x,y,z coordinates to spherical coordinates
  ##---------------------------------------------------------------------------------------------------------------------
  pointcloud<-data.frame(pointcloud,cart2sph(cbind(pointcloud[,"x"],pointcloud[,"y"],pointcloud[,"z"])))
  ##---------------------------------------------------------------------------------------------------------------------
  ## convert to degrees
  ##---------------------------------------------------------------------------------------------------------------------
  pointcloud$theta_grad<-pointcloud$theta*180/pi
  pointcloud$phi_grad<-pointcloud$phi*180/pi
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## determine angle based on laser beam direction, which is more precise than actual coordinates
  ##---------------------------------------------------------------------------------------------------------------------
  range( pointcloud$beam_direction_horizontal)
  pointcloud$pointcloudGradZenith <- pointcloud$beam_direction_vertical * 0.03515625
  pointcloud$angle_horizontal <- pointcloud$beam_direction_horizontal *  0.03515625
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## remove points below scanner height
  ##---------------------------------------------------------------------------------------------------------------------
  pointcloud<-pointcloud[pointcloud[,"z"]>=0,] 
  
  ##---------------------------------------------------------------------------------------------------------------------
  ## Constructing cross-sectional polygons
  ## In total, there are 10236 horizontal beam directions (360?/0.03515625?)
  ## The current version works with only every 4th horizontal and vertical laser beam (as exported from FARO Scene)
  ##---------------------------------------------------------------------------------------------------------------------
  horizontal_1<-seq(0,5116,4)
  horizontal_2<-seq(5120,10236,4)
  opposite<-data.frame(horizontal_1,horizontal_2)
  
  pointcloud<-pointcloud[pointcloud$beam_direction_horizontal<=10240,]
  pointcloud<-pointcloud[pointcloud$pointcloudGradZenith<=90,] # Aenderung Micha
  pointcloud$angle_vertical_mod <- ifelse(pointcloud$angle_horizontal>=180,
                                          180-pointcloud$pointcloudGradZenith,
                                          pointcloud$pointcloudGradZenith)
  pointcloud$angle_horizontal_mod <- 180-pointcloud$angle_horizontal 
  pointcloud<-pointcloud[order(pointcloud$beam_direction_horizontal,
                               pointcloud$angle_vertical_mod),] #Sorting by beam direction and then by corrected zenith angle (angle_vertical_mod)
  
  out<-data.frame(sph2cart(cbind(ifelse(pointcloud$angle_vertical_mod>=90,180,0)/180*pi,
                                 (90-pointcloud$pointcloudGradZenith)/180*pi,pointcloud$r)))
  
  pointcloud$x2D<-out$x
  pointcloud$z2D<-out$z
  pointcloud$y2D<-out$y
  
  rm(out)
  gc()
  
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Computation of canopy openness for 60° opening angle above scanner
  ##---------------------------------------------------------------------------------------------------------------------
  
  hemi<-data.frame(radius_polar=2*tan(abs(deg2rad(90)-pointcloud$phi)/2),
                   winkel_polar=pointcloud$theta,
                   phi_kugel=  pointcloud$phi_grad)
  
  hemi<-hemi[hemi$phi_kugel>=60,] 
  
  if (nrow(hemi)!=0){ 
    
    hemi<-data.frame(hemi,pracma::pol2cart(cbind(hemi$winkel_polar, # phi, r, z
                                                 hemi$radius_polar)))
    
    library(ggplot2)
    ggplot(hemi,aes(x=x,y=y,col=phi_kugel))+
     geom_point( size = 0.01, fill = NA,stroke =0)
    ggplot(hemi,aes(x=x,y=y,col=winkel_polar))+geom_point(size=.0001)
    
    #range(hemi$x)
    #range(hemi$y)
    bereich<-2*tan(abs(deg2rad(90)-deg2rad(60))/2)
    fac=180 # best fit between Mathematica and R
    r <- raster(xmn= -bereich, 
                xmx= bereich,
                ymn= -bereich,
                ymx= bereich,
                nrows=fac,ncols=fac) #
    
    r <- raster::rasterize(cbind(hemi$x, hemi$y), r, hemi$winkel_polar,
                           fun = "count")
    
    A_kreis=  (bereich)^2*pi 
    A_quadrat=(bereich*2)^2 
    area_pixel<-A_quadrat/fac/fac
    A_kreis/A_quadrat 
    
    r[is.na(r[])]<-0 #set NA to zero (white celles outside of hemi-circle)
    r[r>0] <- 1 # set all cells with count >0 to 1
    # plot(r) 
    filled<-cellStats(r,sum) # how many pixels are filled?
    (openness<-1-cellStats(r,sum)*area_pixel/A_kreis) # area of filled pixel as percentage of area circle
  } else{ 
    openness <- 1} # if no there are no scanner points within 60 degrees opening angle (sky is fully visible)
  
  
  ##-------------------------------------------
  # Plot settings
  ##-------------------------------------------
  par(mfrow=c(3,3))
  xylim<-range(c((pointcloud$x),(pointcloud$y),(pointcloud$x2D),(pointcloud$z2D)))
  zlim<-c(0,max(c(pointcloud$z,pointcloud$z2D)))
  
  
  ##-------------------------------------------
  # Loop over all cross-sections
  ##-------------------------------------------
  result<-lapply(seq(1,nrow(opposite)),function(i){
    
    ##-------------------------------------------
    # get front and rear side of cross-sections
    ##-------------------------------------------
    test_left <- pointcloud[pointcloud$beam_direction_horizontal %in% c(opposite$horizontal_2[i]),]
    test_right <- pointcloud[pointcloud$beam_direction_horizontal %in% c(opposite$horizontal_1[i]),]
    
    ##-------------------------------------------
    # adding 0-0-0 (scanner position) as x-y-z coordinate 
    ##-------------------------------------------
    test_middle<-test_left[1,] 
    test_middle[nrow(test_middle),]<-rep(NA,ncol(test_middle))
    test_middle$x2D<-0
    test_middle$z2D<-0
    
    test<-rbind(test_right,test_middle,test_left)
    
    ##-------------------------------------------
    # Connecting points to build polygon
    ##-------------------------------------------
    coordinates(test)<-c("x2D","z2D")
    p = Polygon(test)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    
    #---------------------------------------------------------------------------------------------------------------------
    # Stop if polygon has self-intersection (can happen in very open sites with less than 5 returns)
    #---------------------------------------------------------------------------------------------------------------------
    # stopifnot(gIsValid(sps,byid = TRUE))
    
    #---------------------------------------------------------------------------------------------------------------------
    # Computing area, perimeter and fractal dimension of cross-sectional polygons
    #---------------------------------------------------------------------------------------------------------------------
    area<-(sapply(slot(sps, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area")))
    peri<-(spatialEco::polyPerimeter(sps))
    frac<-2*log(0.25*peri)/log(area)
    
    ##---------------------------------------------------------------------------------------------------------------------
    # Plotting cross-sectional polygons
    ##---------------------------------------------------------------------------------------------------------------------
    plotting <- T  # Plotting cross-sectional polygons (True or False)
    
    mult<- 100 # determines i.th polygon that is plotted for visualization
    
    multiple_of_mult = (i %% mult) == 0
    if(i==1){multiple_of_mult=T}
    
    if(plotting==T & multiple_of_mult ){
      
      plot(sps,axes=T,main=paste0("Cross-sectional Polygon ",i," of ", nrow(opposite)))
      points(coordinates(test),col="red",cex=.5)
      grid()
      gc()
      
    }
    
    test@data$side<-NULL
    out<-apply(test@data,2,mean,na.rm=T)
    out<-data.frame(rbind(out))
    out$area<-area
    out$peri<-peri
    out$frac<-frac
    rownames(out)<-NULL
    
    return(list(out,sps,test))
    
  }) #End of loop over all cross-sections
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Computing mean polygon area, perimeter and fractal dimension (MeanFrac)
  ##---------------------------------------------------------------------------------------------------------------------
  dat<-lapply(result,'[[', 1)  
  dat<-do.call(rbind,dat)
  dat<-as.data.frame(t(colMeans(dat, na.rm=T)))
  dat$name<- gsub(".xyz","",basename(j))
  dat$filename	<-j
  
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Computing the Effective number of layers - ENL (Ehbrecht et al. 2016)
  ##---------------------------------------------------------------------------------------------------------------------
  
  ##-------------------------------------------
  # Subsetting to 20x20m space around scanner
  ##-------------------------------------------
  
  pointcloud_20<-pointcloud[pointcloud$x>-20&pointcloud$x<20 & pointcloud$y>-20 & pointcloud$y<20,]
  pointcloud_20 <- subset(pointcloud_20, pointcloud_20$z >= 0)
  
  ##-------------------------------------------
  ## Point cloud voxelization
  ##-------------------------------------------
  
  pointcloud_20$x<-round(pointcloud_20$x*5)/5
  pointcloud_20$y<-round(pointcloud_20$y*5)/5
  pointcloud_20$z<-round(pointcloud_20$z*5)/5+0.5
  head(pointcloud_20[,c("x","y","z")])
  pointcloud_20<- unique(pointcloud_20[,c("x","y","z")])
  
  pointcloud_20$z_layers<-round(pointcloud_20$z)
  
  
  tab<-table(pointcloud_20$z_layers)
  sum_voxel<-sum(tab)
  
  ##-------------------------------------------
  # Height of highest z-coordinate (top height in m)
  ##-------------------------------------------
  dat$top.height<-max(pointcloud_20$z_layers,na.rm=T) 
  
  ##-------------------------------------------
  # Number of filled layers
  ##-------------------------------------------
  dat$no.filled.layers<-length(tab) 
  
  ##-------------------------------------------
  # Effective number of layers
  ##-------------------------------------------
  dat$enl<-1/sum((tab/sum_voxel)^2) 
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Stand structural complexity index
  ##---------------------------------------------------------------------------------------------------------------------
  dat$ssci <- dat$frac^log(dat$enl) 
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Canopy openness
  ##---------------------------------------------------------------------------------------------------------------------
  
  dat$openness<-openness*100
  rm(openness)
  
  ##---------------------------------------------------------------------------------------------------------------------
  # Exporting Results
  ##---------------------------------------------------------------------------------------------------------------------  
  dput(names(dat))
  results <- dat[,c("filename","area","peri", "frac", "top.height", "no.filled.layers", "enl", 
                    "ssci","openness")]
  colnames(results) <- c("File name", "Mean polygon area (m2)", "Mean polygon perimeter (m)", 
                         "Mean fractal dimension (MeanFrac)",  "Top height (m)", 
                         "Number of filled layers", "Effective number of layers (ENL)",
                         "Stand structural complexity Index (SSCI)","Canopy openness (%)")
  write.csv(results,file= gsub(".xyz$",".csv",j),row.names = F)
  
  
} # End of loop over all files


end_time <- Sys.time()
end_time - start_time


