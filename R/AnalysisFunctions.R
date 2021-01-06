
## Functions used in Hydraulic Geometry script
## Last updated 4 January 2020 by LE Koenig


### ===================================================================== ###
###     Calculate drainage area of each segment by splitting up COMID     ###
### ===================================================================== ###

calc.upstream.area <- function(flines,locs){
  
  # Load flowlines
  flines.sub <- flines[which(flines$COMID == locs$COMID),] %>% st_transform(.,5070)
  
  # Define datum for point locations:
  if(locs$datum == "NAD83"){datum = 4269}
  if(locs$datum == "WGS84"){datum = 4326}
  if(! locs$datum %in% (c("NAD83","WGS84"))){print("Warning: Check that datum of point locations equals NAD83 or WGS84")}
  locs.sp <- st_as_sf(locs,coords=c("LON_DD","LAT_DD"),crs=datum) %>% st_transform(.,5070)
  
  # Cast point locations as sf point object:
  #target.reach.pt <- suppressWarnings(st_cast(locs.sp,"POINT"))
  #target.reach.pt <- target.reach.pt[nrow(target.reach.pt),]
  
  # Cast flowline into points:
  flines.sub2 <- st_cast(st_line_merge(st_union(st_cast(flines.sub, "MULTILINESTRING"))), "LINESTRING")
  flines.sub2.pts <- flines.sub2 %>% st_cast(.,"POINT") %>% st_as_sf(.)
  flines.sub2.pts.coords <- sf::st_coordinates(flines.sub2.pts)
  
  # Match point locations to nearest flowline point segment:
  matched <- st_nearest_feature(locs.sp,flines.sub2.pts)
  empty <- st_as_sfc("POINT(EMPTY)")
  pts <- flines.sub2.pts %>%
    mutate(distance_to_next = as.numeric(sf::st_distance(x, lead(x, default = empty), by_element = TRUE))
    )
  # Calculate the proportion along the reach at which the nearest point location lies:
  reach.proportion <- sum(pts$distance_to_next[c(1:matched)],na.rm=T)/sum(pts$distance_to_next,na.rm=T)
  
  # Calculate upstream area based on proportion of NHD flowline and NHD VAA for that flowline:
  Area_corr <- round(flines.sub$AREASQKM * reach.proportion,5)
  TotArea_corr <- flines.sub$TOTDASQKM - (flines.sub$AREASQKM - Area_corr)
  
  return(TotArea_corr)
  
}


### ===================================================================== ###
###                 Snap point to closest NHDPlusV2 reach                 ###
### ===================================================================== ###

## 5. This function calculates geodesic distances between GRDO points and flowlines to find nearest NHDPlusV2 flowline

link_comid5 <-  function(flowline.data,points,vpu.polygon.data,max.dist){
  
  points.trans <- points %>% st_transform(.,st_crs(vpu.polygon.data))
  vpu <- vpu.polygon.data$VPUID[unlist(st_is_within_distance(points.trans,vpu.polygon.data,sparse=TRUE,dist=max.dist))]
  
  points.transform <- points %>% st_transform(.,st_crs(flowline.data))
  flowline.sub <- flowline.data[which(flowline.data$VPUID %in% vpu), ] 
  
  c <- st_distance(x=points.transform,y=flowline.sub,which="Great Circle",by_element = FALSE)
  join.dat <- bind_cols(points,st_drop_geometry(flowline.sub[which.min(c),c("COMID","GNIS_NAME","REACHCODE","FTYPE","FCODE","STREAMCALC","STREAMORDE","AREASQKM","TOTDASQKM","SLOPE","SLOPELENKM","QE_MA","VE_MA","VPUID")]))
  join.dat$near_dist_m <- as.numeric(c[which.min(c)])
  
  return(st_drop_geometry(join.dat))
}
  
  



