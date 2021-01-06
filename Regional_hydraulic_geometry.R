# Use national river datasets to estimate regional hydraulic geometry relationships  
# L.E. Koenig
# Last modified 4 January 2021


# Load necessary packages:
library(dplyr)         # general data manipulation
library(ggplot2)       # plotting and visualization
library(cowplot)       # plotting and visualization  
library(nhdplusTools)  # interface with NHD in R
library(sf)            # spatial analyses
library(foreach)       # parallel computing in R
library(smatr)         # standardized major axis regression for allometric relationships
library(RColorBrewer)  # create color palette
library(knitr)         # knit Rmarkdown document to HTML
library(patchwork)     # create multi-panel figures
library(grid)          # add annotations to ggplots

options(scipen = 999)

# Call functions defined in AnalysisFunctions.R file:
source("./R/AnalysisFunctions.R")

## Note that user will have to specify local directory where seamless national NHDPlus data is stored ("Load National Hydrography Dataset" section below)


##========================================================##
##                    LOAD RAW DATA                       ##
##========================================================##

# 1. Load raw data - WADEABLE STREAMS ASSESSMENT (2004-2005)
# Downloaded 11/25/2020 from https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys
  wsa.width <- read.csv("./raw_data/WSA_20042005/bankgeometry.csv",header=TRUE,stringsAsFactors = FALSE)
  wsa.sites <- read.csv("./raw_data/WSA_20042005/wsa_siteinfo_ts_final.csv",header=TRUE,stringsAsFactors = FALSE)

# 2. Load raw data - NATIONAL RIVERS AND STREAMS ASSESSMENT (2008-2009)
# Downloaded 11/25/2020 from https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys
  nrsa.width <- read.csv("./raw_data/NRSA_20082009/phablow.csv",header=TRUE,stringsAsFactors = FALSE)
  nrsa.sites <- read.csv("./raw_data/NRSA_20082009/siteinfo_0.csv",header=TRUE,stringsAsFactors = FALSE)

  
##========================================================##
##     Summarize wetted width data from WSA and NRSA      ##
##========================================================##
  
## Summarize measured width by site - WADEABLE STREAMS ASSESSMENT (2004-2005)

  # Calculate mean width and mean bankfull width for each site (average multiple transects into a single site mean):
  wsa.width.mean <- wsa.width %>% 
                    group_by(SITE_ID) %>% 
                    summarize(WettedWidth_m = mean(WT_WID,na.rm=T),
                              BankfullWidth_m = mean(BANKWID,na.rm=T),
                              WaterDepth_cm = mean(DEPTH,na.rm=T))
  # Pull out site info:
  wsa.sites2 <- wsa.sites %>% filter(.,VISIT_NO == 1) %>% select(SITE_ID,YEAR,STATE,LON_DD,LAT_DD,STRAHLER)
  
  # Combine WSA datasets together:
  wsa.dat <- left_join(wsa.sites2,wsa.width.mean,by="SITE_ID")
  
  
## Summarize measured width by site - NATIONAL RIVERS AND STREAMS ASSESSMENT (2008-2009)
  
  # Calculate mean width and mean bankfull width for each site (average multiple transects into one site mean):
  nrsa.width2 <- nrsa.width %>% filter(.,VISIT_NO == 1) %>% group_by(SITE_ID) %>% select(SITE_ID,YEAR,STRAHLERORDER,XWIDTH,XBKF_W,XDEPTH_CM)
  
  # Pull out site info:
  nrsa.sites2 <- nrsa.sites %>% filter(.,VISIT_NO == 1) %>% select(SITE_ID,HUC8,STATE,LON_DD83,LAT_DD83,WSAREA_NARS)
  
## Format NRSA datasets: 
  nrsa.dat <- left_join(nrsa.sites2,nrsa.width2,by="SITE_ID") %>% 
              rename("WettedWidth_m" = "XWIDTH",
                     "BankfullWidth_m" = "XBKF_W",
                     "WaterDepth_cm" = "XDEPTH_CM",
                     "WSAREA" = "WSAREA_NARS",
                     "LON_DD" = "LON_DD83","LAT_DD" = "LAT_DD83") %>%
              mutate(STRAHLER = as.integer(substr(STRAHLERORDER,1,1))) %>%
              select(SITE_ID,YEAR,STATE,LON_DD,LAT_DD,STRAHLER,WettedWidth_m,BankfullWidth_m,WaterDepth_cm)

  # Combine WSA and NRSA sites:
  wsa.dat$assess <- "wsa"
  nrsa.dat$assess <- "nrsa"
  dat <- rbind(wsa.dat,nrsa.dat)
  dat$lon <- dat$LON_DD
  dat$lat <- dat$LAT_DD
  dat$datum <- "NAD83"
  
  # Create spatial object:
  dat.sp <- dat %>% st_as_sf(.,coords=c("lon","lat"),crs=4269)
  
  # Convert data to long format:
  dat.long <- tidyr::gather(dat[,c("SITE_ID","WettedWidth_m","BankfullWidth_m")], width_measure, value, WettedWidth_m:BankfullWidth_m, factor_key=TRUE)
  levels(dat.long$width_measure) <- c("Wetted width (m)", "Bankfull width (m)")
  
  # Plot distribution of wetted width and bankfull width across national datasets:
  dat.long %>% 
    ggplot() + geom_density(aes(x=value,fill=width_measure)) + facet_wrap(~width_measure,nrow=2) + 
    scale_x_log10() + scale_fill_manual(values=c("#8da0cb","#66c2a5"))+
    theme_cowplot() + theme(legend.position="none",axis.text = element_text(size=9),axis.title=element_text(size=10))
  
  
##========================================================##
##               Load NHDPlusV2 flowlines                 ##
##========================================================##
  
## Load National Hydrography Dataset (NHDPlus_V2) 
  # Seamless National-scale data [download here: https://www.epa.gov/waterdata/nhdplus-national-data] (note that this dataset may be memory-intensive)
  
  # Indicate where NHD data is stored locally (using helper functions from nhdplusTools package):
  nhdplusTools::nhdplus_path("/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb/")  
  staged_data <- nhdplusTools::stage_national_data(include="flowline",output_path = "./output")
  
  # Read in flowline data and filter out coastlines (to prevent point locations from being spatially joined to a coastline feature):
  flowline <- readRDS(staged_data$flowline)
  flowline <- flowline[-which(flowline$FTYPE=="Coastline"),]
  rm(staged_data)
  
  # Create consistent (all-caps) column names: 
  names(flowline)[1:(length(names(flowline))-1)] <- toupper(names(flowline)[1:(length(names(flowline))-1)])
  
## Load U.S. VPU Hydroregions:
  nhd.region <- read_sf(dsn="/Users/LKoenig/Documents/SpatialData/NHDPlus/NHDPlusNationalData/",layer="USA_HydroRegions_VPU02")
  nhd.region.Lower48 <- nhd.region[-which(is.na(nhd.region$VPUID)|nhd.region$VPUID=="20"|nhd.region$VPUID=="21"|nhd.region$VPUID=="22"),] %>% st_transform(.,5070)
  
  
##========================================================##
##         Join data to nearest NHDPlus flowline          ##
##========================================================##
  
# Cache results of this code chunk since it takes awhile to run:
  if(file.exists("./output/EPAdat_JoinNHD.rds")){
      dat_JoinNHD <- readRDS("./output/EPAdat_JoinNHD.rds")
  } else {
  
    # Use get_flowline_index from NHDPlusTools package to identify nearest flowline segment:
    dat_index <- get_flowline_index(st_transform(flowline, 5070),
                                      st_transform(dat.sp, 5070),
                                      search_radius = 700)
    # Join flowline index data with spatial points:
    dat.sp2 <- left_join(mutate(dat.sp,id=as.integer(row.names(dat.sp))),dat_index[,c("id","COMID")],by="id") %>%
               rename(.,comid_nhdplustools = COMID)
    
    # Find the geodesic distance between each field site location and the closest NHDPlusV2 flowline:
    
    # define function to monitor progress of foreach:
    progress <- function(n) cat(sprintf("row %d is complete\n", n))
    opts <- list(progress=progress)
    # register parallel backend:
    cl <- parallel::makeCluster(2,outfile="")
    doSNOW::registerDoSNOW(cl)
    # join points to NHD flowlines (minimize geodesic distance):
    dat_JoinNHD_ls <- foreach(i=1:length(dat.sp2$SITE_ID),.packages = c("dplyr","sf"),.options.snow = opts,.errorhandling = 'remove') %dopar% {
      dat <- link_comid5(flowline.data = flowline,points = dat.sp2[i,],vpu.polygon.data = nhd.region.Lower48,max.dist=20000)
      return(dat)
    }
    # close parallel cores:
    parallel::stopCluster(cl)
    
    # combine list into a data frame:
    dat_JoinNHD <- do.call("rbind",dat_JoinNHD_ls)
    rm(dat_JoinNHD_ls)
    
    saveRDS(dat_JoinNHD,"./output/EPAdat_JoinNHD.rds")
    
  }
  
 
##========================================================##
##         Estimate upstream area for each site           ##
##========================================================##
  
# Cache results of this code chunk:
if(file.exists("./output/EPAdat_JoinNHD2.rds")){
  
  dat_JoinNHD2 <- readRDS("./output/EPAdat_JoinNHD2.rds")
  
} else {
    
    # omit lines which near_dist_m > 100 m 
    dat_JoinNHD <- dat_JoinNHD %>% mutate_at(vars(COMID:VE_MA), funs(ifelse(near_dist_m > 100, NA, .)))
    
    # Create new data frame that omits the rows where near_dist_m was greater than 100 m (and that we are therefore not confident that the sample location was matched to the appropriate flowline):
    dat_JoinNHD2 <- dat_JoinNHD[which(!is.na(dat_JoinNHD$COMID)),]
    
    # Calculate upstream area for each segment after adjusting for proportional length along COMID flowline:
    dat_JoinNHD2$Calc_UpstrArea_km2 <- NA
    for(i in 1:length(dat_JoinNHD2$SITE_ID)){
      up.area <- calc.upstream.area(flines=flowline,locs = dat_JoinNHD2[i,])
      dat_JoinNHD2$Calc_UpstrArea_km2[i] <- up.area
      print(i)
    }
    
    saveRDS(dat_JoinNHD2,"./output/EPAdat_JoinNHD2.rds")
    
  }
  
  
##========================================================##
##          Estimate width/depth: national-scale          ##
##========================================================##
  
## 1. Use empirical equations from Raymond et al. 2012 to estimate mean channel width and depth [https://doi.org/10.1215/21573689-1597669]
  dat_JoinNHD2$EstWidth_Raymond <- 12.88*(dat_JoinNHD2$QE_MA*0.0283168)^0.42
  dat_JoinNHD2$EstDepth_Raymond <- 0.408*(dat_JoinNHD2$QE_MA*0.0283168)^0.294

## 2. Estimate coefficients of standardizd major axis regression (SMA) between field-measured width/depth and upstream area:
  # wetted width:
  sma.national <- sma(log10(dat_JoinNHD2$WettedWidth_m) ~ log10(dat_JoinNHD2$Calc_UpstrArea_km2), robust=T)
  print(sma.national)
  
  # bankfull width:
  sma.national.bkfl <- sma(log10(dat_JoinNHD2$BankfullWidth_m[-which(dat_JoinNHD2$BankfullWidth_m==0)]) ~ log10(dat_JoinNHD2$Calc_UpstrArea_km2[-which(dat_JoinNHD2$BankfullWidth_m==0)]), robust=T)
  
  # thalweg depth:
  dat_JoinNHD2$depth_m <- dat_JoinNHD2$WaterDepth_cm/100
  sma.national.depth <- sma(log10(dat_JoinNHD2$depth_m) ~ log10(dat_JoinNHD2$Calc_UpstrArea_km2), robust=T)

## 3. Apply coefficients to estimate width/depth given national-scale SMA:  
  # wetted width:  
  dat_JoinNHD2$EstWidth_NatlSMA <- 10^((log10(dat_JoinNHD2$Calc_UpstrArea_km2) * sma.national$coef[[1]][2,1])+sma.national$coef[[1]][1,1])
  
  # bankfull width:
  dat_JoinNHD2$EstBkflWidth_NatlSMA <- 10^((log10(dat_JoinNHD2$Calc_UpstrArea_km2) * sma.national.bkfl$coef[[1]][2,1])+sma.national.bkfl$coef[[1]][1,1])
  
  # thalweg depth:
  dat_JoinNHD2$EstDepth_NatlSMA <- 10^((log10(dat_JoinNHD2$Calc_UpstrArea_km2) * sma.national.depth$coef[[1]][2,1])+sma.national.depth$coef[[1]][1,1])
  
## 4. Estimate model RMSE:
  # wetted width:  
  sma.fit.lm <- lm(dat_JoinNHD2$WettedWidth_m~dat_JoinNHD2$EstWidth_NatlSMA)
  nat.RMSE.sma <- sqrt(mean(sma.fit.lm$residuals^2))
  
  raymond.fit.lm <- lm(dat_JoinNHD2$WettedWidth_m ~ dat_JoinNHD2$EstWidth_Raymond)
  nat.RMSE.raymond <- sqrt(mean(raymond.fit.lm$residuals^2))
  
  # bankfull width:
  sma.bkfl.fit.lm <- lm(dat_JoinNHD2$BankfullWidth_m~dat_JoinNHD2$EstBkflWidth_NatlSMA)
  nat.bkfl.RMSE.sma <- sqrt(mean(sma.bkfl.fit.lm$residuals^2))
  
  # thalweg depth:
  sma.depth.fit.lm <- lm(dat_JoinNHD2$depth_m~dat_JoinNHD2$EstDepth_NatlSMA)
  nat.depth.RMSE.sma <- sqrt(mean(sma.depth.fit.lm$residuals^2))
  
  raymond.depth.fit.lm <- lm(dat_JoinNHD2$depth_m ~ dat_JoinNHD2$EstDepth_Raymond)
  nat.depth.RMSE.raymond <- sqrt(mean(raymond.depth.fit.lm$residuals^2))
  
  
##========================================================##
##   Plot area-width/depth relationships: national-scale  ##
##========================================================## 
  
## Width/depth-Area plots: 
  grob <- grid::grobTree(textGrob(paste("slope = ",round(sma.national$coef[[1]][2,1],3),"\nint = ",round(sma.national$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.national$r2),2),sep=""), x=0.08,  y=0.85, hjust=0,
                                  gp=gpar(col="black", fontsize=9)))
  
  plot1 <- dat_JoinNHD2 %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(WettedWidth_m)),color="#66c2a5",alpha=.5) + 
    geom_abline(slope = sma.national$coef[[1]][2,1],intercept=sma.national$coef[[1]][1,1],color="royalblue3")+
    labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Wetted~width~(m)))+
    theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10))+
    annotation_custom(grob)
  
  grob2 <- grid::grobTree(textGrob(paste("slope = ",round(sma.national.bkfl$coef[[1]][2,1],3),"\nint = ",round(sma.national.bkfl$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.national.bkfl$r2),2),sep=""), x=0.08,  y=0.85, hjust=0,
                                   gp=gpar(col="black", fontsize=9)))
  
  plot2 <- dat_JoinNHD2 %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(BankfullWidth_m)),color="#8da0cb",alpha=.5) + 
    geom_abline(slope = sma.national.bkfl$coef[[1]][2,1],intercept=sma.national.bkfl$coef[[1]][1,1],color="royalblue3")+
    labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Bankfull~width~(m)))+
    theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10))+
    annotation_custom(grob2)
  
  grob3 <- grid::grobTree(textGrob(paste("slope = ",round(sma.national.depth$coef[[1]][2,1],3),"\nint = ",round(sma.national.depth$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.national.depth$r2),2),sep=""), x=0.08,  y=0.85, hjust=0,
                                   gp=gpar(col="black", fontsize=9)))
  
  plot3 <- dat_JoinNHD2 %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(depth_m)),color="#fc8d62",alpha=.5) + 
    geom_abline(slope = sma.national.depth$coef[[1]][2,1],intercept=sma.national.depth$coef[[1]][1,1],color="royalblue3")+
    labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Thalweg~depth~(m)))+
    theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10))+
    annotation_custom(grob3)
  
  plot1 + plot2 + plot3 + plot_layout(ncol=3)
  
## Density plots of the distribution of estimated width/depth:
  dens1 <- dat_JoinNHD2 %>% ggplot() + geom_density(aes(x=WettedWidth_m,fill="Field"),alpha=.5) + 
    geom_density(aes(x=EstWidth_NatlSMA,fill="Mod-NationalSMA"),alpha=.5) + 
    geom_density(aes(x=EstWidth_Raymond,fill="Mod-Raymond2012"),alpha=.5) +
    #scale_fill_manual(values=c("#33a02c","#1f78b4","#b2df8a"))+
    scale_x_log10()+ theme_cowplot() + labs(x=expression(width~(m)))+
    theme(legend.title = element_blank(),axis.text = element_text(size=9),axis.title=element_text(size=10),legend.position="none")
  dens2 <- dat_JoinNHD2 %>% ggplot() + geom_density(aes(x=BankfullWidth_m,fill="Field"),alpha=.5) + 
    geom_density(aes(x=EstBkflWidth_NatlSMA,fill="Mod-NationalSMA"),alpha=.5) + 
    scale_x_log10()+ theme_cowplot() + labs(x=expression(Bankfull~width~(m)))+
    theme(legend.title = element_blank(),axis.text = element_text(size=9),axis.title=element_text(size=10),legend.position="none")
  dens3 <- dat_JoinNHD2 %>% ggplot() + geom_density(aes(x=depth_m,fill="Field"),alpha=.5) + 
    geom_density(aes(x=EstDepth_NatlSMA,fill="Mod-NationalSMA"),alpha=.5) + 
    geom_density(aes(x=EstDepth_Raymond,fill="Mod-Raymond2012"),alpha=.5) +
    scale_x_log10()+ theme_cowplot() + labs(x=expression(Thalweg~depth~(m)))+
    theme(legend.title = element_blank(),axis.text = element_text(size=9),axis.title=element_text(size=10))
  
  dens1 + dens2 + dens3 + plot_layout(ncol=3)
  
  
##========================================================##
##              Estimate width/depth by region            ##
##========================================================##

# Classify HUC02 as a factor variable:
dat_JoinNHD2$VPUID <- as.factor(dat_JoinNHD2$VPUID)

# Re-assign bankfull width values equal to zero:
length(which(dat_JoinNHD2$BankfullWidth_m==0))
dat_JoinNHD2$BankfullWidth_m <- ifelse(is.na(dat_JoinNHD2$WettedWidth_m),NA,dat_JoinNHD2$BankfullWidth_m)
dat_JoinNHD2$BankfullWidth_m <- ifelse(dat_JoinNHD2$BankfullWidth_m==0,dat_JoinNHD2$BankfullWidth_m+0.0001,dat_JoinNHD2$BankfullWidth_m)
  
# Create empty values that will be filled in within the loop:
dat_JoinNHD2$EstBkflWidth_ReglSMA <- NA
dat_JoinNHD2$EstDepth_ReglSMA <- NA
regional.coef.width <- list()
regional.coef.bkflw <- list()
regional.coef.depth <- list()
regional.RMSE <- list()
plot.list.width <- list()
plot.list.bkfl <- list()
plot.list.depth <- list()
  
# Estimate coefficients of standardized major axis regression (SMA):
dat_JoinNHD2_ls <- split(dat_JoinNHD2,f = dat_JoinNHD2$VPUID)
colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(names(dat_JoinNHD2_ls)))
  
for(i in 1:length(names(dat_JoinNHD2_ls))){
  # Estimate coefficients of standardizd major axis regression (SMA):
  sma.regional <- sma(log10(dat_JoinNHD2_ls[[i]]$WettedWidth_m) ~ log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2), robust=T)
  sma.regional.bkfl <- sma(log10(dat_JoinNHD2_ls[[i]]$BankfullWidth_m) ~ log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2), robust=T)
  sma.regional.depth <- sma(log10(dat_JoinNHD2_ls[[i]]$depth_m) ~ log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2), robust=T)
  
  # Apply coefficients to estimate width given national-scale relationship:  
  dat_JoinNHD2_ls[[i]]$EstWidth_ReglSMA <- 10^((log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2) * sma.regional$coef[[1]][2,1])+sma.regional$coef[[1]][1,1])
  dat_JoinNHD2_ls[[i]]$EstBkflWidth_ReglSMA <- 10^((log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2) *  sma.regional.bkfl$coef[[1]][2,1])+sma.regional.bkfl$coef[[1]][1,1])
  dat_JoinNHD2_ls[[i]]$EstDepth_ReglSMA <- 10^((log10(dat_JoinNHD2_ls[[i]]$Calc_UpstrArea_km2) * sma.regional.depth$coef[[1]][2,1])+sma.regional.depth$coef[[1]][1,1])
    
  # Save coefficients:
  regional.coef.width[[i]] <- sma.regional
  regional.coef.bkflw[[i]] <- sma.regional.bkfl
  regional.coef.depth[[i]] <- sma.regional.depth
  
  # Calculate RMSE- 
   # wetted width:
   rgn.regional.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = WettedWidth_m ~ EstWidth_ReglSMA)
   rgn.RMSE.smaR <- sqrt(mean(rgn.regional.fit.lm$residuals^2))
   
   rgn.national.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = WettedWidth_m ~ EstWidth_NatlSMA)
   rgn.RMSE.smaN <- sqrt(mean(rgn.national.fit.lm$residuals^2))
   
   rgn.raymond.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = WettedWidth_m ~ EstWidth_Raymond)
   rgn.RMSE.raymond <- sqrt(mean(rgn.raymond.fit.lm$residuals^2))
   
   # bankfull width:
   rgn.regional.bkfl.fit.lm <- dat_JoinNHD2_ls[[i]] %>% filter(BankfullWidth_m>0) %>% lm(formula = BankfullWidth_m ~ EstBkflWidth_ReglSMA)
   rgn.RMSE.smaR.bkfl <- sqrt(mean(rgn.regional.bkfl.fit.lm$residuals^2))
   
   rgn.national.bkfl.fit.lm <- dat_JoinNHD2_ls[[i]] %>% filter(BankfullWidth_m>0) %>% lm(formula = BankfullWidth_m ~ EstBkflWidth_NatlSMA)
   rgn.RMSE.smaN.bkfl <- sqrt(mean(rgn.national.bkfl.fit.lm$residuals^2))
   
   # thalweg depth:
   rgn.regional.depth.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = depth_m ~ EstDepth_ReglSMA)
   rgn.RMSE.smaR.depth <- sqrt(mean(rgn.regional.depth.fit.lm$residuals^2))
   
   rgn.national.depth.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = depth_m ~ EstDepth_NatlSMA)
   rgn.RMSE.smaN.depth <- sqrt(mean(rgn.national.depth.fit.lm$residuals^2))
   
   rgn.raymond.depth.fit.lm <- dat_JoinNHD2_ls[[i]] %>% lm(formula = depth_m ~ EstDepth_Raymond)
   rgn.RMSE.raymond.depth <- sqrt(mean(rgn.raymond.depth.fit.lm$residuals^2))
   
   rmse.summary <- data.frame(region = unique(dat_JoinNHD2_ls[[i]]$VPUID),
                              measurement = c(rep("wetted width",3),rep("bankfull width",3),rep("thalweg depth",3)),
                              model=rep(c("nationalSMA","Raymond2012","regionalSMA"),3),
                              model_RMSE = c(rgn.RMSE.smaN,rgn.RMSE.raymond,rgn.RMSE.smaR,
                                             rgn.RMSE.smaN.bkfl,NA,rgn.RMSE.smaR.bkfl,
                                             rgn.RMSE.smaN.depth,rgn.RMSE.raymond.depth,rgn.RMSE.smaR.depth))
   regional.RMSE[[i]] <- rmse.summary
   
  # Create plot:
   name <- case_when(dat_JoinNHD2_ls[[i]]$VPUID[1] == "01" ~ "Northeast",dat_JoinNHD2_ls[[i]]$VPUID[1] == "02" ~ "Mid-Atlantic",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "03N" ~ "South-Atlantic North",dat_JoinNHD2_ls[[i]]$VPUID[1] == "03S" ~ "South-Atlantic South",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "03W" ~ "South-Atlantic West",dat_JoinNHD2_ls[[i]]$VPUID[1] == "04" ~ "Great Lakes",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "05" ~ "Ohio",dat_JoinNHD2_ls[[i]]$VPUID[1] == "06" ~ "Tennessee",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "07" ~ "Upper Mississippi", dat_JoinNHD2_ls[[i]]$VPUID[1] == "08" ~ "Lower Mississippi",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "09" ~ "Souris-Red-Rainy",dat_JoinNHD2_ls[[i]]$VPUID[1] == "10L" ~ "Lower Missouri",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "10U" ~ "Upper Missouri",dat_JoinNHD2_ls[[i]]$VPUID[1] == "11" ~ "Ark-Red-White",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "12" ~ "Texas",dat_JoinNHD2_ls[[i]]$VPUID[1] == "13" ~ "Rio Grande",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "14" ~ "Upper Colorado",dat_JoinNHD2_ls[[i]]$VPUID[1] == "15" ~ "Lower Colorado",
                     dat_JoinNHD2_ls[[i]]$VPUID[1] == "16" ~ "Great Basin",dat_JoinNHD2_ls[[i]]$VPUID[1] == "17" ~ "Pacific Northwest",dat_JoinNHD2_ls[[i]]$VPUID[1] == "18" ~ "California")
   
   grob.width <- grid::grobTree(textGrob(paste("n= ",length(dat_JoinNHD2_ls[[i]]$SITE_ID),"\nslope = ",round(sma.regional$coef[[1]][2,1],3),"\nint = ",round(sma.regional$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.regional$r2),2),sep=""), x=0.08,  y=0.9, hjust=0,
                                         gp=gpar(col="black", fontsize=10)))
   
   p.width <- dat_JoinNHD2_ls[[i]] %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(WettedWidth_m)),shape=21,color="black",fill=colors[i]) + 
     geom_abline(slope = sma.regional$coef[[1]][2,1],intercept=sma.regional$coef[[1]][1,1],color=colors[i],size=1)+   
     ggtitle(label = name)+
     labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Wetted~width~(m)))+
     theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10)) + annotation_custom(grob.width)
   
   plot.list.width[[i]] <- p.width
   
   grob.bkfl <- grid::grobTree(textGrob(paste("n= ",length(dat_JoinNHD2_ls[[i]]$SITE_ID),"\nslope = ",round(sma.regional.bkfl$coef[[1]][2,1],3),"\nint = ",round(sma.regional.bkfl$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.regional.bkfl$r2),2),sep=""), x=0.08,  y=0.9, hjust=0,
                                        gp=gpar(col="black", fontsize=10)))
   
   p.bkfl <- dat_JoinNHD2_ls[[i]] %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(BankfullWidth_m)),shape=21,color="black",fill=colors[i]) + 
     geom_abline(slope = sma.regional.bkfl$coef[[1]][2,1],intercept=sma.regional.bkfl$coef[[1]][1,1],color=colors[i],size=1)+   
     ggtitle(label = name)+
     labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Bankfull~width~(m)))+
     theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10)) + annotation_custom(grob.bkfl)
   
   plot.list.bkfl[[i]] <- p.bkfl
   
   grob.depth <- grid::grobTree(textGrob(paste("n= ",length(dat_JoinNHD2_ls[[i]]$SITE_ID),"\nslope = ",round(sma.regional.depth$coef[[1]][2,1],3),"\nint = ",round(sma.regional.depth$coef[[1]][1,1],3),"\nr2 = ",round(as.numeric(sma.regional.depth$r2),2),sep=""), x=0.08,  y=0.9, hjust=0,
                                         gp=gpar(col="black", fontsize=10)))
   
   p.depth <- dat_JoinNHD2_ls[[i]] %>% ggplot() + geom_point(aes(x=log10(Calc_UpstrArea_km2),y=log10(depth_m)),shape=21,color="black",fill=colors[i]) + 
     geom_abline(slope = sma.regional.depth$coef[[1]][2,1],intercept=sma.regional.depth$coef[[1]][1,1],color=colors[i],size=1)+   
     ggtitle(label = name)+
     labs(x=expression(log[10]~Upstream~drainage~area~(km^2)),y=expression(log[10]~Thalweg~depth~(m)))+
     theme_cowplot() + theme(axis.text = element_text(size=9),axis.title=element_text(size=10)) + annotation_custom(grob.bkfl)
   
   plot.list.depth[[i]] <- p.depth
   
}   

# Print region-by-region table listing SMA RMSE:
print(regional.RMSE)

# Create multi-panel area-width plot across HUC02 regions:
plot.list.width[[1]] + plot.list.width[[2]] + plot.list.width[[3]] + plot.list.width[[4]] + plot.list.width[[5]] + plot.list.width[[6]] + 
  plot.list.width[[7]] + plot.list.width[[8]] + plot.list.width[[9]] + plot.list.width[[10]] + plot.list.width[[11]] + plot.list.width[[12]] + 
  plot.list.width[[13]] + plot.list.width[[14]] + plot.list.width[[15]] + plot.list.width[[16]] + plot.list.width[[17]] + plot.list.width[[18]] + 
  plot.list.width[[19]] + plot.list.width[[20]] + plot.list.width[[21]] + plot_layout(ncol=4)

# Create table of SMA model coefficients:

Table <- data.frame(HUC02 = names(dat_JoinNHD2_ls),HUC_check = NA,name = NA,slope = NA,
                    intercept = NA,slope_lowCI = NA,slope_upperCI = NA,int_lowCI = NA,int_upperCI = NA)

for(i in 1:length(Table$HUC02)){
  Table$HUC_check[i] <- as.character(dat_JoinNHD2_ls[[i]]$VPUID[1])
  
  Table$name[i] <- case_when(Table$HUC02[i] == "01" ~ "Northeast",Table$HUC02[i] == "02" ~ "Mid-Atlantic",
                             Table$HUC02[i] == "03N" ~ "South-Atlantic North",Table$HUC02[i] == "03S" ~ "South-Atlantic South",
                             Table$HUC02[i] == "03W" ~ "South-Atlantic West",Table$HUC02[i] == "04" ~ "Great Lakes",
                             Table$HUC02[i] == "05" ~ "Ohio",Table$HUC02[i] == "06" ~ "Tennessee",
                             Table$HUC02[i] == "07" ~ "Upper Mississippi", Table$HUC02[i] == "08" ~ "Lower Mississippi",
                             Table$HUC02[i] == "09" ~ "Souris-Red-Rainy",Table$HUC02[i] == "10L" ~ "Lower Missouri",
                             Table$HUC02[i] == "10U" ~ "Upper Missouri",Table$HUC02[i] == "11" ~ "Ark-Red-White",
                             Table$HUC02[i] == "12" ~ "Texas",Table$HUC02[i] == "13" ~ "Rio Grande",
                             Table$HUC02[i] == "14" ~ "Upper Colorado",Table$HUC02[i] == "15" ~ "Lower Colorado",
                             Table$HUC02[i] == "16" ~ "Great Basin",Table$HUC02[i] == "17" ~ "Pacific Northwest",Table$HUC02[i] == "18" ~ "California")
  
  Table$slope[i] <- round(regional.coef.width[[i]]$coef[[1]][2,1],3)
  Table$intercept[i] <- round(regional.coef.width[[i]]$coef[[1]][1,1],3)
  Table$slope_lowCI[i] <- round(regional.coef.width[[i]]$coef[[1]][2,2],3)
  Table$slope_upperCI[i] <- round(regional.coef.width[[i]]$coef[[1]][2,3],3)
  Table$int_lowCI[i] <- round(regional.coef.width[[i]]$coef[[1]][1,2],3)
  Table$int_upperCI[i] <- round(regional.coef.width[[i]]$coef[[1]][1,3],3)
}

if(all.equal(as.character(Table$HUC02),Table$HUC_check)){Table <- Table[,-2]}
print(Table)


  