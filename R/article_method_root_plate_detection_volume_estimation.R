
#Root plate detection and root plate volume estimation
#Methods described in the paper entitled
#"Root plates of uprooted trees – automatic detection and biotransport estimation 
#using LiDAR data and field mapping"
# paper link: https://www.sciencedirect.com/science/article/pii/S1569843224003467?via%3Dihub
#author: Janusz Godziek



#read required libraries
library(lidR)
library(sf)
library(terra)
library(stars)
library(dplyr)


#SET DIRECTORY--------------------------------------------------------------------------------
setwd('directory/to/your/folder')
getwd()

#local CRS should be defined by the user depending on the analysis site using EPSG code
EPSG <- 2180 #example EPSG for area located in Poland (ETRF2000-PL / CS92)


#READ DATA - point cloud and area of interest extent-----------------------------------------
pcl <- readLAS('path/to/your/point_cloud.laz')
crs(pcl) <- EPSG #set your local Coordinate Reference System (CRS) using EPSG code
#local CRS should be defined by the user depending on the analysis site


#read study area extent
aoi <- st_read('GPN_test_plot.shp')

#transform aoi to the same CRS as the PCL (if required)
aoi <- st_transform(aoi, EPSG)

aoi_b20 <- st_buffer(aoi, 20) #create buffer 20 m around aoi

pcl <- clip_roi(pcl, aoi_b20) #clip point cloud to the aoi buffer extent
plot(pcl)

#save the clipped point cloud
writeLAS(pcl, 'name_of_your_clipped_pcl.las')



#ROOT PLATE DETECTION - functions--------------------------------------------------------------

#process point cloud data
classify_pcl <- function(pcl, #character; path to point cloud
                         aoi_name, #character; name of the area of interest to be putted in the output_pth filenames
                         crs, #numeric; EPSG code of the coordinate reference system of the point cloud
                         output_pth) #character; path to save the output (processed point cloud)
{
  
  pcl = readLAS(pcl)
  crs(pcl) = crs #set crs
  
  #classify "low points" and "heigh points"
  pcl = classify_noise(pcl, sor())
  unique(pcl$Classification) #check unique classes occurence in the point cloud data
  
  #remove high noise points (18) and low noise points (7)
  pcl = pcl[(pcl$Classification != 18 & pcl$Classification != 7), ]
  
  
  #classify point cloud (own classification, reclassification)
  #parameters for the PMF algorithm (values from lidR documentation)
  ws = seq(3,12, 3)
  th = seq(0.1, 1.5, length.out = length(ws))
  
  #classify ground -> algorithm Progressive Morphological Filter
  pcl = classify_ground(pcl, algorithm = pmf(ws, th))
  
  #normalize height of point cloud (reduce the influence of the terrain)
  npcl = normalize_height(pcl, tin())
  
  #check unique class values
  unique(npcl$Classification)
  
  #classify low vegetation class (0-0.4 m)
  npcl$Classification[(npcl$Z > 0 & npcl$Z < 0.4)] = as.integer(3)
  
  #classify medium vegetation class (0.4-2 m)
  npcl$Classification[(npcl$Z >= 0.4 & npcl$Z < 2)] = as.integer(4)
  
  #classify high vegetation class (above 2 m)
  npcl$Classification[(npcl$Z >= 2)] = as.integer(5)
  
  #classify ground reflections (treat minus Z values as ground -> the lowest min value: -0.6)
  npcl$Classification[(npcl$Z <= 0)] = as.integer(2)
  
  #restore original heights  (unnormalize point cloud)
  pcl = unnormalize_height(npcl)
  
  writeLAS(pcl, paste0(output_pth, '/', aoi_name, '_pcl_rcl.laz'))
  
  return(pcl)
}


#compute elevation models from point cloud
#1) DTM - Digital Terrain Model
#2) DSM - Digital Surface Model
#3) VHM - Vegetation Height Model
#4) GLP - Ground and Low Points Model (last returns & points below 2 m above the ground)
#5) DM - differential model (DM = GLP - DTM)
compute_elev_models <- function(pcl, #LAS object; name of the point cloud LAS object (package lidR, readLAS())
                                aoi_name, #character; name of the area of interest to be putted in the output_pth filenames
                                output_pth) #character; path to save the output_pth (elevation models)
{
  
  #dtm algorithm -> tin()
  dtm = rasterize_terrain(pcl, algorithm = tin(), 
                          res= 0.25, use_class = c(2L, 9L))
  writeRaster(dtm, paste0(output_pth, '/', aoi_name, '_dtm025m.tif'))
  
  dsm = rasterize_canopy(pcl, algorithm = p2r(), 
                         res= 0.25, use_class = c(2L, 3L, 4L, 5L, 9L))
  writeRaster(dsm, paste0(output_pth, '/', aoi_name, '_dsm025m.tif'))
  
  vhm = dsm - dtm
  writeRaster(vhm, paste0(output_pth, '/', aoi_name, '_vhm025m.tif'))
  
  #smooth vhm
  kernel = matrix(1,3,3)
  vhm_smth = terra::focal(vhm, w = kernel, fun = median, na.rm = TRUE)
  writeRaster(vhm_smth, paste0(output_pth, '/', aoi_name, '_vhm025m_smth.tif'))
  
  
  #select only last returns
  pcl_lr = filter_last(pcl)
  unique(pcl_lr$Classification) #check point classes (should be only classes 2,3,4)
  
  #compute GLP model (GLP = ground and low points)
  glp = rasterize_terrain(pcl_lr, res = 0.25, algorithm = tin(), 
                          use_class = c(2L, 3L, 4L))
  
  #compute differential model
  dm = glp - dtm
  plot(dm)
  
  #write
  writeRaster(dm, paste0(output_pth, '/', aoi_name, '_dm025m.tif'))
}


#process differential model - extract polygons using contour lines
extract_plg_rp <- function(aoi_name, #character; name of the area of interest to be putted in the output filenames
                           dm, #character; path to Differential Model
                           vhm, #character; path to Vegetation Height Model
                           crs, #numeric; EPSG code of the coordinate reference system to be used during computations
                           cnt_h, #numeric; height levels of contour lines to be used to "cut off" root plates (to extract root plates polygons)
                           output_pth) #character; path to save the output shapefiles
{
  dm = rast(dm) #read differential model
  
  #1) contour lines computation
  
  #convert to stars object
  dm_stars = st_as_stars(dm)
  
  #calculate contour lines (package stars)
  cnt = st_contour(dm_stars, contour_lines = T, #produces contours as lines (in case of 'F' -> contour lines as multipolygons)
                   breaks = seq(0,5, 0.1)) #seq(start_value, end_value, interval)
  
  cnt = st_transform(cnt, crs) #set crs
  
  #rename column by name
  colnames(cnt)[1] = 'cnt_elev'
  
  #remove additional 0 digits (values like 1.200000000)
  #contour line heights now stored as characters!!!
  cnt$cnt_elev = as.character(cnt$cnt_elev)
  
  #set unique ID of contour lines
  cnt$ID = 1:nrow(cnt)
  
  #save
  st_write(cnt, paste0(output_pth, '/', aoi_name, '_dm025m_cnt01m.shp'))
  
  
  #EXTRACT contour lines by elevation
  cnt_h = cnt_h
  
  #extract contour lines of all given elevations
  cnt_sel = cnt[cnt$cnt_elev %in% cnt_h, ]
  
  #convert contours to polygons
  plg_sel = st_cast(cnt_sel, "POLYGON")
  
  #repair topology errors
  plg_sel = st_make_valid(plg_sel)
  
  #calculate area of polygons
  plg_sel$area = as.vector(st_area(plg_sel)) #area
  
  #remove the smallest artefacts
  plg_sel = plg_sel %>% filter(area > 0.1)
  
  #calculate perimeter of polygons
  plg_sel$perim = as.numeric(lwgeom::st_perimeter(plg_sel))
  
  #create convex hulls
  plg_sel_cvh = st_convex_hull(plg_sel)
  #compute area and perimeter of convex hulls
  plg_sel_cvh$cvh_perim = as.vector(lwgeom::st_perimeter(plg_sel_cvh)) #perimeter
  plg_sel_cvh$cvh_area = as.vector(st_area(plg_sel_cvh)) #area
  
  #add info about convex hull parameters to plg_sel
  plg_sel = cbind(plg_sel, st_drop_geometry(plg_sel_cvh[, c('cvh_perim', 'cvh_area')]))
  
  #compute polygon shape characteristics
  #https://github.com/pondrejk/PolygonComplexity 
  #polygon compactness index (CS) [pComp]
  plg_sel$pComp = as.vector(plg_sel$perim/(3.45*sqrt(plg_sel$area)))
  
  #compute zonal statistics of DM (mean)
  dm_mean = zonal(x = dm, z = vect(plg_sel), fun = 'mean')
  colnames(dm_mean)[1] = 'dm_mean'
  
  #merge columns for all data
  plg_sel = cbind(plg_sel, dm_mean)
  
  st_write(plg_sel, paste0(output_pth, '/', aoi_name, '_plg_sel_ar01_stats.shp'))
  return(plg_sel)
}

#Filter polygons to extract root plates
filter_rp <- function(plg, #root plate polygons (output of the function 'extract_plg_rp)
                      aoi_name, #character; name of the area of interest to be putted in the output filenames
                      ar_min, #border value of area (area greater than ar_min)
                      pc_max, #border value of polygon compactness index (polygon compactness index lower than pc_max)
                      output_pth) #character; path to save the output shapefiles
{
  #filtering
  rp = plg %>% filter (area < 5 & dm_mean > 0.5)
  rp = rp %>% filter (area > ar_min & pComp < pc_max)
  
  #dissolve root plate polygons
  rp = st_union(rp$geometry, by_feature = FALSE) #union
  rp = st_cast(rp, 'POLYGON')
  rp = st_sf(rp) #convert geometry set to simple feature collection (sf)
  
  st_write(rp, paste0(output_pth, '/', aoi_name, '_rp.shp'))
}


#ROOT PLATE VOLUME ESTIMATION - functions----------------------------------------------

#delineate the boundary of root plates
extract_rp_boundary <- function(rp, #sfc object, polygon; root plates - polygons being a result of 'filter_rp' function
                                dm) #terra SpatRaster; raster of the Differential Model (DM)
{
  
  #read Differential Model (DM) data
  dm = dm
  
  #read root plates data (root plates polygons, removed artefacts)
  rp = rp
  
  #delineate root plate boundaries on the basic of extracted closed contour lines and DM values
  #create 1m buffer around selected root plates
  rp_buf = st_buffer(rp, 1)
  
  #reclassify dm (0.1 value is the best threshold)
  dm_rcl = dm >= 0.1 #(create raster with values : 1 - above or equal to 0.1, 0 - less than 0.1)
  dm_rcl[dm_rcl < 1] = NA #set 0 values to NA
  
  dm_plg = as.polygons(dm_rcl) #convert to polygons (class SpatVector)
  dm_plg = st_as_sf(dm_plg) #convert to sf
  
  #smooth the boundaries
  dm_plg = rmapshaper::ms_simplify(dm_plg, keep = 0.2, keep_shapes = T) #simplify polygons
  dm_plg = smoothr::smooth(dm_plg, "chaikin") #smooth polygons
  
  #repair topology
  rp_buf = st_make_valid(rp_buf)
  dm_plg = st_make_valid(dm_plg)
  
  #bnd = boundary
  #clip dm_plg with the buffer of root plates
  rp_bnd = st_intersection(rp_buf, dm_plg)
  
  #convert from multipolygons to polygons
  rp_bnd = st_collection_extract(rp_bnd, 'POLYGON')
  #multipolygons to polygons
  rp_bnd = st_cast(st_cast(rp_bnd, 'MULTIPOLYGON'), 'POLYGON')
  
  #remove additional boundary polygons (1 root plate = 1 boundary)
  rp_centroi = st_centroid(rp) #compute centroids of root plates
  rp_centroi_b = st_buffer(rp_centroi, 0.5) #compute 0.5-m buffer around centroids of root plates
  
  #check if polygon intersects with root plate centroid 0.5-m buffer
  rp_bnd$intersec = st_intersects(rp_bnd, rp_centroi_b) %>% lengths > 0 #obtaining result as boolean variable
  rp_bnd = rp_bnd[rp_bnd$intersec == T, ] #select polygons that intersect with root plate polygons
  rp_bnd = select(rp_bnd, -c(intersec)) #remove 'intersec' column
  
  
  #handle with the overlap
  #split the overlapping parts of the polygons into two neighbouring polygons
  
  st_no_overlap = function(polygons) {
    
    centroids <- st_centroid(polygons)
    
    # Voronoi tesselation
    voronoi <- 
      centroids %>% 
      st_geometry() %>%
      st_union() %>%
      st_voronoi() %>%
      st_collection_extract()
    
    # Put them back in their original order
    voronoi <-
      voronoi[unlist(st_intersects(centroids,voronoi))]
    
    # Keep the attributes
    result <- centroids
    
    # Intersect voronoi zones with buffer zones
    st_geometry(result) <-
      mapply(function(x,y) st_intersection(x,y),
             #st_buffer(st_geometry(centroids),dist), 
             polygons$geometry,
             voronoi,
             SIMPLIFY=FALSE) %>%
      st_sfc(crs=st_crs(centroids))
    
    result
  }
  
  rp_bnd = st_no_overlap(rp_bnd)
  
  #set root plate ID (important for volume computations)
  rp_bnd$rp_ID <- 1:nrow(rp_bnd)
  
  return(rp_bnd)
}


#estimate the root plate volume
rp_volume <- function(rp_bnd, #sfc object; boundaries of root plates (result from function 'extract_rp_boundary')
                      cnt, #sfc object; contour lines of Differential Model (DM)
                      method) #character; one of two available methods (cnt, zstat)
{
  if (method == 'cnt'){
    #Approach 1: estimate volume basing on contour lines
    cnt_v_est = st_cast(cnt, 'POLYGON') #convert contour lines to polygons
    cnt_v_est = st_make_valid(cnt_v_est) #repair topology
    
    
    #for each root plate select polygons refering to contour lines
    l = vector("list", length = nrow(rp_bnd)) #list to append loop results
    
    for (i in 1:nrow(rp_bnd)){
      #clip contour-polygons with root plate boundaries
      clip = st_intersection(x = cnt_v_est, #feature to be clipped
                             y = rp_bnd[i, ]) #clipping feature
      
      #select only polygons intersecting with the selected root plate polygon
      clip$contains_rp = st_intersects(clip, rp) %>% lengths > 0 #obtaining result as boolean variable
      clip = clip[clip$contains_rp == T, ]
      
      l[[i]] = clip #append to list
    }
    
    cnt_v_est_rp = st_sf(do.call(rbind, l)) #convert to sf
    
    #calculate area of each polygon
    cnt_v_est_rp$cnt_area = as.numeric(st_area(cnt_v_est_rp))
    
    #for each root plate sum areas of contours
    cnt_v_est_rp = cnt_v_est_rp%>% group_by(rp_ID) %>%
      summarize(rp_v_area_sum = sum(cnt_area, na.rm = TRUE))
    
    #calculate volume of each root plate - multiply by contour lines interval
    cnt_v_est_rp$rp_volume = cnt_v_est_rp$rp_v_area_sum*0.1
    
    #return volume data as vector
    vol_cnt = cnt_v_est_rp$rp_volume
    
    return(vol_cnt)
  }
  
  if (method == 'zstat'){
    #Approach 2: estimate volume basing on zonal statistics (sum)
    rp_bnd_sv = vect(rp_bnd) #convert to SpatVector
    #perform zonal statistics (sum all values of DM within the rp_bnd polygons)
    dm_sum = terra::extract(dm, rp_bnd_sv, fun = sum)
    #calculate the volume (multiply by squared spatial resolution of DM, i.e 0.25*0.25)
    dm_sum$rp_V_zst = dm_sum$Z*0.25*0.25
    
    #return output as vector
    vol_zstat = dm_sum$rp_V_zst
    
    return(vol_zstat)
  }
}




#APPLY FUNCTIONS---------------------------------------------------------

#root plate detection

setwd('set/your/working/directory')
getwd() #check your working directory

#insert the name of the area of interest (to be putted in the filenames)
aoi_name <- 'name_of_your_aoi'

#process point cloud
pcl <- classify_pcl(pcl = 'path/to/your/point_cloud.laz',
                    crs = EPSG,
                    aoi_name = aoi_name,
                    output_pth = getwd())

#compute elevation models
compute_elev_models(pcl = pcl,
                    aoi_name = aoi_name,
                    output_pth = getwd())


#extract polygons from differential model
plg <- extract_plg_rp(aoi_name = aoi_name,
                      dm = paste0(getwd(), '/', aoi_name, '_dm025m.tif'),
                      vhm = paste0(getwd(), '/', aoi_name, '_vhm025m.tif'),
                      crs = EPSG,
                      cnt_h = c(0.5, 1, 1.5),
                      output_pth = getwd())


#filter polygons to extract root plates
rp <- filter_rp(plg = plg, 
                aoi_name = aoi_name,
                ar_min = 0.9,
                pc_max = 2.2,
                output_pth = getwd())



#root plate volume estimation


#read Differential Model (DM) data
dm <- rast(paste0(getwd(), '/', aoi_name, '_dm025m.tif'))

#read contour lines data
cnt <- st_read(paste0(getwd(), '/', aoi_name, '_dm025m_cnt01m.shp'))

#read root plates polygons
rp <- st_read(paste0(getwd(), '/', aoi_name, '_rp.shp'))

#extract boundaries of root plates
rp_bnd <- extract_rp_boundary(rp, dm)

#calculate volume of root plates
rp_bnd$vol_cnt <- rp_volume(rp_bnd, cnt, method = 'cnt') #contour method
rp_bnd$vol_zstat <- rp_volume(rp_bnd, cnt, method = 'zstat') #zonal statistics method

#save results
st_write(rp_bnd, 'rp_vol_est_results.shp')






