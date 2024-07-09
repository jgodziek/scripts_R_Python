# R SCRIPT

# Pit-mound topography detection
#Method described in the paper entitled
#"Indicators of wind-driven forest disturbances – pit–mound topography, its automatic detection and significance"
# paper link: https://www.sciencedirect.com/science/article/abs/pii/S0341816222007433?via%3Dihub
#author: Janusz Godziek




library(rgdal)
library(raster)
library(terra)
library(lidR)
library(gstat)
library(Rsagacmd)
library(spatialEco)
library(spatstat)
library(sf)
library(units)
library(dplyr)
library(rgugik)
library(smoothr)
library(rmapshaper)
library(exactextractr)
library(rgeos)
library(hablar)
library(ggplot2)
library(sp)
library(stars)
library(stringr)
library(writexl)
library(downloader)


options(scipen = 999) #turn of scientific notatnion

#basic parameters:
site_name = "Your_site_name"
epsg = 2180
dem_res = 0.1
cnt_interv = 0.1
min_cnt_len = 2.5
max_cnt_len = 20



#paths
wd_pth = "Your_path"
CM_pth = "Your_path"

#study area shp path
study_area = "path/to/your/study_area.shp"


##Functions
####################################
##Functions


extract_PCl <- function(study_area, site_name, pth, buffer_size){
        #read shp
        study_area <- read_sf(study_area)
        study_area_buffer <- st_buffer(study_area, buffer_size)
        
        #check available gugik data for study area
        gugik_data <- DEM_request(study_area_buffer)
        #select only point cloud data in LAS format
        gugik_PointCloud <- gugik_data[(gugik_data$product == "PointCloud")&(gugik_data$format=="LAS"), ]
        gugik_PCl2 <- gugik_PointCloud[(gugik_PointCloud$year < 2016), ] #select only data from first LiDAR survery (before 2016)
        
        #download tile
        dir.create(file.path(pth, "las"))
        url <- gugik_PCl2$URL #extract url data
        for (u in url){
                downloader::download(url = u, destfile = paste0(pth, '/las/', str_sub(u, -33, -1 )), mode = 'wb') #download las from url to directory
        }
        
        #merge LAS tiles from folder
        las_catalog <- readLAScatalog(file.path(pth, "las")) #read LAS catalog
        las_merged <- rbind(las_catalog)
        
        #clip las to studyPlot
        las_studyPlot <- clip_roi(las_merged, study_area_buffer)
        
        #check corectness of data
        las_check(las_studyPlot)
        
        #filter duplicates
        las_studyPlot <- filter_duplicates(las_studyPlot)
        
        #saving clipped point cloud to laz file
        writeLAS(las_studyPlot, file.path(pth, "las", gsub(" ", "", paste(site_name, "_KRON86.laz")))
        )
        
}

compute_DEM <- function(site_name, pth, point_cloud, study_area, dem_res, epsg){
        #create DEM from LiDAR Point Cloud
        #read las data
        las_studyPlot <- readLAS(point_cloud)
        
        #DEM computation, k- nearest neighbour interpolation
        DEM_studyPlot <- grid_terrain(las_studyPlot, algorithm = knnidw(k = 10L, p = 2), 
                                      res= dem_res, use_class = c(2L, 9L))
        #set crs
        raster::crs(DEM_studyPlot) <- paste0("EPSG:", c(epsg))
        
        #read input shp data
        studyPlot <- st_read(study_area)
        
        #clip DEM to study area boundary
        DEM_studyPlot2 <- mask(DEM_studyPlot, studyPlot)
        
        dem_resx <- gsub('\\.', '', as.character(dem_res)) #remove dot from DEM resolution
        
        #save DEM files
        writeRaster(DEM_studyPlot, file.path(pth, "/", site_name, "_buffer_DEM", dem_resx,"m.tif", fsep = ""))
        writeRaster(DEM_studyPlot2, file.path(pth, "/", site_name, "_DEM", dem_resx, "m.tif", fsep = ""))
}

DEM_cnt <- function(DEM, pth, site_name, interv) {
        saga = saga_gis()
        cnt_01 <- saga$shapes_grid$contour_lines_from_grid(DEM, zstep = interv)
        #convert MULTILINESTRING TO LINESTRING
        cnt_01 <- st_cast(cnt_01[["contour"]], "LINESTRING")
        #calculate length of each line (add "length" column)
        cnt_01$length <- st_length(cnt_01)
        #convert length column to vector, remove unit '[m]' from column
        cnt_01$length <- as.vector(cnt_01$length)
        #writing sf object to shapefile
        st_write(cnt_01, dsn = file.path(pth, gsub(" ", "", paste(site_name, "_cnt01m.shp"))))
}


#CONTOUR METHOD ALGORITHM

#two versions of function due to data processing issue:
#dissolve polygons: sometimes works st_cast, sometimes st_collection_extract

#1st version -> dissolve polygons with st_collection_extract
#output: polygons_PM, polygons_pm_pairs,
#raster_PM[-1,0,1], raster_pm_pairs[-1,0,1], max_altitude, min_altitude;  
#minimal distance between closed contour and border: 0.5 m
detect_PitMound_cnt <- function(site_name,
                                pth,
                                contours,
                                study_area,
                                DEM,
                                max_cnt_length,
                                min_cnt_length,
                                epsg) {
        #read contour lines shp
        cnt_StudyPlot <- st_read(contours) #
        
        #select contours with given length interval
        cnt_lenInterv <- cnt_StudyPlot[(cnt_StudyPlot$length>min_cnt_length) & (cnt_StudyPlot$length<max_cnt_length), ]
        #transform sf data to given EPSG
        cnt_lenInterv <- st_transform(cnt_lenInterv, epsg)
        
        #read boundary of study plot
        study_area <- st_read(study_area)
        plot(study_area)
        
        #convert sf polygon geometry to sf linestring geometry
        study_area <- st_cast(study_area, "LINESTRING")
        study_area <- st_transform(study_area, epsg)
        
        #add "border" column -> 1 if feature is within a distance of 0.5 m form the border
        cnt_lenInterv$border <- as.numeric(st_is_within_distance(cnt_lenInterv, study_area, 0.5))
        
        #select all rows with 'nan' values in the 'border' column
        cnt_lenInterv <- cnt_lenInterv[is.na(cnt_lenInterv$border), ]
        
        #change contours to polygons
        polygons1 <- st_cast(cnt_lenInterv, "POLYGON") 
        
        #repair topology errors
        polygons1 <- st_make_valid(polygons1)
        
        #dissolve polygons
        polygons2 <- st_union(polygons1$geometry, by_feature = FALSE) #union
        polygons2 <- st_collection_extract(polygons2, type = "POLYGON") #dissolve by extracting polygons from a list of geometries
        
        #simplify polygons
        polygons2 <- ms_simplify(polygons2, keep = 0.2, keep_shapes = T)
        #smooth polygons
        polygons2 <- smooth(polygons2, "chaikin")
        
        polygons2 = polygons2[!st_is_empty(polygons2), drop=FALSE] #remove empty geometries
        #convert sfc polygons to SpatialObject
        polygons2 <- as(polygons2, "Spatial")
        df <- data.frame(ID=character(), stringsAsFactors=FALSE ) #create empty df
        #extract IDs of polygons from SpatialPolygon to df
        for (i in polygons2@polygons) { 
                df <- rbind(df, data.frame(ID=i@ID, stringsAsFactors=FALSE))  }
        df[,1] <- sub("ID", "", df[,1]) #remove 'ID' prefix
        #convert SpatialPolygons to SpatialPolygonsDataFrame
        poly_PitMound <-SpatialPolygonsDataFrame(polygons2, data=df, match.ID = F)
        
        #
        DEM <- raster(DEM)
        #set crs for raster
        crs(DEM) = paste0("EPSG:", c(epsg))
        
        #add raster cell index -> takes quite long time
        points_SVcells <- raster::extract(DEM, polygons2, cellnumbers = TRUE)#values plus information about the raster cell index
        
        #MAX
        #extract max values
        points_SVmax <- t(sapply(points_SVcells, function(i) i[which.max(i[,2]), ] ))
        #add coords, merge vectors
        points_SVmax_coords <- data.frame(points_SVmax, raster::xyFromCell(DEM, points_SVmax[,1]))
        colnames(points_SVmax_coords)[2:4] <- c("maxElv", "XmaxElv", "YmaxElv") #change column names
        
        #add row index
        points_SVmax_coords <- points_SVmax_coords %>%
                tibble::rownames_to_column('ID')
        points_SVmax_coords <- points_SVmax_coords %>%
                convert(dbl(ID))
        
        
        #MIN
        #extract min values
        points_SVmin <- t(sapply(points_SVcells, function(i) i[which.min(i[,2]), ] ))
        points_SVmin_coords <- data.frame(points_SVmin, raster::xyFromCell(DEM, points_SVmin[,1]))
        colnames(points_SVmin_coords)[2:4] <- c("minElv", "XminElv", "YminElv")
        #add row index
        points_SVmin_coords <- points_SVmin_coords %>%
                tibble::rownames_to_column('ID')
        points_SVmin_coords <- points_SVmin_coords %>%
                convert(dbl(ID))
        
        
        poly_PitMound <- st_as_sf(poly_PitMound) #convert to sfc
        poly_PitMound$ID <- as.numeric(poly_PitMound$ID) #convert ID to numeric
        
        #merge SPDF of pit and mound and data of max and min elevation value
        poly_PitMound <- dplyr::left_join(poly_PitMound, points_SVmax_coords, by = c("ID" = "ID"))
        poly_PitMound <- dplyr::left_join(poly_PitMound, points_SVmin_coords, by = c("ID" = "ID"))
        
        
        poly_PitMound <- as_Spatial(poly_PitMound)
        #calculate values basing on existing columns (difference of elevation within polygons)
        poly_PitMound@data <- transform(poly_PitMound@data, ElvDiff=maxElv-minElv)
        poly_PitMound$PPerim <- polyPerimeter(poly_PitMound) #calculate perimeter of each poly
        poly_PitMound$PArea = sapply(slot(poly_PitMound, "polygons"), slot, "area") #extract area of each poly
        
        #convert max&min points to sfc objects
        points_SVmax_coords_sf <- st_as_sf(points_SVmax_coords, coords = c('XmaxElv','YmaxElv'), crs=epsg)
        points_SVmin_coords_sf <- st_as_sf(points_SVmin_coords, coords = c('XminElv','YminElv'), crs=epsg)
        
        poly_PitMound <- st_as_sf(poly_PitMound) #convert PM polygons back to sfc
        line_PM <- st_cast(poly_PitMound, 'LINESTRING') #convert poly_PM to line
        
        #calculate distance from border of polygon to maxElv and minElv
        line_PM$D_mxE_B <- as.numeric(st_distance(points_SVmax_coords_sf, line_PM, which = 'Euclidean', by_element = T))
        line_PM$D_mnE_B <- as.numeric(st_distance(points_SVmin_coords_sf, line_PM, which = 'Euclidean', by_element = T))
        
        poly_PitMound <- st_cast(line_PM, 'POLYGON') #convert sfc linestring to polygon
        
        #calculate difference between 'D_mxE_B' and 'D_mnE_B'
        poly_PitMound$Diff1 <- abs(poly_PitMound$D_mxE_B - poly_PitMound$D_mnE_B)
        
        
        #divide data into 'Pit', 'Mound' and 'Unclassified' basing on the computed distances
        poly_PitMound$Class <- ifelse(poly_PitMound$Diff1 < 0.1, 'Unclassified', 
                                      ifelse(poly_PitMound$D_mnE_B > poly_PitMound$D_mxE_B,
                                             'PIT', 'MOUND'))
        
        #add column with numeric values of classes
        poly_PitMound <- poly_PitMound %>% 
                mutate(Class_num = case_when(Class == "PIT" ~ -1,
                                             Class == "MOUND" ~ 1,
                                             Class == "Unclassified" ~ 0))
        
        #extract bounding box
        bound_box <- st_bbox(study_area)
        
        #rasterize data
        poly_PitMound_rast <- st_rasterize(poly_PitMound[, "Class_num"], driver = "GTiff", 
                                           dx = 0.5, dy = 0.5, crs = c(epsg), #size of cells in x and y directions, set crs
                                           xlim = c(as.numeric(bound_box[1]), as.numeric(bound_box[3])),
                                           ylim = c(as.numeric(bound_box[2]), as.numeric(bound_box[4]))) #set extent of x axis and y axis !extremely important!!!
        
        poly_PitMound_rast[is.na(poly_PitMound_rast[])] <- 0 #replace Na values
        
        pol_p <- poly_PitMound[poly_PitMound$Class == 'PIT', ] #select only 'pit' polygons
        pol_m <- poly_PitMound[poly_PitMound$Class == 'MOUND', ] #select only 'mound' polygons
        
        #extract distance to the nearest feature
        #Pits
        pol_p$nearest_mound <- st_nearest_feature(pol_p, pol_m)
        pol_p$nearest_dist <- as.numeric(st_distance(pol_p, pol_m[pol_p$nearest_mound, ], by_element = T))
        #Mounds
        pol_m$nearest_pit <- st_nearest_feature(pol_m, pol_p)
        pol_m$nearest_dist <- as.numeric(st_distance(pol_m, pol_p[pol_m$nearest_pit, ], by_element = T))
        
        #extract pairs of Pits and Mounds
        pit_pair <- pol_p[(pol_p$nearest_dist <1.5), ] #select objects
        mound_pair <- pol_m[(pol_m$nearest_dist <1.5), ] #select objects
        
        #sort by nearest mound (asc), nearest dist (asc), area(desc)
        pit_pair <- pit_pair[order(pit_pair$nearest_mound, 
                                   pit_pair$nearest_dist, -pit_pair$PArea),]
        mound_pair <- mound_pair[order(mound_pair$nearest_pit, 
                                       mound_pair$nearest_dist, -mound_pair$PArea),]
        
        #add column with info: is nearest mound/pit id duplicated? [T/F]
        pit_pair$duplicated <- duplicated(pit_pair$nearest_mound)
        mound_pair$duplicated <- duplicated(mound_pair$nearest_pit)
        
        #select only not duplicated values
        pit_pair <- pit_pair[pit_pair$duplicated == 'FALSE', ]
        mound_pair <- mound_pair[mound_pair$duplicated == 'FALSE', ]
        
        #add columns to enable concatenation of pit data and mound data
        pit_pair$nearest_pit <- 0
        mound_pair$nearest_mound <- 0
        
        #concatenate sf objects
        pit_mound_pair <- rbind(pit_pair, mound_pair)
        
        #rasterize data
        pit_mound_pair_rast <- st_rasterize(pit_mound_pair[, "Class_num"], driver = "GTiff", 
                                            dx = 0.5, dy = 0.5, crs = c(epsg), #size of cells in x and y directions, set crs
                                            xlim = c(as.numeric(bound_box[1]), as.numeric(bound_box[3])),
                                            ylim = c(as.numeric(bound_box[2]), as.numeric(bound_box[4]))) #set extent of x axis and y axis !extremely important!!!
        
        pit_mound_pair_rast[is.na(pit_mound_pair_rast[])] <- 0 #replace Na values
        
        #save data to raster file (tif)
        write_stars(pit_mound_pair_rast, driver = "GTiff",
                    file.path(pth, '/', site_name, "_pm_pairs.tif", fsep = ''))
        write_stars(poly_PitMound_rast, driver = "GTiff",
                    file.path(pth, '/', site_name, '_rast_PM.tif', fsep = ''))
        
        #WRITE DATA TO SHP
        st_write(poly_PitMound, file.path(pth, '/', site_name, '_poly_PM.shp', fsep = ''))
        st_write(pit_mound_pair, file.path(pth, '/', site_name, '_pm_pairs.shp', fsep = ''))
        st_write(points_SVmax_coords_sf, file.path(pth, '/', site_name, '_pointsMax.shp', fsep = ''))
        st_write(points_SVmin_coords_sf, file.path(pth, '/', site_name, '_pointsMin.shp', fsep = ''))
}


#2nd version -> dissolve polygons with st_cast
#output: polygons_PM, polygons_pm_pairs,
#raster_PM[-1,0,1], raster_pm_pairs[-1,0,1], max_altitude, min_altitude; 
#minimal distance between closed contour and border: 0.5 m
detect_PitMound_cnt2 <- function(site_name,
                                 pth,
                                 contours,
                                 study_area,
                                 DEM,
                                 max_cnt_length,
                                 min_cnt_length,
                                 epsg) {
        #read contour lines shp
        cnt_StudyPlot <- st_read(contours) #
        
        #select contours with given length interval
        cnt_lenInterv <- cnt_StudyPlot[(cnt_StudyPlot$length>min_cnt_length) & (cnt_StudyPlot$length<max_cnt_length), ]
        #transform sf data to given EPSG
        cnt_lenInterv <- st_transform(cnt_lenInterv, epsg)
        
        #read boundary of study plot
        study_area <- st_read(study_area)
        plot(study_area)
        
        #convert sf polygon geometry to sf linestring geometry
        study_area <- st_cast(study_area, "LINESTRING")
        study_area <- st_transform(study_area, epsg)
        
        #add "border" column -> 1 if feature is within a distance of 0.5 m form the border
        cnt_lenInterv$border <- as.numeric(st_is_within_distance(cnt_lenInterv, study_area, 0.5))
        
        #select all rows with 'nan' values in the 'border' column
        cnt_lenInterv <- cnt_lenInterv[is.na(cnt_lenInterv$border), ]
        
        #change contours to polygons
        polygons1 <- st_cast(cnt_lenInterv, "POLYGON") 
        
        #repair topology errors
        polygons1 <- st_make_valid(polygons1)
        
        #dissolve polygons
        polygons2 <- st_union(polygons1$geometry, by_feature = FALSE) #union
        polygons2 <- st_cast(polygons2, 'POLYGON') #dissolve (convert MULTIPOLYGON to POLYGON)
        
        #simplify polygons
        polygons2 <- ms_simplify(polygons2, keep = 0.2, keep_shapes = T)
        #smooth polygons
        polygons2 <- smooth(polygons2, "chaikin")
        
        polygons2 = polygons2[!st_is_empty(polygons2), drop=FALSE] #remove empty geometries
        #convert sfc polygons to SpatialObject
        polygons2 <- as(polygons2, "Spatial")
        df <- data.frame(ID=character(), stringsAsFactors=FALSE ) #create empty df
        #extract IDs of polygons from SpatialPolygon to df
        for (i in polygons2@polygons) { 
                df <- rbind(df, data.frame(ID=i@ID, stringsAsFactors=FALSE))  }
        df[,1] <- sub("ID", "", df[,1]) #remove 'ID' prefix
        #convert SpatialPolygons to SpatialPolygonsDataFrame
        poly_PitMound <-SpatialPolygonsDataFrame(polygons2, data=df, match.ID = F)
        
        #
        DEM <- raster(DEM)
        #set crs for raster
        crs(DEM) = paste0("EPSG:", c(epsg))
        
        #add raster cell index -> takes quite long time
        points_SVcells <- raster::extract(DEM, polygons2, cellnumbers = TRUE)#values plus information about the raster cell index
        
        #MAX
        #extract max values
        points_SVmax <- t(sapply(points_SVcells, function(i) i[which.max(i[,2]), ] ))
        #add coords, merge vectors
        points_SVmax_coords <- data.frame(points_SVmax, raster::xyFromCell(DEM, points_SVmax[,1]))
        colnames(points_SVmax_coords)[2:4] <- c("maxElv", "XmaxElv", "YmaxElv") #change column names
        
        #add row index
        points_SVmax_coords <- points_SVmax_coords %>%
                tibble::rownames_to_column('ID')
        points_SVmax_coords <- points_SVmax_coords %>%
                convert(dbl(ID))
        
        
        #MIN
        #extract min values
        points_SVmin <- t(sapply(points_SVcells, function(i) i[which.min(i[,2]), ] ))
        points_SVmin_coords <- data.frame(points_SVmin, raster::xyFromCell(DEM, points_SVmin[,1]))
        colnames(points_SVmin_coords)[2:4] <- c("minElv", "XminElv", "YminElv")
        #add row index
        points_SVmin_coords <- points_SVmin_coords %>%
                tibble::rownames_to_column('ID')
        points_SVmin_coords <- points_SVmin_coords %>%
                convert(dbl(ID))
        
        
        poly_PitMound <- st_as_sf(poly_PitMound) #convert to sfc
        poly_PitMound$ID <- as.numeric(poly_PitMound$ID) #convert ID to numeric
        
        #merge SPDF of pit and mound and data of max and min elevation value
        poly_PitMound <- dplyr::left_join(poly_PitMound, points_SVmax_coords, by = c("ID" = "ID"))
        poly_PitMound <- dplyr::left_join(poly_PitMound, points_SVmin_coords, by = c("ID" = "ID"))
        
        
        poly_PitMound <- as_Spatial(poly_PitMound)
        #calculate values basing on existing columns (difference of elevation within polygons)
        poly_PitMound@data <- transform(poly_PitMound@data, ElvDiff=maxElv-minElv)
        poly_PitMound$PPerim <- polyPerimeter(poly_PitMound) #calculate perimeter of each poly
        poly_PitMound$PArea = sapply(slot(poly_PitMound, "polygons"), slot, "area") #extract area of each poly
        
        #convert max&min points to sfc objects
        points_SVmax_coords_sf <- st_as_sf(points_SVmax_coords, coords = c('XmaxElv','YmaxElv'), crs=epsg)
        points_SVmin_coords_sf <- st_as_sf(points_SVmin_coords, coords = c('XminElv','YminElv'), crs=epsg)
        
        poly_PitMound <- st_as_sf(poly_PitMound) #convert PM polygons back to sfc
        line_PM <- st_cast(poly_PitMound, 'LINESTRING') #convert poly_PM to line
        
        #calculate distance from border of polygon to maxElv and minElv
        line_PM$D_mxE_B <- as.numeric(st_distance(points_SVmax_coords_sf, line_PM, which = 'Euclidean', by_element = T))
        line_PM$D_mnE_B <- as.numeric(st_distance(points_SVmin_coords_sf, line_PM, which = 'Euclidean', by_element = T))
        
        poly_PitMound <- st_cast(line_PM, 'POLYGON') #convert sfc linestring to polygon
        
        #calculate difference between 'D_mxE_B' and 'D_mnE_B'
        poly_PitMound$Diff1 <- abs(poly_PitMound$D_mxE_B - poly_PitMound$D_mnE_B)
        
        
        #divide data into 'Pit', 'Mound' and 'Unclassified' basing on the computed distances
        poly_PitMound$Class <- ifelse(poly_PitMound$Diff1 < 0.1, 'Unclassified', 
                                      ifelse(poly_PitMound$D_mnE_B > poly_PitMound$D_mxE_B,
                                             'PIT', 'MOUND'))
        
        #add column with numeric values of classes
        poly_PitMound <- poly_PitMound %>% 
                mutate(Class_num = case_when(Class == "PIT" ~ -1,
                                             Class == "MOUND" ~ 1,
                                             Class == "Unclassified" ~ 0))
        
        #extract bounding box
        bound_box <- st_bbox(study_area)
        
        #rasterize data
        poly_PitMound_rast <- st_rasterize(poly_PitMound[, "Class_num"], driver = "GTiff", 
                                           dx = 0.5, dy = 0.5, crs = c(epsg), #size of cells in x and y directions, set crs
                                           xlim = c(as.numeric(bound_box[1]), as.numeric(bound_box[3])),
                                           ylim = c(as.numeric(bound_box[2]), as.numeric(bound_box[4]))) #set extent of x axis and y axis !extremely important!!!
        
        poly_PitMound_rast[is.na(poly_PitMound_rast[])] <- 0 #replace Na values
        
        pol_p <- poly_PitMound[poly_PitMound$Class == 'PIT', ] #select only 'pit' polygons
        pol_m <- poly_PitMound[poly_PitMound$Class == 'MOUND', ] #select only 'mound' polygons
        
        #extract distance to the nearest feature
        #Pits
        pol_p$nearest_mound <- st_nearest_feature(pol_p, pol_m)
        pol_p$nearest_dist <- as.numeric(st_distance(pol_p, pol_m[pol_p$nearest_mound, ], by_element = T))
        #Mounds
        pol_m$nearest_pit <- st_nearest_feature(pol_m, pol_p)
        pol_m$nearest_dist <- as.numeric(st_distance(pol_m, pol_p[pol_m$nearest_pit, ], by_element = T))
        
        #extract pairs of Pits and Mounds
        pit_pair <- pol_p[(pol_p$nearest_dist <1.5), ] #select objects
        mound_pair <- pol_m[(pol_m$nearest_dist <1.5), ] #select objects
        
        #sort by nearest mound (asc), nearest dist (asc), area(desc)
        pit_pair <- pit_pair[order(pit_pair$nearest_mound, 
                                   pit_pair$nearest_dist, -pit_pair$PArea),]
        mound_pair <- mound_pair[order(mound_pair$nearest_pit, 
                                       mound_pair$nearest_dist, -mound_pair$PArea),]
        
        #add column with info: is nearest mound/pit id duplicated? [T/F]
        pit_pair$duplicated <- duplicated(pit_pair$nearest_mound)
        mound_pair$duplicated <- duplicated(mound_pair$nearest_pit)
        
        #select only not duplicated values
        pit_pair <- pit_pair[pit_pair$duplicated == 'FALSE', ]
        mound_pair <- mound_pair[mound_pair$duplicated == 'FALSE', ]
        
        #add columns to enable concatenation of pit data and mound data
        pit_pair$nearest_pit <- 0
        mound_pair$nearest_mound <- 0
        
        #concatenate sf objects
        pit_mound_pair <- rbind(pit_pair, mound_pair)
        
        #rasterize data
        pit_mound_pair_rast <- st_rasterize(pit_mound_pair[, "Class_num"], driver = "GTiff", 
                                            dx = 0.5, dy = 0.5, crs = c(epsg), #size of cells in x and y directions, set crs
                                            xlim = c(as.numeric(bound_box[1]), as.numeric(bound_box[3])),
                                            ylim = c(as.numeric(bound_box[2]), as.numeric(bound_box[4]))) #set extent of x axis and y axis !extremely important!!!
        
        pit_mound_pair_rast[is.na(pit_mound_pair_rast[])] <- 0 #replace Na values
        
        #save data to raster file (tif)
        write_stars(pit_mound_pair_rast, driver = "GTiff",
                    file.path(pth, '/', site_name, "_pm_pairs.tif", fsep = ''))
        write_stars(poly_PitMound_rast, driver = "GTiff",
                    file.path(pth, '/', site_name, '_rast_PM.tif', fsep = ''))
        
        #WRITE DATA TO SHP
        st_write(poly_PitMound, file.path(pth, '/', site_name, '_poly_PM.shp', fsep = ''))
        st_write(pit_mound_pair, file.path(pth, '/', site_name, '_pm_pairs.shp', fsep = ''))
        st_write(points_SVmax_coords_sf, file.path(pth, '/', site_name, '_pointsMax.shp', fsep = ''))
        st_write(points_SVmin_coords_sf, file.path(pth, '/', site_name, '_pointsMin.shp', fsep = ''))
}




#extract LiDAR Point Cloud data (for teritory located in Poland)
extract_PCl(study_area = study_area, #path to shp with the study area boundary (geometry: polygon)
            site_name = site_name,  #the name of study area (or the prefix which we want to add to the name of each output file) 
            pth = wd_pth, #path to save output files; function will create folder 'las' in this directory and will save files in this folder
            buffer_size = 100) #the size of buffer (in meters) -> buffer created in order to prevent DEM' ridge effects to occur on the study area 


compute_DEM(site_name = site_name, #the name of study area (or the prefix which we want to add to the name of each output file)
            pth = wd_pth, #path to save output files; function will create folder 'las' in this directory and will save files in this folder
            point_cloud = file.path(wd_pth, "las", paste0(site_name, "_KRON86.laz")), #path to point cloud
            study_area = study_area, #path to shp with the study area boundary (geometry: polygon)
            dem_res = dem_res) #resolution of computed DEM

#compute 1) contour lines and 2) length of each contour line
DEM_cnt(DEM = file.path(wd_pth, paste0(site_name, "_DEM01m.tif")), #path do DEM 
        pth = wd_pth, #path to save output contour lines file
        site_name = site_name, #the name of study area (or the prefix which we want to add to the name of each output file)
        interv = cnt_interv) #contour lines interval




#RUN the contour method (CM) algorithm

#detect_PitMound_cnt -> uses st_collection_extract as a step in dissolving polygons
#detect_PitMound_cnt2 -> uses st_cast() as a step in dissolving polygons
tryCatch(detect_PitMound_cnt(site_name = site_name, #the name of study area (or the prefix which we want to add to the name of each output file)
                             pth = CM_pth, #path to save output files
                             contours = file.path(wd_pth, paste0(site_name, "_cnt01m.shp")), #path to shp with the contour lines
                             study_area = file.path(wd_pth, paste0(site_name, ".shp")), #path to shp with the study area boundary (geometry: polygon)
                             DEM = file.path(wd_pth, paste0(site_name, "_DEM01m.tif")), #path to DEM
                             max_cnt_length = max_cnt_len, #maximal length of contour line
                             min_cnt_length = min_cnt_len, #minimal length of contour line,
                             epsg = epsg), #epsg to be used in output files (epsg-code as integer)
         error = function(e)
                 detect_PitMound_cnt2(site_name = site_name,
                                      pth = CM_pth,
                                      contours = file.path(wd_pth, paste0(site_name, "_cnt01m.shp")),
                                      study_area = file.path(wd_pth, paste0(site_name, ".shp")),
                                      DEM = file.path(wd_pth, paste0(site_name, "_DEM01m.tif")),
                                      max_cnt_length = max_cnt_len,
                                      min_cnt_length = min_cnt_len,
                                      epsg = epsg))





