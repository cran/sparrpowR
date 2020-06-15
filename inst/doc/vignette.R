## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 7, fig.show = "hold")
options("rgdal_show_exportToProj4_warnings" = "none") # to mute warnings of possible GDAL/OSR exportToProj4() degradation

## ----echo=TRUE, warning=FALSE, message=FALSE, cache=TRUE----------------------
loadedPackages <- c("sparrpowR", "rgdal")
invisible(lapply(loadedPackages, require, character.only = TRUE))
set.seed(1234) # for reproducibility

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
# Washington, D.C. boundary
gis_path1 <- "https://opendata.arcgis.com/datasets/7241f6d500b44288ad983f0942b39663_10.geojson"
dc <- geojsonio::geojson_read(gis_path1,  what = "sp")

# American Community Survey 2018 Census Tracts
gis_path2 <- "https://opendata.arcgis.com/datasets/faea4d66e7134e57bf8566197f25b3a8_0.geojson"
census <- geojsonio::geojson_read(gis_path2,  what = "sp")

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.align='center'----
clipwin <- maptools::unionSpatialPolygons(census, IDs = rep(1, length(census)))
dcc <- rgeos::gIntersection(dc, clipwin, byid = TRUE)
# Plot
plot(dc, main = "DC Boundary")
plot(census,  main = "American Community Survey\n2018")
plot(dcc, main = "Clipped Boundary")

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
dcp <- sp::spTransform(dcc, CRSobj = sp::CRS(projargs = "+init=EPSG:32618"))
dco <- spatstat::as.owin(dcp)

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
navy <- data.frame(lon = 326414.70444451, lat = 4304571.1539442)
spf <- sp::SpatialPoints(coords = navy, proj4string = sp::CRS(projargs = "+init=EPSG:32618"))
# Plot
sp::plot(dcp, main = "Location of Hypothetical\nDisease Cluster")
sp::plot(spf, col = "red", add = T)
legend("bottom", xpd = T, y.intersp = -1.5, legend = c("Navy Yard"), col = "red", pch = 3, bty = "n")

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
start_time <- Sys.time() # record start time
sim_power <- spatial_power(x_case = navy[[1]], y_case = navy[[2]], # center of cluster
                           x_control = navy[[1]], y_control = navy[[2]], # center of cluster
                           n_case = 50, n_control = 950, # sample size of case/control
                           samp_case = "MVN", samp_control = "MVN", # samplers
                           s_case = 1000, s_control = 2000, # approximate size of clusters
                           cascon = FALSE, # power for case cluster(s) only
                           lower_tail = 0.025, upper_tail = 0.975, # two-tailed alpha
                           sim_total = 1, # number of iterations
                           win = dco, # study area
                           resolution = 100, # number gridded knots on x-axis
                           edge = "diggle", # correct for edge effects
                           adapt = FALSE, # fixed-bandwidth
                           h0 = NULL, # automatically select bandwidth for each iteration
                           verbose = FALSE) # no printout
end_time <- Sys.time() # record end time
time_srr <- end_time - start_time # Calculate run time

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.align='center'----
cols <- c("cornflowerblue", "green", "red", "lightgreen", "blue") # colors for plots
chars <- c(4,5) # symbols for point-locations
sizes <- c(0.5,0.5) # size of point-locations
p_thresh <- 0.8 # 80% of iterations with statistically significant results

## Data Visualizaiton of Input and Power
spatial_plots(input = sim_power, # use output of above simulation
              p_thresh = p_thresh, # power cut-off
              plot_pts = TRUE, # display the points in the second plot
              chars = chars, # case, control
              sizes = sizes, # case, control
              cols = cols) # colors of plot

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE, fig.align='center'----
expandbb <- function(bb, f) {
  x <- bb[3] - bb[1] # range of x values
  y <- bb[4] - bb[2] # range of y values
  nb <- bb # make a copy
  nb[1] <- bb[1] - (f * x) # xmin - left
  nb[3] <- bb[3] + (f * x) # xmax - right
  nb[2] <- bb[2] - (f * y) # ymin - bottom
  nb[4] <- bb[4] + (f * y) # ymax - top
  return(nb)}
dcbb <- expandbb(bb = sp::bbox(dc), f = 0.01)
base_map <- ggmap::get_map(location = dcbb, maptype = "toner", source = "stamen")

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
sim_pts <- sim_power$sim  # extract points from first iteration
sim_pts <- maptools::as.SpatialPointsDataFrame.ppp(sim_pts) # convert to spatial data frame
raster::crs(sim_pts) <- sp::proj4string(dcp) # set initial projection
sim_pts_wgs84 <-  sp::spTransform(sim_pts, CRSobj = sp::CRS(projargs = "+init=epsg:4326")) # project to basemap
sim_pts_df <- tibble::tibble(data.frame(sim_pts_wgs84)) # convert to tidy data frame

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
dc_df <- broom::tidy(dcc) # conver to a tidy dataframe
dcc$polyID <- sapply(slot(dcc, "polygons"), function(x) slot(x, "ID")) # preserve polygon id for merge
dc_df <- merge(dc_df, dcc, by.x = "id", by.y="polyID") # merge data

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
pvalprop <- tibble::tibble(x = sim_power$rx, y = sim_power$ry,
                           z = sim_power$pval_prop) # extract proportion significant
lrr_narm <- na.omit(pvalprop) # remove NAs
sp::coordinates(lrr_narm) <- ~ x + y # coordinates
gridded(lrr_narm) <- TRUE # gridded
pvalprop_raster <- raster::raster(lrr_narm) # convert to raster
rm(pvalprop, lrr_narm) # conserve memory
raster::crs(pvalprop_raster) <- raster::crs(dcp) # set output project (UTM 18N)
pvalprop_raster <- raster::projectRaster(pvalprop_raster, crs = raster::crs(dc)) # unproject (WGS84)
rtp <- raster::rasterToPolygons(pvalprop_raster) # convert to polygons
rtp@data$id <- 1:nrow(rtp@data)   # add id column for join
rtpFort <- broom::tidy(rtp, data = rtp@data) # convert to tibble
rtpFortMer <- merge(rtpFort, rtp@data, by.x = 'id', by.y = 'id')  # join data
rampcols <- grDevices::colorRampPalette(colors = c(cols[5], cols[2]), space="Lab")(length(raster::values(pvalprop_raster))) # set colorramp

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
ggmap::ggmap(base_map) + # basemap
  ggplot2::geom_polygon(data = dc_df, # original boundary
               ggplot2::aes(x = long, y = lat, group = group),
               fill = "transparent",
               colour = "black") +
  ggplot2::geom_polygon(data = rtpFortMer, # output raster as polygons
               ggplot2::aes(x = long, y = lat, group = group, fill = z), 
               size = 0, 
               alpha = 0.5) +
  ggplot2::scale_fill_gradientn(colours = rampcols) + # colors for polygons
  ggplot2::geom_point(data = sim_pts_df, # simulated point-locations
             ggplot2::aes(x = mx, y = my, color = marks, shape = marks),
             alpha = 0.8) + 
  ggplot2::scale_color_manual(values = cols[4:5]) + # fill of point-locations
  ggplot2::scale_shape_manual(values = chars) + # shope of point-locations
  ggplot2::labs(x = "", y = "", fill = "Power", color = "", shape = "") # legend labels

## ----echo = TRUE, warning=FALSE, message=FALSE, cache=TRUE--------------------
pvalprop_reclass <- raster::reclassify(pvalprop_raster, c(-Inf, p_thresh-0.0000001, 1,
                                                          p_thresh-0.0000001, Inf, 2))
rtp <- raster::rasterToPolygons(pvalprop_reclass) # convert to polygons
rtp@data$id <- 1:nrow(rtp@data)   # add id column for join
rtpFort <- broom::tidy(rtp, data = rtp@data) # convert to tibble
rtpFortMer <- merge(rtpFort, rtp@data, by.x = 'id', by.y = 'id')  # join data

ggmap::ggmap(base_map) + # basemap 
  ggplot2::geom_polygon(data = dc_df, # original boundary
               ggplot2::aes(x = long, y = lat, group = group),
               fill = "transparent",
               colour = "black") +
  ggplot2::geom_polygon(data = rtpFortMer, # output raster as polygons
               ggplot2::aes(x = long, y = lat, group = group, fill = as.factor(z)), 
               size = 0, 
               alpha = 0.5) +
  ggplot2::scale_fill_manual(values = cols[c(5,2)],
                             labels = c("insufficient", "sufficient")) + # colors for polygons
  ggplot2::labs(x = "", y = "", fill = "Power") # legend labels

