#This code was mainly developed by Moritz Adam (https://github.com/MoritzAdam) for the Climatology and Biosphere (GUZ Tübingen) group, led by Kira Rehfeld. 

library(tidyverse)
library(ggnewscale)
library(rgdal)
library(RColorBrewer)
library(sp)
library(raster)
library(sf)
library(viridisLite)
select <- dplyr::select

# environment for data caching
.STACYmap <- environment()

# style options
GLOBAL_STACY_OPTIONS <- list(
  GLOBAL_CRS = list(equirectangular = 'fixed',
                    robinson = '+proj=robin',
                    wintri = '+proj=wintri', 
                    azequidistant = '+proj=aeqd',
                    area = '+proj=aea +lat_1=29.5 +lat_2=42.5'),
  GLOBAL_GREY_LIGHT = grey.colors(1, start = 0.97, end = 0.99),
  GLOBAL_GREY_LIGHT_ALPHA_LOW = grey.colors(1, start = 0.97, end = 0.99, alpha = 0.45),
  GLOBAL_GREY_DARK = grey.colors(1, start = 0.6),
  GLOBAL_GREY_MEDIUM = grey.colors(1, start = 0.88),
  GLOBAL_GREY_MEDIUM_LIGHT = grey.colors(1, start = 0.95),
  GLOBAL_FONT_FACE_TITLE = 'bold',
  GLOBAL_FONT_FACE_TEXT = 'plain',
  GLOBAL_FONT_SIZE = 9,
  GLOBAL_FONT_FAMILY = 'sans',
  GLOBAL_ARROW_SIZE = 0.03,
  GLOBAL_FIELD_SIZE = 0.2,
  GLOBAL_FIELD_ALPHA = 0.6,
  GLOBAL_POINT_SIZE = 3,
  GLOBAL_POINT_SIZE_DISCR = 3.5,
  GLOBAL_POINT_STROKE = 0.6,
  GLOBAL_LEG_TITLE_VJUST = 1,
  #GLOBAL_BATHYMETRY_COLORS = c('#d9ebf9', '#cae1f4', '#afd3ef', '#aacde9', '#96c1e3', '#83b9df', '#6fadd6', '#5ba2d0', '#589cc9', '#337fb2', '#2a77ac', '#2371a6'),
  GLOBAL_LAND_COLOR = '#f0e6c2',
  GLOBAL_OCEAN_COLOR = '#aacde9',
  GLOBAL_SHAPE_VALUES = c(3, 4, 6, 7, 9, 8, 10, 15, 16, 17, 18)# c('\u25BE', '\u25C6', '\u25A0', '\u25CF', '\u25B2', 'circle cross', 'cross', 'asterisk', 'plus', 'circle plus') # c('\u25BE', '\u25C6', '\u25A0', '\u25CF', '\u25B2', '\u1F7AC', '\u1F7CA', '\u1F7CE', '\u1F7C2', '\u1F7A5')#c(23, 22, 24, 21, 0:20)
)
# helpers
title_and_axis <- function(){
  ax <- theme_bw() + 
    theme(title = element_text(face = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE, size = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE + 1, family = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FAMILY),
          text = element_text(face = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TEXT, size = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE, family = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FAMILY), 
          #line = element_blank(),
          legend.direction = 'vertical', legend.box = 'vertical', legend.box.background = element_blank(), legend.background = element_rect(colour = 'black'), 
          strip.background.x = element_blank(), strip.text.x = element_text(face = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE, hjust = 0), 
          strip.background.y = element_rect(fill = 'white'), strip.text.y = element_text(face = GLOBAL_STACY_OPTIONS$GLOBAL_FONT_FACE_TITLE, vjust = 0))
  return(ax)
}

load_natural_earth_data <- function(file, dir = 'data/shapes/50m_physical/', ...) {
  NAT_EARTH_PATH <- 'data/shapes/50m_physical/'
  if (dir == NAT_EARTH_PATH & str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(NAT_EARTH_PATH, file), verbose = FALSE, ...)
  } else if (str_detect(file, '.shp')) {
    data <- readOGR(dsn = file.path(dir, file), verbose = FALSE, ...)
  } else if (str_detect(file, '.tif')) {
    data <- raster::raster(x = file.path(NAT_EARTH_PATH, file))
  }
  return(data)
}

transform_shapefile <- function(data, transform) {
  if (transform == 'fixed') {
    return(data)
  }
  else {
    return(spTransform(data, CRSobj = CRS(transform)))
  }
}





#' Plot appealing maps of latlon-gridded and point-wise data
#'
#' @param gridlyr latlon-grid to plot of class raster, matrix (both lon columns, lat rows, cell value) or data.frame/tibble (either lon columns, lat rows, cell value or lon, lat, cell columns)
#' @param ptlyr point layer of class matrix or data.frame (lon, lat, cell columns)
#' @param ptlyr_size_var character. optional pointlayer variable to map to point size (given as additional column in ptlyr). values can only be TRUE or FALSE in column!! (e.g. occurences)
#' @param ptlyr_shape_var character. optional pointlayer variable to map to point shape (given as additional column in ptlyr). (e.g. different databases)
#' @param fldlyr field layer of class data.frame (lon,lat, angle(, radius)) containing the angle of the field (mandatory) and the radius (optional)
#' @param coastline logical - plot coastline?
##' @param bathymetry logical - plot ocean bathymetry?, option can be used independently of coastlines and will underlay a climate field
#' @param projection object of class CRS or character (is converted to CRS), see package rgdal
#' @param graticules logical, should a projected panel.grid be plotted?
#' @param zoom numeric, c(west long, south lat, east long, north lat) to be cut off
#' @param legend_names list(grid = nm1, pt = nm2)
#' @param colorscheme color palette(s) or character from 'temp', 'prcp_div', 'prcp_grd', 'spec' - pass one color palette to be respected by grid and points or list(grid = pal1, pt = pal2) for independent color schemes
#' @param legend_inside logical - should the legend be placed inside the map plot?
#' @param legend_cb logical - should the color/fill legend be a continous colorbar? if FALSE: legend as a discrete legend
#' @param legend_num_breaks numeric - number of legend breaks
#' @param revcolorscheme logical - has effect only on built-in schemes
#' @param flipy logical - flip the grid input in y direction, try option if map orientation is wrong
#' @param transpose logical - transpose the grid input, try option if map orientation is wrong
#' @param rotate logical - rotate the map by 180 degrees
#' @param box logical - should the entire plot be bounded by a box?
#' @param filledbg logical - should land and ocean surface be coloured?, option can be used independently of coastlines and will underlay a climate field
#' @param centercolor numeric - or list(grid = val1, pt = val2) of numeric to center colorscheme to a specific value 
#' @param limits  c(-1.5, 3)
#'
#' @return
#' @export
#'
#' @examples 
#' # plotting a field and points
#' nc <- nc_open(PLASIM_TESTSLICE_PI)
#' testfld <- ncvar_get(nc, varid = "ts")
#' lat <- ncvar_get(nc, varid = "lat")
#' lon <- ncvar_get(nc, varid = "lon")
#' nc_close(nc)
#' testfld <- testfld[,,1]
#' ptlyr <- tibble(long = c(65, 30), lat = c(45, -60), value = c(300, 260))
#' STACYmap(gridlyr = testfld, ptlyr = ptlyr, colorscheme = 'temp',
#'          graticules = T, zoom = c(-100, -45, 100, 45), 
#'          legend_names = list(grid = 'temp', pt = 'temp'))
#' 
#' # plotting points on a world map
#' STACYmap(ptlyr = ptlyr, colorscheme = 'temp', filledbg = T, bathymetry = T)
#' 
#' # using other projections than the default robinson (here: Lambert azimuthal equal-area)
#' # disclaimer: some combinations of projection and raster plots are not working yet
#' plt <- STACYmap(ptlyr = ptlyr, colorscheme = 'temp', filledbg = T, projection = '+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000')
STACYmap <- function(gridlyr = NULL, 
                     ptlyr = NULL, 
                     ptlyr_size_var = NULL,
                     ptlyr_shape_var = NULL,
                     fldlyr = NULL,
                     coastline = TRUE,
                     filledbg = FALSE,
                    # bathymetry = FALSE,
                     projection = as.character('+proj=robin +datum=WGS84'), #+proj=longlat +datum=WGS84
                     graticules = TRUE, 
                     zoom = NULL, 
                     legend_names = list(grid = 'grid', pt = 'points'), 
                     legend_inside = FALSE,
                     legend_cb = TRUE,
                     legend_num_breaks = 5,
                     colorscheme = rev(RColorBrewer::brewer.pal(11, 'RdBu'))[2:10],
                     centercolor = NULL,
                     revcolorscheme = F,
                     flipy = FALSE, 
                     transpose = TRUE, 
                     rotate = FALSE,
                     box = TRUE,
                     limits = c(-1.5, 3)){

  # projections
  ptlyro <- ptlyr; gridlyro <- gridlyr
  #if (class(projection) == 'character') {
  #  projection <- CRS(projection)
  #}
  splcol <- FALSE; if (all(class(colorscheme) == 'list')) {splcol <- TRUE}
  if ((class(colorscheme) == 'list' & all(lapply(colorscheme[1:length(colorscheme)], length) == 1)) | (class(colorscheme) == 'character' & length(colorscheme) == 1)) {
    colorschemeo <- colorscheme
    colorscheme <- lapply(colorscheme, 
                          function (x) {
                            if (is.character(x) & length(x) == 1) {
                              x <- list(
                                temp = rev(RColorBrewer::brewer.pal(11, 'RdBu')),#[2:10],
                                prcp_grd = RColorBrewer::brewer.pal(9, 'YlGnBu'), 
                                prcp_div = RColorBrewer::brewer.pal(11, 'BrBG'),
                                spec = rev(RColorBrewer::brewer.pal(11, 'RdBu'))
                              )[[x]]
                              
                              if (is.null(x)) {stop('unknwon colorscheme given to STACYmap')}
                              if (revcolorscheme) {x <- rev(x)}
                              return(x)
                            } else {
                              return(x)
                            }
                          })
    if (class(colorschemeo) == 'character' & length(colorschemeo) == 1) {
      colorscheme <- unlist(colorscheme)
    }
  } #else if (any(lapply(colorscheme[1:length(colorscheme)], length) != 1)) {
  #print(lapply(colorscheme[1:length(colorscheme)], length))
  #stop('mixing pre-defined and built-in colorschemes is not supported, please choose one option for all layers')
  #}
  
  
  ## checks and data preparation/projection
  #if (length(colorscheme) == 1) {colorscheme <- unlist(colorscheme)}
  
  xlim <- NULL; ylim <- NULL
  if (!is.null(zoom)) {
    zoom <- matrix(c(zoom[1], 0, zoom[3], 0, 0, zoom[2], 0, zoom[4]), nrow = 4)
    zoom <- project(zoom,
                    proj = as.character(projection))
    xlim <- c(zoom[1,1], zoom[3,1])
    ylim <- c(zoom[2,2], zoom[4,2])
  }
  
  if (!is.null(gridlyr)) {
    if (!any(class(gridlyr) %in% c('raster', 'matrix', 'data.frame'))) {
      stop('class(gridlyr) has to be one of \'raster\', \'matrix\', \'data.frame\'')
    } else {
      prepg <- list(raster = function(x) x, 
                    matrix = function(x) {
                      if (!rotate) {
                        o <- raster::raster({if (transpose) {t(x)} else {x}}, 
                                            crs = "+proj=longlat +datum=WGS84", xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% 
                                            {if (flipy) {raster::flip(., 'y')} else {.}}
                      } else {
                        o <- raster::raster({if (transpose) {t(x)} else {x}}, 
                                            crs = "+proj=longlat +datum=WGS84", xmn = 0, xmx = 360, ymn = -90, ymx = 90) %>% 
                          raster::rotate() %>% 
                          {if (flipy) {raster::flip(., 'y')} else {.}}
                      }
                      o
                    },
                    data.frame = function(x) as.matrix(x) %>% prepg$raster(.)
      )
      
      gridlyr <- gridlyr %>% 
        prepg[[class(gridlyr)]](.) %>% 
        raster::rasterToPolygons() %>% 
        spTransform(CRS(projection)) %>%
        as('SpatialPolygonsDataFrame') %>% 
        as('sf')
    }
  }
  if (!is.null(ptlyr)) {
    if (!any(class(ptlyr) %in% c('matrix', 'data.frame'))) {
      stop('class(ptlyr) has to be one of \'matrix\', \'data.frame\'')
    } else {
      prep <- list(matrix = function(x) as_tibble(as.data.frame(x)), 
                   data.frame = function(x) as_tibble(x)
      )
      ptlyr <- ptlyr %>% 
        prep[[class(ptlyr)[length(class(ptlyr))]]](.) %>% 
        rename(long = 1, lat = 2, layer = 3)
      
      ptlyr <- ptlyr %>% 
        rownames_to_column('ptid')
      ptlyr_trf <- project(cbind(ptlyr$long, ptlyr$lat),
                           proj = as.character(projection)) %>% 
        as_tibble() %>% 
        rename(long = V1, lat = V2) %>% 
        bind_cols(select(ptlyr, ptid)) %>% 
        inner_join(select(ptlyr, -lat, -long), by = 'ptid')
      ptlyr <- ptlyr_trf
      rm(ptlyr_trf)
    }
  }
  #fieldlyr projection
  if (!is.null(fldlyr)) {
    if (!any(class(fldlyr) %in% c('matrix', 'data.frame'))) {
      stop('class(fldlyr) has to be one of \'matrix\', \'data.frame\'')
    } else {
      prep <- list(matrix = function(x) as_tibble(as.data.frame(x)), 
                   data.frame = function(x) as_tibble(x)
      )
      fldlyr <- fldlyr %>% 
        prep[[class(fldlyr)[length(class(fldlyr))]]](.) %>% 
        rename(long = 1, lat = 2, angle = 3, radius = 4)
      
      fldlyr <- fldlyr %>% 
        rownames_to_column('fldid')
      fldlyr_trf <- project(cbind(fldlyr$long, fldlyr$lat),
                            proj = as.character(projection)) %>% 
        as_tibble() %>% 
        rename(long = V1, lat = V2) %>% 
        bind_cols(select(fldlyr, fldid, angle, radius))
      fldlyr <- fldlyr_trf
      rm(fldlyr_trf)
    }
  }
  
  
  #lims <- NULL
  if (!is.null(centercolor)) {
    if ((class(centercolor) != 'list' | length(centercolor) == 1) & splcol) {
      centercolor <- list(
        grid = centercolor[[1]], 
        pt = centercolor[[1]]
      )
    }
    
    if (!splcol & !is.null(gridlyro) & !is.null(ptlyro)) {
      if (class(gridlyro) == 'raster') {
        allmx <- max(max(abs(raster::values(gridlyro)), na.rm = T), 
                     max(abs(ptlyro[3]), na.rm = T))
      } else {
        allmx <- max(max(abs(gridlyro), na.rm = T), 
                     max(abs(ptlyro[3]), na.rm = T))
      }
      lims <- list(
        grid = c(-allmx, allmx), pt = c(-allmx, allmx)
      )
    } else if(!is.null(gridlyro) & !is.null(ptlyro)) {
      if (class(gridlyro) == 'raster') {
        allmxg <- max(abs(raster::values(gridlyro)), na.rm = T)
      } else {
        allmxg <- max(abs(gridlyro), na.rm = T)
      }
      allmxp <- max(abs(ptlyro[3]), na.rm = T)
      lims <- list(
        grid = c(-allmxg, allmxg), pt = c(-allmxp, allmxp)
      )
    } else if (!is.null(ptlyro)) {
      allmx <- max(abs(ptlyro[3]), na.rm = T)
      lims <- list(
        pt = c(-allmx, allmx)
      )
    } else {
      if (class(gridlyro) == 'raster') {
        allmx <- max(abs(raster::values(gridlyro)), na.rm = T)
      } else {
        allmx <- max(abs(gridlyro), na.rm = T)
      }
      lims <- list(
        grid = c(-allmx, allmx)
      )
    }
  } else if (!is.null(gridlyr) & !is.null(ptlyr)) {
    if (splcol) {
      lims <- list(
        grid = c(min(gridlyr$layer), max(gridlyr$layer)),
        pt = c(min(ptlyr$layer), max(ptlyr$layer))
      )
    } else {
      lims <- list (
        grid = c(min(c(min(ptlyr$layer, na.rm = T), min(gridlyr$layer))), max(c(max(ptlyr$layer, na.rm = T), max(gridlyr$layer))))
      )
    }
  } else if (!is.null(gridlyr)) {
    lims <- list(
      grid = c(min(gridlyr$layer), max(gridlyr$layer))
    )
  } else if (!is.null(ptlyr)) {
    lims <- list(
      grid = c(min(ptlyr$layer), max(ptlyr$layer))
    )
  }
  
  ## plotting
  map_plot <- ggplot()
  
  # bathymetry
#  load_bathy <- FALSE
#  if (bathymetry | filledbg) {
#    if (exists('bathy_data', envir = .STACYmap)) {
#      if (!is.null(.STACYmap$bathy_data[[as.character(projection)]])) {
#        bathy_data <- .STACYmap$bathy_data[[as.character(projection)]]#lapply(as.vector(names(.STACYmap$bathy_data[[projection]])), function(d, x){as_tibble(x[[d]]) %>% mutate(depth = d)},
#        #       x = .STACYmap$bathy_data[[projection]]) %>% setNames(., sort(c(200, seq(0, 10000, 1000))))
#      } else {
#        load_bathy <- TRUE
#      }
#    } else {
#      load_bathy <- TRUE
#    }
#  }
  
  #if (load_bathy) {
  #  bathy_data <- lapply(sort(c(200, seq(0, 10000, 1000))), function (x, L) {
  #    load_natural_earth_data(file = paste0('ne_10m_bathymetry_', L[[as.character(x)]], '_', x, '.shp')) %>% 
  #      spTransform(., CRSobj = CRS(projection)) %>% 
  #      fortify(.)}, L = LETTERS[seq(1, 12, 1) %>% rev()] %>% 
  #      setNames(., c(200, seq(0, 10000, 1000)) %>% sort())) %>% 
  #    setNames(., sort(c(200, seq(0, 10000, 1000))))
  #  .STACYmap$bathy_data[[as.character(projection)]] <- bathy_data
  #}
  
  if (filledbg) {
    map_data <- load_natural_earth_data(file = 'ne_50m_land.shp') %>% 
      spTransform(., CRSobj = CRS(projection)) %>% 
      fortify(.)
    
   # if (bathymetry) {
   #   deps <- c(0, 200, seq(1e3, 1e4, 1e3)) %>% as.character()
   #   cols <- GLOBAL_STACY_OPTIONS$GLOBAL_BATHYMETRY_COLORS
   #   for (i in 1:12) {
   #     map_plot <- map_plot + 
   #       geom_polygon(data = bathy_data[[deps[i]]],
   #                    mapping = aes(x = long, y = lat, group = group),
   #                    fill = cols[i], 
   #                    na.rm = T)
   #   }
   #   map_plot <- map_plot + 
   #     geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole)) + 
   #     scale_fill_manual(values = c(GLOBAL_STACY_OPTIONS$GLOBAL_LAND_COLOR, NA), guide = FALSE) #+ 
      #new_scale_fill() + new_scale_color()
   # } else {
      # TODO fix Caspian sea hole
      map_plot <- map_plot + 
        geom_polygon(data = bathy_data[['0']], mapping = aes(x = long, y = lat, group = group,  fill = TRUE, color = TRUE)) + 
        geom_polygon(data = map_data, mapping = aes(x = long, y = lat, group = group, fill = hole, color = hole)) + 
        scale_fill_manual(values = c(GLOBAL_STACY_OPTIONS$GLOBAL_LAND_COLOR, GLOBAL_STACY_OPTIONS$GLOBAL_OCEAN_COLOR), guide = FALSE) + 
        scale_color_manual(values = c(GLOBAL_STACY_OPTIONS$GLOBAL_LAND_COLOR, GLOBAL_STACY_OPTIONS$GLOBAL_OCEAN_COLOR), guide = FALSE) #+ 
      #new_scale_fill() + new_scale_color()
    #}
    
    if (!is.null(gridlyr) | !is.null(ptlyr)) {
      map_plot <- map_plot + 
        new_scale_fill() + new_scale_color()
    }
  }
  
  # define legend types
  if (legend_cb) {
    clrgd <- function(title, ...) {
      return(
        guide_colorbar(title,
                       direction = 'horizontal', 
                       title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST,
                       barwidth = 10, barheight = 0.3)
      )
    }
  } else {
    clrgd <- function(title, breaks) {
      return(
        guide_legend(title, 
                     direction = 'horizontal',
                     label.position = 'bottom',
                     title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST,
                     breaks = breaks)
      )
    }
  }
  
  # gridlyr
  # TODO: option fopr multiple gridlyrs
  if (!is.null(gridlyr)) {
    map_plot <- map_plot + 
      geom_sf(data = gridlyr,
              mapping = aes(color = layer, fill = layer))
    if(!legend_cb){
      map_plot <- map_plot + scale_fill_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}},
                                                  limits = limits,
                                                  breaks = signif(seq(-1., 3., length.out = legend_num_breaks), digits = 3),
                                                  labels = signif(seq(-1., 3., length.out = legend_num_breaks), digits = 3),
                                                  guide = guide_legend(title = legend_names$grid, title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST,
                                                                       ncol = legend_num_breaks, direction = 'horizontal', label.position = 'bottom')) +
        scale_color_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}}, 
                              limits = limits, guide = FALSE) 
      
      
    } else {
      map_plot <- map_plot + scale_fill_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}},
                                                  limits = limits,
                                                  guide = clrgd(legend_names$grid,
                                                                breaks = signif(seq(-1., 3, length.out = legend_num_breaks), digits = 3))) +
        scale_color_gradientn(colors = {if (splcol) {colorscheme$grid} else {colorscheme}}, 
                              limits = limits, guide = FALSE)
    } 
    
    if (splcol) {
      map_plot <- map_plot + 
        new_scale_color() + 
        new_scale_fill()
    }
  }
  
  # coastline contour
  if (coastline) {
    load_map <- FALSE
    if (exists('map_data', envir = .STACYmap)) {
      if (!is.null(.STACYmap$map_data[[as.character(projection)]])) {
        map_data <- .STACYmap$map_data[[as.character(projection)]]
      } else {
        load_map <- TRUE
      }
    } else {
      load_map <- TRUE
    }
    
    if (load_map) {
      map_data <- load_natural_earth_data(file = 'ne_50m_land.shp') %>% 
        spTransform(., CRSobj = CRS(projection)) %>% 
        fortify(.)
      .STACYmap$map_data[[as.character(projection)]] <- map_data
    }
    
    map_data$id <- as.numeric(map_data$id)
    map_data[map_data$id <= 0, 'id'] <- NA
    
    map_plot <- map_plot + 
      geom_polygon(data = map_data,
                   mapping = aes(x = long, y = lat, group = group),
                   fill = NA,
                   colour = 'black',
                   size = 0.25)
  }
  
  #fieldlyr
  if(!is.null(fldlyr) & !is.null(fldlyr$radius)){
    map_plot <- map_plot +
      geom_spoke(data = fldlyr,
                 mapping = aes(x=long, y= lat, angle = angle, radius = scales::rescale(radius, c(1e4, 8e5))),
                 arrow = arrow(length = unit(GLOBAL_STACY_OPTIONS$GLOBAL_ARROW_SIZE, 'inches')),
                 size = GLOBAL_STACY_OPTIONS$GLOBAL_FIELD_SIZE,
                 alpha = GLOBAL_STACY_OPTIONS$GLOBAL_FIELD_ALPHA)
  } else if(!is.null(fldlyr)){
    map_plot <- map_plot +
      geom_spoke(data = fldlyr,
                 mapping = aes(x=long, y= lat, angle = angle, radius = 100000),
                 arrow = arrow(length = unit(GLOBAL_STACY_OPTIONS$GLOBAL_ARROW_SIZE, 'inches')),
                 size = GLOBAL_STACY_OPTIONS$GLOBAL_FIELD_SIZE,
                 alpha = GLOBAL_STACY_OPTIONS$GLOBAL_FIELD_ALPHA)
  }
  
  # pointlyr
  if (!is.null(ptlyr)) {
    if (!is.null(ptlyr_size_var) & !is.null(ptlyr_shape_var)) {
      ptlyr[[ptlyr_size_var]] <- factor(as.character(ptlyr[[ptlyr_size_var]]), levels = c('TRUE', 'FALSE'))
      map_plot <- map_plot + 
        geom_point(data = ptlyr,
                   mapping = aes(x = long, y = lat, fill = layer, size = !!sym(ptlyr_size_var), shape = !!sym(ptlyr_shape_var)),  
                   alpha = I(0.7), stroke = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_STROKE,
                   show.legend = c(fill = {if (is.null(gridlyr) | splcol) {TRUE} else {FALSE}}, size = TRUE)) + 
        scale_shape_manual(values = c(23, 22, 24, 21, 0:20), 
                           guide = guide_legend(nrow = 1, direction = 'horizontal', title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST,
                                                override.aes = list(size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, fill = GLOBAL_STACY_OPTIONS$GLOBAL_GREY_DARK, alpha = I(0.7)), order = 4)) + 
        scale_size_manual(values = c(GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE / 2), 
                          guide = guide_legend(ptlyr_size_var, nrow = 1, direction = 'horizontal', title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST, order = 3))
    } else if (!is.null(ptlyr_size_var)) {
      ptlyr[[ptlyr_size_var]] <- factor(as.character(ptlyr[[ptlyr_size_var]]), levels = c('TRUE', 'FALSE'))
      map_plot <- map_plot + 
        geom_point(data = ptlyr,
                   mapping = aes(x = long, y = lat, fill = layer, size = !!sym(ptlyr_size_var)),  
                   alpha = I(0.7),  shape = 21, stroke = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_STROKE,
                   show.legend = c(fill = {if (is.null(gridlyr) | splcol) {TRUE} else {FALSE}}, size = TRUE)) + 
        scale_size_manual(values = c(GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE / 2), 
                          guide = guide_legend(ptlyr_size_var, nrow = 1, direction = 'horizontal', title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST, order = 3))
    } else if (!is.null(ptlyr_shape_var)) {
      map_plot <- map_plot + 
        geom_point(data = ptlyr,
                   mapping = aes(x = long, y = lat, fill = layer, shape = !!sym(ptlyr_shape_var)),  
                   alpha = I(0.7), size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, stroke = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_STROKE,
                   show.legend = c(fill = {if (is.null(gridlyr) | splcol) {TRUE} else {FALSE}}, size = TRUE)) + 
        scale_shape_manual(values = c(23, 22, 24, 21, 0:20), 
                           guide = guide_legend(nrow = 1, direction = 'horizontal', 
                                                override.aes = list(size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, fill = GLOBAL_STACY_OPTIONS$GLOBAL_GREY_DARK, alpha = I(0.7)), title.vjust = GLOBAL_STACY_OPTIONS$GLOBAL_LEG_TITLE_VJUST, order = 3))
    } else {
      map_plot <- map_plot + 
        geom_point(data = ptlyr,
                   mapping = aes(x = long, y = lat, fill = layer),  
                   alpha = I(0.7), shape = 21, size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE, stroke = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_STROKE,
                   show.legend = c(fill = {if (is.null(gridlyr) | splcol) {TRUE} else {FALSE}}, size = TRUE))
    }
    
    
    
    if (splcol | is.null(gridlyr)) {
      map_plot <- map_plot + 
        scale_fill_gradientn(colors = {if (splcol) {colorscheme$pt} else {colorscheme}}, 
                             limits = limits, 
                             guide = clrgd(legend_names$pt,
                                           breaks = signif(seq(-1, 3, length.out = legend_num_breaks), digits = 3)))
    }
  }
  
  map_plot <- map_plot + 
    title_and_axis() + 
    theme(panel.ontop = T,
          panel.background = element_blank(),
          axis.title = element_blank(), 
          legend.position = 'bottom',
          legend.box = 'horizontal',
          panel.grid = element_line(colour = "darkgrey", size = 0.05),#element_blank(),, 
          plot.background = element_blank(), 
          axis.ticks = element_blank()) # element_line(linetype = 'af', colour = 'grey50'))
  
  if (!graticules) {
    map_plot <- map_plot + 
      coord_sf(label_graticule = 'SW', label_axes = '--EN',
               crs = CRS(projection),
               datum = NA,
               xlim = {if (!is.null(xlim)) {xlim}},
               ylim = {if (!is.null(ylim)) {ylim}},#c(-0.63*1e7, 0.77*1e7),
               expand = F) + 
      theme(axis.text = element_blank())
  } else {
    map_plot <- map_plot  + 
      coord_sf(label_graticule = 'SW', label_axes = '--EN',
               crs = CRS(projection),
               xlim = xlim,#{if (!is.null(xlim)) {xlim}},
               ylim = ylim,#{if (!is.null(ylim)) {ylim}},#c(-0.63*1e7, 0.77*1e7),
               expand = F)
  }
  
  map_plot <- map_plot + 
    scale_x_continuous() +#seq(-90, 90, 90)) + 
    scale_y_continuous() #seq(-60, 60, 30))
  
  if (!box) {
    map_plot <- map_plot + 
      theme(panel.border = element_blank())
  }
  
  if (legend_inside) {map_plot <- map_plot + 
    theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0), legend.box = 'vertical')
  }
  
  return(map_plot)
}