# Download daily MODIS-Aqua Level-3 fields from NOAA CoastWatch ERDDAP (no account needed)
# and convert them to one projected GeoTIFF per variable and day.

#' Download a daily ERDDAP grid subset as NetCDF
#' @param dataset  ERDDAP dataset id, e.g. "erdMH1chla1day" (chlorophyll-a) or "erdMH1sstd1day" (SST)
#' @param variable variable name inside the dataset, e.g. "chlorophyll" or "sst"
#' @param bbox     c(xmin, ymin, xmax, ymax) in degrees
#' @param start,end dates ("YYYY-MM-DD")
#' @param dest     output .nc path
erddap_download <- function(dataset, variable, bbox, start, end, dest,
                            server = "https://coastwatch.pfeg.noaa.gov/erddap") {
  q <- sprintf("%s[(%s):1:(%s)][(%s):1:(%s)][(%s):1:(%s)]", variable, start, end,
               bbox[4], bbox[2], bbox[1], bbox[3])        # latitude is stored north -> south
  url <- sprintf("%s/griddap/%s.nc?%s", server, dataset, utils::URLencode(q, reserved = TRUE))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  invisible(dest)
}

#' Split a multi-day NetCDF into one projected GeoTIFF per day
#' @param nc      path to NetCDF downloaded with erddap_download()
#' @param var     short name used in output file names, e.g. "chla"
#' @param out_dir output folder
#' @param crs     target CRS (default UTM zone 40N, as used for the Persian Gulf / Gulf of Oman)
#' @param res     target resolution in metres
#' @param log10   log-transform (recommended for chlorophyll-a, which is log-normally distributed)
#' @param template optional SpatRaster grid to write onto; pass the same template for every variable so
#'                 that all daily layers share exactly one grid
nc_to_daily_tifs <- function(nc, var, out_dir, crs = "EPSG:32640", res = 1000, log10 = FALSE,
                             template = NULL) {
  r <- terra::rast(nc)
  dates <- as.Date(terra::time(r))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(template)) template <- grid_template(nc, crs, res)
  files <- character(0)
  for (i in seq_len(terra::nlyr(r))) {
    layer <- r[[i]]
    if (log10) layer <- log10(layer)
    out <- terra::project(layer, template, method = "bilinear")
    names(out) <- var
    f <- file.path(out_dir, sprintf("%s_%s.tif", var, format(dates[i], "%Y-%m-%d")))
    terra::writeRaster(out, f, overwrite = TRUE)
    files <- c(files, f)
  }
  invisible(files)
}

#' Projected grid covering a NetCDF file (use one template for all variables)
grid_template <- function(nc, crs = "EPSG:32640", res = 1000) {
  r <- terra::rast(nc)
  terra::rast(terra::project(terra::ext(r), terra::crs(r), crs), crs = crs, res = res)
}
