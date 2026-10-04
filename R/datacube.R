# Event-based satellite data cubes.
#
# For every point event (e.g. a field-observed algal bloom) a daily stack of satellite variables is
# assembled for each day in a window around the event date, and predictor values are extracted at the
# event location. The window is explicit, because it decides what the model can claim:
#
#   window = -10:-1   days BEFORE the event only  -> a genuine forecast with 1-10 days lead time
#   window = 0        the event day only          -> same-day detection / nowcasting
#   window = 0:10     the event day and AFTER     -> assumes the event persists; detection, not forecasting
#
# Only days on which every requested variable exists are used ("complete cubes"), and pixel values that
# are missing (clouds, glint) are flagged rather than silently filled.

#' Index daily GeoTIFFs written by nc_to_daily_tifs() ("<var>_<YYYY-MM-DD>.tif")
index_daily_tifs <- function(dir) {
  f <- list.files(dir, pattern = "_\\d{4}-\\d{2}-\\d{2}\\.tif$", full.names = TRUE, recursive = TRUE)
  b <- basename(f)
  data.frame(var  = sub("_\\d{4}-\\d{2}-\\d{2}\\.tif$", "", b),
             date = as.Date(regmatches(b, regexpr("\\d{4}-\\d{2}-\\d{2}", b))),
             file = f, stringsAsFactors = FALSE)
}

#' Stack all requested variables for one date; NULL if any variable is missing that day
daily_cube <- function(index, date, vars) {
  sel <- index[index$date == as.Date(date) & index$var %in% vars, ]
  if (!all(vars %in% sel$var)) return(NULL)
  r <- terra::rast(sel$file[match(vars, sel$var)])
  names(r) <- vars
  r
}

#' One row per event x day offset
event_windows <- function(events, window) {
  out <- merge(events, data.frame(offset = window), by = NULL)
  out$cube_date <- out$date + out$offset
  out[order(out$id, out$offset), ]
}

#' Extract predictor values at each event for every day in the window
#' @param events data.frame with columns id, date (Date), lon, lat (WGS84)
#' @return data.frame: id, date, offset, cube_date, lon, lat, x, y, <vars>, complete
extract_event_cubes <- function(events, index, vars, window = -10:-1, crs = "EPSG:32640") {
  ew  <- event_windows(events, window)
  pts <- terra::project(terra::vect(ew, geom = c("lon", "lat"), crs = "EPSG:4326"), crs)
  xy  <- terra::crds(pts)
  ew$x <- xy[, 1]; ew$y <- xy[, 2]
  days <- sort(unique(ew$cube_date))
  res <- list()
  for (k in seq_along(days)) {
    cube <- daily_cube(index, days[k], vars)
    if (is.null(cube)) next                               # skip days without a complete cube
    rows <- which(ew$cube_date == days[k])
    v <- terra::extract(cube, pts[rows], ID = FALSE)
    res[[length(res) + 1]] <- cbind(ew[rows, c("id", "date", "offset", "cube_date", "lon", "lat", "x", "y")], v)
  }
  out <- do.call(rbind, res)
  out$complete <- stats::complete.cases(out[, vars])      # FALSE where a pixel is cloudy / missing
  out[order(out$id, out$offset), ]
}

#' Random background (pseudo-absence) points from the complete cubes of the given dates
sample_background <- function(index, dates, vars, n_per_date = 500, seed = 1) {
  set.seed(seed)
  res <- list()
  for (k in seq_along(dates)) {
    cube <- daily_cube(index, dates[k], vars)
    if (is.null(cube)) next
    # cloudy days can have fewer valid cells than requested (or none): that is expected, so the
    # corresponding terra warnings are suppressed and the day simply contributes fewer points
    s <- suppressWarnings(
      terra::spatSample(cube, size = n_per_date, method = "random", na.rm = TRUE, xy = TRUE))
    s <- s[stats::complete.cases(s), , drop = FALSE]
    if (nrow(s) == 0) next
    s$cube_date <- as.Date(dates[k])
    res[[length(res) + 1]] <- s
  }
  do.call(rbind, res)
}

#' Presence (label 1) + background (label 0) table ready for e.g. the sdm package
modelling_table <- function(presence, background, vars) {
  p <- presence[presence$complete, c("cube_date", "x", "y", vars)]
  p$label <- 1
  b <- background[, c("cube_date", "x", "y", vars)]
  b$label <- 0
  rbind(p, b)
}
