# End-to-end example: MODIS-Aqua chlorophyll-a and SST (public, via NOAA ERDDAP) -> daily GeoTIFFs
# -> event-based data cubes -> presence / background modelling table.
#
# The EVENTS in this example are SYNTHETIC (random points and dates), used only to demonstrate the
# workflow. They are not observations.
#
# Run from the repository root:   Rscript examples/run_example.R

suppressPackageStartupMessages(library(terra))
source("R/download.R")
source("R/datacube.R")
set.seed(2008)
dir.create("outputs", showWarnings = FALSE)

bbox   <- c(56.5, 23.5, 59.5, 26.0)          # Gulf of Oman (xmin, ymin, xmax, ymax)
start  <- "2008-10-01"; end <- "2008-11-30"
window <- -10:-1                              # 1-10 days BEFORE each event: forecasting set-up
vars   <- c("chla", "sst")

## 1. Download and convert (skipped if already present) -----------------------
if (!file.exists("data/chla.nc")) erddap_download("erdMH1chla1day", "chlorophyll", bbox, start, end, "data/chla.nc")
if (!file.exists("data/sst.nc"))  erddap_download("erdMH1sstd1day", "sst",         bbox, start, end, "data/sst.nc")
grid <- grid_template("data/chla.nc")                     # one shared 1 km UTM-40N grid
nc_to_daily_tifs("data/chla.nc", "chla", "data/daily", log10 = TRUE, template = grid)
nc_to_daily_tifs("data/sst.nc",  "sst",  "data/daily", template = grid)
idx <- index_daily_tifs("data/daily")
message("Daily GeoTIFFs: ", nrow(idx), " (", length(unique(idx$date)), " days)")

## 2. Synthetic events on water pixels ----------------------------------------
water <- terra::app(terra::rast(idx$file[idx$var == "chla"]), fun = mean, na.rm = TRUE)
pts   <- terra::spatSample(water, 25, method = "random", na.rm = TRUE, as.points = TRUE)
ll    <- terra::crds(terra::project(pts, "EPSG:4326"))
events <- data.frame(id = seq_len(nrow(ll)), lon = ll[, 1], lat = ll[, 2],
                     date = as.Date("2008-10-15") + sample(0:46, nrow(ll), replace = TRUE))

## 3. Event cubes and predictor extraction ------------------------------------
pres <- extract_event_cubes(events, idx, vars, window = window)
write.csv(pres, "outputs/event_cube_values.csv", row.names = FALSE)

comp <- aggregate(complete ~ offset, data = pres, FUN = function(z) round(100 * mean(z), 1))
names(comp)[2] <- "percent_cloud_free"
write.csv(comp, "outputs/completeness_by_offset.csv", row.names = FALSE)
print(comp)

## 4. Background points from the same cubes and the modelling table -----------
bg  <- sample_background(idx, sort(unique(pres$cube_date)), vars, n_per_date = 100)
tab <- modelling_table(pres, bg, vars)
write.csv(tab, "outputs/modelling_table.csv", row.names = FALSE)
message("Modelling table: ", sum(tab$label == 1), " presence rows, ", sum(tab$label == 0), " background rows")

## 5. Figure: one complete daily cube with the events whose window includes it --
d   <- as.Date(names(sort(table(pres$cube_date[pres$complete]), decreasing = TRUE))[1])
cube <- daily_cube(idx, d, vars)
on   <- terra::project(terra::vect(pres[pres$cube_date == d, ], geom = c("lon", "lat"), crs = "EPSG:4326"), cube)
png("outputs/example_cube.png", width = 1400, height = 600, res = 130)
par(mfrow = c(1, 2))
terra::plot(cube$chla, main = paste("log10 chlorophyll-a,", d)); terra::points(on, pch = 21, bg = "red")
terra::plot(cube$sst,  main = paste("SST (deg C),", d));         terra::points(on, pch = 21, bg = "red")
dev.off()
message("Done. Outputs in outputs/")
