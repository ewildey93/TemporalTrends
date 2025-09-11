library(sswids)
library(sf)
library(readr)
library(dplyr)
library(leaflet)


sswids::connect_to_sswidb(db_version = 'PROD')
classprec <- sswidb::sswidb_species_precision(conn)

grid <- "SSWI"

prec <- 0.95
daterange <- 
  create_season_dates(
    min_date = "-01-01",
    max_date = "-12-31",
    years = c(2019,2024)
  )
species <- c("Bear",
             "Deer",
             "Turkey")

Q <- query_effort(conn = conn, prec = prec, grid = grid, daterange = daterange, remove0Timelapse = TRUE)
#removes location in UP, location in middle of lake and locations without needed precision
Q2 <- rm_bad_locations(locationeffort = Q, coordinate_precision = 4)
#assigns camsite id, creates average location coordinates, removes overlapping effort of more than 1 day
Q3 <- merge_nearby_cameras(locationeffort = Q2, cam_distance = 100)
#pull in detections
detections <- query_detections(conn = conn, species = species, grid = "SSWI", daterange = daterange, prec = prec)
saveRDS(detections, "./rawdetections.rds")
detectionslist <- split(detections, detections$species)
#summarize effort by sampling occasions make prop_classified column
Q4 <-create_sampling_occasions(daterange = daterange, locationeffort = Q3, num_occasions = 52,class_threshold = 0.95)
#join effort and detections data, summarize detections by occasion
Q5.max <- summarize_detections(detections = detections, locationeffort = Q4, summary_value = "max count")
Q5.ct <- lapply(detectionslist, function (x) summarize_detections(detections = x, locationeffort = Q4, summary_value = "count triggers"))


sswids::list_spatial_layers()
# merge spatial information with effort wide form
effort_by_occ_df_maxcount_sf = Q5.max %>%
  sf::st_as_sf(coords=c("lon", "lat"), crs=4326)%>%
  sf::st_join(
    .,
    sf::st_transform(get_spatial_data("furbearer_zones"), 4326),
    join = sf::st_within
  )

effort_by_occ_df_counttriggers_sf = mapply(function(x,y) x %>%
  sf::st_as_sf(coords=c("lon", "lat"), crs=4326)%>%
  sf::st_join(
    .,
    sf::st_transform(get_spatial_data(y), 4326),
    join = sf::st_within
  ), x=Q5.ct, y=c("bear_zones", "dmus", "turkey_mgt_zones"))


write_rds(effort_by_occ_df_counttriggers_sf, "./BearDeerTurkeyList.rds")
