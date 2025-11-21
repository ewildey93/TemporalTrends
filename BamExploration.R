library(sf)
library(mgcv)
library(dplyr)
library(gratia)
library(ggplot2)
library(modelbased)
library(insight)
library(performance)
library(purrr)
library(landscapemetrics)
library(sswids)
library(terra)
library(readxl)
library(tidyr)
library(DBI)

sswids::connect_to_sswidb(db_version = 'PROD')

Turkey <- readRDS("./TurkeyCPUE_CT.rds")
BearBobTurkeyList <- readRDS("./BearBobTurkeyList.rds")
Turkey2week <- BearBobTurkeyList[["Turkey"]]

##############################################################################################################
####                                        cam site level                                                ####
##############################################################################################################
Turkey <- Turkey%>%filter(!grepl(x = camera_version, pattern = ","))
Turkey2week <- Turkey2week%>%filter(!grepl(x = camera_version, pattern = ","))
camsites <- Turkey%>%select(cam_site_id)%>%distinct()%>%st_transform(., 3071)
camsites2week <- Turkey2week%>%select(cam_site_id)%>%distinct()%>%st_transform(., 3071)

Wiscland3 <- get_spatial_data(layer_name = 'wiscland2', level = 3)
Wiscland3.3071 <- project(Wiscland3, "EPSG:3071")
WisclandGuide <- read_xlsx("C:/Users/wildeefb/Documents/GeoSpatial/wiscland2/user_guide/Wiscland2 Color Scheme.xlsx",
                           range = "B2:C70", .name_repair = make.names)

newocc <- data.frame(occ=1:52, occ2week=rep(1:26, each=2))

# Wiscland2 has 30m resolution

buffers <- c(5000, 8000)

camsites <- Forest

trailtype <- dbGetQuery(conn, "SELECT
g83100.sswi_metadata.metadata_name,
g83100.sswi_metadata.metadata_text,
g83100.sswi_camera_location.camera_location_seq_no,
g83100.sswi_event.event_seq_no
FROM
g83100.sswi_metadata
INNER JOIN g83100.sswi_event ON g83100.sswi_metadata.event_seq_no = g83100.sswi_event.event_seq_no
INNER JOIN g83100.sswi_camera_location ON g83100.sswi_event.event_seq_no = g83100.sswi_camera_location.event_seq_no
WHERE
g83100.sswi_metadata.metadata_name = 'TRAVEL_CORRIDOR_TYPE_CODE'")

#occasions are 7 days long
days_active_threshold <- 8

ppn_class_threshold <-  0.95

Turkey2 <- Turkey2week %>%
  filter(days_active >= days_active_threshold) %>%
  filter(prop_classified >= ppn_class_threshold) %>%
  mutate(across(matches("[A-Z]*_AMT", ignore.case = FALSE), 
                ~ifelse(.>0,1,0), .names = "{sub('_AMT', '_binary',col)}"),
         lat=scale(st_coordinates(.)[,2]),
         lon=scale(st_coordinates(.)[,1]),
         camera_version2=as.factor(ifelse(camera_version == "V4", 2, 1)))%>%
  rename(zone=turkey_mgmt_unit_id)%>%st_drop_geometry()
Turkey2$zone <- as.factor(Turkey2$zone)
Turkey2$camera_version <- as.factor(Turkey2$camera_version)
levels(Turkey2$camera_version)
Turkey2$cam_site_id <- as.factor(Turkey2$cam_site_id)
Turkey2$year <- Turkey2$season+2018
Turkey3 <- left_join(Turkey2, camsites)
Turkey3 <- Turkey3%>%mutate(across(matches("Forest"), scale))


nocc <- length(unique(Turkey2$occ))
knots <- list(occ = c(0.5, nocc+0.5))
nyears <- length(unique(Turkey2$season))

start <- Sys.time()
TurkeyBam <- mgcv::bam(TURKEY_AMT ~ zone + camera_version2 + s(lat, lon, k=60) + s(Forest_5000, k=12) +
                         s(season, k=nyears, by=zone) + s(occ, bs = "cc", k=nocc, by=zone) +
                         ti(season, occ, bs = c("tp", "cc")),
                       data = Turkey3,
                       family = ziP(),
                       knots = knots) #9:51
TurkeyBam2 <- mgcv::bam(TURKEY_AMT ~ zone + camera_version2  + 
                          s(season, k=nyears, by=zone) + s(occ, bs = "cc", k=nocc, by=zone) +
                          ti(season, occ, bs = c("tp", "cc")),
                        data = Turkey3,
                        family = nb(),
                        knots = knots) #9:51


end <- Sys.time()
end-start
##############################################################################################################################
####                                         aggregated to zone level                                                     ####
##############################################################################################################################
daterange <- Turkey2%>%group_by(season)%>% #recreate date ranges from data frame
  dplyr::summarise(start_date=as.Date(min(start_date)), end_date=as.Date(max(end_date)))%>%
  sf::st_drop_geometry()

df.byocc <- Turkey2 %>%
  filter(days_active >= 4) %>%
  filter(prop_classified >= 0.95) %>%
  mutate(across(matches("[A-Z]*_AMT", ignore.case = FALSE), ~ifelse(.>0,1,0), .names = "{sub('_AMT', '_binary',col)}")) %>%
  group_by(season,occ,zone) %>%
  dplyr::summarise(across(matches("[A-Z]*_AMT", ignore.case = FALSE), ~sum(.),.names = "{sub('_AMT', '_sum',col)}"),
                   across(matches("[A-Z]*_binary", ignore.case = FALSE), ~sum(.),.names = "{sub('_binary', '_occ',col)}"),
                   num.sites = dplyr::n(),
                   num.days = sum(days_active)) %>%
  mutate(across(matches("[A-Z]*_occ", ignore.case = FALSE), ~./num.sites,.names = "{sub('_occ', '_propocc',col)}"),
         across(matches("[A-Z]*_sum", ignore.case = FALSE), ~./num.days,
                .names = "{sub('_sum','_trigsperday',col)}"))%>%
  mutate(yearocc = paste0(lubridate::year(min(daterange$start_date))+(season-1),stringr::str_pad(occ, width=2, side="left", pad="0"))) %>%
  dplyr::arrange(yearocc) %>%
  group_by(yearocc) %>%
  mutate(time = dplyr::cur_group_id())

df.byocc.long = df.byocc %>%
  select(time, season,occ,zone,num.sites,
         matches("_occ")) %>%
  tidyr::pivot_longer(cols=matches("_occ"), names_pattern  = "(.*)_occ", names_to = "Spp")%>%
  mutate(year=season+2018) #this may need to be modified

df.byocc.long$zone <- as.factor(df.byocc.long$zone)
df.byocc.long$binomresponse <- with(df.byocc.long, cbind(value, num.sites - value))

nocc <- length(unique(df.byocc.long$occ))
knots <- list(occ = c(0.5, nocc+0.5))
nyears <- length(unique(df.byocc.long$season))


#model with year x occ interaction as well as occ x zone interaction
ZoneGam <- mgcv::gam(binomresponse ~ zone + s(season, k=nyears, by=zone) + s(occ, bs = "cc", k=nocc, by=zone) +
                   ti(season, occ, bs = c("tp", "cc")),
                 data = df.byocc.long,
                 family = binomial,
                 knots = knots)





summary(TurkeyBam)
summary(ZoneGam)
k.check(TurkeyBam)
r2(TurkeyBam)
r2(ZoneGam)
check_overdispersion(TurkeyBam)
check_zeroinflation(TurkeyBam)
check_model(TurkeyBam, check="pp_check", residual_type = "normal")
resids <- simulate_residuals(TurkeyBam)
model_performance(TurkeyBam)
gam.check(TurkeyBam)

gratia::draw(TurkeyBam)
variance_comp(TurkeyBam)
meansTurkeyBAM <- estimate_means(TurkeyBam, by = c("zone", "season"))%>%mutate(time=(nocc/2)+nocc*(season-1))
avgTurkeyBAM <- estimate_means(TurkeyBam, by = c("zone", "season"), estimate = "average")%>%mutate(time=(nocc/2)+nocc*(season-1))#takes a while
occTurkeyBAM <- estimate_means(TurkeyBam, by = c("zone", "season", "occ"))%>%mutate(time=occ+nocc*(season-1))
occTurkeyBAM2 <- estimate_means(TurkeyBam, by = c("zone", "season", "occ"), estimate = "average")#takes a while
occTurkeyBAM3 <- occTurkeyBAM2%>%mutate(time=occ+nocc*(season-1))
datagrid <- get_datagrid(TurkeyBam, by = c("zone", "season", "occ"), factors = "all", numerics = "all")
datagrid <- get_datagrid(TurkeyBam, by = c("occ", "zone", "season"), factors = "all", numerics = "all")

newdatayr <- expand.grid(season=unique(Turkey2$season), zone=unique(Turkey2$zone), camera_version=levels(Turkey2$camera_version), occ=nocc/2)%>%
            arrange(season, zone,camera_version)

occeffects <- paste("s(occ)", paste0("zone",unique(Turkey2$zone)), sep=":")
#predict just year trend
yrtrendm2y <- gratia::fitted_values(TurkeyBam, data=newdatayr, exclude=c(occeffects,"ti(season,occ)"), scale="response")%>%
  mutate(time=(nocc/2)+nocc*(season-1))


#predict whole model
newdataocc <- expand.grid(season=unique(Turkey2$season), zone=unique(Turkey2$zone), camera_version=levels(Turkey2$camera_version), occ=1:nocc)%>%
  arrange(season, zone,camera_version)
occtrendm2y <- gratia::fitted_values(TurkeyBam, data=newdataocc, scale = "response")%>%mutate(time=occ+nocc*(season-1))


speciesframe <- Turkey2
plottemp <-  ggplot() +
  facet_wrap(~zone)+
  geom_line(data=occtrendm2y, aes(x = time, y = .fitted, color=camera_version), lwd=0.5) +
  geom_line(data=yrtrendm2y, aes(x = time, y = .fitted, color=camera_version), lwd=2) +
  geom_pointrange(data=yrtrendm2y, aes(x= time, y= .fitted,ymin = .lower_ci, ymax = .upper_ci, color=camera_version), size=1, lwd=1) +
  labs(title=stringr::str_wrap("Weekly Turkey Detection Rate Turkey (triggers/week)",75),
       y = "Proportion of sites",
       x = "Time",
       subtitle = sprintf("Year Round, %s - %s", min(speciesframe$year), max(speciesframe$year))) +
  geom_vline(xintercept=seq(1,(nocc+1)*nyears,nocc)) +
  scale_x_continuous(labels = seq(min(speciesframe$year),max(speciesframe$year),1), breaks = seq(nocc/2,nocc*nyears,nocc)) +
  scale_color_brewer(palette = "Set2",
                     name = "Camera Version",
                     labels = levels(speciesframe$camera_version)) +
  scale_fill_brewer(palette = "Set2",
                    name = "Camera Version",
                    labels = levels(speciesframe$camera_version))
plottemp



ggplot() +
  geom_line(data=occTurkeyBAM, aes(x = time, y = Mean, color=zone), lwd=0.5) +
  geom_line(data=avgTurkeyBAM, aes(x = time, y = Mean, color=zone), lwd=2) +
  geom_pointrange(data=avgTurkeyBAM, aes(x= time, y= Mean ,ymin = CI_low, ymax = CI_high, color=zone), size=1, lwd=1) +
  labs(title=stringr::str_wrap("Weekly Turkey Detection Rate Turkey (triggers/week)",75),
       y = "Avg. Turkey Triggers/week",
       x = "Time",
       subtitle = sprintf("Year Round, %s - %s", 2019, 2024)) +
  geom_vline(xintercept=seq(1,(nocc+1)*nyears,nocc)) +
  scale_x_continuous(labels = seq(2019,2024,1), breaks = seq(nocc/2,nocc*nyears,nocc)) +
  scale_color_brewer(palette = "Set2",
                     name = "Zone",
                     labels = levels(Turkey2$zone)) +
  scale_fill_brewer(palette = "Set2",
                    name = "Zone",
                    labels = levels(Turkey2$zone))

ggplot() +
  facet_wrap(~zone)+
  geom_line(data=occTurkeyBAM, aes(x = time, y = Mean, color=zone), lwd=0.5) +
  geom_line(data=avgTurkeyBAM, aes(x = time, y = Mean, color=zone), lwd=2) +
  geom_pointrange(data=avgTurkeyBAM, aes(x= time, y= Mean ,ymin = CI_low, ymax = CI_high, color=zone), size=1, lwd=1) +
  labs(title=stringr::str_wrap("Weekly Turkey Detection Rate Turkey (triggers/week)",75),
       y = "Avg. Turkey Triggers/week",
       x = "Time",
       subtitle = sprintf("Year Round, %s - %s", 2019, 2024)) +
  geom_vline(xintercept=seq(1,(nocc+1)*nyears,nocc)) +
  scale_x_continuous(labels = seq(2019,2024,1), breaks = seq(nocc/2,nocc*nyears,nocc)) +
  scale_color_brewer(palette = "Set2",
                     name = "Zone",
                     labels = levels(Turkey2$zone)) +
  scale_fill_brewer(palette = "Set2",
                    name = "Zone",
                    labels = levels(Turkey2$zone))
######################################################################################################
####                             functions                                                        ####
######################################################################################################
ForestProp <- function(camsites, buffers, layer){
  # Wiscland 3 prop land cover
  lm_output <- 
    buffers %>% 
    set_names() %>% 
    # produce a dataframe after this is all done
    map_dfr( 
      ~sample_lsm(
        # raster layer
        landscape = layer,
        # camera locations
        y = camsites,
        # get landcover class level metrics
        level = "class",
        # return NA values for classes not in buffer
        # all_classes = TRUE, 
        # camera site IDs here
        plot_id = camsites$cam_site_id,
        # can do multiple metrics at once
        what = 'lsm_c_pland',
        # buffer sizes to use
        size = ., 
        # default is square buffer
        shape = "circle", 
        # turn warnings on or off
        verbose = FALSE 
      ), 
      # get buffer size column in the output
      .id = "buffer_size"
    )
  
  lm_output <- left_join(lm_output, WisclandGuide, by= join_by(class == dn.label))
  lm_output$label <- gsub(pattern = "\\W", replacement = "", x = lm_output$label)
  
  # in this data frame plot_id = camera ID
  # class = landcover type
  # value = % of that landcover type in the buffer
  # make each landcover type x buffer into a column
  lm_output <- 
    lm_output %>%
    # MAY NEED TO ADD distinct() HERE???
    #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
    pivot_wider(
      id_cols = plot_id,
      names_from = c(label, buffer_size),
      values_from = c(value),
      # give class 0 if it doesn't exist in buffer
      values_fill = 0
    ) %>%
    # clean up names
    rename(cam_site_id = plot_id)
  
  
  # join pland back to camera sf object
  forestLCs <- c("AspenPaperBirch", "RedMaple", "Oak", "CentralHardwoods", "NorthernHardwoods","AspenForestedWetland", "BottomlandHardwoods", "SwampHardwoods",
                 "MixedDeciduousConiferousForest", "MixedDeciduousConiferousForestedWetland")
  
  lm_output2 <- cbind(lm_output, lapply(buffers, function(x) rowSums(lm_output[,paste0(forestLCs, "_", x)])))
  forestcols <- grep(pattern = "c\\(", x = colnames(lm_output2))
  colnames(lm_output2)[forestcols] <- paste0("Forest", "_",buffers)
  
  ForestProp <- lm_output2[,c(1,forestcols)]
  camsites <- left_join(camsites, ForestProp, by="cam_site_id")
}
