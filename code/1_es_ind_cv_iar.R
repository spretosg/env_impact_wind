### Local stakeholders rated ecosystem capacity for each land use land cover (LULC) class
### This script transforms the ratings to a  ES capacity map, using the modified LULC map (see step 0)

### input: es_rating_participants.xls participants (1-14) ratings, LULC_pred (adjusted landcover map) 
### output: individual ES capacity maps, CV and IAR raster maps for each ES stored in "/data/2_ecosystem_service_modelling/..."

packages <- c(
  "terra",
  "sf",
  "ggplot2",
  "dplyr",
  "readxl",
  "purrr"
)

# Install missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org/")
}

# Load all packages
lapply(packages, library, character.only = TRUE)



# 0.global parameters
resolution<-250 ## resolution of the es mapping (250m cell size)

## reclassification of the LULC types to match IMAGINE questionnaire (xls data)
transf<-c("Dyrket mark (fulldyrka + overflatedyrka)",6,
          "Innmarksbeite",11,
          "Mindre grøntområder, lekeplasser, rekreasjonsområder",1,
          "Grupper av trær (f.eks. alléer, treklynger, åkerøyer)",12,
          "Løvskog",3,
          "Skrenter",2,
          "Innsjø",4,
          "Dam",4,
          "Bekk",4,
          "Elv",4,
          "Kystsone",10,
          "Elvesone",13,
          "Gran- og blandingsskog",3,
          "Åpen mark",8,
          "Myr",9
)


transf<-matrix(transf, ncol=2, byrow=TRUE)
transf<-as.data.frame(transf)
colnames(transf)<-c("Landskapstype","rast_val")

id<-c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13","ID14")

lulc_rast<-terra::rast(here::here("data/02_ecosystem_service_modelling/LULC_pred.tif"))
stud_site<-here::here("data/01_base_data/stud_area_large.gpkg")


stud_site<-st_read(stud_site)
#vector format for clipping raster
p <- vect(stud_site)

bbox<-as.data.frame(as.matrix(sf::st_bbox(stud_site)))

## standard raster high resolution
r1<-rast(res = resolution, xmin = bbox$V1[1], xmax = bbox$V1[3], ymin = bbox$V1[2], ymax = bbox$V1[4],crs=st_crs(stud_site)$wkt)
out_ind_path<-here::here("data/02_ecosystem_service_modelling/ind_maps")

for(r in 1:length(id)){
  print(id[r])
  ## attach coefficient of variation values from IMAGINE project to further calculate IAR of ES
  es_lt<-read_excel(here::here("data/02_ecosystem_service_modelling/es_rating_participants.xlsx"), sheet = id[r],skip = 2)
  es_lt<-es_lt %>%  filter(!row_number() %in% c(16:nrow(es_lt)))
  es_lt<-es_lt[1:(length(es_lt)-5)]
  lt_types<-unique(es_lt$Landskapstype)
  es<-colnames(es_lt)[-1]
  es<-es[1:length(es)-1]
  
  es_lt<-left_join(es_lt,transf,by="Landskapstype")
  es_lt<-as.data.frame(es_lt)

  for(t in 1:length(es)){
    es_tmp<-es[t]
    print(paste0("Calculate next es: ", es_tmp))
    a<-as.numeric(unlist(es_lt%>%dplyr::select(all_of(es_tmp))))
    df_val<-as.data.frame(cbind(es_lt$rast_val,a))
    colnames(df_val)<-c("lt_class","es_val_ind")
    
    df_val$es_val_ind<-as.numeric(df_val$es_val_ind)
    # df_val$mean_val<-as.numeric(df_val$mean_val)
    df_val<-df_val%>%group_by(lt_class)%>%summarise(es_val = mean(es_val_ind))
    all<-as.data.frame(as.character(c(1:16)))
    colnames(all)<-"lt_class"
    all<-left_join(all,df_val,by="lt_class")
    
    # Build a matrix "from-to" this comes per ES from the IMAGINE project manual
    rcl_es <- as.matrix(data.frame(from = as.integer(all$lt_class) , to = all$es_val))
    
    
    es_rast<-terra::classify(lulc_rast,rcl_es)
    es_rast<-terra::project(es_rast,st_crs(stud_site)$wkt)
    
    #resampling to standard raster
    es_rast <- resample(es_rast, r1)
    es_rast[is.na(es_rast)] <- 0
    
    out_folder<-paste0(out_ind_path,"/",es_tmp)
    if (!dir.exists(out_folder)) {
      dir.create(out_folder, recursive = TRUE)
    }
    writeRaster(es_rast,paste0(out_folder,"/",id[r],".tif"),overwrite=TRUE)

  }
  
  
}


### calculate pixel wise IAR for each es across raters
## read just for es matrix
es_lt<-read_excel(here::here("data/02_ecosystem_service_modelling/es_rating_participants.xlsx"), sheet = id[1],skip = 2)
es_lt<-es_lt %>%  filter(!row_number() %in% c(16:nrow(es_lt)))
es_lt<-es_lt[1:(length(es_lt)-5)]
lt_types<-unique(es_lt$Landskapstype)
es<-colnames(es_lt)[-1]
es<-es[1:length(es)-1]
rm(es_lt)


for (n  in 1:length(es)) {
  folder_path<-paste0("data/02_ecosystem_service_modelling/ind_maps/",es[n])
  # List all raster files (e.g., .tif)
  raster_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

  # Load all rasters as a SpatRaster list and then stack
  rstack <- rast(raster_files)
  mean_es<-mean(rstack,na.rm=T)
  sds_es <- app(rstack, sd, na.rm = TRUE)
  cov_es <- sds_es/mean_es
  IAR_es <- 1/(cov_es^2)
  plot(cov_es)
  writeRaster(IAR_es,paste0("data/02_ecosystem_service_modelling/iar/",es[n],".tif"))
  writeRaster(mean_es,paste0("data/02_ecosystem_service_modelling/mean/",es[n],".tif"))
  writeRaster(mean_es,paste0("data/02_ecosystem_service_modelling/cv/",es[n],".tif"))
}
