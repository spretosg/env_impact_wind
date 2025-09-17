### CODE TO CALCULATE PDF FOR IMPACT PATHWAYS D/LO OF WIND ENERGY ON ES
### in a final step the PDF maps are cropped to the final analysis area to avoid border effects for analysis

### input: windturbines, z_coef.csv (see 2_z_coef.R), CV and IAR rasters for each ES (see 1_es_ind_cv_iar.R)
### output: impact PDF maps for each ES for land occupation (LO) and if available for disturbance (D)


packages <- c(
  "terra",
  "sf",
  "ggplot2",
  "dplyr",
  "readxl",
  "stars"
)

# Install missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org/")
}

# Load all packages
lapply(packages, library, character.only = TRUE)


source("code/functions.R")

# load wind turbines
wt<-st_read("data/07_base_data/wind_turbines.gpkg")

# crop for final study area to avoid border effects
fin_area<-st_read("data/0_base/stud_area_fin.gpkg")

# load z
z_val<-as.data.frame(read.csv("data/02_ecosystem_service_modelling/z_coef_sens.csv"))
z_val<-z_val%>%group_by(es)%>%summarise(z = mean(z,na.rm=T))


##0. wt parameters
#aEP is area lost per km2 and MW
aEP<-0.01 # from Denholm et al (0.003 for permanent plus 0.007 for temporary losses)


##1. Area lost land occupation LO
## EP is capacity of WT (e.g. 3MW) according to project description
EP <-3.075 

## area lost due to occupation in km2
A_lost_LO<-aEP*EP

#param list new
param <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(param)<-colnames(z_val)


IAR_cube<-read_stars(list.files(here::here("data/02_ecosystem_service_modelling/iar"),full.names = T, pattern = "\\.tif$"))


PDF_cube<-IAR_cube

PDF<-list()

### calculate scaling factors per ES for the later summing of pdfs
nlayers <- length(IAR_cube)

# Compute sum for each raster layer
# raster_sums <- sapply(1:nlayers, function(i) {
#   sum(IAR_cube[[i]],na.rm=T)
# })
# 
# raster_sums<-raster_sums/sum(raster_sums)

## ecosystem services that in addition follow a LCA disturbance pathway
disturbance_es<-c("CULT_RECR","CULT_IDEN", "CULT_PERC")


for(t in 1: length(IAR_cube)){
  tmp_IAR<-IAR_cube[t]
  tmpes<-gsub("\\..*","",names(tmp_IAR))
  tmp_z<-as.numeric(z_val%>%filter(es == tmpes)%>%dplyr::select(z))
  print(paste0("Calculate next es: ", tmpes))
  
  #calculate the wind farm area, given the wind turbines wt and the area lost land occupation
  tmp_A_LO<-wind_farm_polygon(wt,A_lost_LO)
  
  if(tmpes %in% disturbance_es){
    tmp_A_lost_D<-A_lost_D(tmpes = tmpes, iar_rast = tmp_IAR, q_min = 0, q_max = 0.8, store_values = T, resolution = 1000)
    tmp_A_dist<-wind_farm_polygon(wt,tmp_A_lost_D)
    # moving window disturbance
    PDF_D<-pdf_sliding_window(tmp_IAR,tmp_A_dist,tmp_z,F)
    # crop to final area
    PDF_D<-st_crop(PDF_D, fin_area)
    #PDF_D<-nor(PDF_D,na.rm=T)
    file<- here::here("data", "03_output_pdf/PDF_ES/disturbance", paste0(tmpes,".tif"))
    write_stars(PDF_D,file)
  }
  
  
  ## moving window PDF land occupation
  PDF_LO<-pdf_sliding_window(tmp_IAR,tmp_A_LO,tmp_z,F)
  # crop to final area
  PDF_LO<-st_crop(PDF_LO, fin_area)
  file_LO<- here::here("output", "03_output_pdf/PDF_ES/land_occupation", paste0(tmpes,".tif"))
  write_stars(PDF_LO,file_LO)

}


