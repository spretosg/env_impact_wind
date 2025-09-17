## Analysis script
packages <- c(
  "terra",
  "dplyr",
  "ggplot2",
  "sf",
  "rstatix",
  "cowplot",
  "viridis",
  "stringr",
  "forcats",
  "lme4",
  "tidyr",
  "classInt",
  "car"
  
)

# Install missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org/")
}

# Load all packages
lapply(packages, library, character.only = TRUE)

## Source helper functions
source("analysis/analysis_fcts.R")


### color schema for publication
cols_bio <- colorRampPalette(c("#F2CFEE", "#A02B93"))(5)
## es cols
cols_es <- colorRampPalette(c("#F6C6AD", "#E97132"))(5)

## load basic data
path_base<-"data/01_base_data"
stud_site<-st_read(paste0(path_base,"/stud_area_fin.gpkg"))
wt<-st_read(paste0(path_base,"/wind_turbines.gpkg"))

wp<-st_convex_hull(st_union(wt))
#wp cent
wp_cent<-st_centroid(wp)

## convert sf
wp<-st_sf(wp)

######### input PDF files
folder_path <- here::here("data/03_output_pdf")
# List all .tif files in the folder
PDF_D_files <- list.files(paste0(folder_path,"/PDF_ES/disturbance"), pattern = "\\.tif$", full.names = TRUE)
PDF_LO_files <- list.files(paste0(folder_path,"/PDF_ES/land_occupation"), pattern = "\\.tif$", full.names = TRUE)

## biodiversity PDF
PDF_BIO_D_files <- list.files(paste0(folder_path,"/PDF_BIO/disturbance"), pattern = "\\.tif$", full.names = TRUE)
PDF_BIO_HL_files <- list.files(paste0(folder_path,"/PDF_BIO/habitat_loss"), pattern = "\\.tif$", full.names = TRUE)
PDF_BIO_B_files <- list.files(paste0(folder_path,"/PDF_BIO/barrier"), pattern = "\\.tif$", full.names = TRUE)
PDF_BIO_C_files <- list.files(paste0(folder_path,"/PDF_BIO/collision"), pattern = "\\.tif$", full.names = TRUE)
bio_pdf<-terra::rast(here::here("data/03_output_pdf/PDF_BIO/lca_region_cumulative_scaled_total.tif"))


pdf_bio_d<-rast(PDF_BIO_D_files)
pdf_bio_hl<-rast(PDF_BIO_HL_files)
pdf_bio_b<-rast(PDF_BIO_B_files)
pdf_bio_c<-rast(PDF_BIO_C_files)

## just keep species group as names
names(pdf_bio_d) <- gsub("^lca_region_pdf_scaled_disturbance_", "", names(pdf_bio_d))
names(pdf_bio_hl) <- gsub("^lca_region_pdf_scaled_habitat_loss_", "", names(pdf_bio_hl))
names(pdf_bio_b) <- gsub("^lca_region_pdf_scaled_barrier_", "", names(pdf_bio_b))
names(pdf_bio_c) <- gsub("^lca_region_pdf_scaled_collision_", "", names(pdf_bio_c))

pdf_D<-rast(PDF_D_files)
pdf_LO<-rast(PDF_LO_files)
pdf_es_all<-c(pdf_D,pdf_LO)

## cummulative impact
sum_D<-sum(pdf_D)
sum_LO<-sum(pdf_LO)

w_lo<-(nlyr(pdf_LO)/(nlyr(pdf_D)+nlyr(pdf_LO)))
w_d<-(nlyr(pdf_D)/(nlyr(pdf_D)+nlyr(pdf_LO)))
#weighted sum across pathways and es
sum_all_pdf_es<-sum_D*w_d+sum_LO*w_lo

#resample to fit bio
sum_all_pdf_es <- resample(sum_all_pdf_es, bio_pdf, method = "bilinear") 
sum_all_pdf_es<-terra::mask(sum_all_pdf_es,bio_pdf)

#normalized ES pdf sum to compare with bio for NEP
nor_sum_pdf_es<-nor(sum_all_pdf_es)


##### 1. EIA statistics for drafted wind park
### 1.1 Avian diversity
facet_labels <- c(
  "Habitat loss (H)" = "2a Habitat loss (H)",
  "Disturbance (D)" = "2b Disturbance (D)",
  "Barrier (B)"     = "2c Barrier (B)",
  "Collision (C)"   = "2d Collision (C)"
)

pdf_bio_all<-read.csv(here::here("data/03_output_pdf/PDF_BIO/pdf_bio_all.csv"))
pdf_bio_all <- pdf_bio_all %>%
  mutate(
    pdf_scaled = (pdf_wp - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp)),
    Q5_scaled   = (wp_q5 - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp)),
    Q95_scaled  = (wp_q95 - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp)),
    nat_val = rep("Biodiversity",nrow(pdf_bio_all))
  )

pdf_bio_all <- pdf_bio_all %>%
  mutate(LCA_path = fct_recode(LCA_path,
                               "Habitat loss (H)" = "LO",
                               "Disturbance (D)" = "D",
                               "Barrier (B)" = "B",
                               "Collision (C)" = "C"
  ))

pdf_bio_all$LCA_path <- factor(
  pdf_bio_all$LCA_path,
  levels = c("Habitat loss (H)", "Disturbance (D)", "Barrier (B)", "Collision (C)")
)

## mean per lca path
facet_means <- pdf_bio_all %>%
  group_by(LCA_path) %>%
  summarise(mean = mean(pdf_wp, na.rm = TRUE)) %>%
  arrange(mean)

ggplot(pdf_bio_all, aes(x = pdf_wp, y = ID)) +
  geom_point(position = position_dodge(width = 0.9)) +
  # horizontal error bars
  geom_errorbarh(
    aes(xmin = wp_q5, xmax = wp_q95),
    position = position_dodge(width = 0.9),
    height = 0.25
  ) +
  # add mean line per facet
  geom_vline(
    data = facet_means,
    aes(xintercept = mean),
    linetype = "dashed",
    color = "blue"
  ) +
  xlab(bquote(PDF[Bio])) +
  labs(fill = "LCA impact pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4)
  ) +
  facet_wrap(~LCA_path, scales = "free_x", labeller = labeller(LCA_path = facet_labels))


## some statistics
leveneTest(pdf_wp ~ LCA_path, data = pdf_bio_all)

kruskal.test(pdf_wp ~ LCA_path, data = pdf_bio_all)
pairwise.wilcox.test(pdf_bio_all$pdf_wp, pdf_bio_all$LCA_path, p.adjust.method = "BH")
pdf_bio_all%>% wilcox_effsize(pdf_wp ~ LCA_path)

## differences in barrier effects
b<-pdf_bio_all%>%filter(LCA_path =="Barrier (B)")
leveneTest(pdf_wp ~ ID, data = b)

kruskal.test(pdf_wp ~ ID, data = b)
pairwise.wilcox.test(b$pdf_wp, b$LCA_path, p.adjust.method = "BH")
pdf_bio_all%>% wilcox_effsize(pdf_wp ~ LCA_path)

### 1.2 Ecosystem services

results_list <- list()

for (i in 1:nlyr(pdf_es_all)) {
  layer <- pdf_es_all[[i]]
  layer_name <- names(pdf_es_all)[i]
  layer_name <- sub("_D$", "", layer_name)
  layer_name <- sub("_LO$", "", layer_name)
  
  res <- terra::extract(layer, wp, fun = layer_stats)
  res_df <- as.data.frame(res)
  if(i %in% c(1,2,3)){
    res_df$LCA_path <- "D"
  }
  else{
    res_df$LCA_path <- "L"
  }
  # Rename columns with layer-specific names
  colnames(res_df) <- c("ID","_mean", "_q5","_q25","_q75", "_q95","LCA_path")
  res_df[1]<-layer_name
  results_list[[i]] <- unlist(res_df)
}

# Combine all into one data.frame (by column)
pdf_es_res <- do.call(rbind, results_list)
pdf_es_res<-as.data.frame(pdf_es_res)
colnames(pdf_es_res) <- c("ID","pdf_wp", "wp_q5","wp_q25","wp_q75", "wp_q95","LCA_path")
#pdf_es_all<-pdf_es_all%>%filter(pdf_wp>0)
pdf_es_res <- pdf_es_res %>%
  mutate(across(c(pdf_wp, wp_q5,wp_q25,wp_q75, wp_q95), as.numeric))

## rescale to min-max to compare it to bio
pdf_es_res <- pdf_es_res %>%
  mutate(
    pdf_scaled = (pdf_wp - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp)),
    Q5_scaled   = (wp_q5 - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp)),
    Q95_scaled  = (wp_q95 - min(pdf_wp)) / (max(pdf_wp) - min(pdf_wp))
  )
pdf_es_res$nat_val <- rep("Ecosystem services",nrow(pdf_es_res))


pdf_es_res <- pdf_es_res %>%
  mutate(ID = fct_recode(ID,
                         "Existence value of biodiversity" = "CULT_BIO",
                         "Place identity" = "CULT_IDEN",
                         "Spiritual experiences" = "CULT_PERC",
                         "Recreation" = "CULT_RECR",
                         "Education" = "CULT_SCIEN",
                         "Wild plants" = "PROV_PLANTS",
                         "Agricultural products" = "PROV_PROD",
                         "Wild animals" = "PROV_WILD",
                         "Climate regulation" = "REG_ATMO",
                         "Erosion control" = "REG_EROSION",
                         "Pest control" = "REG_PEST",
                         "Pollination" = "REG_POL",
                         "Habitat quality" = "REG_POP",
                         "Water quality" = "REG_WATQUAL"
  ))


pdf_es_res <- pdf_es_res %>%
  mutate(LCA_path = fct_recode(LCA_path,
                               "Land occupation (L)" = "L",
                               "Disturbance (D)" = "D"
  ))


facet_means <- pdf_es_res %>%filter(pdf_wp>0)%>%
  group_by(LCA_path) %>%
  summarise(mean = mean(pdf_wp), na.rm = TRUE)

ggplot(pdf_es_res, aes(x = pdf_wp, y = ID)) +
  geom_point(position = position_dodge(width = 0.9)) +
  # horizontal error bars
  geom_errorbarh(
    aes(xmin = wp_q5, xmax = wp_q95),
    position = position_dodge(width = 0.9),
    height = 0.25
  ) +
  # add mean line per facet
  geom_vline(
    data = facet_means,
    aes(xintercept = mean),
    linetype = "dashed",
    color = "blue"
  ) +
  xlab(bquote(PDF[Es])) +
  labs(fill = "LCA impact pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4)
  ) +
  facet_wrap(~LCA_path, scale = "free_x")


leveneTest(pdf_wp ~ LCA_path, data = pdf_bio_all)

kruskal.test(pdf_wp ~ ID, data = pdf_es_res%>%filter(LCA_path == "Disturbance (D)"))
pairwise.wilcox.test(pdf_es_res$pdf_wp, pdf_es_res$ID, p.adjust.method = "BH")
pdf_es_res%>% wilcox_effsize(pdf_wp ~ ID)


### 2. Jenks classified EIA
## 2.1 Avian diversity
stacks <- list(
  d  = pdf_bio_d,
  c  = pdf_bio_c,
  hl = pdf_bio_hl,
  b  = pdf_bio_b
)

# apply get_jenks_breaks to each layer of each stack
jenks_df_all <- bind_rows(
  lapply(names(stacks), function(stack_name) {
    s <- stacks[[stack_name]]
    bind_rows(lapply(1:nlyr(s), function(i) {
      get_jenks_breaks(
        r     = s[[i]],
        n     = 5,
        name  = names(s)[i],     # layer name
        group = stack_name,      # stack name passed explicitly
        wp    = wp
      )
    }))
  })
)


df_long <- jenks_df_all %>%
  pivot_longer(cols = starts_with("class"),
               names_to = "class",
               values_to = "range") %>%
  separate(range, into = c("min", "max"), sep = " - ") %>%
  mutate(min = as.numeric(min),
         max = as.numeric(max)) %>%
  mutate(class = fct_recode(class, "very low impact"="Class1",
                            "low impact"="Class2",
                            "medium impact"="Class3",
                            "high impact"="Class4",
                            "very high impact"="Class5"))%>%
  mutate(group = fct_recode(group, "Disturbance (D)"="d",
                            "Habitat loss (H)"="hl",
                            "Collison (C)"="c",
                            "Barrier (B)"="b"))

bio_stats<-df_long%>%
  # assign wt_class based on interval
  mutate(wt_class = ifelse(between(wt_mean, min, max), class, NA_character_)) %>%
  # keep only the matching row
  filter(!is.na(wt_class))

write.csv(bio_stats,"classif_bio.csv")


df_long$group <- factor(
  df_long$group,
  levels = c("Habitat loss (H)", "Disturbance (D)", "Barrier (B)", "Collison (C)")
)

ggplot(df_long, aes(y = factor(ID))) +
  # plot class intervals
  geom_segment(aes(x = min, xend = max, yend = factor(ID), color= factor(class)),
               size = 3) +
  geom_errorbarh(aes(xmin = wtQ25, xmax = wtQ75),
                 height = 0.2, color = "blue") +
  geom_point(aes(x = wt_mean), color = "blue", size = 2) +
  labs(x = "PDF", y = "") +
  scale_color_manual(values = cols_bio, name = "EIA impact categories") +
  theme_minimal(base_size = 14)+
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",         # Move legend to bottom
    legend.direction = "horizontal",    # Make legend horizontal
    legend.box = "horizontal"           # Spread legend items horizontally
  )+
  facet_wrap(~group,scales = "free_x")



## 2.2 Ecosystem services
pdf_LO <- resample(pdf_LO, bio_pdf, method = "bilinear")  # or "near" for categorical

# Clip to reference extent (if you want strict clipping)
pdf_LO <- crop(pdf_LO, bio_pdf)

# Optionally mask (to apply ref NA mask)
pdf_LO <- mask(pdf_LO, bio_pdf)


pdf_D <- resample(pdf_D, bio_pdf, method = "bilinear")  # or "near" for categorical

# Clip to reference extent (if you want strict clipping)
pdf_D <- crop(pdf_D, bio_pdf)

# Optionally mask (to apply ref NA mask)
pdf_D <- mask(pdf_D, bio_pdf)


rasters_es<-list(
  d=pdf_D,
  lo = pdf_LO
)

jenks_df_all_es <- bind_rows(
  lapply(names(rasters_es), function(stack_name) {
    s <- rasters_es[[stack_name]]
    bind_rows(lapply(1:nlyr(s), function(i) {
      get_jenks_breaks(
        r     = s[[i]],
        n     = 5,
        name  = names(s)[i],     # layer name
        group = stack_name,      # stack name passed explicitly
        wp    = wp
      )
    }))
  })
)


df_long_es <- jenks_df_all_es %>%
  pivot_longer(cols = starts_with("class"),
               names_to = "class",
               values_to = "range") %>%
  separate(range, into = c("min", "max"), sep = " - ") %>%
  mutate(min = as.numeric(min),
         max = as.numeric(max)) %>%
  mutate(ID = fct_recode(ID,
                         "Existence value of biodiversity" = "CULT_BIO",
                         "Place identity" = "CULT_IDEN",
                         "Spiritual experiences" = "CULT_PERC",
                         "Recreation" = "CULT_RECR",
                         "Education" = "CULT_SCIEN",
                         "Wild plants" = "PROV_PLANTS",
                         "Agricultural products" = "PROV_PROD",
                         "Wild animals" = "PROV_WILD",
                         "Climate regulation" = "REG_ATMO",
                         "Erosion control" = "REG_EROSION",
                         "Pest control" = "REG_PEST",
                         "Pollination" = "REG_POL",
                         "Habitat quality" = "REG_POP",
                         "Water quality" = "REG_WATQUAL"
  ))%>%mutate(group = fct_recode(group, "Disturbance (D)"="d",
                                 "Land occupation (L)"="lo"))%>%
  mutate(class = fct_recode(class, "very low impact"="Class1",
                            "low impact"="Class2",
                            "medium impact"="Class3",
                            "high impact"="Class4",
                            "very high impact"="Class5"))

es_stats<-df_long_es%>%
  # assign wt_class based on interval
  mutate(wt_class = ifelse(between(wt_mean, min, max), class, NA_character_)) %>%
  # keep only the matching row
  filter(!is.na(wt_class))

write.csv(es_stats,"classif_es.csv")




ggplot(df_long_es, aes(y = factor(ID))) +
  # plot class intervals
  geom_segment(aes(x = min, xend = max, yend = factor(ID), color= factor(class)),
               size = 3) +
  geom_errorbarh(aes(xmin = wtQ25, xmax = wtQ75),
                 height = 0.2, color = "blue") +
  geom_point(aes(x = wt_mean), color = "blue", size = 2) +
  labs(x = "PDF", y = "") +
  scale_color_manual(values = cols_es, name = "EIA impact categories") +
  theme_minimal(base_size = 14)+
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",         # Move legend to bottom
    legend.direction = "horizontal",    # Make legend horizontal
    legend.box = "horizontal"           # Spread legend items horizontally
  )+facet_wrap(~group, scales = "free_x")

######## 3. Regional PDF
bio_pdf_nor<-nor(bio_pdf,na.rm = T)

# Extract values
vals_bio <- values(bio_pdf_nor, mat=FALSE)

# Compute Jenks natural breaks (5 classes)
jenks_bio <- classIntervals(vals_bio, n=5, style="jenks")

# Reclassify raster based on Jenks breaks
r_class_bio <- classify(bio_pdf_nor, rbind(
  c(-Inf, jenks_bio$brks[2], 1),
  c(jenks_bio$brks[2], jenks_bio$brks[3], 2),
  c(jenks_bio$brks[3], jenks_bio$brks[4], 3),
  c(jenks_bio$brks[4], jenks_bio$brks[5], 4),
  c(jenks_bio$brks[5], Inf, 5)
))


vals_es <- values(nor_sum_pdf_es, mat=FALSE)

# Compute Jenks natural breaks (5 classes)
jenks_es <- classIntervals(vals_es, n=5, style="jenks")

# Reclassify raster based on Jenks breaks
r_class_es <- classify(nor_sum_pdf_es, rbind(
  c(-Inf, jenks_es$brks[2], 1),
  c(jenks_es$brks[2], jenks_es$brks[3], 2),
  c(jenks_es$brks[3], jenks_es$brks[4], 3),
  c(jenks_es$brks[4], jenks_es$brks[5], 4),
  c(jenks_es$brks[5], Inf, 5)
))

## municipalities of study area
com<-st_read(here::here("data/01_base_data/municipalities.gpkg"))
com<-sf::st_transform(com,crs(wp))
com<-com%>%st_intersection(stud_site)

r_section_df_bio <- as.data.frame(r_class_bio, xy = TRUE, na.rm = TRUE)


r_sec_long_bio <- tidyr::pivot_longer(r_section_df_bio, cols = all_of(names(r_class_bio)),
                                      names_to = "Layer", values_to = "bio_impact")


r_section_df_es <- as.data.frame(r_class_es, xy = TRUE, na.rm = TRUE)


r_sec_long_es <- tidyr::pivot_longer(r_section_df_es, cols = all_of(names(r_class_es)),
                                     names_to = "Layer", values_to = "es_impact")



## EIA categories ES
labels <- paste0(
  c("very low","low","moderate","high","very high"), ": ",
  sprintf("%.2f", jenks_es$brks[-length(jenks_es$brks)]), " – ",
  sprintf("%.2f", jenks_es$brks[-1])
)

ggplot() +
  geom_raster(
    data = r_sec_long_es,
    aes(x = x, y = y, fill = factor(es_impact))
  ) +
  scale_fill_manual(
    values = cols_es,
    na.value = "transparent",
    name = bquote(paste("Cumulative PDF","*")[ES], 12),
    labels = labels
  ) +
  geom_sf(data = wp, color = "blue", size = 10, linetype = "dashed", alpha = 0.2) +
  geom_sf(data = com, fill = NA, color = "black", size = 4) +
  #geom_sf(data = stud_small, fill = NA, color= "red")+
  coord_sf()+
  geom_sf_text(data = com, aes(label = kommunenavn), color = "black", size = 4,
               nudge_x = 0.2, nudge_y = 0.1) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())

## EIA categories Bio
labels <- paste0(
  c("very low","low","moderate","high","very high"), ": ",
  sprintf("%.2f", jenks_bio$brks[-length(jenks_bio$brks)]), " – ",
  sprintf("%.2f", jenks_bio$brks[-1])
)

ggplot() +
  geom_raster(
    data = r_sec_long_bio,
    aes(x = x, y = y, fill = factor(bio_impact))
  ) +
  scale_fill_manual(
    values = cols_bio,
    na.value = "transparent",
    name = bquote(paste("Cumulative PDF","*")[Bio], 12),
    labels = labels
  ) + 
  geom_sf(data = wp,  color = "blue", size = 10, linetype = "dashed", alpha = 0.2) +
  geom_sf(data = com, fill = NA, color = "black", size = 4) +
  #geom_sf(data = stud_small, fill = NA, color= "red")+
  coord_sf()+
  geom_sf_text(data = com, aes(label = kommunenavn), color = "black", size = 4,
               nudge_x = 0.2, nudge_y = 0.1) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())


## some area statistics
r_vals<-values(nor_sum_pdf_es)
r_vals <- r_vals[!is.na(r_vals)]  # remove NA if any

# total area covered
tot_area <-res(nor_sum_pdf_es)^2 * length(r_vals)/10^6

area_Q50 <-res(nor_sum_pdf_es)^2 * length(r_vals[r_vals<=.5])/10^6
area_Q50_p <-res(nor_sum_pdf_es)^2 * length(r_vals[r_vals>.5])/10^6

# Create data frame
df <- data.frame(value = r_vals)

# Plot histogram with vertical quantile lines
ggplot(df, aes(x = value)) +
  geom_histogram(binwidth = 0.1, fill = "grey", color = "black", alpha = 0.7) +
  xlab(bquote(paste(PDF,"*")[Es]))+
  ylab(bquote(paste(Area, " km")^2))+
  theme_minimal()

##### 4. NEP
nep <- 1-(bio_pdf_nor+nor_sum_pdf_es)/2

## 1. sample nep in wp
res_nep <- terra::extract(nep, wp, fun = layer_stats)



r_section_df <- as.data.frame(nep, xy = TRUE, na.rm = TRUE)


r_sec_long <- tidyr::pivot_longer(r_section_df, cols = all_of(names(nep)),
                                  names_to = "Layer", values_to = "env_imp")

nep_vals<-values(nep)


## jenks classes in low to high env perf...
jenks_nep <- classIntervals(nep_vals, n=5, style="jenks")

# Reclassify raster based on Jenks breaks
r_class_nep <- classify(nep, rbind(
  c(-Inf, jenks_nep$brks[2], 1),
  c(jenks_nep$brks[2], jenks_nep$brks[3], 2),
  c(jenks_nep$brks[3], jenks_nep$brks[4], 3),
  c(jenks_nep$brks[4], jenks_nep$brks[5], 4),
  c(jenks_nep$brks[5], Inf, 5)
))



r_section_df_nep <- as.data.frame(r_class_nep, xy = TRUE, na.rm = TRUE)


r_sec_long_nep <- tidyr::pivot_longer(r_section_df_nep, cols = all_of(names(r_class_nep)),
                                      names_to = "Layer", values_to = "nep_class")



## NEP categories
labels <- paste0(
  c("very low","low","moderate","high","very high"), ": ",
  sprintf("%.2f", jenks_nep$brks[-length(jenks_nep$brks)]), " – ",
  sprintf("%.2f", jenks_nep$brks[-1])
)

greens <- c("#E5F5E0", "#A1D99B", "#41AB5D", "#238B45", "#005A32")
ggplot() +
  geom_raster(
    data = r_sec_long_nep,
    aes(x = x, y = y, fill = factor(nep_class))
  ) +
  scale_fill_manual(
    values = greens,
    na.value = "transparent",
    name = "Environmental performance",
    labels = labels
  ) +
  geom_sf(data = wp, color = "blue", size = 10, linetype = "dashed", alpha = 0.2) +
  geom_sf(data = com, fill = NA, color = "black", size = 4) +
  coord_sf()+
  geom_sf_text(data = com, aes(label = kommunenavn), color = "black", size = 4,
               nudge_x = 0.2, nudge_y = 0.1) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())


######## 5. LCOE and NEP
lcoe<-rast(here::here("data/04_wind_data/LCOE.tif"))

lcoe<-terra::project(lcoe,crs(bio_pdf_nor))
lcoe <- resample(lcoe, bio_pdf_nor, method = "bilinear") 
lcoe<-terra::mask(lcoe,bio_pdf_nor)


r_section_df <- as.data.frame(lcoe, xy = TRUE, na.rm = TRUE)


r_sec_long <- tidyr::pivot_longer(r_section_df, cols = all_of(names(lcoe)),
                                  names_to = "Layer", values_to = "lcoe")






### trade-off map wind NEP
r_stack <- c(nep,lcoe)

# Convert to data frame
r_tradeoff <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)
colnames(r_tradeoff)[3:4] <- c("NEP", "LCOE")

# Define quantile-based categories (or use custom thresholds)
r_tradeoff <- r_tradeoff %>%
  mutate(
    nep_cat = case_when(
      NEP >= quantile(NEP, .5, na.rm = TRUE) ~ "High NEP",
      NEP < quantile(NEP, .5, na.rm = TRUE) ~ "Low NEP",
      TRUE ~ "Mid NEP"
    ),
    lcoe_cat = case_when(
      LCOE >= quantile(LCOE, .5, na.rm = TRUE) ~ "High LCOE",
      LCOE < quantile(LCOE, .5, na.rm = TRUE) ~ "Low LCOE",
      TRUE ~ "Mid LCOE"
    ),
    bivar_class = paste(nep_cat, lcoe_cat, sep = " / ")
  )

# Define bivariate color scale
bivar_colors <- c(
  "Low NEP / High LCOE" = "#f46d43",  # medium red (bad env, expensive)
  "Low NEP / Low LCOE"  = "#fde0dd",  # very light red (bad env, cheap)
  "High NEP / High LCOE"= "#d9f0d3",  # very light green (good env, expensive)
  "High NEP / Low LCOE" = "#66c2a4"   # medium green (good env, cheap)
)



# Plot
ggplot() +
  geom_raster(data = r_tradeoff, aes(x = x, y = y, fill = bivar_class)) +
  scale_fill_manual(values = bivar_colors,
                    na.value = "transparent",
                    name = "Trade-off LCOE and NEP") +
  geom_sf(data = wp, , color = "blue", size = 8, linetype = "dashed", alpha = 0.2) +
  geom_sf(data = com, fill = NA, color = "black", size = 4) +
  #geom_sf(data = stud_small, fill = NA, color= "red")+
  coord_sf()+
  geom_sf_text(data = com, aes(label = kommunenavn), color = "black", size = 4,
               nudge_x = 0.2, nudge_y = 0.1) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())+
  theme(legend.position = "right")



########### Appendix plots
### A3 Mean Ecosystem service capacity maps
es<-z_coef%>%distinct(es)

## es capacity maps
raster_list <- list()

for(i in 1:14) {
  es_name <- es[[i, 1]]  # extract character value from tibble
  ES_cap_files <- list.files(here::here(paste0("data/02_ecosystem_service_modelling/ind_maps/", es_name)), 
                             pattern = "\\.tif$", full.names = TRUE)
  
  rast_es <- rast(ES_cap_files)
  raster_list[[es_name]] <- mean(rast_es)  # name the list entry by ES
}

# Combine into a SpatRaster stack
raster_stack <- rast(raster_list)
raster_stack <- resample(raster_stack, bio_mask, method = "bilinear") 
raster_stack<-terra::mask(raster_stack,bio_mask)
raster_stack <- project(raster_stack, "EPSG:4326")  # reproject to lon/lat

es_labels <- c(
  "CULT_BIO"   = "Existence value of biodiversity",
  "CULT_IDEN"  = "Place identity",
  "CULT_PERC"  = "Spiritual experiences",
  "CULT_RECR"  = "Recreation",
  "CULT_SCIEN" = "Education",
  "PROV_PLANTS"= "Wild plants",
  "PROV_PROD"  = "Agricultural products",
  "PROV_WILD"  = "Wild animals",
  "REG_ATMO"   = "Climate regulation",
  "REG_EROSION"= "Erosion control",
  "REG_PEST"   = "Pest control",
  "REG_POL"    = "Pollination",
  "REG_POP"    = "Habitat quality",
  "REG_WATQUAL"= "Water quality"
)

# Rename raster layers based on the lookup
names(raster_stack) <- es_labels[names(raster_stack)]

raster_df <- as.data.frame(raster_stack, xy = TRUE, na.rm = FALSE) %>%
  tidyr::pivot_longer(-c(x, y), names_to = "Service", values_to = "Value")

# Plot facetted maps
ggplot(raster_df) +
  geom_raster(aes(x = x, y = y, fill = Value)) +
  coord_equal() +
  scale_fill_viridis(option = "magma", na.value = "transparent") +
  facet_wrap(~ Service, ncol = 3) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())+
  theme(legend.position = "bottom")+
  
  guides(fill = guide_colorbar(title = "Mean ecosystem service capacity"))

### A4 sensitivity z coefficient
z_coef<-z_coef %>%
  mutate(ID = fct_recode(es,
                         "Existence value of biodiversity" = "CULT_BIO",
                         "Place identity" = "CULT_IDEN",
                         "Spiritual experiences" = "CULT_PERC",
                         "Recreation" = "CULT_RECR",
                         "Education" = "CULT_SCIEN",
                         "Wild plants" = "PROV_PLANTS",
                         "Agricultural products" = "PROV_PROD",
                         "Wild animals" = "PROV_WILD",
                         "Climate regulation" = "REG_ATMO",
                         "Erosion control" = "REG_EROSION",
                         "Pest control" = "REG_PEST",
                         "Pollination" = "REG_POL",
                         "Habitat quality" = "REG_POP",
                         "Water quality" = "REG_WATQUAL"
  ))



summary_df <- z_coef %>%
  group_by(ID, area, logA) %>%
  summarise(
    mean_logI = mean(logI, na.rm = TRUE),
    sd_logI = sd(logI, na.rm = TRUE),
    n = n(),
    se_logI = sd_logI / sqrt(n),
    .groups = "drop"
  )

z_summary <- z_coef %>%
  group_by(ID) %>%
  summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

summary_df <- summary_df %>%
  left_join(z_summary, by = "ID")




ggplot(summary_df, aes(x = logA, y = mean_logI)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean_logI - se_logI, ymax = mean_logI + se_logI),
              alpha = 0.2, linetype = 0) +
  
  labs(
    x = "log(Area)",
    y = "Mean log(invariability ES)",
    title = "",
  ) +
  theme_minimal(base_size = 16) +
  facet_wrap(~ID)


## stats on z
model <- lmer(z ~ es + (1 | run), data = z_coef)

summary(model)


## A8 LCOE map
ggplot() +
  geom_raster(
    data = r_sec_long,
    aes(x = x, y = y, fill = lcoe)
  ) +
  scale_fill_continuous(
    na.value = "transparent",
    name = "LCOE [ø/kwh]",
    labels = labels
  ) +
  geom_sf(data = wp, color = "blue", size = 10, linetype = "dashed", alpha = 0.2) +
  geom_sf(data = com, fill = NA, color = "black", size = 4) +
  coord_sf()+
  geom_sf_text(data = com, aes(label = kommunenavn), color = "white", size = 4,
               nudge_x = 0.2, nudge_y = 0.1) +
  xlab("") +
  ylab("")+
  theme_cowplot(font_size = 14) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.background = element_blank())

## A9 NEP - LCOE relation
# Correlation (quantitative measure of synergy/tradeoff)
cor_val <- cor(r_tradeoff$NEP, r_tradeoff$LCOE, method="spearman") 
cor.test(r_tradeoff$NEP, r_tradeoff$LCOE, method="spearman")
print(cor_val)

ggplot(r_tradeoff, aes(x = NEP, y = LCOE)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = paste("Raster Tradeoff/Synergy (Spearman r =", round(cor_val, 2), ")"),
       x = "NEP", y = "LCOE")+
  theme_minimal()


nor_lcoe<-nor(lcoe)

tradeoff <- (nor_lcoe - nep) / (nor_lcoe + nep)

# Handle division by zero (optional)
tradeoff[is.nan(tradeoff)] <- NA
