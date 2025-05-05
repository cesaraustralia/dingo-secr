library(tidyverse)
library(secr)
library(seizer)
library(sf)
library(terra)
library(tidyterra)

source("code/functions.R")

covars <-
  rast("spatialData/covs_stack.tif") %>%
  setNames(c("rds_dist", "water_dist", "late_GS", "shdi"))

hab <- rast("C:/GIS/NVIS_V7_0_AUST_RASTERS_EXT_ALL/NVIS_V7_0_AUST_EXT.gdb")
bd <-
  read_csv("data/dingo_lookup.csv") %>%
  filter(region == "Mallee") %>%
  mutate(Session = year(dmy(collection_date)),
         M = month(dmy(collection_date))) %>%
  filter((Session == 2023 & M %in% c(10, 11)) |
           (Session == 2024 & M %in% c(3, 5))) %>%
  st_as_sf(coords = c("dec_long", "dec_lat"), crs = "WGS84") %>%
  st_transform(crs(hab)) %>%
  vect()

hab <- (hab %>% crop(covars %>% project(crs(hab))))[[1]]

hab_low <-
  hab %>%
  aggregate(fact = 25, fun = "modal")

covars <- covars %>% project(hab) %>% resample(hab)
covars[[1]] <- log(covars[[1]])
covars[[2]] <- log(covars[[2]])

roads <-
  read_sf("spatialData/roads_BD.shp") %>%
  st_transform(crs(hab)) %>%
  vect() %>%
  crop(hab) %>%
  st_as_sf()

bigdesert <-
  read_sf("spatialData/Big_Desert.shp") %>%
  filter(NAME == "Wyperfeld") %>%
  st_transform(crs(hab))

# bd_mask <-
#   hab %>%
#   as.polygons() %>%
#   filter(!NVIS7_0_AUST_EXT_MVG_ALB %in%
#            c("Cleared, non-native vegetation, buildings",
#              "Inland aquatic - freshwater, salt lakes, lagoons")) %>%
#   st_as_sf() %>%
#   summarise() %>%
#   st_cast("POLYGON")%>%
#   mutate(area = st_area(geometry)) %>%
#   arrange(desc(area)) %>%
#   slice(1L) %>%
#   st_buffer(100) %>%
#   st_make_valid() %>%
#   vect()

bd_mask <-
  covars[[1]] %>%
  as.polygons() %>%
  st_as_sf() %>%
  summarise() %>%
  st_cast("POLYGON")%>%
  mutate(area = st_area(geometry)) %>%
  vect()

ind <-
  covars[[1]]

values(ind) <-
  1:(dim(ind)[1]*dim(ind)[2])

ind <-
  ind %>%
  mask(bd_mask)

tdf_full <-
  read_csv("data/dingo_lookup.csv") %>%
  filter(region == "Mallee") %>%
  mutate(Session = year(dmy(collection_date)),
         M = month(dmy(collection_date))) %>%
  filter(Session == 2023 & M %in% c(10, 11) |
           Session == 2024 & M %in% c(3, 5)) %>%
  dplyr::select(site, dec_lat, dec_long, Session, M) %>%
  st_as_sf(coords = c("dec_long", "dec_lat"), crs = "WGS84") %>%
  st_transform(crs(hab)) %>%
  st_intersection(., st_as_sf(bd_mask)) %>%
  dplyr::select(site, Session, M) %>%
  mutate(Detector = extract(ind, vect(.))[,2],
         X = extract(ind, vect(.), xy = T)[,3],
         Y = extract(ind, vect(.), xy = T)[,4],
         rds_dist = extract(covars, vect(.))[,2],
         water_dist = extract(covars, vect(.))[,3],
         late_GS = extract(covars, vect(.))[,4],
         shdi = extract(covars, vect(.))[,5],
         long = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2]) %>%
  as_tibble() %>%
  mutate(sep = "/") %>%
  dplyr::select(Detector, X, Y, Session, sep, rds_dist, water_dist, late_GS, shdi, M, site, long, lat) %>%
  as.data.frame()

tdf_full %>%
  ggplot(aes(x = long, y = lat, colour = factor(Session, labels = c("Spring 2023", "Autumn 2024")))) +
  geom_spatraster(data = subst(hab, "Cleared, non-native vegetation, buildings", NA), inherit.aes = F) +
  geom_sf(data = roads, inherit.aes = F, colour = warm_grey) +
  geom_sf(data = bigdesert, colour = galliano, fill = NA, inherit.aes = F, linewidth = 1) +
  geom_point() +
  scale_colour_cesar_d(name = "Session") +
  scale_fill_cesar_d(palette = "cesar_light", name = "Habitat", na.translate = F) +
  theme_cesar() +
  drop_titles("both") +
  facet_wrap(~Session) +
  coord_sf(xlim = (ext(bd) * 1.2)[c(1,2)],
           ylim = (ext(bd) * 1.2)[c(3,4)])

cesar_save("plots/sampling.png",
           scale = 2,
           height = 5.625,
           width = 10,
           dpi = 300)

all_dat <-
  dir("data", full.names = T)[str_detect(dir("data"), "capture_history")]

all_dat <-
  tibble(path = all_dat,
         h = as.numeric(str_remove(sapply(str_split(str_remove(all_dat, ".csv"), "_"), "[[", 3), "h")),
         at = as.numeric(sapply(str_split(str_remove(all_dat, ".csv"), "_"), "[[", 4)),
         m = as.numeric(sapply(str_split(str_remove(all_dat, ".csv"), "_"), "[[", 5))) %>%
  left_join(read_csv("data/ScatResultsComp.csv"))

models <- list()
pred.aic <- list()
pred.r <- list()
pred.param <- list()
pred.tot <- list()
for (i in 1:nrow(all_dat)){
  edf <-
    read_csv(as.character(all_dat[i, 1])) %>%
    dplyr::select(individual, site, `Oct 2023`, `Nov 2023`, `Mar 2024`, `May 2024`) %>%
    pivot_longer(3:6) %>%
    drop_na() %>%
    mutate(Session = year(my(name))) %>%
    mutate(site = as.character(site)) %>%
    filter(site %in% tdf_full$site) %>%
    left_join(read_csv("data/dingo_lookup.csv") %>%
                dplyr::select(site, collection_date)) %>%
    arrange(dmy(collection_date)) %>%
    group_by(Session) %>%
    mutate(Occasion = as.numeric(factor(collection_date, levels = unique(collection_date), ordered = T))) %>%
    left_join(read_csv("data/sex.csv") %>%
                mutate(site = as.character(sample))) %>%
    mutate(sex = case_when(sex == 1 ~ "F",
                           sex == 2 ~ "M")) %>%
    rename(ID = individual,
           Sex = sex,
           Date = collection_date) %>%
    left_join(tdf_full %>%
                dplyr::select(site, Detector)) %>%
    ungroup() %>%
    dplyr::select(Session, ID, Occasion, Detector, Sex, Date, site) %>%
    mutate(Session = factor(Session),
           ID = as.numeric(ID),
           Occasion = as.numeric(Occasion),
           Sex = factor(Sex),
           Date = factor(Date, levels = unique(Date))) %>%
    dplyr::select(-site) %>%
    update_sex_values() %>%
    as.data.frame()
  
  all_dat[i, 6] <- nrow(edf)
  all_dat[i, 7] <- nrow(edf %>% count(ID))
  
  sexes <- edf %>% group_by(ID) %>% count(Sex) %>% count(ID) %>% pull(n)
  
  if(any(sexes > 1)) {
    sex <- F
  } else {
    sex <- T
  }
  
  tdf_full %>%
    mutate(rds_dist = paste0("/", rds_dist)) %>%
    dplyr::select(Detector, X, Y, rds_dist, water_dist, late_GS, shdi) %>%
    unique() %>%
    as.data.frame() %>%
    write.table("temp/tdf.txt", row.names = F, col.names = F, quote = F)
  
  if(sex){
    edf %>%
      dplyr::select(Session, ID, Occasion, Detector, Sex) %>%
      as.data.frame() %>%
      write.table("temp/edf.txt", row.names = F, col.names = F, quote = F)
  } else {
    edf %>%
      dplyr::select(Session, ID, Occasion, Detector) %>%
      as.data.frame() %>%
      write.table("temp/edf.txt", row.names = F, col.names = F, quote = F)
  }
  
  captfile <- "temp/edf.txt"
  trapfile <- "temp/tdf.txt"
  if(sex){
    BDW_CH <- read.capthist(captfile,
                            trapfile,
                            detector = "count",
                            covnames = "Sex",
                            trapcovnames = c("rds_dist", "water_dist", "late_GS", "shdi"))
  } else {
    BDW_CH <- read.capthist(captfile,
                            trapfile,
                            detector = "count",
                            trapcovnames = c("rds_dist", "water_dist", "late_GS", "shdi"))
  }
  
  summary(BDW_CH)
  
  par(mfrow=c(1,2))
  plot(BDW_CH, tracks = TRUE)
  
  m <- unlist(moves(BDW_CH)) 
  par(mfrow=c(1,1))
  hist(m, xlab = "Movement m", main = "")
  
  initialsigma <- RPSV(BDW_CH, CC = TRUE)
  cat("Quick and biased estimate of sigma =", initialsigma[[1]], "m\n")
  
  mask <-
    make.mask(traps(BDW_CH),
              buffer = 20000,
              spacing = 2500,
              poly = bd_mask %>% st_as_sf())
  
  mask <- addCovariates(mask, covars[[1]])
  mask <- addCovariates(mask, covars[[2]])
  mask <- addCovariates(mask, covars[[3]])
  mask <- addCovariates(mask, covars[[4]])
  
  par(mfrow=c(1,2))
  plot(mask)
  
  if(sex){
    base.args <-
      list(capthist = BDW_CH,
           mask = mask,
           trace = TRUE,
           ncores = 10,
           hcov = "Sex")
  } else {
    base.args <-
      list(capthist = BDW_CH,
           mask = mask,
           ncores = 10,
           trace = TRUE)
  }
  
  # use covariates for g0?
  
  varyingmodels <- list(g0 ~ rds_dist,
                        list(D ~ water_dist, g0 ~ rds_dist),
                        list(D ~ late_GS, g0 ~ rds_dist),
                        list(D ~ shdi, g0 ~ rds_dist),
                        list(D ~ water_dist + late_GS, g0 ~ rds_dist),
                        list(D ~ water_dist + shdi, g0 ~ rds_dist),
                        list(D ~ late_GS + shdi, g0 ~ rds_dist),
                        list(D ~ water_dist + late_GS + shdi, g0 ~ rds_dist))
  
  fits <-
    list.secr.fit(model = varyingmodels,
                  constant = base.args)
  
  pred.aic[[i]] <- AIC(fits, criterion = "AICc")[,-2]
  
  bestmod = which(names(fits) == rownames(AIC(fits, criterion = "AICc")[1,]))
  
  models[[i]] <- fits[[bestmod]]
  
  par(mfrow=c(1,1))
  esaPlot(fits[[bestmod]])
  abline(v = 4 * initialsigma[[1]], lty = 2, col = 'red')
  abline(v = 4 * initialsigma[[2]], lty = 2, col = 'blue')
  
  surfaceD <- predictDsurface(fits[[bestmod]], se.D = T, cl.D = T)
  
  pred.r[[i]] <-
    lapply(1:2, function(x)
      surfaceD[[x]] %>%
        as_tibble() %>%
        bind_cols(covariates(surfaceD[[x]])) %>%
        mutate(Session = c(2023, 2024)[x]) %>%
        convert_df_to_raster(resolution = c(2500, 2500), crs = crs(hab)) %>%
        setNames(c(2023, 2024)[x])) %>%
    do.call(c, .)
  
  pred.tot[[i]] <-
    bind_rows(lapply(1:2,
                     function(y)
                       region.N(fits[[bestmod]])[[y]] %>%
                       mutate(Session = c(2023, 2024)[y]))) %>%
    rownames_to_column("type") %>%
    mutate(type = substr(type, start = 1, stop = 1),
           h = as.numeric(all_dat[i, 2]),
           at = as.numeric(all_dat[i, 3]),
           m = as.numeric(all_dat[i, 4]),
           scats = as.numeric(all_dat[i, 6]),
           genotypes = as.numeric(all_dat[i, 7]))
  
  if(sex){
    pred.param[[i]] <-
      secr:::predict.secr(fits[[bestmod]])[[2]] %>%
      rownames_to_column("parameter") %>%
      mutate(estimate = round(case_when(parameter == "D" ~ estimate*100,
                                        T ~ estimate), 4),
             SE.estimate = round(case_when(parameter == "D" ~ SE.estimate*100,
                                           T ~ SE.estimate), 4),
             lcl = round(case_when(parameter == "D" ~ lcl*100,
                                   T ~ lcl), 4),
             ucl = round(case_when(parameter == "D" ~ ucl*100,
                                   T ~ ucl), 4),
             sex = "male",
             h = as.numeric(all_dat[i, 2]),
             at = as.numeric(all_dat[i, 3]),
             m = as.numeric(all_dat[i, 4]),
             scats = as.numeric(all_dat[i, 6]),
             genotypes = as.numeric(all_dat[i, 7]))
  } else {
    pred.param[[i]] <-
      secr:::predict.secr(fits[[bestmod]])[[1]] %>%
      rownames_to_column("parameter") %>%
      mutate(estimate = round(case_when(parameter == "D" ~ estimate*100,
                                        T ~ estimate), 4),
             SE.estimate = round(case_when(parameter == "D" ~ SE.estimate*100,
                                           T ~ SE.estimate), 4),
             lcl = round(case_when(parameter == "D" ~ lcl*100,
                                   T ~ lcl), 4),
             ucl = round(case_when(parameter == "D" ~ ucl*100,
                                   T ~ ucl), 4),
             sex = "total",
             h = as.numeric(all_dat[i, 2]),
             at = as.numeric(all_dat[i, 3]),
             m = as.numeric(all_dat[i, 4]),
             scats = as.numeric(all_dat[i, 6]),
             genotypes = as.numeric(all_dat[i, 7]))
  }
}

## plots

bind_rows(lapply(pred.param,
                 function(p)
                   if(any(p$parameter == "pmix")){
                     tibble(name = calculate_sex_specific_val(p[1,3], p[1,4], p[4,3], p[4,4]) %>% unlist() %>% names(),
                            value = calculate_sex_specific_val(p[1,3], p[1,4], p[4,3], p[4,4]) %>% unlist()) %>%
                       mutate(h = unique(p$h),
                              at = unique(p$at),
                              m = unique(p$m),
                              scats = unique(p$scats),
                              genotypes = unique(p$genotypes)) 
                   } else {
                     tibble(name = c("total.estimate", "total.se"),
                            value = c(p[1,3], p[1,4])) %>%
                       mutate(h = unique(p$h),
                              at = unique(p$at),
                              m = unique(p$m),
                              scats = unique(p$scats),
                              genotypes = unique(p$genotypes)) 
                   })) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m),
         sex = sapply(str_split(name, "\\."), "[[", 1),
         name = sapply(str_split(name, "\\."), "[[", 2)) %>%
  pivot_wider() %>%
  mutate(ucl = estimate + 1.96*se,
         lcl = estimate - 1.96*se) %>%
  ggplot(aes(x = Dataset, y = estimate, fill = sex)) +
  labs(y = "Density (km2)") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, size = 1, position = position_dodge(width = .5)) +
  geom_point(shape=21, size=6, aes(fill=sex), position = position_dodge(width = .5)) +
  geom_text(aes(label = paste("s:", scats),
                y = .05)) +
  geom_text(aes(label = paste("g:", genotypes),
                y = .055)) +
  theme_cesar() +
  scale_colour_cesar_d() +
  scale_fill_cesar_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_density.png",
           preset = "print")

bind_rows(pred.param) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  filter(parameter == "g0") %>%
  dplyr::select(Dataset, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = Dataset, y = estimate, fill = sex)) +
  labs(y = "Detection probability") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, size = 1, position = position_dodge(width = .5)) +
  geom_point(shape=21, size=6, fill = galliano, position = position_dodge(width = .5)) +
  geom_text(aes(label = paste("s:", scats),
                y = .03)) +
  geom_text(aes(label = paste("g:", genotypes),
                y = .032)) +
  theme_cesar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_prob.png",
           preset = "print")

bind_rows(pred.param) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  filter(parameter == "sigma") %>%
  dplyr::select(Dataset, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = Dataset, y = estimate, fill = sex)) +
  labs(y = "Sigma (m)") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, size = 1, position = position_dodge(width = .5)) +
  geom_point(shape=21, size=6, fill = galliano, position = position_dodge(width = .5)) +
  geom_text(aes(label = paste("s:", scats),
                y = 15000)) +
  geom_text(aes(label = paste("g:", genotypes),
                y = 16000)) +
  theme_cesar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_sigma.png",
           preset = "print")

bind_rows(lapply(1:length(pred.tot),
                 function(p)
                   if(any(pred.param[[p]]$parameter == "pmix")){
                     tibble(name = calculate_sex_specific_val(pred.tot[[p]][,2], pred.tot[[p]][,3], pred.param[[p]][4,3], pred.param[[p]][4,4]) %>% unlist() %>% names(),
                            value = calculate_sex_specific_val(pred.tot[[p]][,2], pred.tot[[p]][,3], pred.param[[p]][4,3], pred.param[[p]][4,4]) %>% unlist()) %>%
                       mutate(h = unique(pred.tot[[p]]$h),
                              at = unique(pred.tot[[p]]$at),
                              m = unique(pred.tot[[p]]$m),
                              scats = unique(pred.tot[[p]]$scats),
                              genotypes = unique(pred.tot[[p]]$genotypes))
                   } else {
                     tibble(name = c("total.estimate1", "total.estimate2", "total.estimate3", "total.estimate4",
                                     "total.se1", "total.se2", "total.se3", "total.se4"),
                            value = c(pred.tot[[p]][,2], pred.tot[[p]][,3])) %>%
                       mutate(h = unique(pred.tot[[p]]$h),
                              at = unique(pred.tot[[p]]$at),
                              m = unique(pred.tot[[p]]$m),
                              scats = unique(pred.tot[[p]]$scats),
                              genotypes = unique(pred.tot[[p]]$genotypes))
                   })) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m),
         sex = sapply(str_split(name, "\\."), "[[", 1),
         name = sapply(str_split(name, "\\."), "[[", 2)) %>%
  mutate(Session = case_when(str_detect(name, "1|2") ~ 2023,
                             str_detect(name, "3|4") ~ 2024),
         type = case_when(str_detect(name, "1|3") ~ "E",
                          str_detect(name, "2|4") ~ "R"),
         name = str_remove(name, "1|2|3|4")) %>%
  pivot_wider() %>%
  mutate(ucl = estimate + 1.96*se,
         lcl = estimate - 1.96*se,
         Session = factor(Session, labels = c("Spring 2023", "Autumn 2024"))) %>%
  filter(type == "E", Session == "Spring 2023") %>%
  unique() %>%
  ggplot(aes(x = Dataset, y = estimate, fill = sex)) +
  geom_errorbar(aes(ymax = ucl, ymin = lcl), width = 0, size = 1, position = position_dodge(width = .5)) +
  geom_point(shape=21, size=6, position = position_dodge(width = .5)) +
  geom_text(aes(label = paste("s:", scats),
                y = 120)) +
  geom_text(aes(label = paste("g:", genotypes),
                y = 130)) +
  theme_cesar() +
  scale_colour_cesar_d() +
  scale_fill_cesar_d() +
  labs(y = "Population size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_popsize.png",
           preset = "print")

bind_rows(lapply(1:length(models),
                 function(p)
                   plot(models[[p]], xval = seq(0, 20000, by = 1000), limits = T)[[1]] %>%
                   mutate(h = unique(pred.tot[[p]]$h),
                          at = unique(pred.tot[[p]]$at),
                          m = unique(pred.tot[[p]]$m),
                          scats = unique(pred.tot[[p]]$scats),
                          genotypes = unique(pred.tot[[p]]$genotypes)))) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  ggplot(aes(x = x, y = y, colour = Dataset, fill = Dataset)) +
  geom_ribbon(aes(ymax = ucl, ymin = lcl), colour = NA, alpha = .05) +
  geom_line() +
  labs(x = "Distance (m)",
       y = "Detection probability") +
  theme_cesar() +
  scale_colour_cesar_d() +
  scale_fill_cesar_d()

cesar_save("plots/BDW_pdetection.png",
           preset = "print")


preds <- lapply(models,
                function(m)
                  predict(m, newdata = data.frame(water_dist = seq(7, 11, by = .5),
                                                  h2 = factor("M", levels = c("F", "M")),
                                                  rds_dist = 0,
                                                  late_GS = 100)))
bind_rows(
  lapply(1:length(models),
         function(m)
           bind_rows(
             lapply(1:length(seq(7, 11, by = .5)),
                    function(p) preds[[m]][[p]][1,] %>%
                      mutate(water_dist = seq(7, 11, by = .5)[p]))) %>%
           mutate(scats = unique(pred.tot[[m]]$scats),
                  genotypes = unique(pred.tot[[m]]$genotypes),
                  h = unique(pred.tot[[m]]$h),
                  at = unique(pred.tot[[m]]$at),
                  m = unique(pred.tot[[m]]$m))
  )
) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  ggplot(aes(x = exp(water_dist), y = estimate*100, colour = Dataset, fill = Dataset)) +
  geom_ribbon(aes(ymax = ucl*100, ymin = lcl*100), colour = NA, alpha = .05) +
  geom_line() +
  labs(x = "Distance from water (m)",
       y = "Density (km2)") +
  theme_cesar() +
  scale_colour_cesar_d() +
  scale_fill_cesar_d()

cesar_save("plots/BDW_dwater.png",
           preset = "print")

preds <- lapply(models,
                function(m)
                  predict(m, newdata = data.frame(water_dist = 7,
                                                  h2 = factor("M", levels = c("F", "M")),
                                                  rds_dist = 0,
                                                  late_GS = seq(0, 100))))
bind_rows(
  lapply(1:length(models),
         function(m)
           bind_rows(
             lapply(1:100,
                    function(p) preds[[m]][[p]][1,] %>%
                      mutate(late_GS = p))) %>%
           mutate(scats = unique(pred.tot[[m]]$scats),
                  genotypes = unique(pred.tot[[m]]$genotypes),
                  h = unique(pred.tot[[m]]$h),
                  at = unique(pred.tot[[m]]$at),
                  m = unique(pred.tot[[m]]$m))
  )
) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  ggplot(aes(x = late_GS, y = estimate*100, colour = Dataset, fill = Dataset)) +
  geom_ribbon(aes(ymax = ucl*100, ymin = lcl*100), colour = NA, alpha = .05) +
  geom_line() +
  labs(x = "Proportion of plants in late vegetation growth",
       y = "Density (km2)") +
  theme_cesar() +
  scale_colour_cesar_d() +
  scale_fill_cesar_d()

cesar_save("plots/BDW_dfire.png",
           preset = "print")

preds <- lapply(models,
                function(m)
                  predict(m, newdata = data.frame(water_dist = 7,
                                                  h2 = factor("M", levels = c("F", "M")),
                                                  rds_dist = log(seq(200, 11000, by = 500)),
                                                  late_GS = 100)))

lapply(1:length(models),
       function(m){
         lapply(1:length(seq(200, 11000, by = 500)),
                function(r)
                {
                  tibble(p = sapply(seq(0, 20000, by = 1000),
                                    function(x) preds[[m]][[r]][2,2]*exp(-x^2/(2*(preds[[m]][[r]][3,2])^2))),
                         dist = seq(0, 20000, by = 1000),
                         rds_dist = seq(200, 11000, by = 500)[[r]])  
                }
         ) %>%
           bind_rows() %>%
           mutate(scats = unique(pred.tot[[m]]$scats),
                  genotypes = unique(pred.tot[[m]]$genotypes),
                  h = unique(pred.tot[[m]]$h),
                  at = unique(pred.tot[[m]]$at),
                  m = unique(pred.tot[[m]]$m))
       }) %>%
  bind_rows() %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  ggplot(aes(x = dist, y = p, colour = rds_dist, group = factor(rds_dist))) +
  geom_line() +
  theme_cesar() +
  scale_colour_cesar_g(name = "Distance from road (m)", palette = "collin_d", mid = (200 + 11000)/2) +
  labs(x = "Distance (m)",
       y = "Detection probability") +
  facet_wrap(~Dataset) +
  theme(legend.text = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_dprob.png",
           preset = "print")

## tables

bind_rows(lapply(1:length(pred.tot),
                 function(p)
                   if(any(pred.param[[p]]$parameter == "pmix")){
                     tibble(name = calculate_sex_specific_val(pred.tot[[p]][,2], pred.tot[[p]][,3], pred.param[[p]][4,3], pred.param[[p]][4,4]) %>% unlist() %>% names(),
                            value = calculate_sex_specific_val(pred.tot[[p]][,2], pred.tot[[p]][,3], pred.param[[p]][4,3], pred.param[[p]][4,4]) %>% unlist()) %>%
                       mutate(h = unique(pred.tot[[p]]$h),
                              at = unique(pred.tot[[p]]$at),
                              m = unique(pred.tot[[p]]$m),
                              scats = unique(pred.tot[[p]]$scats),
                              genotypes = unique(pred.tot[[p]]$genotypes))
                   } else {
                     tibble(name = c("total.estimate1", "total.estimate2", "total.estimate3", "total.estimate4",
                                     "total.se1", "total.se2", "total.se3", "total.se4"),
                            value = c(pred.tot[[p]][,2], pred.tot[[p]][,3])) %>%
                       mutate(h = unique(pred.tot[[p]]$h),
                              at = unique(pred.tot[[p]]$at),
                              m = unique(pred.tot[[p]]$m),
                              scats = unique(pred.tot[[p]]$scats),
                              genotypes = unique(pred.tot[[p]]$genotypes))
                   })) %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m),
         sex = sapply(str_split(name, "\\."), "[[", 1),
         name = sapply(str_split(name, "\\."), "[[", 2)) %>%
  mutate(Session = case_when(str_detect(name, "1|2") ~ 2023,
                             str_detect(name, "3|4") ~ 2024),
         type = case_when(str_detect(name, "1|3") ~ "E",
                          str_detect(name, "2|4") ~ "R"),
         name = str_remove(name, "1|2|3|4")) %>%
  pivot_wider() %>%
  mutate(ucl = estimate + 1.96*se,
         lcl = estimate - 1.96*se,
         Session = factor(Session, labels = c("Spring 2023", "Autumn 2024"))) %>%
  filter(type == "E", Session == "Spring 2023") %>%
  dplyr::select(Dataset, sex, estimate, se, lcl, ucl) %>%
  write_csv("output/BDW_popsize.csv")

lapply(1:length(pred.aic),
       function(p)
         pred.aic[[p]] %>%
         as_tibble() %>%
         mutate(h = unique(pred.tot[[p]]$h),
                at = unique(pred.tot[[p]]$at),
                m = unique(pred.tot[[p]]$m),
                scats = unique(pred.tot[[p]]$scats),
                genotypes = unique(pred.tot[[p]]$genotypes))) %>%
  bind_rows() %>%
  mutate(Dataset = paste0("h", h, "_", at, "_", m)) %>%
  dplyr::select(Dataset, model, npar, logLik, AIC, AICc, dAICc, AICcwt) %>%
  write_csv("output/BDW_modelcomp.csv")

write_rds(models,
          "output/best_models.rds")

datasets <-
  read_csv("output/BDW_modelcomp.csv") %>%
  pull(Dataset) %>%
  unique()

lapply(1:length(datasets),
       function(r){
         writeRaster(pred.r[[r]],
                     paste0("output/rasters/", datasets[[r]], ".tiff"))
       })

## more plots!

read_csv("output/BDW_popsize_total.csv") %>%
  mutate(un_rat = round(genotypes/scats, 2),
         estimate = round(estimate, 2),
         ucl = round(ucl, 2),
         lcl = round(lcl, 2)) %>%
  dplyr::select(un_rat, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = scats, y = estimate, colour = un_rat, group = factor(un_rat))) +
  geom_errorbar(aes(ymax = ucl, ymin = lcl), width = 0, size = 1, position = position_dodge(width = 1)) +
  geom_point(shape=21, size=6, position = position_dodge(width = 1), fill = warm_grey) +
  geom_smooth(aes(x = scats, y = estimate), colour = "black", inherit.aes = F, method = "lm") +
  theme_cesar() +
  labs(y = "Population size",
       x = "No. of scats") +
  scale_colour_cesar_g(palette = "gold_teal_d", name = "Proportion of unique individuals", mid = 0.5, limits = c(0.3, 0.7))

cesar_save("plots/BDW_scats.png",
           preset = "print")

read_csv("output/BDW_popsize_total.csv") %>%
  mutate(un_rat = round(genotypes/scats, 2),
         estimate = round(estimate, 2),
         ucl = round(ucl, 2),
         lcl = round(lcl, 2)) %>%
  dplyr::select(un_rat, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = genotypes, y = estimate, colour = un_rat, group = factor(un_rat))) +
  geom_errorbar(aes(ymax = ucl, ymin = lcl), width = 0, size = 1, position = position_dodge(width = 1)) +
  geom_point(shape=21, size=6, position = position_dodge(width = 1), fill = warm_grey) +
  geom_smooth(aes(x = genotypes, y = estimate), colour = "black", inherit.aes = F, method = "lm") +
  theme_cesar() +
  labs(y = "Population size",
       x = "No. of unique individuals") +
  scale_colour_cesar_g(palette = "gold_teal_d", name = "Proportion of unique individuals", mid = 0.5, limits = c(0.3, 0.7))

cesar_save("plots/BDW_genotypes.png",
           preset = "print")

read_csv("output/BDW_popsize_total.csv") %>%
  mutate(un_rat = round(genotypes/scats, 2),
         estimate = round(estimate, 2),
         ucl = round(ucl, 2),
         lcl = round(lcl, 2)) %>%
  dplyr::select(un_rat, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = un_rat, y = estimate, colour = scats, group = factor(scats))) +
  geom_errorbar(aes(ymax = ucl, ymin = lcl), width = 0, size = 1, position = position_dodge(width = 1)) +
  geom_point(shape=21, size=6, position = position_dodge(width = 1), fill = warm_grey) +
  geom_smooth(aes(x = un_rat, y = estimate), colour = "black", inherit.aes = F, method = "lm") +
  theme_cesar() +
  labs(y = "Population size",
       x = "Proportion of unique individuals") +
  scale_colour_cesar_g(palette = "gold_teal_d", name = "No. of scats", mid = 75) +
  xlim(0.45, 0.65)

cesar_save("plots/BDW_unrat.png",
           preset = "print")

read_csv("output/BDW_popsize_total.csv") %>%
  mutate(un_rat = round(genotypes/scats, 2),
         estimate = round(estimate, 2),
         ucl = round(ucl, 2),
         lcl = round(lcl, 2)) %>%
  dplyr::select(un_rat, estimate, lcl, ucl, scats, genotypes) %>%
  unique() %>%
  ggplot(aes(x = un_rat, y = estimate, colour = genotypes, group = factor(genotypes))) +
  geom_errorbar(aes(ymax = ucl, ymin = lcl), width = 0, size = 1, position = position_dodge(width = 1)) +
  geom_point(shape=21, size=6, position = position_dodge(width = 1), fill = warm_grey) +
  geom_smooth(aes(x = un_rat, y = estimate), colour = "black", inherit.aes = F, method = "lm") +
  theme_cesar() +
  labs(y = "Population size",
       x = "Proportion of unique individuals") +
  scale_colour_cesar_g(palette = "gold_teal_d", name = "No. of unique individuals", mid = 40) +
  xlim(0.45, 0.65)

cesar_save("plots/BDW_unrat2.png",
           preset = "print")

tdf_full %>%
  ggplot(aes(x = long, y = lat)) +
  geom_spatraster(data = (lapply(pred.r, "[[", 1) %>% do.call(c, .) %>%
                            setNames(datasets))*100, inherit.aes = F, alpha = 0.5) +
  geom_sf(data = roads, inherit.aes = F, colour = warm_grey) +
  geom_sf(data = bigdesert, colour = galliano, fill = NA, inherit.aes = F, linewidth = 1) +
  geom_point(size = .2) +
  scale_fill_cesar_c(palette = "rufous_c", reverse = T, name = "Density (km2)", na.value = NA) +
  theme_cesar() +
  drop_titles("both") +
  coord_sf(xlim = (ext(bd) * 1.2)[c(1,2)],
           ylim = (ext(bd) * 1.2)[c(3,4)]) +
  facet_wrap(~lyr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(angle = 45, hjust = 1))

cesar_save("plots/BDW_densityrasts.png",
           scale = 1.5,
           height = 10,
           width = 10,
           dpi = 300)
