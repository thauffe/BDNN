library(sf)
library(terra)
library(rgplates)
library(viridisLite)
library(rnaturalearth)
# remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(countrycode)
library(foreach)
library(doParallel)


##########################################
# Setup rgplates with the paleomap model #
##########################################

archive <- file.path(system.file("extdata", package = "rgplates"), 
                     "paleomap_model_v19o_r1c.zip")
# extract to temporary directory
unzip(archive, exdir = tempdir())
# path to the combined model/rotation file
path <- file.path(tempdir(), 
                  "paleomap_model_v19o_r1c/paleomap_model_v19o_r1c.mod")
# register in R - to be used in reconstruct()
GplatesModel <- platemodel(path)
# Example usage
# reconstruct(matrix(c(95, 54), nrow=1), 140.5, model = GplatesModel)


#####################
# Utility functions #
#####################

getCountry <- function(Name) {
  CountSf <- ne_states(country = Name, returnclass = "sf")
  if (Name == 'Chile') {
    Clip <- cbind(c(-78, -60, -60, -78, -78),
                  c(-15, -15, -60, -60, -15))
  }
  else if (Name == 'Ecuador') {
    Clip <- cbind(c(-85, -75, -75, -85, -85),
                  c(  3,   3,  -7,  -7,   3))
  }
  else if (Name == 'France') {
    Clip <- cbind(c(-20, 20, 20, -20, -20),
                  c( 60, 60, 30,  30,  60))
  }
  else if (Name == 'India') {
    Clip <- cbind(c(60, 100, 80, 60, 60),
                  c(40,  40,  5,  5, 40))
  }
  else if (Name == 'Netherlands') {
    Clip <- cbind(c(0,  10, 10, 0,   0),
                  c(60, 60, 40, 40, 60))
  }
  else if (Name == 'Portugal') {
    Clip <- cbind(c(-13, 5, 5, -13, -13),
                  c( 45, 45, 35, 35, 45))
  }
  else if (Name == 'South Africa') {
    Clip <- cbind(c( 10,  40,  40,  10,  10),
                  c(-20, -20, -40, -40, -20))
  }
  else if (Name == 'Spain') {
    Clip <- cbind(c(-15,  5,  5, -15, -15),
                  c( 45, 45, 32,  32, 45))
  }
  else if (Name == 'United States of America') {
    Df <- st_drop_geometry(CountSf)
    Hawaii <- which(Df$name == 'Hawaii')
    CountSf <- CountSf[-Hawaii, ]
  }
  else if (Name == 'Russia') {
    Df <- st_drop_geometry(CountSf)
    Problem <- which(Df$name %in% c('Karelia', 'Chukchi Autonomous Okrug'))
    CountSf <- CountSf[-Problem, ]
  }
  if (Name %in% c('Chile', 'Ecuador', 'France', 'India',
                  'Netherlands', 'Portugal', 'South Africa',
                  'Spain')) {
    Clip <- st_polygon(list(Clip))
    Clip <- st_sfc(Clip, crs = "+proj=longlat +datum=WGS84")
    CountSf <- suppressWarnings(st_intersection(CountSf, Clip))
    CountSf <- st_cast(CountSf, to = 'MULTIPOLYGON')
  }
  return(CountSf)
}



#########################
# Set working directory #
#########################

setwd("/.../BDNN/Proboscidea")


#########################################################
# Load paleo temperatures from Hagen et al. 2021 (PNAS) #
# Only needed once if Ages.txt will be loaded instead   #
#########################################################

FirstTime <- FALSE

if (FirstTime) {
  Env <- read.table(file.path(getwd(), "PaleoTemp", "dataset_s3_trimmed.csv"),
                    header = TRUE, sep = ",", check.names = FALSE) 
} else {
  Ages <- read.table(file.path(getwd(), "PaleoTemp", "Ages.txt"),
                     sep = "\t", header = TRUE)
  Ages <- Ages$Name
}



###################################################################
# Correct names of time bins by re-ordering and make better names #
#                e.g. 00_00Ma, 17_50Ma, 24_17Ma                   #
#                  ! Only needed once !                           #
###################################################################

if (FirstTime) {
  # Correct names of ages 
  # (only needed for the original file from Oscar Hagen. 
  # dataset_s3_trimmed.csv is truncated at 83 Ma 
  # for not exceeding the maximum file size on GitHub)
  # colnames(Env)[3:ncol(Env)] <- colnames(Env)[ncol(Env):3]
  # Make usefull names for writting files (no. ".", trailing 0 etc.)
  for (i in 3:ncol(Env)) {
    Age <- colnames(Env)[i]
    if (grepl(".", Age, fixed = TRUE)) {
      Age <- strsplit(Age, ".", fixed = TRUE)[[1]]
      Age[1] <- sprintf("%02d", as.numeric(Age[1]))
      if (nchar(Age[2]) == 3) {
        Age[2] <- paste0(strsplit(Age[2], "Ma", fixed = TRUE)[[1]], "0", "Ma")
      }
    }
    else {
      Age <- strsplit(Age, "Ma", fixed = TRUE)[[1]]
      Age <- c(sprintf("%02d", as.numeric(Age[1])), "00Ma")
    }
    Age <- paste0(Age[1], "_", Age[2])
    colnames(Env)[i] <- Age
  }
  # Write to hard disk to avoid loading the entire heavy data set
  Ages <- data.frame(Name = colnames(Env))
  write.table(Ages, file.path(getwd(), "PaleoTemp", "Ages.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  Ages <- Ages$Name
}


########################################## 
# Write rasters until 45 Ma to hard disk #
#          ! Only needed once !          #
##########################################

if (FirstTime) {
  for (i in 3:423) {
    R <- rast(Env[, c(1:2, i)], type = "xyz", crs = "+proj=longlat +datum=WGS84")
    Age <- Ages[i]
    FileName <- file.path(getwd(), "PaleoTemp", "TempRasters",
                          paste0(Age, ".tif"))
    writeRaster(R, FileName, overwrite = TRUE)
  }
}


##############################################
# Copy shapefiles of continents through time #
#          ! Only needed once !              #
##############################################

Continents <- vector(mode = "list", length = 5)
if (FirstTime) {
  Continents[[1]] <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       "00_00Ma_Africa.shp"),
                             quiet = TRUE)
  Continents[[2]] <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       "00_00Ma_SouthAmerica.shp"),
                             quiet = TRUE)
  Continents[[3]] <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       "00_00Ma_NorthAmerica.shp"),
                             quiet = TRUE)
  Continents[[4]] <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       "00_00Ma_Europe.shp"),
                             quiet = TRUE)
  Continents[[5]] <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       "00_00Ma_Asia.shp"),
                             quiet = TRUE)
  names(Continents) <- c("Africa", "SouthAmerica", "NorthAmerica",
                         "Europe", "Asia")
  for (j in 1:length(Continents)) {
    ContShp <- Continents[[j]]
    if (j <= 4) {
      ContVertices <- st_coordinates(ContShp)[, 1:2]
    }
    else {
      ContVertices <- rbind(st_coordinates(ContShp$geometry[1])[, 1:2],
                            st_coordinates(ContShp$geometry[2])[, 1:2])
      RussiaFarEast <- ContVertices[, 1] < 0
      SetPlus180 <- ContVertices[, 1] <= 180 & ContVertices[, 1] >= 179
      SetMinus180 <- ContVertices[, 1] >= -180 & ContVertices[, 1] <= -179
    }
    print(j)
    for (i in 4:423) {
      # print(i)
      AgeName <- Ages[i]
      Age <- gsub("Ma", "", AgeName)
      Age <- as.numeric(gsub("_", ".", Age))
      PaleoVertices <- NA
      GplatesOK <- FALSE
      while (!GplatesOK) {
        PaleoVertices <- tryCatch(reconstruct(ContVertices,
                                              age = Age,
                                              model = GplatesModel),
                                  error = function(e) NA)
        if (!all(is.na(PaleoVertices))) {
          GplatesOK <- TRUE
        }
        else {
          print("Repeating Gplates")
        }
      }
      # PaleoVertices <- reconstruct(ContVertices, age = Age,
      #                              model = GplatesModel)
      if (j <= 4) {
        PaleoPoly <- st_polygon(list(PaleoVertices))
        # Would be nicer if these were elements of a single shapefile
        PaleoCont <- st_sfc(PaleoPoly, crs = "+proj=longlat +datum=WGS84")
      }
      else {
        # Part of Asia next to the Bering strait
        PaleoVertices[SetPlus180 , 1] <- 180
        PaleoVertices[SetMinus180 , 1] <- -180
        PaleoVertices[!RussiaFarEast, 1] <- abs(PaleoVertices[!RussiaFarEast, 1])
        Rfe <- PaleoVertices[RussiaFarEast, ]
        RfeWrong <- Rfe[, 1] > 0
        if (any(RfeWrong)) {
          Rfe <- Rfe[!RfeWrong, ]
          if (all(Rfe[1, ] != Rfe[nrow(Rfe), ])) {
            Rfe <- rbind(Rfe, Rfe[1, ])
          }
        }
        PaleoPoly <- st_polygon(list(PaleoVertices[!RussiaFarEast, ], Rfe))
        # Would be nicer if these were elements of a single shapefile
        PaleoCont <- st_sfc(PaleoPoly, crs = "+proj=longlat +datum=WGS84")
      }
      st_write(PaleoCont,
               file.path(getwd(), "PaleoTemp", "ContinentsThroughTime",
                         paste0(AgeName, "_", names(Continents)[j], ".shp")),
               delete_layer = TRUE, quiet = TRUE)
    }
  }
}




if (FirstTime) {
  Bio01Present <- rast(file.path(getwd(), "PaleoTemp", "Bioclim_WC21",
                                 "wc2.1_30s_bio", "wc2.1_30s_bio_1.tif"))
  OskarPresent <- rast(file.path(getwd(), "PaleoTemp",
                                 "TempRasters", "00_00Ma.tif"))
  MidHol <- rast(file.path(getwd(), "PaleoTemp", "MidHolocene_CCM4",
                           "ccmidbi_30s", "ccmidbi1.tif")) # Res 0.0083
  Lgm <- rast(file.path(getwd(), "PaleoTemp", "LGM_CCM4",
                        "cclgmbi_2-5m", "cclgmbi1.tif")) # Res 0.0416
  Lig <- rast(file.path(getwd(), "PaleoTemp", "LastInterglacial",
                        "lig_30s_bio", "lig_30s_bio_1.bil")) # Res 0.0083
  
  # Decrease spatial resolution
  MidHol <- aggregate(MidHol, fact = 4)
  
  # Same extent for all
  MidHol <- extend(MidHol, ext(-180, 180, -90, 90))
  Lig <- extend(Lig, ext(-180, 180, -90, 90))
  
  writeRaster(MidHol,
              file.path(getwd(), "PaleoTemp",
                        "MidHol.tif"), overwrite = TRUE)
  
  # Set to the same coordinate reference system
  crs(MidHol) <- crs(Bio01Present)
  crs(Lgm) <- crs(Bio01Present)
  crs(Lig) <- crs(Bio01Present)
  
  # Same resolution for all
  Bio01Present <- resample(Bio01Present, MidHol,
                           filename = file.path(getwd(),
                                                "PaleoTemp",
                                                "Bio01.tif"),
                           overwrite = TRUE)
  OskarPresent <- resample(OskarPresent, MidHol,
                           filename = file.path(getwd(),
                                                "PaleoTemp",
                                                "OskarPresent.tif"),
                           overwrite = TRUE)
  Lgm <- resample(Lgm, MidHol,
                  filename = file.path(getwd(),
                                       "PaleoTemp",
                                       "Lgm.tif"),
                  overwrite = TRUE)
  Lig <- resample(Lig, MidHol,
                  filename = file.path(getwd(),
                                       "PaleoTemp",
                                       "Lig.tif"),
                  overwrite = TRUE)

  # Difference between Oskar's estimate and Bioclim
  DiffPres <- OskarPresent - Bio01Present
  writeRaster(DiffPres,
              file.path(getwd(), "PaleoTemp",
                        "DiffPres.tif"), overwrite = TRUE)
  DiffPres[is.na(DiffPres)] <- 0
  writeRaster(DiffPres,
              file.path(getwd(), "PaleoTemp",
                        "DiffPres.tif"), overwrite = TRUE)
  
  # Adjust Mid Holocene, LGM and Last Interglacial by this offset
  MidHolAdj <- (MidHol / 10) + DiffPres
  writeRaster(MidHolAdj,
              file.path(getwd(), "PaleoTemp",
                        "MidHolocene_LGM_LastInterglacial", 
                        "MidHolocene.tif"), overwrite = TRUE)
  
  LgmAdj <- (Lgm / 10) + DiffPres
  writeRaster(LgmAdj,
              file.path(getwd(), "PaleoTemp",
                        "MidHolocene_LGM_LastInterglacial",
                        "LastGlacialMaximum.tif"), overwrite = TRUE)
  
  LigAdj <- (Lig / 10) + DiffPres
  writeRaster(LigAdj,
              file.path(getwd(), "PaleoTemp",
                        "MidHolocene_LGM_LastInterglacial", 
                        "LastInterglacial.tif"), overwrite = TRUE)
}
  
###############################################
# Continent specific temperature through time #
###############################################

WgsRas <- rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90, resolution = 0.01,
               crs = "+proj=longlat +datum=WGS84")
# Mollweide projection
RasProj <- project(WgsRas, y = "ESRI:54009", res = c(100000, 100000)) 
# Equal Earth projection, odd duplication of Alaska
# RasProj <- project(WgsRas, y = "EPSG:8857", res = c(100000, 100000))
Ccm <- vector(mode = "list", length = 3)
Ccm[[1]] <- rast(file.path(getwd(), "PaleoTemp",
                           "MidHolocene_LGM_LastInterglacial", 
                           "MidHolocene.tif"))
Ccm[[2]] <- rast(file.path(getwd(), "PaleoTemp",
                           "MidHolocene_LGM_LastInterglacial", 
                           "LastGlacialMaximum.tif"))
Ccm[[3]] <- rast(file.path(getwd(), "PaleoTemp",
                           "MidHolocene_LGM_LastInterglacial", 
                           "LastInterglacial.tif"))
LenI <- 423
ContNames <- c("Africa", "SouthAmerica", "NorthAmerica",
               "Europe", "Asia", "Eurasia",
               "NAmericaEurasia", "Americas")
TempThroughTime <- matrix(NA_real_,
                          ncol = length(ContNames) + 1,
                          nrow = LenI + 1)
colnames(TempThroughTime) <- c("Time", ContNames)


for (i in 3:LenI) {
  print(i)
  AgeName <- Ages[i]
  Age <- gsub("Ma", "", AgeName)
  Age <- as.numeric(gsub("_", ".", Age))
  TempThroughTime[i - 2, 1] <- Age
  WorldTemp <- rast(file.path(getwd(), "PaleoTemp",
                              "TempRasters", paste0(AgeName, ".tif")))
  WorldTempHighRes <- disagg(WorldTemp, fact = 20, method = 'bilinear')
  for (j in 1:length(ContNames)) {
    if (ContNames[j] %in% c("Eurasia", "NAmericaEurasia")) {
      AsiaShp <- st_read(file.path(getwd(), "PaleoTemp",
                                   "ContinentsThroughTime",
                                   paste0(AgeName, "_Asia.shp")),
                         quiet = TRUE)
      EuropeShp <- st_read(file.path(getwd(), "PaleoTemp",
                                     "ContinentsThroughTime",
                                     paste0(AgeName, "_Europe.shp")),
                           quiet = TRUE)
      if (ContNames[j] %in% c("NAmericaEurasia")) {
        NAmericaShp <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       paste0(AgeName, "_NorthAmerica.shp")),
                             quiet = TRUE)
      }
      AsiaShp <- suppressWarnings(st_cast(AsiaShp, "POLYGON"))
      L <- list(st_coordinates(AsiaShp$geometry[[1]])[, 1:2],
                st_coordinates(AsiaShp$geometry[[2]])[, 1:2],
                st_coordinates(EuropeShp$geometry[[1]])[, 1:2])
      if (ContNames[j] %in% c("NAmericaEurasia"))  {
        L[[4]] <- st_coordinates(NAmericaShp$geometry[[1]])[, 1:2]
      }
      ContShp <- st_polygon(L)
      ContShp <- st_sfc(ContShp, crs = "+proj=longlat +datum=WGS84")
      ContShp <- st_sf(data.frame(A = 1), geom = ContShp)
    }
    else if (ContNames[j] %in% c("Americas")) {
      NAmericaShp <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       paste0(AgeName, "_NorthAmerica.shp")),
                             quiet = TRUE)
      SAmericaShp <- st_read(file.path(getwd(), "PaleoTemp",
                                       "ContinentsThroughTime",
                                       paste0(AgeName, "_SouthAmerica.shp")),
                             quiet = TRUE)
      L <- list(st_coordinates(NAmericaShp$geometry[[1]])[, 1:2],
                st_coordinates(SAmericaShp$geometry[[1]])[, 1:2])
      ContShp <- st_polygon(L)
      ContShp <- st_sfc(ContShp, crs = "+proj=longlat +datum=WGS84")
      ContShp <- st_sf(data.frame(A = 1), geom = ContShp)
    }
    else {
      ContShp <- st_read(file.path(getwd(), "PaleoTemp",
                                   "ContinentsThroughTime",
                                   paste0(AgeName, "_", ContNames[j], ".shp")),
                         quiet = TRUE)
    }
    ContMask <- rasterize(ContShp, WorldTempHighRes)
    ContTemp <- mask(WorldTempHighRes, ContMask)
    ContTempProj <- project(ContTemp, RasProj)
    ContTempProj <- trim(ContTempProj, padding = 0)
    FileName <- paste0(AgeName, '_', ContNames[j], '.tif')
    writeRaster(ContTempProj,
                file.path(getwd(), 'PaleoTemp', 'TempRasters', FileName),
                overwrite = TRUE, filetype = 'GTiff')
    M <- summary(ContTempProj, size = 500000)
    MeanTemp <- as.numeric(strsplit(M[4], ":")[[1]][2])
    TempThroughTime[i - 2, j + 1] <- MeanTemp
    if (i == 3) {
      AgesCCM <- c("00_006Ma", '00_02Ma', '00_12Ma')
      for (k in 1:3) {
        CcmK <- resample(Ccm[[k]], WorldTempHighRes)
        ContTemp <- mask(CcmK, ContMask)
        ContTempProj <- project(ContTemp, RasProj)
        ContTempProj <- trim(ContTempProj, padding = 0)
        FileName <- paste0(AgesCCM[k], '_', ContNames[j], '.tif')
        writeRaster(ContTempProj,
                    file.path(getwd(), 'PaleoTemp', 'TempRasters', FileName),
                    overwrite = TRUE, filetype = 'GTiff')
        M <- summary(ContTempProj, size = 500000)
        MeanTemp <- as.numeric(strsplit(M[4], ":")[[1]][2])
        TempThroughTime[(LenI-2) + k, j + 1] <- MeanTemp
      }
    }
  }
}
TempThroughTime[(LenI-1):(LenI+1), 1] <- c(0.006, 0.020, 0.12)
TempThroughTime <- TempThroughTime[order(TempThroughTime[, 1]), ]

write.table(TempThroughTime,
            file.path(getwd(),
                      'PaleoTemp',
                      'Paleo_temp_per_continent_Hagen2021.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)




#########################################
# Check if we can omit Russian Far East #
#########################################

# Problems due to crossing the -180 meridian

if (FirstTime) {
  Occ <- read.table(file.path(getwd(), "Cantalapiedra2021",
                              "Proboscidea_occurrences_def_Nov_2020_for_Torsten.csv"),
                    sep = ",", header = TRUE, fill = NA) 
  
  # Check occurrence coordinates in QGIS
  Keep <- !is.na(Occ$LONG)
  Points <- st_as_sf(Occ[Keep, ], coords = c('LONG', 'LAT'), dim = "XY",
                     crs = "+proj=longlat +datum=WGS84")
  Points <- st_multipoint(as.matrix(Points[Keep, ]), dim = "XY")
  st_write(Points, file.path(getwd(), 'PaleoTemp', "Occurrences.shp"),
           delete_layer = TRUE, quiet = TRUE, driver = 'ESRI Shapefile')
  
  # Check if any occurrence overlaps with Russian Far East
  Russia <- getCountry('Russia')
  st_write(Russia, file.path(getwd(), 'PaleoTemp', "Russia.shp"),
           delete_layer = TRUE, quiet = TRUE, driver = 'ESRI Shapefile')
}



#################################
# Species specific temperatures # 
#################################

# Load all GIS files to avoid loading them a million times
WgsRas <- rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90, resolution = 0.01,
               crs = "+proj=longlat +datum=WGS84")
# Mollweide projection
RasProj <- project(WgsRas, y = "ESRI:54009", res = c(100000, 100000))

TempList <- vector(mode = 'list', length = 271)
ContNames <- c("Africa", "SouthAmerica", "NorthAmerica", "Europe", "Asia")
for (i in 1:274) {
  TmpShp <- vector(mode = 'list', length = 5)
  TmpRas <- vector(mode = 'list', length = 5)
  # This is terrible ...
  if (i == 1) {
    # Present
    AgeName <- "00_00Ma"
  }
  else if (i == 2) {
    # Mid Holocene
    AgeName <- "00_006Ma"
  }
  else if (i == 3) {
    # LGM
    AgeName <- "00_02Ma"
  }
  else if (i == 4) {
    # Last Interglacial
    AgeName <- "00_12Ma"
  }
  else {
    AgeName <- Ages[i - 1]
  }
  for (j in 1:5) {
    TmpRas[[j]] <- rast(file.path(getwd(), "PaleoTemp",
                                         "TempRasters",
                                  paste0(AgeName, "_", ContNames[j], ".tif")))
  }
  TempList[[i]] <- TmpRas
}


Occ <- read.table(file.path(getwd(), "Cantalapiedra2021",
                            "Proboscidea_occurrences_def_Nov_2020_for_Torsten.csv"),
                  sep = ",", header = TRUE, fill = NA) 


# Country names for rnaturalearth
Occ[Occ$COUNTRY == 'Burma', 'COUNTRY'] <- 'Myanmar'
Occ[Occ$COUNTRY == 'Congo - Kinshasa', 'COUNTRY'] <- 'Democratic Republic of the Congo'
Occ[Occ$COUNTRY == 'Congo: Democratic republic of (prev. Zaire)', 'COUNTRY'] <- 'Democratic Republic of the Congo'
Occ[Occ$COUNTRY == 'Côte d’Ivoire', 'COUNTRY'] <- 'Ivory coast'
Occ[Occ$COUNTRY == 'Czechia', 'COUNTRY'] <- 'Czech Republic'
Occ[Occ$COUNTRY == 'Macedonia: The Former Yugoslav Republic of', 'COUNTRY'] <- 'Macedonia'
Occ[Occ$COUNTRY == 'Serbia', 'COUNTRY'] <- 'Republic of Serbia'
Occ[Occ$COUNTRY == 'Perú', 'COUNTRY'] <- 'Peru'
Occ[Occ$COUNTRY == 'Tanzania', 'COUNTRY'] <- 'United Republic of Tanzania'
Occ[Occ$COUNTRY == 'United States', 'COUNTRY'] <- 'United States of America'

Occ[Occ$continent == 'Central America', 'continent'] <- 'NorthAmerica'
Occ[Occ$continent == 'Northern America', 'continent'] <- 'NorthAmerica'
Occ[Occ$continent == 'Central Asia', 'continent'] <- 'Asia'
Occ[Occ$continent == 'Eastern Africa', 'continent'] <- 'Africa'
Occ[Occ$continent == 'Eastern Asia', 'continent'] <- 'Asia'
Occ[Occ$continent == 'Eastern Europe', 'continent'] <- 'Europe'
Occ[Occ$continent == 'Middle Africa', 'continent'] <- 'Africa'
Occ[Occ$continent == 'Northern Africa', 'continent'] <- 'Africa'
Occ[Occ$continent == 'Northern Europe', 'continent'] <- 'Europe'
Occ[Occ$continent == 'South-Eastern Asia', 'continent'] <- 'Asia'
Occ[Occ$continent == 'Southern Africa', 'continent'] <- 'Africa'
Occ[Occ$continent == 'Southern Asia', 'continent'] <- 'Asia'
Occ[Occ$continent == 'Southern Europe', 'continent'] <- 'Europe'
Occ[Occ$continent == 'Western Africa', 'continent'] <- 'Africa'
Occ[Occ$continent == 'Western Asia', 'continent'] <- 'Asia'
Occ[Occ$continent == 'Western Europe', 'continent'] <- 'Europe'
Occ[Occ$continent == 'South America', 'continent'] <- 'SouthAmerica'


# Add underscore between genus and species
Occ$species_complete_corrected <- gsub(" ", "_",
                                       Occ$species_complete_corrected,
                                       fixed = TRUE)

# Age range for which we want to get the temperature
SpAges <- aggregate(MAX_AGE_homog ~ species_complete_corrected, 
                  data = Occ, FUN = 'max')
SpAges$MIN_AGE_homog <- aggregate(MIN_AGE_homog ~ Occ$species_complete_corrected,
                                data = Occ, FUN = 'min')[, 2]


# Species specific temperature trajectory
Nspecies <- nrow(SpAges)
SpTemp <- as.data.frame(matrix(NA_real_, ncol = Nspecies, nrow = 271 + 3))
colnames(SpTemp) <- SpAges$species_complete_corrected



AgesTemp <- rep(NA_real_, 274)
for (i in 4:273) {
  AgeName <- Ages[i]
  Age <- gsub("Ma", "", AgeName)
  AgesTemp[i + 1] <- as.numeric(gsub("_", ".", Age))
}
AgesTemp[1:4] <- c(0, 0.006, 0.02, 0.12)


# Build for each species a convex hull around its countries on each continent
ConvHullVertices <- matrix(NA_real_, nrow = 0, ncol = 2)
SpeciesConvHull <- c()
ContinentConvHull <- c()
for (s in 1:Nspecies) {
  S <- SpAges$species_complete_corrected[s]
  SpOcc <- Occ[Occ$species_complete_corrected == S, ]
  SpCont <- unique(SpOcc$continent)
  for (i in 1:length(SpCont)) {
    SpOccCont <- SpOcc[SpOcc$continent == SpCont[i], ]
    Countries <- unique(SpOccCont$COUNTRY)
    Coords <- matrix(NA_real_, nrow = 0, ncol = 2)
    for (j in 1:length(Countries)) {
      CountrySf <- getCountry(Countries[j])
      if (Countries[j] %in% c('United States of America', 'Russia')) {
        # Get states by occurrence coordinates
        OccCoord <- SpOcc[SpOcc$COUNTRY == Countries[j], c('LONG', 'LAT')]
        OccCoord <- OccCoord[!is.na(OccCoord$LONG), ]
        Points <- st_multipoint(as.matrix(OccCoord), dim = "XY")
        Points <- st_sfc(Points, crs = "+proj=longlat +datum=WGS84")
        Keep <- st_intersects(Points, CountrySf, sparse = FALSE)
        CountrySf <- CountrySf[Keep, ]
      }
      CoordsTmp <- st_coordinates(CountrySf$geometry)[, 1:2]
      Coords <- rbind(Coords, CoordsTmp)
    }
    Points <- st_multipoint(Coords, dim = "XY")
    Points <- st_sfc(Points, crs = "+proj=longlat +datum=WGS84")
    ConvHull <- st_convex_hull(Points)
    ConvHull <- st_coordinates(ConvHull)[, 1:2]
    ConvHullVertices <- rbind(ConvHullVertices, ConvHull)
    SpeciesConvHull <- c(SpeciesConvHull, rep(S, nrow(ConvHull)))
    ContinentConvHull <- c(ContinentConvHull, rep(SpCont[i], nrow(ConvHull)))
  }
}

# Get species specific temperatures
Nspecies <- nrow(SpAges)
SpTemp <- as.data.frame(matrix(NA_real_,
                               ncol = Nspecies,
                               nrow = length(AgesTemp)))
colnames(SpTemp) <- SpAges$species_complete_corrected
for (tt in 1:length(AgesTemp)) {
  print(tt)
  Age <- 0.0
  ConvHulltt <- ConvHullVertices
  # Get Paleocoordiantes for the convex hulls
  if (tt > 3) {
    Ages <- AgesTemp[tt]
    GplatesOK <- FALSE
    while (!GplatesOK) {
      PaleoVertices <- tryCatch(reconstruct(ConvHullVertices,
                                            age = Age,
                                            model = GplatesModel),
                                error = function(e) NA)
      if (!all(is.na(PaleoVertices))) {
        GplatesOK <- TRUE
      }
    }
    ConvHulltt <- PaleoVertices
  }
  for (s in 1:Nspecies) {
    S <- SpAges$species_complete_corrected[s]
    SpOcc <- Occ[Occ$species_complete_corrected == S, ]
    SpCont <- unique(SpOcc$continent)
    TempS <- rep(NA_real_, ncol = length(SpCont))
    for (i in 1:length(SpCont)) {
      Keep <- ContinentConvHull == SpCont[i] & SpeciesConvHull == S
      ConvHullPoints <- st_multipoint(ConvHulltt[Keep, ], dim = "XY")
      ConvHullPoints <- st_sfc(ConvHullPoints,
                               crs = "+proj=longlat +datum=WGS84")
      ConvHullMollweide <- st_transform(ConvHullPoints, crs = "ESRI:54009")
      LatRange <- range(st_coordinates(ConvHullMollweide)[, 2])
      ContIdx <- SpCont[i] == c("Africa", "SouthAmerica", "NorthAmerica",
                                "Europe", "Asia")
      ContTemp <- TempList[[tt]][ContIdx][[1]]
      ContTemp <- crop(ContTemp,
                       ext(ext(ContTemp)[1:2], LatRange[1], LatRange[2])) 
      M <- summary(ContTemp, size = 500000)
      MeanTemp <- as.numeric(strsplit(M[4], ":")[[1]][2])
      TempS[i] <- MeanTemp
    }
    SpTemp[tt, s] <- mean(TempS)
  }
}
rownames(SpTemp) <- AgesTemp

write.table(SpTemp,
            file.path(getwd(), "PaleoTemp", "SpeciesTempThroughTime.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)


par(las = 1, mar = c(4, 4, 0.5, 0.5))
plot(1, 1, type = 'n', xlim = c(35, 0), ylim = c(-20, 30))
for (i in 2:ncol(SpTemp)) {
  lines(AgesTemp, SpTemp[, i])
}


####################################
# Bin temperature for each species #
####################################

getBinnedTemp <- function(x, Intervals) {
  Bt <- aggregate(x, by = list(Intervals), FUN = mean)[, 2]
  return(Bt)
}

approxTemp <- function(x, Time, Breaks) {
  Bt <- approx(x = Time, y = x, xout = Breaks, method = 'constant')$y
  return(Bt)
}


SpTemp <- read.table(file.path(getwd(), "PaleoTemp",
                               "SpeciesTempThroughTime.txt"),
                     header = TRUE, sep = "\t")

AgesTemp <- as.numeric(rownames(SpTemp))



# Using few bins in BDNN inference
Breaks <- c(0.0, 0.0117, 0.126, 0.781, 1.8, 2.58, 3.6, 5.333, 7.246,
            11.63, 13.82, 15.97, 20.44, 23.03, 28.1, 33.9)
Intervals <- findInterval(AgesTemp, Breaks)

BinnedTemp <- apply(SpTemp, 2, function(x) getBinnedTemp(x, Intervals))
rownames(BinnedTemp) <- paste0("Bin", 0:15)
# Reorder from past to present as in BDNN
BinnedTemp <- BinnedTemp[nrow(BinnedTemp):1, ]


# Taxa older than 40 Ma and removed from the BDNN analysis
RemoveTaxa <- c("Daouitherium_rebouli", "Eritherium_azzouzorum",
                "Phosphatherium_escuilliei", "Khamsaconus_bulbosus",
                "Numidotherium_koholense",
                # These are species without traits
                "Serbelodon_burnhami", "Stegodon_mindanensis",
                "Stegodon_namadicus", "Stegodon_zhaotongensis",
                "Stegoloxodon_indonesicus")
BinnedTemp <- BinnedTemp[, !(colnames(BinnedTemp) %in% RemoveTaxa)]
MeanTemp <- mean(BinnedTemp) # 14.49948
SdTemp <- sd(BinnedTemp) # 4.939096
BinnedTemp <- (BinnedTemp - MeanTemp) / SdTemp

write.table(BinnedTemp,
            file.path(getwd(),
                      "Species_specific_temperature",
                      "BinnedSpeciesTempThroughTime.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

