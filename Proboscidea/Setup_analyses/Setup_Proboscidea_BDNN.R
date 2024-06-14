library(ape)
library(PVR)


# Auxilliary function
#####################
# Get phylogenetic eigenvectors from tree
get_eigenvec <- function(tree, variance_fraction) {
  
  decomp <- PVRdecomp(tree, type = "newick") # produces object of class 'PVR'
  eigvec <- as.data.frame(decomp@Eigen$vectors) # extract eigenvectors
  
  egval <- decomp@Eigen$values # extract eigenvalues
  eigPerc <- egval/sum(egval) # calculate % of variance
  eigPercCum <- t(cumsum(eigPerc)) # cumulated variance
  
  # eigenvectors representing more than X% variance
  numEigen <- sum(eigPercCum < variance_fraction)
  
  eigOK <- eigvec[, 1:numEigen, drop = TRUE] 
  
  if(numEigen == 0){
    print(paste("The variance_fraction should at leat be equal to ", 
                eigPercCum[1]))
    eigOK <- eigvec[1]
  }
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  # Eigenvectors generated in object 'eigenTobind'
  # rename eigenTobind species column so it matches trait dataset species column
  eigOK <- as.data.frame(eigOK)
  names(eigOK)[1] <- "c1"
  row.names(eigOK) <- decomp@phylo$tip.label
  
  return(eigOK)
}


# Set seed for reproducible input generation
############################################
set.seed(73)


# Set paths
###########
PathData <- '/home/torsten/Work/BDNN/DataGitHub/BDNN/Proboscidea'


# Number of replicates for BDNN analyses
########################################
NumReplicates <- 10


# Load data
###########
# Posterior distribution of proboscidean pyhlogenetic analysis (Cantalapiedra et al. 2021)
ProbTrees <- read.nexus(file.path(PathData,
                                  'Cantalapiedra2021',
                                  'proboscidea_trees_zero_free.tre'))

ProbSize <- read.table(file.path(PathData,
                                 'Cantalapiedra2021',
                                 'Extended_Data_1_Sheet_Traits_and_PFTs.csv'),
                       sep = '\t', header = TRUE, row.names = 1)

NMDS <- read.table(file.path(PathData,
                             'Cantalapiedra2021',
                             'NMDS_data.txt'),
                   sep = '\t', header = TRUE, row.names = 1)

GeoDistr <- read.table(file.path(PathData,
                                 'Cantalapiedra2021',
                                 'Geographic_distribution.txt'),
                       sep = '\t', header = TRUE, row.names = 1)

# Species specific temperature curves binned into epochs
# (top row: most recent)
Temperature <- read.table(file.path(PathData,
                                    'Species_specific_temperature',
                                    'BinnedSpeciesTempThroughTime.txt'),
                          sep = '\t', header = TRUE)

FilterTaxa <- read.table(file.path(PathData,
                                   'PyRate_analyses',
                                   'Body_size',
                                   'FilterTaxa.txt'),
                         sep = '\t', header = FALSE)[, 1]


# Get NumReplicates random trees from posterior phylogenies
NumTrees <- length(ProbTrees)
ProbTrees <- ProbTrees[sample(NumTrees, NumReplicates)]

# Drop taxa from phylogeny that are not in ProbSize
ProbTrees <- drop.tip(ProbTrees,
                      setdiff(ProbTrees[[1]]$tip.label, rownames(ProbSize)))

# Order traits alphabetically
ProbSize <- ProbSize[order(rownames(ProbSize)), ]


# Island taxa
IslandTaxa <- c('Stegoloxodon_celebensis', 'Palaeoloxodon_falconeri',
                'Stegodon_florensis', 'Stegodon_hypsilophus', 
                'Stegodon_luzonensis', 'Stegodon_mindanensis',
                'Stegodon_sompoensis', 'Stegodon_sondaari',
                'Stegodon_sumbaensis', 'Stegodon_timorensis',
                'Stegodon_trigonocephalus', 'Stegolophodon_pseudolatidens',
                'Anancus_perimensis', 'Palaeoloxodon_mnaidriensis',
                'Mammuthus_lamarmorai', 'Palaeoloxodon_cypriotes', 
                'Mammuthus_exilis', 'Mammuthus_creticus',
                'Palaeoloxodon_tiliensis', 'Palaeoloxodon_lomolinoi',
                'Palaeoloxodon_creutzburgi')


# Create initial trait tables  B O D Y  S I Z E
###############################################
TraitsMaster <- data.frame(Species = rownames(ProbSize),
                           # Center body size to it's median value
                           Body_size = ProbSize$body_size_cat - 1 - 4,
                           Island = 0,
                           Temperature = 0,
                           Humans = 0)

for (i in 1:NumReplicates) {
  Traits <- TraitsMaster
  # Assign island taxa
  Traits$Island[rownames(ProbSize) %in% IslandTaxa] <- 1
  # Calculate phylogenetic eigenvectors
  PhyloEigVec <- get_eigenvec(ProbTrees[[i]], 0.95)
  # Order as body size
  PhyloEigVec <- PhyloEigVec[match(Traits$Species, rownames(PhyloEigVec)), 1:2]
  PVR1 <- PhyloEigVec[, 1]
  PVR2 <- PhyloEigVec[, 2]
  Traits$PVR1 <- scale(PVR1)
  Traits$PVR2 <- scale(PVR2)
  # Get mean and standard deviation of the features.
  # We need them later for plotting the effect of the feature 
  # with partial dependence plots.
  BackscaleStats <- data.frame(PVR1 = c(mean(PVR1), sd(PVR1)),
                               PVR2 = c(mean(PVR2), sd(PVR2)),
                               Body_size = c(4, 1),
                               Temperature = c(14.49948, 4.939096),
                               row.names = c('Mean', 'SD'))
  write.table(BackscaleStats,
              file.path(PathData,
                        'PyRate_analyses',
                        'Body_size',
                        paste0('Backscale_stats_', i, '.txt')),
              sep = '\t', row.names = FALSE, quote = FALSE)
  # Create folder where we will place the predictor files.
  dir.create(file.path(PathData,
                       'PyRate_analyses',
                       'Body_size',
                       'BDNN_predictors', i,
                       'speciation'),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PathData,
                       'PyRate_analyses',
                       'Body_size',
                       'BDNN_predictors', i,
                       'extinction'),
             showWarnings = FALSE, recursive = TRUE)
  IdxTempTraitsBin <- match(colnames(Temperature), Traits$Species)
  for (j in 1:16) {
    # There are ten time bins
    TraitsBin <- Traits
    # Assign temperature
    TraitsBin[IdxTempTraitsBin, 'Temperature'] <- unlist(Temperature[j, ])
    # Write to speciation predictors after filtering for targeted taxa
    TraitsBinSp <- TraitsBin[TraitsBin$Species %in% FilterTaxa, -1]
    write.table(TraitsBinSp,
                file.path(PathData,
                          'PyRate_analyses',
                          'Body_size',
                          'BDNN_predictors', i,
                          'speciation',
                          paste0(sprintf("%02d", j), '.txt')),
                sep = '\t', row.names = FALSE, quote = FALSE)
    # Write to extinction predictors after adding humans 
    # and filtering for targeted taxa
    if (j == 1) {
      # 0.00-0.0117
      TraitsBin$Humans <- 2
    }
    else if (j == 2) {
      # 0.0117-0.126
      # Africa and/or Eurasia, but not in America
      TraitsBin$Humans[GeoDistr$America == 0] <- 2
      # Mammuthus_meridionalis and Mammuthus_primigenius last record in Eurasia
      TraitsBin$Humans[TraitsBin$Species %in%
                         c('Mammuthus_meridionalis',
                           'Mammuthus_primigenius')] <- 2
    }
    else if (j %in% 3:4) {
      # 0.126-1.8
      # Africa and/or Eurasia, but not in America
      TraitsBin$Humans[GeoDistr$America == 0] <- 1
      TraitsBin$Humans[TraitsBin$Species %in%
                         c('Mammuthus_meridionalis',
                           'Mammuthus_primigenius')] <- 1
    }
    TraitsBinEx <- TraitsBin[TraitsBin$Species %in% FilterTaxa, -1]
    write.table(TraitsBinEx,
                file.path(PathData,
                          'PyRate_analyses',
                          'Body_size',
                          'BDNN_predictors', i,
                          'extinction',
                          paste0(sprintf("%02d", j), '.txt')),
                sep = '\t', row.names = FALSE, quote = FALSE)
  }
}


# Create initial trait tables  N M D S
######################################
TraitsMaster <- data.frame(Species = rownames(ProbSize),
                           NMDS,
                           GeoDistr,
                           Island = 0,
                           Grassland = 0,
                           Temperature = 0,
                           Humans = 0)
SpIndia <- TraitsMaster$Species %in% c("Elephas_hysudricus",
                                       "Gomphotherium_browni",
                                       "Stegodon_insignis")

for (i in 1:NumReplicates) {
  Traits <- TraitsMaster
  # Assign island taxa
  Traits$Island[rownames(ProbSize) %in% IslandTaxa] <- 1
  Traits$Africa[rownames(ProbSize) %in% IslandTaxa] <- 0
  Traits$America[rownames(ProbSize) %in% IslandTaxa] <- 0
  Traits$Eurasia[rownames(ProbSize) %in% IslandTaxa] <- 0
  # Calculate phylogenetic eigenvectors
  PhyloEigVec <- get_eigenvec(ProbTrees[[i]], 0.95)
  # Order as body size
  PhyloEigVec <- PhyloEigVec[match(Traits$Species, rownames(PhyloEigVec)), 1:2]
  PVR1 <- PhyloEigVec[, 1]
  PVR2 <- PhyloEigVec[, 2]
  Traits$PVR1 <- scale(PVR1)
  Traits$PVR2 <- scale(PVR2)
  # Get mean and standard deviation of the features.
  # We need them later for plotting the effect of the feature 
  # with partial dependence plots.
  BackscaleStats <- data.frame(PVR1 = c(mean(PVR1), sd(PVR1)),
                               PVR2 = c(mean(PVR2), sd(PVR2)),
                               Body_size = c(4, 1),
                               Temperature = c(14.49948, 4.939096),
                               row.names = c('Mean', 'SD'))
  write.table(BackscaleStats,
              file.path(PathData,
                        'PyRate_analyses',
                        'NMDS',
                        paste0('Backscale_stats_', i, '.txt')),
              sep = '\t', row.names = FALSE, quote = FALSE)
  # Create folder where we will place the predictor files.
  dir.create(file.path(PathData,
                       'PyRate_analyses',
                       'NMDS',
                       'BDNN_predictors', i,
                       'speciation'),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PathData,
                       'PyRate_analyses',
                       'NMDS',
                       'BDNN_predictors', i,
                       'extinction'),
             showWarnings = FALSE, recursive = TRUE)
  IdxTempTraitsBin <- match(colnames(Temperature), Traits$Species)

  for (j in 1:16) {
    # There are ten time bins
    TraitsBin <- Traits
    # Assign temperature
    TraitsBin[IdxTempTraitsBin, 'Temperature'] <- unlist(Temperature[j, ])
    # Grassland
    if (j <= 9) {
      # India
      TraitsBin$Grassland[SpIndia] <- 1
    }
    if (j <= 11) {
      # Africa
      TraitsBin$Grassland[GeoDistr$Africa == 1] <- 1
    }
    if (j <= 12) {
      # Eurasia, except India
      TraitsBin$Grassland[GeoDistr$Eurasia == 1 & !SpIndia] <- 1
    }
    if (j <= 13) {
      # America
      TraitsBin$Grassland[GeoDistr$America == 1] <- 1
    }
    # Write to speciation predictors after filtering for targeted taxa
    TraitsBinSp <- TraitsBin[TraitsBin$Species %in% FilterTaxa, -1]
    write.table(TraitsBinSp,
                file.path(PathData,
                          'PyRate_analyses',
                          'NMDS',
                          'BDNN_predictors', i,
                          'speciation',
                          paste0(sprintf("%02d", j), '.txt')),
                sep = '\t', row.names = FALSE, quote = FALSE)
    # Write to extinction predictors after adding humans 
    # and filtering for targeted taxa
    if (j == 1) {
      # 0.00-0.0117
      TraitsBin$Humans <- 2
    }
    else if (j == 2) {
      # 0.0117-0.126
      # Africa and/or Eurasia, but not in America
      TraitsBin$Humans[GeoDistr$America == 0] <- 2
      # Mammuthus_meridionalis and Mammuthus_primigenius last record in Eurasia
      TraitsBin$Humans[TraitsBin$Species %in%
                         c('Mammuthus_meridionalis',
                           'Mammuthus_primigenius')] <- 2
    }
    else if (j %in% 3:4) {
      # 0.126-1.8
      # Africa and/or Eurasia, but not in America
      TraitsBin$Humans[GeoDistr$America == 0] <- 1
      TraitsBin$Humans[TraitsBin$Species %in%
                         c('Mammuthus_meridionalis',
                           'Mammuthus_primigenius')] <- 1
    }
    TraitsBinEx <- TraitsBin[TraitsBin$Species %in% FilterTaxa, -1]
    write.table(TraitsBinEx,
                file.path(PathData,
                          'PyRate_analyses',
                          'NMDS',
                          'BDNN_predictors', i,
                          'extinction',
                          paste0(sprintf("%02d", j), '.txt')),
                sep = '\t', row.names = FALSE, quote = FALSE)
  }
}
