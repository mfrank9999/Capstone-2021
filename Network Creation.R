install.packages(c("stringr", "jsonlite", "sjmisc", "MASS", "UpSetR", "ggplot2", "googledrive","dplyr","varhandle", "plyr", "ggpmisc", "spatstat.utils", "umap", "ggthemes","RColorBrewer","infotheo","BiocManager", "gprofiler2", "devtools"))
BiocManager::install(c("arrow", "Seurat", "RcisTarget","zoo", "mixtools", "rbokeh","DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne","doMC", "doRNG","BiocParallel"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
devtools::install_github("aertslab/AUCell")

library(remotes)
library(Seurat)
library(stringr)
library(sjmisc)
library(MASS)
library(SCENIC)
library(gprofiler2)

remotes::install_version("spatstat", version = "1.64-1")

rm(list = ls()) # clean workspace

#===============================================PARAMETERS===================================================

#SCENIC Options Parameters and Initialization
cpuCores <-24
dbDir <- "databases"
dataLoc <- "CM_Data"
organism <- "hgnc"
dbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-5kb-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions <- SCENIC::initializeScenic(org = organism, dbDir = dbDir, dbs = dbs)  
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- cpuCores


#Seurat Preprocessing Parameters
numGenes = 1000 #number of variable features to keep

#Preprocessing Parameters
maxCellsPerCondition <- 0

#GENIE3 Parameters
targetsPerPart <- 1000
resumePreviousRun <- FALSE

#Network Postprocessing Parameters
edgeWeightThreshold <- 0.001
correlationThreshold <- 0.03


#IMPORTANT: Currently data is expected to be in the following format:
#file type = .csv file
#Rows =  cells
#Columns = genes
#Gene notation = any

#===============================================CODE TO RUN===================================================
files <- list.files(dataLoc, pattern = ".\\.csv" )

transitions <- data.frame(name = c("Sequential 2-3", "Sequential 3-4", "Sequential 4-5", "Combined 2-5"), start = c(2,3,4,2), end = c(3,4,5,5), data = vector(mode = "list", length = 4))
for (i in 1:length(transitions$name)) {
  print(paste0("Creating data for transition ", transitions$name[i]))
  data <- readFiles(files[transitions$start[i]:transitions$end[i]])
  transitions$data[[i]] <- list(mergeData(data))
}

network <- runNetworkCreation(transitions$data[[4]], scenicOptions = scenicOptions, maxCellsPerCondition = maxCellsPerCondition, targetsPerPart = targetsPerPart, edgeWeightThreshold = edgeWeightThreshold, correlationThreshold = correlationThreshold, outputLocation = "CM_Network")




#=================================================FUNCTIONS===================================================


#-----------------------------------------------PREPROCESSING-------------------------------------------------
readFiles <- function(files, maxCellsPerCondition = 0) {
  data <- vector(mode = "list", length = length(files))
  for(i in seq_along(files)){
    #Format & retrieve data
    filePath <- file.path(dataLoc,files[i])
    expressionMatrix <- read.csv(filePath)
    expressionMatrix <- t(expressionMatrix)
    #seuratObj <- Seurat::CreateSeuratObject(expressionMatrix, project = tools::file_path_sans_ext(files[i]), assay = "RNA", min.cells = 0, min.features = 0)
    #seuratObj <- RenameCells(seuratObj, add.cell.id = tools::file_path_sans_ext(files[i]))
    colnames(expressionMatrix) <- paste0(tools::file_path_sans_ext(files[i]),"_",colnames(expressionMatrix))
    data[[i]] <- expressionMatrix
  }
  return(data)
}
mergeData <- function(data){
  mergedData <- NULL
  for(j in 1:length(data)){
    if(is.null(mergedData)){
      mergedData <- data[[j]]
    } else {
      mergedData <- merge(mergedData, data[[j]], by = "row.names", all = TRUE)
      #fixes bug where merge creates a new column called "Row.names" with the row names from the two matrices and assigning the row names of the new merged matrix as integers 
      if("Row.names" %in% colnames(mergedData)){
        rownames(mergedData) <- t(mergedData["Row.names"])
        mergedData <- mergedData[,!(colnames(mergedData) %in% c("Row.names"))]
      }
      mergedData[is.na(mergedData)] <- 0
    }
  }
  mergedData <- as.matrix(mergedData)
  return(mergedData)
}
convertGeneNames <- function(data, target = "HGNC", verbose = FALSE){
  rows <- rownames(data)
  conversion <- gconvert(query = rows, organism = "hsapiens", target = target, mthreshold = 1, filter_na = FALSE)
  geneNames <-conversion[,4]
  if(verbose){
    genesLost <- 0
    for (j in 1:length(geneNames)){
      if(is.na(geneNames[j])){
        genesLost <- genesLost+1
      }
    }
    print(paste0(genesLost," genes lost during gene name conversion"))
  }
  .rowNamesDF(data, TRUE) <- geneNames
  return(data)
}
filterData <- function(data, scenicOptions = NULL){
  if(is.null(scenicOptions)){
    print("No scenic options provide for scenic filtering")
    scenicOptions <- reinitializeScenicOptions()
  }
  data <- as.matrix(data)
  genesKept <- geneFiltering(data, scenicOptions=scenicOptions)
  data <- data[genesKept, ]
  
  return(data)
}
seuratPreprocess <- function(mergedData, minCells = 3, minGenes = 200, numGenes = 0, projName = "NetworkCreation"){
  if(typeof(mergedData) == "list"){
    mergedSeurat <- Seurat::CreateSeuratObject(mergedData[[1]], project ="Network Creation", assay = "RNA", min.cells = minCells, min.features = minGenes)
    for(i in 2:length(mergedData)){
      tmp <- Seurat::CreateSeuratObject(mergedData[[i]], project ="Network Creation", assay = "RNA", min.cells = minCells, min.features = minGenes)
      mergedSeurat <- merge(mergedSeurat, tmp, project = projName, merge.data = TRUE )
    }
  } else {
    mergedSeurat <- Seurat::CreateSeuratObject(mergedData, project =projName, assay = "RNA", min.cells = minCells, min.features = minGenes)
  }
  mergedSeurat <- NormalizeData(mergedSeurat, normalization.method = "LogNormalize")
  dataReturn <- as.matrix(Seurat::GetAssayData(mergedSeurat, slot = "data"))
  
  if(numGenes > 0){
    mergedSeurat <- FindVariableFeatures(mergedSeurat, selection.method = "vst", nfeatures = numGenes)
    print(paste0("Keeping ", length(VariableFeatures(mergedSeurat)), " variable features"))
    dataReturn <- dataReturn[rownames(dataReturn) %in% VariableFeatures(mergedSeurat),]
    if(length(rownames(dataReturn))!=length(VariableFeatures(mergedSeurat))){
      print(paste0("Only ", length(rownames(dataReturn)), "remain"))
    }
  }
  
  return(dataReturn)
}
cutMergedData <- function(mergedData, maxCellsPerCondition = 0, maxGenes = 0, verbose = FALSE){
  cutMat<-mergedData
  cutMat <- as.matrix(cutMat)
  #cut cells
  if(maxCellsPerCondition>0){
    groupNames <- gsub(pattern = "_(?:.(?!_))+$", "", colnames(cutMat), perl = TRUE)
    groupNames <- groupNames[!duplicated(groupNames)]
    for(i in 1:length(groupNames)){
      print(paste0("cutting group: ", groupNames[i]))
      columnsInGroup <- grep(groupNames[i], colnames(cutMat))
      if(length(columnsInGroup)>maxCellsPerCondition){
        toRemove <- columnsInGroup[(maxCellsPerCondition+1):length(columnsInGroup)]
        cutMat <- cutMat[,-toRemove]
      }
      
    }
  }
  
  if(maxGenes>0){
    cutMat <- cutMat[1:maxGenes,]
  }
  return(cutMat)
}

#---------------------------------------------NETWORK CREATION------------------------------------------------

createNetwork <- function(edgeWeightThreshold = 0.001, correlationThreshold = 0.03, scenicOptions = scenicOptions, coexMethodName){
  #Process SCENIC output
  regulonTargetsInfo <- readRDS(getIntName(scenicOptions, "regulonTargetsInfo"))
  #keep coexmodule containing all spearman correlations
  regulonTargetsInfo <- regulonTargetsInfo[which(regulonTargetsInfo$coexModule == paste0(coexMethodName,"IgnCorr"))]
  #remove all rows that have genes that aren't TFs
  regulonTargetsInfo <- regulonTargetsInfo[regulonTargetsInfo$gene %in% regulonTargetsInfo$TF,]
  #remove all self regulators (should these be removed?)
  regulonTargetsInfo <- regulonTargetsInfo[!which(regulonTargetsInfo$TF == regulonTargetsInfo$gene),]
  #remove all non-correlated relationships
  regulonTargetsInfo <- regulonTargetsInfo[which(abs(regulonTargetsInfo$spearCor) > correlationThreshold), ]
  #trim further if needed
  regulonTargetsInfo <- regulonTargetsInfo[which(regulonTargetsInfo$CoexWeight > 2*edgeWeightThreshold), ]
  regulons <- createRegulons(regulonTargetsInfo = regulonTargetsInfo)
  network <- regulonTargetsInfo[,c(1,2)]
  network$edge <- as.numeric(regulonTargetsInfo$spearCor > correlationThreshold) - as.numeric(regulonTargetsInfo$spearCor < -correlationThreshold)
  colnames(network) <- c("source","target", "type")
  write.table(network,row.names = FALSE,col.names = TRUE, sep = ",", file=file.path(dirname(getOutName(scenicOptions, "s2_motifEnrichment")), "networkForRACIPE.csv"),quote = FALSE)
  
  
  
  return(network)
}
createRegulons <- function(regulonTargetsInfo){
  #create and save regulons
  regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$highConfAnnot)
  regulons <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["TRUE"]]))
  {
    regulons <- lapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
  }
  regulons_extended <- NULL
  if(!is.null(regulonTargetsInfo_splitByAnnot[["FALSE"]]))
  {
    regulons_extended <- lapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(unlist(x[,"gene"])))
    regulons_extended <- setNames(lapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], unlist(regulons_extended[[tf]]))))), names(regulons_extended))
    names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
  }
  regulons <- c(regulons, regulons_extended)
  saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))
  return(regulons)
}
runNetworkCreation <-function(mergedData, scenicOptions, numGenes = 10000, maxCellsPerCondition = 0, targetsPerPart = 1000, edgeWeightThreshold, correlationThreshold, outputLocation = NA){
  print(paste0("Beginning network creation at ",Sys.time()))
  #create new save location and reset scenic options
  print("Initializing Scenic Options")
  scenicOptions <- reinitializeScenicOptions(scenicOptions)
  scenicOptions <- setOutputDirectory(scenicOptions, outputDirectory = outputLocation)
  
  print(paste0("Preprocessing data ",Sys.time()))
  #preprocess data
  mergedData <- convertGeneNames(mergedData)
  mergedData <- filterData(mergedData, scenicOptions)
  mergedData <- seuratPreprocess(mergedData, minCells = 3, minGenes = 200, numGenes = numGenes)
  expMat <- as.matrix(cutMergedData(mergedData, maxCellsPerCondition = maxCellsPerCondition))
  write.table(expMat,row.names = TRUE,col.names = NA, sep = ",", file=file.path(dirname(getOutName(scenicOptions, "s2_motifEnrichment")), "processesedData.csv"),quote = FALSE)
  
  #numParts
  if(targetsPerPart>0){
    nParts <- length(rownames(expMat))/targetsPerPart
  } else {
    nParts <- 10
  }
  
  
  print(paste0("Creating Network (Using ", nParts, " parts) ", Sys.time()))
  #create network
  runGenie3(expMat, scenicOptions, nParts = nParts, resumePreviousRun = resumePreviousRun)
  runCorrelation(expMat, scenicOptions)
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions, topThr = c(edgeWeightThreshold)) #saving to scenicOptions just saves current progress
  coexMethodName <- paste("w", format(edgeWeightThreshold, scientific=FALSE), sep="")
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c(coexMethodName), onlyPositiveCorr = FALSE)
  network <- createNetwork(scenicOptions = scenicOptions, edgeWeightThreshold = edgeWeightThreshold, correlationThreshold = correlationThreshold, coexMethodName = coexMethodName)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expMat, skipBinaryThresholds = TRUE, skipHeatmap = TRUE , skipTsne = TRUE)
  regulonAUC <- getAUC(readRDS(getIntName(scenicOptions, "aucell_regulonAUC")))
  
  print(paste0("Network Complete! Nodes: ", length(unique(network$source)), " , Edges: ", length(newtork$source)))
  return(network)
}

#----------------------------------------AUTOMATION AND ORGANIZATION------------------------------------------
runTime <- function(mergedData, lengthPoints, scenicOptions){
  tsData <- data.frame(timeDif = rep(0,length(lengthPoints)), cells = rep(0, length(lengthPoints)), genes = rep(0, length(lengthPoints)))
  for (i in 1:length(lengthPoints)){
    scenicOptions <- reinitializeScenicOptions(scenicOptions)
    print(paste0("running length ", lengthPoints[i]))
    expMat <- as.matrix(cutMergedData(mergedData, maxCellsPerCondition = lengthPoints[i], filter = TRUE, scenicOptions = scenicOptions))
    start = Sys.time()
    runGenie3(expMat, scenicOptions)
    end = Sys.time()
    timeDif <- as.numeric(difftime(end,start, units = c("mins")))
    tsData$timeDif[i] <- timeDif
    tsData$cells[i] <- length(colnames(expMat))
    tsData$genes[i] <- length(rownames(expMat))
    print(paste0("# of cells:", tsData$cells[i], " # of genes:", tsData$genes[i]," time (mins):", tsData$timeDif[i]))
    
  }
  
  return(tsData)
  
}
createNetworkParameters <- function(data, numGenes = 10000, maxCellsPerCondition = 0, edgeWeightThreshold, correlationThreshold){
  data = mergedData
  maxCellsPerCondition = c(200,0)
  numGenes = 1000
  edgeWeightThreshold = c(0.0005, 0.001)
  correlationThreshold = c(0.01, 0.03)
  if(typeof(data) != "list"){
    data <- list(data)
  }
  nwConditions <- expand.grid(maxCellsPerCondition,numGenes,edgeWeightThreshold,correlationThreshold, 1:length(data))
  
  colnames(nwConditions) <- c("maxCellsPerCondition", "numGenes", "edgeWeightThreshold", "correlationThreshold", "dataPosition")
  networks <- vector(mode = "list", length = length(rownames(nwConditions)))
  
  for(i in 1:length(rownames(nwConditions))){
    networks[[i]] <- list(parameters = nwConditions[i,], data = data[[nwConditions[i,5]]])
  }
  return(networks)
}
reinitializeScenicOptions <- function(scenicOptions = NULL){
  if(is.null(scenicOptions)){
    print(paste0("No scenic options provided, initializing with default parameters"))
    organism <- "hgnc"
    databases <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-5kb-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
    scenicOptions <- SCENIC::initializeScenic(org = organism, dbDir = dbDir, dbs = databases)  
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- cpuCores
  } else {
    cores <- scenicOptions@settings$nCores
    seed <- scenicOptions@settings$seed
    scenicOptions <- initializeScenic(org = scenicOptions@inputDatasetInfo$org, dbDir = scenicOptions@settings$dbDir, dbs = scenicOptions@settings$dbs)
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- cores
    scenicOptions@settings$seed <- seed
  }
  
  return(scenicOptions)
}
setOutputDirectory <- function(scenicOptions, outputDirectory = NA){
  if(!(is.na(outputDirectory) | outputDirectory=="")){
    if(!dir.exists(outputDirectory)){
      dir.create(outputDirectory)
    }
    if(!dir.exists(file.path(outputDirectory,"int"))){
      dir.create(file.path(outputDirectory,"int"))
    }
    if(!dir.exists(file.path(outputDirectory,"output"))){
      dir.create(file.path(outputDirectory,"output"))
    }
  }
  
  intFnames <-c(
    "int/1.1_genesKept.Rds",                        
    "int/1.2_corrMat.Rds",                          
    "int/1.3_GENIE3_weightMatrix.Rds",              
    "int/1.4_GENIE3_linkList.Rds",                  
    "int/1.5_weightPlot",                           
    "int/1.6_tfModules_asDF.Rds",                   
    "int/2.1_tfModules_forMotifEnrichmet.Rds",      
    "int/2.2_motifs_AUC.Rds",                       
    "int/2.3_motifEnrichment.Rds",                  
    "int/2.4_motifEnrichment_selfMotifs_wGenes.Rds",
    "int/2.5_regulonTargetsInfo.Rds",               
    "int/2.6_regulons_asGeneSet.Rds",               
    "int/2.6_regulons_asIncidMat.Rds",              
    "int/3.1_regulons_forAUCell.Rds",               
    "int/3.2_aucellGenesStats",                     
    "int/3.3_aucellRankings.Rds",                   
    "int/3.4_regulonAUC.Rds",                       
    "int/3.5_AUCellThresholds.Rds",                 
    "int/3.5_AUCellThresholds_Info.tsv",            
    "int/4.1_binaryRegulonActivity.Rds",            
    "int/4.2_binaryRegulonActivity_nonDupl.Rds",    
    "int/4.3_regulonSelections.Rds",                
    "int/4.4_binaryRegulonOrder.Rds")
  
  outFnames=c(
    "output/Step2_MotifEnrichment.tsv", 
    "output/Step2_MotifEnrichment_preview.html", 
    "output/Step2_regulonTargetsInfo.tsv", 
    "output/Step3_RegulonActivity_heatmap",
    "output/Step3_RegulonActivity_tSNE_colByActivity",
    "output/Step3_RegulonActivity_tSNE_colByCellProps",
    "output/Step4_BoxplotActiveCellsRegulon",
    "output/Step4_BinaryRegulonActivity_Heatmap_",
    "output/Step4_BinaryRegulonActivity_tSNE_colByActivity",
    "output/Step4_BinaryRegulonActivity_tSNE_colByCellProps",
    "output/scenic.loom")
  
  for(i in 1:length(intFnames)){
    if(!(is.na(outputDirectory) | outputDirectory=="")){
      scenicOptions@fileNames$int[i] <- file.path(outputDirectory,intFnames[i])
    } else {
      scenicOptions@fileNames$int[i] <- intFnames[i]
    }
  }
  for(i in 1:length(outFnames)){
    if(!(is.na(outputDirectory) | outputDirectory=="")){
      scenicOptions@fileNames$output[i] <- file.path(outputDirectory,outFnames[i])
    } else {
      scenicOptions@fileNames$output[i] <- outFnames[i]
    }
  }
  return(scenicOptions)
}
