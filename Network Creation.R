#===============================================INITIALIZATION==============================================
install.packages(c("stringr", "jsonlite", "sjmisc", "MASS", "UpSetR", "ggplot2", "googledrive","dplyr"))
install.packages(c("varhandle", "plyr", "ggpmisc", "spatstat.utils", "umap", "ggthemes","RColorBrewer","infotheo")) # rowr?
install.packages("BiocManager")
#need for database reading (feather file)
#needed to convert between gene ID types
install.packages("gprofiler2")
BiocManager::install("arrow")
BiocManager::install("Seurat")
BiocManager::install(c("RcisTarget"))
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
BiocManager::install("BiocParallel")
install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")
devtools::install_github("aertslab/AUCell")
packageVersion("SCENIC")

library(remotes)
remotes::install_version("spatstat", version = "1.64-1")
rm(list = ls()) # clean workspace

library(Seurat)
library(stringr)
library(sjmisc)
library(MASS)
library(SCENIC)
library(gprofiler2)
#===============================================PARAMETERS===================================================
cpuCores <-24 #number of cpu cores genie3 will use (should be 4 if using personal device or equal to the # of cores used in OOD)
targetsPerPart <- 1000 #targets per part. Genie3 splits its run by targets, with weight matrices being saved after each part in a seperate RDS file for later use
#setting targets per part to 0 will default to 10 parts, regardless of # of genes
numGenes = 10000
resumePreviousRun <- FALSE #set to true if your previous run was interrupted. Will resume from the last completed part
maxCellsPerCondition <- 200 #maximum number of cells used for each file provided.
#setting cells per condition to 0 will default to no maximum (all cells will be used)
dbDir <- "databases" #directory for RcisTarget DBs (these must be downloaded seperately) Note all directory paths should be in relation to the working directory
dataLoc <- "CM_Data" #directory for expression data Note all directory paths should be in relation to the working directory

genesPerRegulon <- 50
edgeWeightThreshold <- 0.001
correlationThreshold <- 0.07
#IMPORTANT: Currently data is expected to be in the following format:
#file type = .csv file
#Rows =  cells
#Columns = genes
#Gene notation = any
#================================================PREPROCESSING================================================
files <- list.files(dataLoc, pattern = ".\\.csv" )
data <- readFiles(files[2:5])
mergedData <- mergeData(data)
mergedData <- seuratPreprocess(mergedData, minCells = 3, minGenes = 200, numGenes = numGenes)
mergedData <- convertGeneNames(mergedData)
#===============================================NETWORK CREATION==============================================

organism <- "hgnc"
dbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-5kb-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
scenicOptions <- SCENIC::initializeScenic(org = organism, dbDir = dbDir, dbs = dbs)  
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- cpuCores

#tsData <- runTime(mergedData, lengthPoints = c(50,100,200,500), scenicOptions)

expMat <- as.matrix(filterMergedData(mergedData, maxCellsPerCondition = maxCellsPerCondition, filtered = TRUE, scenicOptions = scenicOptions))
#Save preprocessed data
write.table(expMat,row.names = TRUE,col.names = NA, sep = ",", file=file.path(dirname(getOutName(scenicOptions, "s2_motifEnrichment")), "processesedData.csv"),quote = FALSE)

if(targetsPerPart>0){
  nParts <- length(rownames(expMat))/targetsPerPart
} else {
  nParts <- 10
}
print(Sys.time())
runGenie3(expMat, scenicOptions, nParts = nParts, resumePreviousRun = resumePreviousRun)
runCorrelation(expMat, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions, topThr = c(edgeWeightThreshold)) #saving to scenicOptions just saves current progress
coexMethodName <- paste("w", format(edgeWeightThreshold, scientific=FALSE), sep="")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c(coexMethodName), onlyPositiveCorr = FALSE)

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
#create regulons
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
length(unique(names(regulons)))
saveRDS(regulons, file=getIntName(scenicOptions, "regulons"))

incidList <- reshape2::melt(regulons)
incidMat <- table(incidList[,2], incidList[,1])
saveRDS(incidMat, file=getIntName(scenicOptions, "regulons_incidMat"))

network <- regulonTargetsInfo[,c(1,2)]
network$edge <- as.numeric(regulonTargetsInfo$spearCor > correlationThreshold) - as.numeric(regulonTargetsInfo$spearCor < -correlationThreshold)
saveRDS(network, file=file.path(dirname(getOutName(scenicOptions, "s2_motifEnrichment")), "networkForRACIPE"))
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expMat, skipBinaryThresholds = TRUE, skipHeatmap = TRUE , skipTsne = TRUE)
write.table(network,row.names = FALSE,col.names = TRUE, sep = ",", file=file.path(dirname(getOutName(scenicOptions, "s2_motifEnrichment")), "networkForRACIPE.csv"),quote = FALSE)

# genie3ll <- readRDS(getIntName(scenicOptions, "genie3ll"))
# corrMat <- readRDS(getIntName(scenicOptions, "corrMat"))
# links <- genie3ll[which(genie3ll$weight>edgeWeightThreshold),]
# onlyTfs <- links[links$Target %in% links$TF,]
# modules <- split(onlyTfs, as.factor(onlyTfs$TF))
# modules <- modules[which(!lapply(modules, function(x) length(rownames((x)))) == 0)]
# modulesWithCorr <- lapply(modules, function(geneSet) 
#   {
#     cutGeneSet <- geneSet
#     cutGeneSet <- cutGeneSet[order(cutGeneSet$weight, decreasing = TRUE),]
#     cutGeneSet <- cutGeneSet[1:min(genesPerRegulon,length(rownames(cutGeneSet))),]
#     tf <- as.character(unique(cutGeneSet$TF))
#     targets <- as.character(cutGeneSet$Target)
#     #print(paste0("retrieving TF: ", tf, " and targets: ", targets))
#     cbind(cutGeneSet, correlation=c(as.numeric(corrMat[tf,targets] > correlationThreshold) - as.numeric(corrMat[tf,targets] < -correlationThreshold)))
#   })
# meltedModules<-reshape2::melt(modulesWithCorr, id.vars = c("TF", "Target", "weight", "correlation"))[1:4]
# network <-meltedModules[which(meltedModules$correlation!=0),]
# length(unique(network$TF))
#=================================================FUNCTIONS==================================================

#takes list of file names and outputs list of data (cells are annotated based on which file they are from)
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

#takes a list of data and merges it into one matrix, keeping all genes and filling in 0 for gene expression that are missing from data sets
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
seuratPreprocess <- function(mergedData, minCells = 3, minGenes = 200, numGenes = 2000, projName = "NetworkCreation"){
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
  mergedSeurat <- FindVariableFeatures(mergedSeurat, selection.method = "vst", nFeatures = numGenes)
  dataReturn <- as.matrix(Seurat::GetAssayData(mergedSeurat, slot = "data"))
}
#filters merged data using the file name in the cell annotation to seperate data. Can be used to limit the number of cells per condition, cut the number of genes, or filter genes using SCENICs filter method
filterMergedData <- function(mergedData, maxCellsPerCondition = 0, maxGenes = 0, scenicOptions = NULL, filtered = FALSE, verbose = FALSE){
  cutMat<-mergedData
  cutMat <- as.matrix(cutMat)
  #cut cells
  if(maxCellsPerCondition>0){
    groupNames <- gsub(pattern = "_(?:.(?!_))+$", "", colnames(cutMat), perl = TRUE)
    groupNames <- groupNames[!duplicated(groupNames)]
    for(i in 1:length(groupNames)){
      if(TRUE){
        print(paste0("cutting group: ", groupNames[i]))
      }
      columnsInGroup <- grep(groupNames[i], colnames(cutMat))
      if(length(columnsInGroup)>maxCellsPerCondition){
        toRemove <- columnsInGroup[(maxCellsPerCondition+1):length(columnsInGroup)]
        cutMat <- cutMat[,-toRemove]
      }
      
    }
  }
  
  if(filtered){
    if(is.null(scenicOptions)){
      print("No scenic options provide for scenic filtering")
      scenicOptions <- reinitializeScenicOptions()
    }
    genesKept <- geneFiltering(cutMat, scenicOptions=scenicOptions)
    cutMat <- cutMat[genesKept, ]
  }
  
  if(maxGenes>0){
    cutMat <- cutMat[1:maxGenes,]
  }
  return(cutMat)
}


#converts gene names from any gene naming method to the specified target method
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

#used to easily evaluate run time vs cell# 
runTime <- function(mergedData, lengthPoints, scenicOptions){
  tsData <- data.frame(timeDif = rep(0,length(lengthPoints)), cells = rep(0, length(lengthPoints)), genes = rep(0, length(lengthPoints)))
  for (i in 1:length(lengthPoints)){
    scenicOptions <- reinitializeScenicOptions(scenicOptions)
    print(paste0("running length ", lengthPoints[i]))
    expMat <- as.matrix(filterMergedData(mergedData, maxCellsPerCondition = lengthPoints[i], filter = TRUE, scenicOptions = scenicOptions))
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
#resets scenic options before a run (clean data)
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

#used to change the target directory for SCENIC output
setOutputDirectory <- function(scenicOptions, outputDirectory = NA){
  if(!(is.na(outputDirectory) | outputDirectory=="")){
    if(!dir.exists(outputDirectory)){
      dir.create(outputDirectory)
    }
    if(!dir.exists(file.path(outputDirectory,"int"))){
      dir.create(file.path(outputDirectory,"int"))
    }
    if(!dir.exists(file.path(outputDirectory,"out"))){
      dir.create(file.path(outputDirectory,"out"))
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
