################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### April 20th, 2015
### designed to be executed with SARTools 1.2.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "analysis/projects/NeuroDiff-RNA-seq/30_SARTools_Ctrl"      # working directory for the R session

projectName <- "NeuroDiff-RNA-seq"                         # name of the project
author <- "Arvind Sundaram"                                # author of the statistical analysis/report

targetFile <- "design.txt"                           # path to the design/target file
rawDir <- "analysis/projects/NeuroDiff-RNA-seq/20_featureCounts"                                      # path to the directory containing raw counts files
featuresToRemove <- NULL # NULL if no feature to remove

varInt <- "groups"                                    # factor of interest
condRef <- "Ctrl0"                                      # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library(SARTools)

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- round(loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove))

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)