#SingleR Analysis for LaMannoBrainData, Human Primary Cell Atlas

library(SingleR)
library(celldex)
library(Seurat)

rnaobjects.qcFilt<–ReadRDS(Diff_Merged_Final_EndofPipe_Final_rnaobjects.qcFilt.rds)

rnaobjects.anno <- rnaobjects.qcFilt

# Put normalized counts in an object for singleR
singler.rna <- GetAssayData(rnaobjects.qcFilt, assay = "RNA", slot = "data")

# Datasets
hpca.se <- HumanPrimaryCellAtlasData()

# REF HumanPrimaryCellAtlas
common <- intersect(rownames(singler.rna), rownames(hpca.se))
hpca.se <- hpca.se[common,]
singler.rna <- singler.rna[common,]

pred.hpca <- SingleR(test = singler.rna, ref = hpca.se, 
                     labels = hpca.se$label.main, assay.type.ref = "logcounts")
table(pred.hpca$labels)

# Labels and plot
rnaobjects.anno <- rnaobjects.qcFilt
rnaobjects.anno[["SingleR.labels"]] <- pred.hpca$labels
DimPlot(
  rnaobjects.anno,
  pt.size = 0.5,
  reduction = "umap",
  group.by = "SingleR.labels"
)


# Plots
at1 <- DimPlot(rnaobjects.anno,pt.size = 0.5,reduction = "umap",group.by = "SingleR.labels") + ggtitle(label = "Human Primary Cell Atlas")+ theme(text=element_text(size=16,  family="ArialMT")) + scale_color_brewer(palette = "Set3")
at1
dev.copy(pdf,paste0(currentjob,"_17.FILT.singleR-HPCA_UMAP.pdf"), width=8, height=9, paper='special')
dev.off()


# Heatmaps
plotScoreHeatmap(pred.hpca, max.labels = 5, cellwidth=0.1,main="Human Primary cell Atlas")+ theme(text=element_text(size=5,  family="Arial"))
dev.copy(pdf,paste0(currentjob,"_16.2.FILT.singleR-HCA-labels-heatmap.pdf"), width=20, height=15, paper='special')
dev.off()

table(pred.hpca$labels)

knitr::kable(table(pred.hpca$labels), col.names = c("Cell type: HumanPrimaryCellAtlasData", "Freq"))

#Trying different references: 
library(scRNAseq)


# REF LaMannoBrainData  humn es
hESCs = LaMannoBrainData("human-es", ensembl = FALSE)
hESCs <- logNormCounts(hESCs)

singler.hESCs <- GetAssayData(rnaobjects.qcFilt, assay = "RNA", slot = "data")
common <- intersect(rownames(singler.hESCs), rownames(hESCs))
hESCs <- hESCs[common,]
singler.hESCs <- singler.hESCs[common,]
pred.hESCs <- SingleR(test = singler.hESCs, ref = hESCs, 
                      labels = hESCs$Cell_type, assay.type.ref = "logcounts")
table(pred.hESCs$labels)
rnaobjects.anno[["SingleR.hESCs"]] <- pred.hESCs$labels

at2 <- DimPlot(rnaobjects.anno,pt.size = 0.5,reduction = "umap",group.by = "SingleR.hESCs") + ggtitle(label = "LaMannoBrainData - Human ES")+ theme(text=element_text(size=16,  family="ArialMT")) + scale_color_manual(values = c('#D37295','#FF9888','#CC6666','#996666','#C3CE3D','#8CC2CA','#6666FF','#59A14F','#FFC966', '#FFFFCC','#993399','#996699','#663366','
                                                                                                                                                                                                                                   #A3ACB9','#F28E2B'))

# REF LaMannoBrainData human embryo
hEMBs = LaMannoBrainData("human-embryo", ensembl = FALSE)
hEMBs <- logNormCounts(hEMBs)

singler.hEMBs <- GetAssayData(rnaobjects.qcFilt, assay = "RNA", slot = "data")
common <- intersect(rownames(singler.hEMBs), rownames(hEMBs))
hEMBs <- hEMBs[common,]
singler.hEMBs <- singler.hEMBs[common,]
pred.hEMBs <- SingleR(test = singler.hEMBs, ref = hEMBs, 
                      labels = hEMBs$Cell_type, assay.type.ref = "logcounts")
table(pred.hEMBs$labels)
rnaobjects.anno[["SingleR.hEMBs"]] <- pred.hEMBs$labels

at3 <- DimPlot(rnaobjects.anno,pt.size = 0.5,reduction = "umap",group.by = "SingleR.hEMBs") + ggtitle(label = "LaMannoBrainData - Human embryo")+ theme(text=element_text(size=16,  family="ArialMT")) + scale_color_manual(values = c('#FFC966','#8CC2CA','#B15928','#59A14F','#D37295','#FF9888','#DAB6AF','#8175AA','#023858','#C3CE3D','#008080','#C0D6E4','#FFD700','#A3ACB9','#F28E2B'))


at2 / at3
dev.copy(pdf,paste0(currentjob,"_18.FILT.singleR-LMBD_UMAP.pdf"), width=7, height=9, paper='special')
dev.off()

#Merge stem cells with embryonal data
singler.hESBD = GetAssayData(rnaobjects.qcFilt, assay = "RNA", slot = "data")

common <- intersect(intersect(rownames(hEMBs), rownames(hESCs)),rownames(singler.hESBD))
hEMBs <- hEMBs[common,]
hESCs <- hESCs[common,]
singler.hESBD <- singler.hESBD[common,]

hESBD = cbind(hEMBs,hESCs )
pred.hESBD <- SingleR(test = singler.hESBD, ref = hESBD, 
                      labels = hESBD$Cell_type, assay.type.ref = "logcounts")
table(pred.hESBD$labels)
rnaobjects.anno[["SingleR.hESBD"]] <- pred.hESBD$labels

at4 <- DimPlot(rnaobjects.anno,pt.size = 0.5,reduction = "umap",group.by = "SingleR.hESBD") + ggtitle(label = "LaMannoBrainData - Human embryo")+ theme(text=element_text(size=16,  family="ArialMT")) + scale_color_manual(values = c('#D37295','#FF9888','#CC6666','#996666','#C3CE3D','#8CC2CA','#6666FF','#59A14F','#FFC966', '#FFFFCC','#993399','#996699','#663366','#FFC966','#8CC2CA','#B15928','#59A14F','#A3ACB9','#F28E2B'))
at4
dev.copy(pdf,paste0(currentjob,"_19.FILT.singleR-merged_LMBD_UMAP.pdf"), width=6, height=4.5, paper='special')
dev.off()

plotScoreHeatmap(pred.hESBD, max.labels = 40, main="LaMannoBrainData")+ theme(text=element_text(size=5,  family="Arial"))
dev.copy(pdf,paste0(currentjob,"_19.FILT.singleR-LMBD-labels-heatmap.pdf"), width=7, height=5, paper='special')
dev.off()

knitr::kable(table(pred.hESBD$labels), col.names = c("Cell type: LaMannoBrainData", "Freq»))
