

#Halbrook 1229 New analysis April 2022

Epi_1229 <- subset(PNM_Epi, idents = "1229", invert=F)

#Re normalize and re scale 1238
Epi_1229 <- NormalizeData(object = Epi_1229, normalization.method = "LogNormalize", scale.factor = 10000)
Epi_1229 <- FindVariableFeatures(object = Epi_1229, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
Epi_1229 <- ScaleData(object = Epi_1229, vars.to.regress = "nCount_RNA", features = rownames(Epi_1229))
Epi_1229 <- RunPCA(object = Epi_1229, pc.genes = Epi_1229@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#Find number of PCs that gives 90% variance
st_dev <- Epi_1229@reductions$pca@stdev
var <- st_dev^2
sum(var[1:36])/sum(var)
#Find clusters
Epi_1229 <- FindNeighbors(object = Epi_1229, dims = 1:36, save.SNN = TRUE, force.recalc = T)
Epi_1229 <- FindClusters(object = Epi_1229, resolution = 1.2, verbose = F)
#Run UMAP
Epi_1229 <- RunUMAP(object = Epi_1229, dims = 1:36)
DimPlot(object = Epi_1229, reduction = "umap", pt.size = 1.7, label = T)
Epi_1229[["Epi_1229_clusters"]]<- Idents(object = Epi_1229)
FeaturePlot(Epi_1229, features = c("MUC1","KRT18",'KRT19','PRSS1',"CD2",'CD14',"PTPRC"))
VlnPlot(Epi_1229, features = c("PTPRC"))
save(Epi_1229, file = 'Epi_1229.RData')
#LABEL THE CLUSTERS round 1 for collapsed populations
Idents(object = Epi_1229) <- 'seurat_clusters'
levels(Epi_1229)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids <- c("E","E","E","E","E","J","E","J",'J','E')
names(x = new.cluster.ids) <- levels(x = Epi_1229)
Epi_1229 <- RenameIdents(object = Epi_1229, new.cluster.ids)
#Put these new cluster labels as metadata
Epi_1229[["Epi_labels"]] <- Idents(object = Epi_1229)
levels(Epi_1229)
#__________________________________ROUND1_______________#
#subset out Immune cells
Epi_1229_fil <- subset(Epi_1229, idents = "E")
#Re normalize and re scale 1238
Epi_1229_fil <- NormalizeData(object = Epi_1229_fil, normalization.method = "LogNormalize", scale.factor = 10000)
Epi_1229_fil <- FindVariableFeatures(object = Epi_1229_fil, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, y.cutoff = 0.5)
Epi_1229_fil <- ScaleData(object = Epi_1229_fil, vars.to.regress = "nCount_RNA", features = rownames(Epi_1229_fil))
Epi_1229_fil <- RunPCA(object = Epi_1229_fil, pc.genes = Epi_1229_fil@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#Find number of PCs that gives 90% variance
st_dev <- Epi_1229_fil@reductions$pca@stdev
var <- st_dev^2
sum(var[1:36])/sum(var)
#Find clusters
Epi_1229_fil <- FindNeighbors(object = Epi_1229_fil, dims = 1:36, save.SNN = TRUE, force.recalc = T)
Epi_1229_fil <- FindClusters(object = Epi_1229_fil, resolution = 1.2, verbose = F)
#Run UMAP
Epi_1229_fil <- RunUMAP(object = Epi_1229_fil, dims = 1:36)
DimPlot(object = Epi_1229_fil, reduction = "umap", pt.size = 1.7, label = F)
Epi_1229_fil[["Epi_1229_fil_clusters"]]<- Idents(object = Epi_1229_fil)
FeaturePlot(Epi_1229_fil, features = c("MUC1","KRT18",'KRT19','KRT8'), cols = c("grey80", "darkred"))
save(Epi_1229_fil, file = 'Epi_1229_fil.RData')
#run top ten to determine what is what
Epi_1229_fil.markers <- FindAllMarkers(Epi_1229_fil)
write.csv(Epi_1229_fil.markers, file = "~/Desktop/1238.csv")
top_10 <- Epi_1229_fil.markers%>%group_by(cluster)%>%top_n(n=10, wt = avg_logFC)
DoHeatmap(Epi_1229_fil, features = top_10$gene, size = 5.5)


#Re-run with less clusters__________________________________
Epi_1229_fil <- FindClusters(object = Epi_1229_fil, resolution = 0.25, verbose = F)
#Run UMAP
Epi_1229_fil <- RunUMAP(object = Epi_1229_fil, dims = 1:36)
DimPlot(object = Epi_1229_fil, reduction = "umap", pt.size = 1.7, label = F, cols = c("Cluster_1"='red',
                                                                                      "Cluster_2"='darkblue',
                                                                                      "Cluster_3"='lightblue'))
Epi_1229_fil[["Epi_1229_fil_clusters_3"]]<- Idents(object = Epi_1229_fil)
Idents(object = Epi_1229_fil) <- 'Epi_1229_fil_clusters_3'
levels(Epi_1229_fil)
current.cluster.ids <- c(0,1,2)
new.cluster.ids <- c("Cluster_1","Cluster_2","Cluster_3")
names(x = new.cluster.ids) <- levels(x = Epi_1229_fil)
Epi_1229_fil <- RenameIdents(object = Epi_1229_fil, new.cluster.ids)
#Put these new cluster labels as metadata
Epi_1229_fil[["Epi_labels"]] <- Idents(object = Epi_1229_fil)
levels(Epi_1229_fil)
save(Epi_1229_fil, file = 'Epi_1229_fil_4_clusters.RData')
#run top ten to determine what is what
Epi_1229_fil.markers <- FindAllMarkers(Epi_1229_fil)
write.csv(Epi_1229_fil.markers, file = "~/Desktop/1238.csv")
top_10 <- Epi_1229_fil.markers%>%group_by(cluster)%>%top_n(n=10, wt = avg_logFC)
DoHeatmap(Epi_1229_fil, features = top_10$gene, size = 5.5)
FeaturePlot(Epi_1229_fil, features = c("IL8",'CXCL1'), cols = c("grey80",'red4'), pt.size = 1)
VlnPlot(Epi_1229_fil, features = c("IL8",'CXCL1'), pt.size = 0)


#simple plot for heatmap analysis-v2
levels(Epi_1229_fil)
gene_list <- c("CHAC1",'ASNS','PSAT1','XPOT','DDIT3','TRIB3','ATF3','MTHFD2','SLC1A5','ATF4','SESN2', 'SLC7A11','ATF5','PHGDH','ALDH1L2','SLC7A5','SHMT2','GADD45A','PSPH','GADD45G','SREBF1','DDIT4')
t_cell_data <- FetchData(Epi_1229_fil, vars = c(gene_list, 'Epi_labels'))
#The following requested variables were not found: GARS1, IARS1, NARS1, WARS1

patient_avg <- data.frame()
n <- 1 # number of metadata columns, if this is 1, you need to be careful, to line#19, add as.data.frame() around the code

for (id in levels(factor(t_cell_data$Epi_labels))) {
  data_subset <- t_cell_data %>% filter(Epi_labels == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-n)], 2, mean)
  patient_avg <- rbind(patient_avg, data_subset_avg)
}

colnames(patient_avg) <- colnames(t_cell_data)[1:(ncol(t_cell_data)-n)]
rownames(patient_avg) <- levels(factor(t_cell_data$Epi_labels))
# adding Simple annotations 

metadata <- as.data.frame(unique(t_cell_data %>% select("Epi_labels")))
metadata <- metadata[order(metadata$Epi_labels), , drop = F]
rownames(metadata) <- metadata$Epi_labels

pheatmap(as.matrix(patient_avg), fontsize = 14, color = colorRampPalette(colors = c('#0000FF','#FFFFFF','#FF0000'))(250), annotation_row = metadata,
         border_color = 'black', cellwidth = 20, cellheight = 20,scale = "column", cluster_rows = T, cluster_cols = T)


