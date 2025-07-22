setwd("/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Analysis")

library(dplyr)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(Azimuth)
library(EnhancedVolcano)
library(SingleCellExperiment)
library(tidyr)
library(pheatmap)


tam.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.1_250109/FilteredMatrix/filtered_matrix_tam/force_5000")
bex.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.1_250109/FilteredMatrix/filtered_matrix_bex/force_5000")
tambex.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.1_250109/FilteredMatrix/filtered_matrix_tambex/force_5000")
vehicle.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.1_250109/FilteredMatrix/filtered_matrix_vehicle/force_5000")

naive1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_naive/force_10000")
cpz1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_cpz/force_10000")
tam1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_tam/force_10000")
bex1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_bex/force_10000")
tambex1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_tambex/force_10000")
vehicle1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix3/filtered_matrix_vehicle/force_10000")

naive2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_naive/force_10000")
cpz2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_cpz/force_10000")
tam2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_tam/force_10000")
bex2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_bex/force_10000")
tambex2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_tambex/force_10000")
vehicle2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P3.2_250209/FilteredMatrix4/filtered_matrix_vehicle/force_10000")

naive3.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_Naive/force_5000")
cpz3.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_CPZ/force_5000")
onewkveh1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_1wkVeh/force_5000")
onewktambex1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_1wkTamBex/force_5000")
twowkveh1.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_2wkVeh/force_5000")
twowktambex1.data <- Read10X(data.dir = "//Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix/filtered_matrix_2wkTamBex/force_5000")

naive4.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_Naive/force_5000")
cpz4.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_CPZ/force_5000")
onewkveh2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_1wkVeh/force_5000")
onewktambex2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_1wkTamBex/force_5000")
twowkveh2.data <- Read10X(data.dir = "/Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_2wkVeh/force_5000")
twowktambex2.data <- Read10X(data.dir = "//Users/Scripps/Desktop/Project_Folder/Project_2/Single_Cell/PIPseq_P4.1_250311/Filtered_Matrix_2/filtered_matrix_2wkTamBex/force_5000")



#Merge dataset run on different machines with same sample
naiveb.data <- naive1.data + naive2.data
cpzb.data <- cpz1.data + cpz2.data
tamb.data <- tam1.data + tam2.data
bexb.data <- bex1.data + bex2.data
tambexb.data <- tambex1.data + tambex2.data
vehb.data <- vehicle1.data + vehicle2.data

naivec.data <- naive3.data + naive4.data
cpzc.data <- cpz3.data + cpz4.data
onewkVeh.data <- onewkveh1.data + onewkveh2.data
onewkTamBex.data <- onewktambex1.data + onewktambex2.data
twowkVeh.data <- twowkveh1.data + twowkveh2.data
twowkTamBex.data <- twowktambex1.data + twowktambex2.data

tamA <- CreateSeuratObject(counts = tam.data, project = "tamA")
bexA <- CreateSeuratObject(counts = bex.data, project = "bexA")
tambexA <- CreateSeuratObject(counts = tambex.data, project = "tambexA")
vehA <- CreateSeuratObject(counts = vehicle.data, project = "vehA")

naiveB <- CreateSeuratObject(counts = naiveb.data, project = "naiveB")
cpzB <- CreateSeuratObject(counts = cpzb.data, project = "cpzB")
tamB <- CreateSeuratObject(counts = tamb.data, project = "tamB")
bexB <- CreateSeuratObject(counts = bexb.data, project = "bexB")
tambexB <- CreateSeuratObject(counts = tambexb.data, project = "tambexB")
vehB <- CreateSeuratObject(counts = vehb.data, project = "vehB")

naiveC <- CreateSeuratObject(counts = naivec.data, project = "naiveC")
cpzC <- CreateSeuratObject(counts = cpzc.data, project = "cpzC")
veh1wk <- CreateSeuratObject(counts = onewkVeh.data, project = "veh1wk")
tambex1wk <- CreateSeuratObject(counts = onewkTamBex.data, project = "tambex1wk")
vehC <- CreateSeuratObject(counts = twowkVeh.data, project = "vehC")
tambexC <- CreateSeuratObject(counts = twowkTamBex.data, project = "tambexC")

naiveB$condition <- "naive"
naiveC$condition <- "naive"
cpzB$condition <- "cpz"
cpzC$condition <- "cpz"
tamA$condition <- "tam"
tamB$condition <- "tam"
bexA$condition <- "bex"
bexB$condition <- "bex"
tambexA$condition <- "tambex"
tambexB$condition <- "tambex"
tambexC$condition <- "tambex"
tambex1wk$condition <- "tambex1wk"
vehA$condition <- "veh"
vehB$condition <- "veh"
vehC$condition <- "veh"
veh1wk$condition <- "veh1wk"

tamA[["percent.mito"]] <- PercentageFeatureSet(tamA, pattern = "^mt-")
bexA[["percent.mito"]] <- PercentageFeatureSet(bexA, pattern = "^mt-")
tambexA[["percent.mito"]] <- PercentageFeatureSet(tambexA, pattern = "^mt-")
vehA[["percent.mito"]] <- PercentageFeatureSet(vehA, pattern = "^mt-")
naiveB[["percent.mito"]] <- PercentageFeatureSet(naiveB, pattern = "^mt-")
cpzB[["percent.mito"]] <- PercentageFeatureSet(cpzB, pattern = "^mt-")
tamB[["percent.mito"]] <- PercentageFeatureSet(tamB, pattern = "^mt-")
bexB[["percent.mito"]] <- PercentageFeatureSet(bexB, pattern = "^mt-")
tambexB[["percent.mito"]] <- PercentageFeatureSet(tambexB, pattern = "^mt-")
vehB[["percent.mito"]] <- PercentageFeatureSet(vehB, pattern = "^mt-")
naiveC[["percent.mito"]] <- PercentageFeatureSet(naiveC, pattern = "^mt-")
cpzC[["percent.mito"]] <- PercentageFeatureSet(cpzC, pattern = "^mt-")
tambexC[["percent.mito"]] <- PercentageFeatureSet(tambexC, pattern = "^mt-")
vehC[["percent.mito"]] <- PercentageFeatureSet(vehC, pattern = "^mt-")
tambex1wk[["percent.mito"]] <- PercentageFeatureSet(tambex1wk, pattern = "^mt-")
veh1wk[["percent.mito"]] <- PercentageFeatureSet(veh1wk, pattern = "^mt-")

all.comb <- merge(
  tamA,
  y = c(bexA, tambexA, vehA, tamB, bexB, tambexB, vehB, naiveB, cpzB, naiveC, cpzC, tambexC, vehC, tambex1wk, veh1wk),
  add.cell.ids = c("tamA", "bexA", "tambexA", "vehA", "tamB", "bexB", "tambexB", "vehB", "naiveB", "cpzB", "naiveC", "cpzC", "tambexC", "vehC", "tambex1wk", "veh1wk"),
  project = "allmerged"
)
all.comb
#An object of class Seurat 
#32285 features across 110000 samples within 1 assay 
#Active assay: RNA (32285 features, 0 variable features)
#16 layers present: counts.tamA, counts.bexA, counts.tambexA, counts.vehA, counts.tamB, counts.bexB, counts.tambexB, counts.vehB, counts.naiveB, counts.cpzB, counts.naiveC, counts.cpzC, counts.tambexC, counts.vehC, counts.tambex1wk, counts.veh1wk
table(all.comb$orig.ident)
# bexA      bexB      cpzB      cpzC    naiveB    naiveC      tamA      tamB tambex1wk   tambexA   tambexB 
#5000     10000     10000      5000     10000      5000      5000     10000      5000      5000     10000 
#tambexC    veh1wk      vehA      vehB      vehC 
#5000      5000      5000     10000      5000 

table(all.comb$condition)
#bex       cpz     naive       tam    tambex tambex1wk       veh    veh1wk 
#15000     15000     15000     15000     20000      5000     20000      5000 


VlnPlot(
  all.comb, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
  ncol = 3, 
  raster = FALSE
)
v1 <- VlnPlot(all.comb, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + coord_cartesian(ylim = c(0, 1000))


ggsave("Raw_nfeat.pdf")
save(all.comb, file = "allcomb.RData")
# If the file is saved with the 'save()' function, use 'load()'
load("allcomb.Rdata")

#Subset to get rid of trash cells
allcomb_400qc <- subset(all.comb, subset = nFeature_RNA > 400 & percent.mito < 25)

allcomb_400qc <- NormalizeData(allcomb_400qc)
allcomb_400qc <- FindVariableFeatures(allcomb_400qc)
allcomb_400qc <- ScaleData(allcomb_400qc)
allcomb_400qc <- RunPCA(allcomb_400qc)
ElbowPlot(object = allcomb_400qc, ndims = 30)

allcomb_400qc <- FindNeighbors(allcomb_400qc, dims = 1:30, reduction = "pca")
allcomb_400qc <- FindClusters(allcomb_400qc, resolution = 0.3, cluster.name = "unintegrated_clusters")
allcomb_400qc <- RunUMAP(allcomb_400qc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

load("allcomb400.Rdata")

#Identify low quality clusters
DimPlot(
  allcomb_400qc,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

VlnPlot(
  allcomb_400qc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
  ncol = 3, 
  pt.size = 0,
  raster = FALSE
)
VlnPlot(
  allcomb_400qc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
  group.by = "orig.ident", 
  ncol = 3, 
  pt.size = 0,
  raster = FALSE
)

v2 <- VlnPlot(all.comb, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + coord_cartesian(ylim = c(0, 5000))
v2


#Subset to get rid of trash clusters
obj <- subset(allcomb_400qc, idents = c(0), invert = TRUE)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(object = obj, ndims = 30)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "unintegrated_clusters")
obj <- JoinLayers(obj)
obj <- RunAzimuth(obj, reference = "mousecortexref")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
save(obj, file = "obj.RData")
# If the file is saved with the 'save()' function, use 'load()'
load("obj.Rdata")


DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "predicted.subclass",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

VlnPlot(
  obj, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
  ncol = 3, 
  raster = FALSE
)

DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "orig.ident",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "condition",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

save(obj, file = "obj.RData")
load("obj.Rdata")

# Create a table of cluster counts by 'orig.ident'
table(obj$seurat_clusters, obj$orig.ident)


#Use this to test clusters
FeaturePlot(
  object = obj, 
  features = c('Aif1','Cx3cr1','P2ry12','Tmem119','Ctss','C1qa'),
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)
VlnPlot(obj, features = c('Aif1','Cx3cr1','P2ry12','Tmem119','Ctss','C1qa'), group.by = "seurat_clusters")


UMAP_cluster_ids <- c(
  "0" = "Oligodendrocytes",
  "1" = "L2/3IT",
  "2" = "L5IT",
  "3" = "L6IT",
  "4" = "L5ET",
  "5" = "Microglia",
  "6" = "Astrocytes", 
  "7" = "L6CT",
  "8" = "VIP/Lamp5IN",
  "9" = "L2/3IT",
  "10" = "SST/PVIN",
  "11" = "Oligodendrocytes",
  "12" = "OPC",
  "13" = "L5/6NP",
  "14" = "Microglia",
  "15" = "Astrocytes", 
  "16" = "VLMC",
  "17" = "Endo",
  "18" = "VLMC",
  "19" = "Pericytes", 
  "20" = "L2/3IT",
  "21" = "Microglia")

names(UMAP_cluster_ids) <- levels(obj)
obj <- RenameIdents(obj, UMAP_cluster_ids)

DimPlot(
  obj,
  reduction = "umap.unintegrated",
  group.by = "condition",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)


##SUBSET OLOPC
OLOPC <- subset(obj, seurat_clusters %in% c("0", "11", "12", "14"))
VlnPlot(OLOPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)  
#If you need to see closer
v4 <- VlnPlot(OLOPC, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) + coord_cartesian(ylim = c(0, 3000))
v5 <- VlnPlot(OLOPC, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0) + coord_cartesian(ylim = c(0, 5000))
v6 <- VlnPlot(OLOPC, features = "percent.mito", group.by = "orig.ident", pt.size = 0)


OLOPC <- subset(OLOPC, nFeature_RNA < 3000)
table(OLOPC$orig.ident)


#Normalize
OLOPC <- NormalizeData(OLOPC)
OLOPC <- FindVariableFeatures(OLOPC)
OLOPC <- ScaleData(OLOPC)
OLOPC <- RunPCA(OLOPC)
ElbowPlot(object = OLOPC, ndims = 30)
ggsave("elbow_OLOPC.pdf")
OLOPC <- FindNeighbors(OLOPC, dims = 1:15, reduction = "pca")
OLOPC <- FindClusters(OLOPC, resolution = 0.3, cluster.name = "unintegrated_clusters")
OLOPC <- RunUMAP(OLOPC, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(
  OLOPC,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC,
  reduction = "umap.unintegrated",
  group.by = "orig.ident",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC,
  reduction = "umap.unintegrated",
  group.by = "condition",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

FeaturePlot(
  object = OLOPC, 
  features = c("Mbp", "Pdgfra"),
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)
VlnPlot(OLOPC, features = c("Mbp", "Pdgfra"), group.by = "seurat_clusters")


OLOPC_int <- IntegrateLayers(
  object = OLOPC, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

VlnPlot(OLOPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
        ncol = 3, group.by = "seurat_clusters", pt.size = 0) 
table(OLOPC$seurat_clusters)
#Remove bad seurat clusters and renormalize if needed


OLOPC_19 <- JoinLayers(OLOPC)

save(OLOPC_19, file = "OLOPC_19.RData")
##March 19,2025
load("OLOPC_19.Rdata")

table(OLOPC_19$seurat_clusters, OLOPC_19$orig.ident)


#Identifying what the subclusters are for OLOPC
OLOPC_19_fm <- FindAllMarkers(
  object = OLOPC_19,  # Your Seurat object
  only.pos = TRUE,  # Include both upregulated and downregulated genes
  min.pct = 0.1,  # Only consider genes expressed in at least 10% of cells
  logfc.threshold = 0.25  # Minimum log fold change threshold
)
write.csv(OLOPC_19_fm, "OLOPC_19_findmarkers.csv", row.names = FALSE)

johkietal_genes <- read.csv("johkietal_genes.csv", header = TRUE)
johkietal_genes_list <- johkietal_genes$Gene_list
valid_genes_johki <- johkietal_genes_list[johkietal_genes_list %in% Features(OLOPC_19, assay = "RNA")]
OLOPC_19 <- ScaleData(OLOPC_19, features = valid_genes_johki)
DoHeatmap(OLOPC_19, features = valid_genes_johki, group.by = "seurat_clusters")

FeaturePlot(
  object = OLOPC_19, 
  features = c('Aif1','Cx3cr1','P2ry12','Tmem119','Ctss','C1qa'),
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)

markers_FinalOPC <- c('Pdgfra', 'Cspg5', 'Kcnip3', 'Matn4', 'Serpine2')
markers_FinalCOP <- c('Gpr17', 'Bcas1', 'Fyn', 'S100b', 'Bmp4')
markers_FinalNFOL <- c('9630013A20Rik', 'Sema6a', 'Eml1', 'Ttyh1', 'Tspan2')
markers_FinalNFOL2 <- c('Tmem2', 'Cnksr3', 'Rras2', 'Itpr2', 'Sema4d')
markers_FinalMFOL <- c('Ctps', 'Opalin', 'Ttyh1', 'Mcam', 'Clic4')
markers_FinalMFOL2 <- c('Ccp110', 'Hhip', 'Tmem88b', 'Arap2', 'Efhd1')
markers_FinalMOL <- c('Mbp', 'Ptgds', 'Il33', 'Plp1', 'Mag')

VlnPlot(OLOPC_19, features = markers_FinalMOL, group.by = "seurat_clusters", ncol = 5, pt.size = 0)

FeaturePlot(
  object = OLOPC_19, 
  features = markers_FinalOPC,
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)

##Clusters 3,6 look more like microglia

OLOPC_19_fix <- subset(OLOPC_19, idents = c(3, 6), invert = TRUE)
OLOPC_19_fix <- NormalizeData(OLOPC_19_fix)
OLOPC_19_fix <- FindVariableFeatures(OLOPC_19_fix)
OLOPC_19_fix <- ScaleData(OLOPC_19_fix)
OLOPC_19_fix <- RunPCA(OLOPC_19_fix)
ElbowPlot(object = OLOPC_19_fix, ndims = 30)
ggsave("elbow_OLOPC_19_fix.pdf")
OLOPC_19_fix <- FindNeighbors(OLOPC_19_fix, dims = 1:15, reduction = "pca")
OLOPC_19_fix <- FindClusters(OLOPC_19_fix, resolution = 0.3, cluster.name = "unintegrated_clusters")
OLOPC_19_fix <- RunUMAP(OLOPC_19_fix, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(
  OLOPC_19_fix,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC_19_fix,
  reduction = "umap.unintegrated",
  group.by = "orig.ident",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC_19_fix,
  reduction = "umap.unintegrated",
  group.by = "condition",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

FeaturePlot(
  object = OLOPC_19_fix, 
  features = c("Mbp", "Pdgfra"),
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)

VlnPlot(OLOPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
        ncol = 3, group.by = "seurat_clusters", pt.size = 0) 
table(OLOPC$seurat_clusters)
#Remove bad seurat clusters and renormalize if needed


OLOPC_19 <- JoinLayers(OLOPC)

save(OLOPC_19_fix, file = "OLOPC_19_fix.RData")
##March 19,2025
load("OLOPC_19_fix.Rdata")

table(OLOPC_19_fix$seurat_clusters, OLOPC_19_fix$orig.ident)


#Identifying what the subclusters are for OLOPC
OLOPC_19_fix_fm <- FindAllMarkers(
  object = OLOPC_19_fix,  # Your Seurat object
  only.pos = TRUE,  # Include both upregulated and downregulated genes
  min.pct = 0.1,  # Only consider genes expressed in at least 10% of cells
  logfc.threshold = 0.25  # Minimum log fold change threshold
)
write.csv(OLOPC_19_fix_fm, "OLOPC_20_findmarkers.csv", row.names = FALSE)

johkietal_genes <- read.csv("johkietal_genes.csv", header = TRUE)
johkietal_genes_list <- johkietal_genes$Gene_list
valid_genes_johki <- johkietal_genes_list[johkietal_genes_list %in% Features(OLOPC_19_fix, assay = "RNA")]
OLOPC_19_fix <- ScaleData(OLOPC_19_fix, features = valid_genes_johki)
DoHeatmap(OLOPC_19_fix, features = valid_genes_johki, group.by = "ident")

VlnPlot(OLOPC_19_fix, features = markers_FinalMOL, group.by = "seurat_clusters", ncol = 5, pt.size = 0)

table(OLOPC_19_fix$seurat_clusters, OLOPC_19_fix$orig.ident)

OL_cluster_ids <- c(
  "0" = "MOL",
  "1" = "MOL",
  "2" = "OPC",
  "3" = "MOL_DA",
  "4" = "MOL_DA",
  "5" = "COP",
  "6" = "NFOL",
  "7" = "COP")

names(OL_cluster_ids) <- levels(OLOPC_19_fix)
OLOPC_19_fix <- RenameIdents(object = OLOPC_19_fix, OL_cluster_ids)


DimPlot(OLOPC_19_fix, label = TRUE, label.size = 6, repel = TRUE) + NoLegend()

Idents(OLOPC_19_fix) <- factor(Idents(OLOPC_19_fix), levels = order_cluster)
DoHeatmap(OLOPC_19_fix, features = valid_genes_johki, group.by = "ident") +
  theme(axis.text.x = element_text(size = 4))


##write gene list with biological relevance for each subcluster in csv called "Pettyetal_genes"
pettyetal_genes <- read.csv("Pettyetal_genes.csv", header = TRUE)
pettyetal_genes_list <- pettyetal_genes$Gene_list
valid_genes_petty <- pettyetal_genes_list[pettyetal_genes_list %in% Features(OLOPC_19_fix, assay = "RNA")]
OLOPC_19_fix <- ScaleData(OLOPC_19_fix, features = valid_genes_petty)
DoHeatmap(OLOPC_19_fix, features = valid_genes_petty, group.by = "seurat_clusters")
DoHeatmap(OLOPC_19_fix, features = valid_genes_petty, group.by = "ident") +
  theme(axis.text.x = element_text(size = 4))

#Bar graph for percentages
order_origident <- c("cpzB", "cpzC", "veh1wk", "tambex1wk", "vehA","vehB","vehC","tamA", "tamB", "bexA", "bexB", "tambexA", "tambexB", "tambexC", "naiveB", "naiveC")
order_cluster <- c("OPC", "COP", "NFOL", "MOL_DA", "MOL")  # Replace with your desired order for seurat_clusters


OLOPC_19_fix$orig.ident <- factor(OLOPC_19_fix$orig.ident, levels = order_origident)
Idents(OLOPC_19_fix) <- factor(Idents(OLOPC_19_fix), levels = order_cluster)
OLOPC_19_fix$cluster <- Idents(OLOPC_19_fix)  # Store identities in meta.data


# Step 3: Create a data frame of the count of subclusters per `orig.ident`
subcluster_data <- OLOPC_19_fix@meta.data %>%
  group_by(orig.ident, cluster) %>%
  tally() %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()


# Step 4: Create the stacked bar plot
ggplot(subcluster_data, aes(x = orig.ident, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  # Create the bar plot
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "Original Identity (orig.ident)", 
       y = "Percentage of Cells (%)", 
       fill = "Subclusters") + 
  scale_fill_manual(values = c("violetred1", "pink", "orange", "yellow", "lightskyblue", "lightblue", "purple", "black")) +  # Customize the colors of the bars
  ggtitle("Distribution of Subclusters by Original Identity (as Percentage)")  # Title for the plot

##Also do this by condition
order_condition <- c("cpz", "veh1wk", "tambex1wk", "veh","tam", "bex", "tambex", "naive")
order_cluster <- c("OPC", "COP", "NFOL", "MOL_DA", "MOL")  # Replace with your desired order for seurat_clusters


OLOPC_19_fix$condition <- factor(OLOPC_19_fix$condition, levels = order_condition)
Idents(OLOPC_19_fix) <- factor(Idents(OLOPC_19_fix), levels = order_cluster)
OLOPC_19_fix$cluster <- Idents(OLOPC_19_fix)  # Store identities in meta.data


# Step 3: Create a data frame of the count of subclusters per `orig.ident`
subclustercon_data <- OLOPC_19_fix@meta.data %>%
  group_by(condition, cluster) %>%
  tally() %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()


# Step 4: Create the stacked bar plot
ggplot(subclustercon_data, aes(x = condition, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  # Create the bar plot
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Increase x-axis text size and rotate
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) +
  labs(y = "Percentage of Cells (%)", 
       fill = "Subclusters") + 
  scale_fill_manual(values = c("#bf5186", "#dfa8c3", "#bdbf51", "#51bf8a", "#5351bf")) +  # Better color palette
  ggtitle("Distribution of Subclusters by Condition (as Percentage)")


##Dotplot for enrichement

# Example enrichment data frame
enrichment_data <- data.frame(
  Pathway = c("Axon ensheathment", "Myelination", "Neurogenesis", "Oligodendrocyte development", "Oligodendrocyte differentiation", "CNS Myelination"),
  GeneRatio = c(0.2, 0.2, 0.06, 0.24, 0.17, 0.33),
  Pvalue = c(.0000000000173, 0.0000000000403, 0.000000000000000362, 0.00000162, 0.00000091, 0.00000763),
  Count = c(26, 25, 99, 13, 17, 10)
)

# Create dot plot
ggplot(enrichment_data, aes(x = GeneRatio, y = reorder(Pathway, GeneRatio), size = Count, color = Pvalue)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +  # P-value scale
  theme_minimal() +
  xlab("Gene Ratio") +
  ylab("Enriched Pathways") +
  ggtitle("Enrichment Dot Plot") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )



##### MOLS DEGS
TheMOLS <- subset(OLOPC_19_fix, idents = c("MOL", "MOL_DA"))
Idents(TheMOLS) <- "condition"

TheOPCs <- subset(OLOPC_19_fix, idents = c("OPC"))
Idents(TheOPCs) <- "condition"

TheInts <- subset(OLOPC_19_fix, idents = c("COP", "NFOL"))
Idents(TheInts) <- "condition"
#Nothing significant


#Tambex vs vhe
degs1 <- FindMarkers(TheMOLS, 
                     ident.1 = "tambex",  # group 1
                     ident.2 = "veh",  # group 2
                     test.use = "wilcox",  # statistical test (Wilcoxon by default)
                     min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                     logfc.threshold = 0.25)  # Log fold change threshold
head(degs1)
degs1$gene <- rownames(degs1)

degs1vol <- EnhancedVolcano(degs1,
                            lab = rownames(degs1),           # Gene labels
                            x = 'avg_log2FC',                    # Log2 fold change
                            y = 'p_val_adj',                     # Adjusted p-value
                            xlab = bquote(~Log[2]~ "Fold Change"),
                            ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                            pCutoff = 0.05,                      # Significance threshold
                            FCcutoff = 0.5,                      # Fold-change threshold
                            pointSize = 2.0,
                            labSize = 3.0,
                            col = c("grey30", "forestgreen", "royalblue", "red2"),
                            title = "OPCs TamBex (-) vs Veh",
                            subtitle = "Differentially Expressed Genes",
                            caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs1vol + scale_x_continuous(limits = c(-1, 3))

# Save the differential expression results to a CSV file
write.csv(degs1O, file = "deg1O_tambex_veh.csv", row.names = TRUE)

#naive vs cpz
degs2 <- FindMarkers(TheMOLS, 
                     ident.1 = "naive",  # group 1
                     ident.2 = "cpz",  # group 2
                     test.use = "wilcox",  # statistical test (Wilcoxon by default)
                     min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                     logfc.threshold = 0.25)  # Log fold change threshold
head(degs2)
degs2$gene <- rownames(degs2)

degs2vol <- EnhancedVolcano(degs2,
                            lab = rownames(degs2),           # Gene labels
                            x = 'avg_log2FC',                    # Log2 fold change
                            y = 'p_val_adj',                     # Adjusted p-value
                            xlab = bquote(~Log[2]~ "Fold Change"),
                            ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                            pCutoff = 0.05,                      # Significance threshold
                            FCcutoff = 0.5,                      # Fold-change threshold
                            pointSize = 2.0,
                            labSize = 3.0,
                            col = c("grey30", "forestgreen", "royalblue", "red2"),
                            title = "OPCs Naive vs CPZ",
                            subtitle = "Differentially Expressed Genes",
                            caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs2vol + scale_x_continuous(limits = c(-2.5, 2.5))

# Save the differential expression results to a CSV file
write.csv(degs2O, file = "deg2O_naive_cpz.csv", row.names = TRUE)

#naive vs TamBex
degs3O <- FindMarkers(TheOPCs, 
                     ident.1 = "tambex",  # group 1
                     ident.2 = "naive",  # group 2
                     test.use = "wilcox",  # statistical test (Wilcoxon by default)
                     min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                     logfc.threshold = 0.25)  # Log fold change threshold
head(degs3O)
degs3O$gene <- rownames(degs3O)

degs3Ovol <- EnhancedVolcano(degs3O,
                            lab = rownames(degs3O),           # Gene labels
                            x = 'avg_log2FC',                    # Log2 fold change
                            y = 'p_val_adj',                     # Adjusted p-value
                            xlab = bquote(~Log[2]~ "Fold Change"),
                            ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                            pCutoff = 0.05,                      # Significance threshold
                            FCcutoff = 0.5,                      # Fold-change threshold
                            pointSize = 2.0,
                            labSize = 3.0,
                            col = c("grey30", "forestgreen", "royalblue", "red2"),
                            title = "OPCs Naive (-) vs TamBex",
                            subtitle = "Differentially Expressed Genes",
                            caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs3Ovol + scale_x_continuous(limits = c(-4.5, 2.5))

# Save the differential expression results to a CSV file
write.csv(degs3, file = "deg3O_naive_tambex.csv", row.names = TRUE)

#Vhe1wk vs TamBex1sk
degs4O <- FindMarkers(TheOPCs, 
                     ident.1 = "veh1wk",  # group 1
                     ident.2 = "tambex1wk",  # group 2
                     test.use = "wilcox",  # statistical test (Wilcoxon by default)
                     min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                     logfc.threshold = 0.25)  # Log fold change threshold
head(degs4O)
degs4O$gene <- rownames(degs4O)

degs4Ovol <- EnhancedVolcano(degs4O,
                            lab = rownames(degs4O),           # Gene labels
                            x = 'avg_log2FC',                    # Log2 fold change
                            y = 'p_val_adj',                     # Adjusted p-value
                            xlab = bquote(~Log[2]~ "Fold Change"),
                            ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                            pCutoff = 0.05,                      # Significance threshold
                            FCcutoff = 0.5,                      # Fold-change threshold
                            pointSize = 2.0,
                            labSize = 3.0,
                            col = c("grey30", "forestgreen", "royalblue", "red2"),
                            title = "OPCs TamBex1wk (-) vs Veh1wk",
                            subtitle = "Differentially Expressed Genes",
                            caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs4Ovol + scale_x_continuous(limits = c(-2, 1.5))

# Save the differential expression results to a CSV file
write.csv(degs4, file = "deg4_Tambex1wk_Veh1wk.csv", row.names = TRUE)

##Consider Pseudobulk

#For Publication
deg_genes_pub <- read.csv("interestinggenes1.csv", header = TRUE)
head(deg_genes_pub)
gene_listpub <- deg_genes_pub$gene_list
head(gene_listpub)
expression_data_pub <- FetchData(TheMOLS, vars = c(gene_listpub, "orig.ident"))
head(expression_data_pub)
expression_data_long_pub <- gather(expression_data_pub, key = "gene", value = "expression_level", -orig.ident)
head(expression_data_long_pub)
average_expression_by_ident_pub <- expression_data_long_pub %>%
  group_by(orig.ident, gene) %>%
  summarise(average_expression = mean(expression_level, na.rm = TRUE))
head(average_expression_by_ident_pub)
selected_conditions <- c("cpzB", "cpzC", "tambexA", "tambexB", "tambexC", "naiveB", "naiveC")
average_expression_selected_conditions <- average_expression_by_ident_pub %>%
  filter(orig.ident %in% selected_conditions)
head(average_expression_selected_conditions)
write.csv(average_expression_selected_conditions, file = "interestinggenesPubOI.csv", row.names = FALSE)
dg5pub <- ggplot(average_expression_selected_conditions, aes(x = gene, y = average_expression, color = orig.ident, group = orig.ident)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Average Gene Expression by Orig.Ident-Genes of Interest", 
       x = "Gene", 
       y = "Average Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

heatmap_data_5pub <- average_expression_selected_conditions %>%
  pivot_wider(names_from = orig.ident, values_from = average_expression, values_fill = list(average_expression = 0))

desired_order_pub <- c("cpzB", "cpzC", "tambexA", "tambexB", "tambexC", "naiveB", "naiveC")

heatmap_data_5pub <- heatmap_data_5pub %>%
  select(gene, all_of(desired_order_pub))

heatmap_matrix_5_pub <- as.matrix(heatmap_data_5pub[, -1])  # Remove 'gene' column for matrix conversion
rownames(heatmap_matrix_5_pub) <- heatmap_data_5pub$gene  # Set gene names as row names
# Rename the columns (conditions)

# Desired order of columns
desired_order_pub <- c("CPZ", "Naive", "Tam+Bex")

# Reorder the columns in the matrix
heatmap_matrix_5_pub <- heatmap_matrix_5_pub[, match(desired_order_pub, colnames(heatmap_matrix_5_pub))]

# Verify the column order
colnames(heatmap_matrix_5_pub)


g <- pheatmap(
  heatmap_matrix_5_pub,
  cluster_rows = TRUE,  # Cluster rows (genes)
  cluster_cols = FALSE,  # Do not cluster columns (orig.ident)
  scale = "row",        # Scale by row to normalize gene expression
  show_rownames = TRUE,  # Show gene names
  show_colnames = TRUE,  # Show orig.ident names
  fontsize_row = 16,      # Adjust font size for gene names
  fontsize_col = 20,     # Adjust font size for orig.ident names
  main = "Heatmap of MOL expression",
  column_order = match(desired_order_pub, colnames(heatmap_matrix_5_pub)),  # Explicitly set the column order
  angle_col = 45  # Rotate column names to 45 degrees
)
ggsave("try.eps",g)

#Make Gene list and graph heatmap
deg_genes_int <- read.csv("interesting_genes.csv", header = TRUE)
head(deg_genes_int)
gene_listint <- deg_genes_int$gene_list
head(gene_listint)
expression_data_con <- FetchData(TheMOLS, vars = c(gene_listint, "condition"))
head(expression_data_con)
expression_data_long_con <- gather(expression_data_con, key = "gene", value = "expression_level", -condition)
head(expression_data_long_con)
average_expression_by_ident_con <- expression_data_long_con %>%
  group_by(condition, gene) %>%
  summarise(average_expression = mean(expression_level, na.rm = TRUE))
head(average_expression_by_ident_con)
write.csv(average_expression_by_ident_con, file = "interestinggenesexp.csv", row.names = FALSE)
dg5con <- ggplot(average_expression_by_ident_con, aes(x = gene, y = average_expression, color = orig.ident, group = orig.ident)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Average Gene Expression by Orig.Ident-Genes of Interest", 
       x = "Gene", 
       y = "Average Expression Level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
heatmap_data_5con <- average_expression_by_ident_con %>%
  pivot_wider(names_from = condition, values_from = average_expression, values_fill = list(average_expression = 0))

desired_order_con <- c("cpz",
                       "veh1wk",
                       "tambex1wk",
                       "veh", 
                       "bex", 
                       "tam",
                       "tambex", 
                       "naive")

heatmap_data_5con <- heatmap_data_5con %>%
  select(gene, all_of(desired_order_con))

heatmap_matrix_5_con <- as.matrix(heatmap_data_5con[, -1])  # Remove 'gene' column for matrix conversion
rownames(heatmap_matrix_5_con) <- heatmap_data_5con$gene  # Set gene names as row names


pheatmap(
  heatmap_matrix_5_con,
  cluster_rows = TRUE,  # Cluster rows (genes)
  cluster_cols = FALSE,  # Cluster columns (orig.ident)
  scale = "row",        # Scale by row to normalize gene expression
  show_rownames = TRUE,  # Show gene names
  show_colnames = TRUE,  # Show orig.ident names
  fontsize_row = 8,      # Adjust font size for gene names
  fontsize_col = 10,     # Adjust font size for orig.ident names
  main = "Heatmap of MOL expression",
  column_order = match(desired_order_con, colnames(heatmap_matrix_5_con))  # Explicitly set the column order
)


##Neurons
TheMicroglia <- subset(obj, idents = c("Microglia"))
Idents(TheMicroglia) <- "condition"

#Tambex vs vhe
degs1Microglia <- FindMarkers(TheMicroglia, 
                      ident.1 = "veh",  # group 1
                      ident.2 = "tambex",  # group 2
                      test.use = "wilcox",  # statistical test (Wilcoxon by default)
                      min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                      logfc.threshold = 0.25)  # Log fold change threshold
head(degs1Microglia)
degs1Microglia$gene <- rownames(degs1Microglia)

degs1Microgliavol <- EnhancedVolcano(degs1Microglia,
                             lab = rownames(degs1Microglia),           # Gene labels
                             x = 'avg_log2FC',                    # Log2 fold change
                             y = 'p_val_adj',                     # Adjusted p-value
                             xlab = bquote(~Log[2]~ "Fold Change"),
                             ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                             pCutoff = 0.05,                      # Significance threshold
                             FCcutoff = 0.5,                      # Fold-change threshold
                             pointSize = 2.0,
                             labSize = 3.0,
                             col = c("grey30", "forestgreen", "royalblue", "red2"),
                             title = "Microglias TamBex (-) vs Veh",
                             subtitle = "Differentially Expressed Genes",
                             caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs1Microgliavol + scale_x_continuous(limits = c(-5, 5))

# Save the differential expression results to a CSV file
write.csv(degs1Microglia, file = "deg1Microglia_tambex_veh.csv", row.names = TRUE)

#naive vs cpz
degs2Microglia <- FindMarkers(TheMicroglia, 
                      ident.1 = "cpz",  # group 1
                      ident.2 = "naive",  # group 2
                      test.use = "wilcox",  # statistical test (Wilcoxon by default)
                      min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                      logfc.threshold = 0.25)  # Log fold change threshold
head(degs2Microglia)
degs2Microglia$gene <- rownames(degs2Microglia)

degs2Microgliavol <- EnhancedVolcano(degs2Microglia,
                             lab = rownames(degs2Microglia),           # Gene labels
                             x = 'avg_log2FC',                    # Log2 fold change
                             y = 'p_val_adj',                     # Adjusted p-value
                             xlab = bquote(~Log[2]~ "Fold Change"),
                             ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                             pCutoff = 0.05,                      # Significance threshold
                             FCcutoff = 0.5,                      # Fold-change threshold
                             pointSize = 2.0,
                             labSize = 3.0,
                             col = c("grey30", "forestgreen", "royalblue", "red2"),
                             title = "Microglia Naive (-) vs CPZ",
                             subtitle = "Differentially Expressed Genes",
                             caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs2Microgliavol + scale_x_continuous(limits = c(-4, 2))

# Save the differential expression results to a CSV file
write.csv(degs2Microglia, file = "deg2Microglia_naive_cpz.csv", row.names = TRUE)

#naive vs TamBex
degs3Microglia <- FindMarkers(TheMicroglia, 
                      ident.1 = "tambex",  # group 1
                      ident.2 = "naive",  # group 2
                      test.use = "wilcox",  # statistical test (Wilcoxon by default)
                      min.pct = 0.1,        # Only consider genes expressed in at least 10% of cells in either group
                      logfc.threshold = 0.25)  # Log fold change threshold
head(degs3Microglia)
degs3Microglia$gene <- rownames(degs3Microglia)

degs3Microgliavol <- EnhancedVolcano(degs3Microglia,
                             lab = rownames(degs3Microglia),           # Gene labels
                             x = 'avg_log2FC',                    # Log2 fold change
                             y = 'p_val_adj',                     # Adjusted p-value
                             xlab = bquote(~Log[2]~ "Fold Change"),
                             ylab = bquote(~-Log[10]~ "Adjusted P-value"),
                             pCutoff = 0.05,                      # Significance threshold
                             FCcutoff = 0.5,                      # Fold-change threshold
                             pointSize = 2.0,
                             labSize = 3.0,
                             col = c("grey30", "forestgreen", "royalblue", "red2"),
                             title = "Microglias Naive (-) vs TamBex",
                             subtitle = "Differentially Expressed Genes",
                             caption = "Thresholds: Adj p-value < 0.05, |log2FC| > 0.5")
degs3Microgliavol + scale_x_continuous(limits = c(-3, 2.5))

# Save the differential expression results to a CSV file
write.csv(degs3Microglia, file = "deg3Microglia_naive_tambex.csv", row.names = TRUE)


#Remove exp B to look at heterogeneity 
# Remove cells belonging to the specified groups
OLOPC_19_subs <- OLOPC_19_fix[, !(OLOPC_19_fix$orig.ident %in% orig_to_remove)]
table(OLOPC_19_subs$orig.ident)

OLOPC_19_subs <- NormalizeData(OLOPC_19_subs)
OLOPC_19_subs <- FindVariableFeatures(OLOPC_19_subs)
OLOPC_19_subs <- ScaleData(OLOPC_19_subs)
OLOPC_19_subs <- RunPCA(OLOPC_19_subs)
ElbowPlot(object = OLOPC_19_subs, ndims = 30)
ggsave("elbow_OLOPC_19_subs.pdf")
OLOPC_19_subs <- FindNeighbors(OLOPC_19_subs, dims = 1:15, reduction = "pca")
OLOPC_19_subs <- FindClusters(OLOPC_19_subs, resolution = 0.3, cluster.name = "unintegrated_clusters")
OLOPC_19_subs <- RunUMAP(OLOPC_19_subs, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(
  OLOPC_19_subs,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC_19_subs,
  reduction = "umap.unintegrated",
  group.by = "orig.ident",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)
DimPlot(
  OLOPC_19_subs,
  reduction = "umap.unintegrated",
  group.by = "condition",
  label = FALSE, # Adds the subclass name to the plot
  label.size = 4 # Adjusts the font size of the labels
)

FeaturePlot(
  object = OLOPC_19_subs, 
  features = c("Mbp", "Pdgfra"),
  reduction = "umap.unintegrated", 
  cols = c("lightgrey", "red"),
  combine = TRUE
)

VlnPlot(OLOPC_19_subs, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
        ncol = 3, group.by = "seurat_clusters", pt.size = 0) 
table(OLOPC_19_subs$seurat_clusters)
#Remove bad seurat clusters and renormalize if needed


OLOPC_19_subs <- JoinLayers(OLOPC_19_subs)

save(OLOPC_19_subs, file = "OLOPC_19_subs.RData")
##March 19,2025
load("OLOPC_19_subs.Rdata")

table(OLOPC_19_subs$seurat_clusters, OLOPC_19_subs$orig.ident)

OLOPC_19_subs_fm <- FindAllMarkers(
  object = OLOPC_19_subs,  # Your Seurat object
  only.pos = TRUE,  # Include both upregulated and downregulated genes
  min.pct = 0.1,  # Only consider genes expressed in at least 10% of cells
  logfc.threshold = 0.25  # Minimum log fold change threshold
)
write.csv(OLOPC_19_subs_fm, "OLOPC_21_findmarkers.csv", row.names = FALSE)


###Identifiers for MOL subtypes
MOL1_CB <- c("Fos", "Egr1", "Fosb", "Jun", "Junb", "Dnajb1", "Btg2", "Dusp1", "Nr4a1", "Hspa1a")
MOL3_CB <- c("Klk6", "S100b", "Anxa5", "Enpp6", "Sepp1", "Pmp22", "Ldhb", "Pdgfa", "Hapln2", "Anln")
MOL4_CB <- c("Gsn", "Trf", "Apod", "Prdx1", "Polr3e", "Gapdh", "Serpinb1a", "Hapln2", "Pmp22", "Sepp1")
MOL5_CB <- c("Ptgds", "Apoe", "Cd63", "Sgk1", "Htra1")
MOL6_CB <- c("Il33", "Erbb2ip", "AY036118", "Gm42418", "Ttll7", "Dst", "Ccp110", "Ptprd", "Cdc42bpa", "Cntn2")

VlnPlot(OLOPC_19_subs_named, features = c("Zbtb16", "Fkbp5", "Rhoj", "Zfand4", "Mobp"), 
        ncol = 3, group.by = "cluster", pt.size = 0) 
# Assuming OLOPC_19_subs is your Seurat object

# Create violin plots for each set of genes
VlnPlot(OLOPC_19_subs, features = MOL1_CB, ncol = length(MOL1_CB), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(OLOPC_19_subs, features = MOL3_CB, ncol = length(MOL3_CB), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(OLOPC_19_subs, features = MOL4_CB, ncol = length(MOL4_CB), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(OLOPC_19_subs, features = MOL5_CB, ncol = length(MOL5_CB), group.by = "seurat_clusters", pt.size = 0)
VlnPlot(OLOPC_19_subs, features = MOL6_CB, ncol = length(MOL6_CB), group.by = "seurat_clusters", pt.size = 0)

johkietal_genes <- read.csv("johkietal_genes.csv", header = TRUE)
johkietal_genes_list <- johkietal_genes$Gene_list
valid_genes_johki <- johkietal_genes_list[johkietal_genes_list %in% Features(OLOPC_19_subs, assay = "RNA")]
OLOPC_19_subs <- ScaleData(OLOPC_19_subs, features = valid_genes_johki)
DoHeatmap(OLOPC_19_subs, features = valid_genes_johki, group.by = "ident")


OL_subcluster_ids <- c(
  "0" = "MOL1/2",
  "1" = "MOL5/6",
  "2" = "OPC",
  "3" = "MOL3/4",
  "4" = "MOL_DA",
  "5" = "COP",
  "6" = "COP",
  "7" = "NFOL")

names(OL_subcluster_ids) <- levels(OLOPC_19_subs)
OLOPC_19_subs_named <- RenameIdents(object = OLOPC_19_subs, OL_subcluster_ids)


DimPlot(OLOPC_19_subs_named, label = FALSE, label.size = 5, repel = TRUE)

DimPlot(OLOPC_19_subs_named, 
        group.by = "condition", 
        label = FALSE, 
        label.size = 5, 
        repel = TRUE)

#Bar graph for percentages
order_origident <- c("cpzC", "veh1wk", "tambex1wk", "vehA", "vehC","tamA", "bexA", "tambexA", "tambexC", "naiveC")
order_cluster <- c("OPC", "COP", "NFOL", "MOL_DA", "MOL1/2", "MOL3/4", "MOL5/6")  # Replace with your desired order for seurat_clusters


OLOPC_19_subs_named$orig.ident <- factor(OLOPC_19_subs_named$orig.ident, levels = order_origident)
Idents(OLOPC_19_subs_named) <- factor(Idents(OLOPC_19_subs_named), levels = order_cluster)
OLOPC_19_subs_named$cluster <- Idents(OLOPC_19_subs_named)  # Store identities in meta.data


# Step 3: Create a data frame of the count of subclusters per `orig.ident`
subcluster_data <- OLOPC_19_subs_named@meta.data %>%
  group_by(orig.ident, cluster) %>%
  tally() %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()


# Step 4: Create the stacked bar plot
ggplot(subcluster_data, aes(x = orig.ident, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  # Create the bar plot
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "Original Identity (orig.ident)", 
       y = "Percentage of Cells (%)", 
       fill = "Subclusters") + 
  scale_fill_manual(values = c("violetred1", "pink", "orange", "yellow", "lightskyblue", "lightblue", "purple", "black")) +  # Customize the colors of the bars
  ggtitle("Distribution of Subclusters by Original Identity (as Percentage)")  # Title for the plot
unique(OLOPC_19_subs_named@meta.data$condition)  # Check existing values
##Also do this by condition
OLOPC_19_subs_named@meta.data$condition <- as.character(OLOPC_19_subs_named@meta.data$condition)  # Convert to character

OLOPC_19_subs_named@meta.data$condition <- dplyr::case_when(
  OLOPC_19_subs_named@meta.data$condition == "cpz" ~ "CPZ",
  OLOPC_19_subs_named@meta.data$condition == "veh1wk" ~ "Vehicle (1 Week)",
  OLOPC_19_subs_named@meta.data$condition == "tambex1wk" ~ "Tam+Bex (1 Week)",
  OLOPC_19_subs_named@meta.data$condition == "veh" ~ "Vehicle",
  OLOPC_19_subs_named@meta.data$condition == "tam" ~ "Tamoxifen",
  OLOPC_19_subs_named@meta.data$condition == "bex" ~ "Bexarotene",
  OLOPC_19_subs_named@meta.data$condition == "tambex" ~ "Tam+Bex",
  OLOPC_19_subs_named@meta.data$condition == "naive" ~ "Naive",
  TRUE ~ OLOPC_19_subs_named@meta.data$condition  # Keep existing values if not in the list
)

OLOPC_19_subs_named@meta.data$condition <- as.factor(OLOPC_19_subs_named@meta.data$condition)  # Convert back to factor if needed

order_condition <- c("cpz", "veh1wk", "tambex1wk", "veh","tam", "bex", "tambex", "naive")
order_condition_update <- c("CPZ", "Vehicle (1 Week)", "Tamoxifen+Bexarotene (1 Week)", "Vehicle","Tamoxifen", "Bexarotene", "Tamoxifen+Bexarotene", "Naive")


OLOPC_19_subs_named$condition <- factor(OLOPC_19_subs_named$condition, levels = order_condition_update)
Idents(OLOPC_19_subs_named) <- factor(Idents(OLOPC_19_subs_named), levels = order_cluster)
OLOPC_19_subs_named$cluster <- Idents(OLOPC_19_subs_named)  # Store identities in meta.data


# Step 3: Create a data frame of the count of subclusters per `orig.ident`
subclustercon_data <- OLOPC_19_subs_named@meta.data %>%
  group_by(condition, cluster) %>%
  tally() %>%
  ungroup() %>%
  group_by(condition) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()



# Step 4: Create the stacked bar plot
ggplot(subclustercon_data, aes(x = condition, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  # Create the bar plot
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),  # Increase x-axis text size and rotate
    axis.text.y = element_text(size = 9),  # Increase y-axis text size
    axis.title.x = element_text(size = 0),  # Increase x-axis title size
    axis.title.y = element_text(size = 9),  # Increase y-axis title size
    legend.text = element_text(size = 10),  # Increase legend text size
    legend.title = element_text(size = 0)  # Increase legend title size
  ) +
  labs(y = "Percentage of Cells (%)", 
       fill = "Subclusters") + 
  # Extend the color palette with two more colors
  scale_fill_manual(values = c("#bf5186", "#dfa8c3", "#bdbf51", "#51bf8a", "#5351bf", "#f4a261", "#2d87a8")) +  # Better color palette
  ggtitle("Distribution of Subclusters by Condition (as Percentage)")


# Define the conditions you want to display
selected_conditions <- c("Naive", "CPZ", "Tam+Bex")

# Filter the data
filtered_data <- subclustercon_data %>%
  filter(condition %in% selected_conditions)

# Generate the plot with only selected conditions
g <- ggplot(filtered_data, aes(x = condition, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  
  theme_minimal() +  
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 20),  
    axis.title.x = element_text(size = 0),  
    axis.title.y = element_text(size = 20),  
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 0)  
  ) +
  labs(y = "Percentage of Cells (%)", fill = "Subclusters") + 
  scale_fill_manual(values = c("#e37c72", "#bc9c44", "#6cb356", "#5bb694", "#51b4e4", "#938abd", "#c873a8")) +  
  ggtitle("Distribution of Selected Subclusters by Condition (as Percentage)")


selected_OI <- c("cpzC", "naiveC", "tambexA", "tambexC")

# Filter the data
filtered_data_OI <- subcluster_data %>%
  filter(orig.ident %in% selected_OI)

# Generate the plot with only selected conditions
ggplot(filtered_data_OI, aes(x = orig.ident, y = percentage, fill = cluster)) +
  geom_bar(stat = "identity") +  # Create the bar plot
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "Original Identity (orig.ident)", 
       y = "Percentage of Cells (%)", 
       fill = "Subclusters") + 
  scale_fill_manual(values = c("#e37c72", "#bc9c44", "#6cb356", "#5bb694", "#51b4e4", "#938abd", "#c873a8")) +  # Customize the colors of the bars
  ggtitle("Distribution of Subclusters by Original Identity (as Percentage)")  # Title for the plot
unique(OLOPC_19_subs_named@meta.data$condition)  # Check existing values

ThenewMOLS <- subset(OLOPC_19_subs_named, idents = c("MOL1/2", "MOL_DA", "MOL3/4", "MOL5/6"))
# Add the current active identity to metadata as a new column
ThenewMOLS$active_ident <- Idents(ThenewMOLS)


ThenewMOLS_fm <- FindAllMarkers(
  object = ThenewMOLS,  # Your Seurat object
  only.pos = TRUE,  # Include both upregulated and downregulated genes
  min.pct = 0.1,  # Only consider genes expressed in at least 10% of cells
  logfc.threshold = 0.25  # Minimum log fold change threshold
)
write.csv(ThenewMOLS_fm, "ThenewMOLs_fm.csv", row.names = FALSE)






MOL_subtype <- read.csv("MOL_subtype.csv", header = TRUE)
MOL_subtype_list <- MOL_subtype$Gene_list
valid_genes_MOLsubtype <- MOL_subtype_list[MOL_subtype_list %in% Features(ThenewMOLS, assay = "RNA")]
ThenewMOLS <- ScaleData(ThenewMOLS, features = valid_genes_MOLsubtype)
DoHeatmap(ThenewMOLS, features = valid_genes_MOLsubtype, group.by = "ident")

vln_plot <- VlnPlot(
  ThenewMOLS,
  features = valid_genes_MOLsubtype,
  group.by = "active_ident",
  pt.size = 0,  # Hide points for cleaner look
  stack = TRUE,
  flip = TRUE  # Puts violins horizontal (optional)
)
vln_plot

# Set the order of the factor levels for active identities
# Reorder the active_ident column in metadata
ThenewMOLS$active_ident <- factor(ThenewMOLS$active_ident, levels = c("MOL_DA", "MOL1/2", "MOL3/4", "MOL5/6"))





expression_data <- GetAssayData(ThenewMOLS, assay = "RNA", slot = "scale.data")

z_scores <- t(scale(t(expression_data)))
# Add Z-scores for each gene as metadata (optional)
ThenewMOLS$z_scores <- z_scores
# Example: Plot Z-scores for a specific gene
VlnPlot(ThenewMOLS, features = c(valid_genes_MOLsubtype), group.by = "active_ident")

z_scores_by_ident <- t(scale(t(expression_data)))  # Z-scores by row

# Subset the Seurat object to include only the 'active_ident' in the metadata
z_scores_by_ident_df <- as.data.frame(z_scores_by_ident)
z_scores_by_ident_df$active_ident <- ThenewMOLS$active_ident
agg_z_scores <- aggregate(. ~ active_ident, data = z_scores_by_ident_df, FUN = mean)

ThenewMOLS$active_ident <- factor(ThenewMOLS$active_ident, levels = c("MOL_DA", "MOL1/2", "MOL3/4", "MOL5/6"))

agg_z_scores_matrix <- agg_z_scores[, -1]  # Remove 'active_ident' column

# Generate the heatmap of Z-scores by 'active_ident'
pheatmap(
  agg_z_scores_matrix,                    # Aggregated Z-scores matrix
  cluster_rows = TRUE,                    # Cluster genes (rows)
  cluster_cols = TRUE,                    # Cluster identities (columns)
  scale = "none",                         # No scaling needed, Z-scores already calculated
  show_rownames = TRUE,                   # Show gene names
  show_colnames = TRUE,                   # Show active_ident group names
  fontsize_row = 8,                       # Adjust font size for gene names
  fontsize_col = 10,                      # Adjust font size for identity names
  main = "Heatmap of Gene Expression Z-scores by Active Identity",
  color = colorRampPalette(c("blue", "white", "red"))(100)  # Color scale (blue to red)
)


# Calculate the mean Z-scores for each gene across each active.ident
agg_z_scores <- aggregate(. ~ active_ident, data = z_scores_by_ident_df, FUN = mean)



pheatmap(
  z_scores,                         # Z-scores matrix
  cluster_rows = TRUE,               # Cluster genes (rows)
  cluster_cols = TRUE,               # Cluster cells (columns)
  scale = "none",                    # Don't scale again, since we already have Z-scores
  show_rownames = TRUE,              # Show gene names
  show_colnames = FALSE,             # Hide cell names (optional, can display if you like)
  fontsize_row = 8,                  # Adjust font size for gene names
  fontsize_col = 6,                  # Adjust font size for cell names
  main = "Heatmap of Gene Expression Z-scores",
  color = colorRampPalette(c("blue", "white", "red"))(100)  # Custom color scale (blue to red)
)


# Aggregate gene expression by "ident" (e.g., clusters)
agg_data <- AggregateExpression(ThenewMOLS, group.by = "orig.ident", return.seurat = FALSE)

# Extract the expression matrix
agg_matrix <- agg_data$RNA  # Assuming you are using the "RNA" assay

pheatmap(
  agg_matrix,
  cluster_rows = TRUE,        # Cluster genes (rows)
  cluster_cols = TRUE,        # Cluster identities (columns)
  scale = "row",              # Scale by rows (genes)
  show_rownames = TRUE,       # Show gene names
  show_colnames = TRUE,       # Show cluster/condition names
  fontsize_row = 8,           # Adjust font size for gene names
  fontsize_col = 10,          # Adjust font size for identity names
  main = "Aggregated Gene Expression by Identity",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)  # Customize color scale
)
