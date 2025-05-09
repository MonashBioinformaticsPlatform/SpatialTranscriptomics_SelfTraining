library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(scDblFinder)
library(BiocParallel)
library(SingleR)
library(future)

# For quicly loading datasets
plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)

# Seurat functions mostly taken from here:
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

outdir <- "/oldvol/apattison/Data/Spatial/Comparing_technologies/output/"

# Load the CosMx cell profile. Already looks bad. 'Stroma' is not a cell type
colon_6k_profile <- read_csv("/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx/ColonCRC_6k.profiles.csv")

# Dataset is from this paper
# https://doi.org/10.1186/s13059-025-03554-1

# A 13 mm ×12 mm×5 μm section of colorectal adenocarcinoma was profiled using the 
# standard CosMx RNA protocol and a 6000-plex, pre-commercial version of the CosMx 6 K 
# Discovery RNA panel. Seventy-three fields of view were placed according to markup 
# of a serial hematoxylin and eosin (H&E) stain, focusing on normal intestinal mucosa, 
# lymphoid aggregates, and cancer. The slide was imaged with a 5-channel morphology panel 
# (PanCK, CD68, CD298/B2M, CD45, DAPI).

# Try and read in the CosMx 6k CRC data
six_k <- load("/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx/colon cancer dataset.RData") 

unique(annot$tissue)
unique(annot$slide_ID_numeric)

# Drop NEAT1 as it trends to skew results
raw <- raw[,colnames(raw)!= "NEAT1"]
raw[1:5,1:5]

# Make a 6k seurat object
# Transpose as it is the wrong way around
CosMx.obj_seu = CreateSeuratObject(counts = t(raw), assay="Spatial")

xy[1:5,1:2]
xy_df <- data.frame(xy)

coord.df = data.frame(x=xy_df$x_slide_mm, y=xy_df$y_slide_mm, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = rownames(xy_df)

CosMx.obj_seu@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)

CosMx.obj_seu@assays$Spatial$counts[1:5,1:5]

head(CosMx.obj_seu@meta.data)

# Set up the annotation to join onto the seurat object
annot_bc <- annot%>%
  rownames_to_column("Barcode")

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- CosMx.obj_seu@meta.data%>%
  mutate(Barcode = colnames(CosMx.obj_seu))%>%
  left_join(annot_bc)%>%
  mutate(m = median(Area))%>%
  mutate(s = mad(Area))%>%
  mutate(robzscore_Area = abs((Area - m) / (s)))%>%
  mutate(m = median(nFeature_Spatial))%>%
  mutate(s = mad(nFeature_Spatial))%>%
  mutate(robzscore_nFeature_Spatial = abs((nFeature_Spatial - m) / (s)))%>%
  mutate(m = median(nCount_Spatial))%>%
  mutate(s = mad(nCount_Spatial))%>%
  mutate(robzscore_nCount_Spatial = abs((nCount_Spatial - m) / (s)))%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
CosMx.obj_seu@meta.data <- md

VlnPlot(CosMx.obj_seu, features = c("nCount_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("nFeature_Spatial", "robzscore_nFeature_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("nFeature_Spatial", "robzscore_nFeature_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("Area", "robzscore_Area"), pt.size = -1)

dim(CosMx.obj_seu)

min_QC_robz <- 2

# Subset down based on robust Z score cutoffs for each sample
CosMx.obj_seu <- subset(CosMx.obj_seu, 
                        subset = robzscore_nFeature_Spatial < min_QC_robz & 
                          robzscore_nCount_Spatial < min_QC_robz & 
                          robzscore_Area < min_QC_robz & 
                          nCount_Spatial > 200)

dim(CosMx.obj_seu)

VlnPlot(CosMx.obj_seu, features = c("nCount_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("nFeature_Spatial", "robzscore_nFeature_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("nFeature_Spatial", "robzscore_nFeature_Spatial"), pt.size = -1)
VlnPlot(CosMx.obj_seu, features = c("Area", "robzscore_Area"), pt.size = -1)

# Run without cluster info
# scDblFinder is for scRNA-Seq but might remove the cells with 
# the worst mixing
dbs <- scDblFinder(CosMx.obj_seu@assays$Spatial$counts,
                   BPPARAM=MulticoreParam(20))

# Look at the doublet rate
table(dbs$scDblFinder.class)

CosMx.obj_seu$scDblFinder.class <- dbs$scDblFinder.class

# Keep score for inspection
CosMx.obj_seu$scDblFinder.score <- dbs$scDblFinder.score

# Drop doublets
#CosMx.obj_seu <- CosMx.obj_seu[,CosMx.obj_seu$scDblFinder.class == "singlet"]

# Visualize score distribution and potential cutoff
#VlnPlot(CosMx.obj_seu, features = "scDblFinder.score", pt.size = 0.1) + NoLegend()

# Log normalsie
CosMx.obj_seu <- NormalizeData(CosMx.obj_seu)
CosMx.obj_seu <- FindVariableFeatures(CosMx.obj_seu, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(CosMx.obj_seu), 20) 

# Scale data and run PCA
CosMx.obj_seu <- ScaleData(CosMx.obj_seu)
CosMx.obj_seu <- RunPCA(CosMx.obj_seu, features = VariableFeatures(object = CosMx.obj_seu))

print(CosMx.obj_seu[["pca"]], dims = 1:5, nfeatures = 5)

# Find neighbour cells (in PCA space, not real space)
CosMx.obj_seu <- FindNeighbors(CosMx.obj_seu, dims = 1:30)

# Run with default cluster params
CosMx.obj_seu <- FindClusters(CosMx.obj_seu)

# Run UMAP
CosMx.obj_seu <- RunUMAP(CosMx.obj_seu, dims = 1:30)

# Load up the MMR atlas as a singleR reference
# Rerence dataset is from this paper: https://doi.org/10.1016/j.cell.2021.08.003
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

# Add on the annotations for the MMR atlas
predictions <- SingleR(test=CosMx.obj_seu@assays$Spatial$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       num.threads = 30, aggr.ref = T)

# https://docs.omicsbox.biobam.com/latest/Single-Cell-RNA-Seq-Cell-Type-Prediction-SingleR/

# delta.next, a numeric vector containing the difference between the best and next-best score.

CosMx.obj_seu$SingleR_pred <- predictions$labels 
CosMx.obj_seu$SingleR_delta.next <- predictions$delta.next 
CosMx.obj_seu$SingleR_pred_pruned <- predictions$pruned.labels 

DimPlot(CosMx.obj_seu, label = T)+ NoLegend()

ggplot(data = CosMx.obj_seu@meta.data, aes(x = scDblFinder.class, y = nCount_RNA))+
  geom_boxplot()

ggplot(data = CosMx.obj_seu@meta.data, aes(x = scDblFinder.class, y = SingleR_delta.next))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(0, 0.3))

ggplot(data = CosMx.obj_seu@meta.data, aes(x = SingleR_pred_pruned, y = SingleR_delta.next))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(0, 0.3))

ggplot(data = CosMx.obj_seu@meta.data, aes(x = SingleR_pred_pruned, y = nCount_RNA))+
  geom_boxplot(outlier.shape = NA)

ggplot(data = CosMx.obj_seu@meta.data, aes(x = SingleR_pred_pruned, y = Area))+
  geom_boxplot(outlier.shape = NA)

colnames(CosMx.obj_seu@meta.data)

FeaturePlot(CosMx.obj_seu, features = c("nCount_RNA"))
FeaturePlot(CosMx.obj_seu, features = c("scDblFinder.score"))
FeaturePlot(CosMx.obj_seu, features = c("delta.next"))
FeaturePlot(CosMx.obj_seu, features = c("Area"), order = T)
DimPlot(CosMx.obj_seu, group.by = "scDblFinder.class", label = T)
DimPlot(CosMx.obj_seu, group.by = "SingleR_pred_pruned", label = T, cols = cell_type_colors)
FeaturePlot(CosMx.obj_seu, features = "Mean.PanCK", label = T)
FeaturePlot(CosMx.obj_seu, features = "Mean.PanCK", label = T)
FeaturePlot(CosMx.obj_seu, features = "Mean.PanCK", label = T)
FeaturePlot(CosMx.obj_seu, features = "Mean.DAPI", label = T)

pred_cluster <- DimPlot(CosMx.obj_seu, group.by = "SingleR_pred", label = T, cols = cell_type_colors)
ggsave(plot = pred_cluster, filename = paste0(outdir, "/plots/CosMx_ct_predictions_vs_clusters.pdf"), width = 10, height = 8)

SpatialDimPlot(CosMx.obj_seu,cols = cell_type_colors,
                    group.by = "SingleR_pred")

SpatialDimPlot(CosMx.obj_seu,
               group.by = "seurat_clusters")
  
fovs <- DimPlot(CosMx.obj_seu, group.by = "fov", label = T)+ NoLegend()

ggsave(plot = fovs, filename = paste0(outdir, "/plots/CosMx_FOVs.pdf"), width = 10, height = 8)

FeaturePlot(CosMx.obj_seu, features = "EPCAM")
FeaturePlot(CosMx.obj_seu, features = "FN1")
FeaturePlot(CosMx.obj_seu, features = "LGR5")
FeaturePlot(CosMx.obj_seu, features = "PIGR")

marks <- FindAllMarkers(CosMx.obj_seu)



# Save the object
#qsave(CosMx.obj_seu, paste0(outdir_sample, slide, "_", sample,"_seurat_processed.qs"))


# 10X Xenium dataset ----
# Xenium vs Xenium prime
# https://x.com/NeBanovich/status/1824207483902234995

# *TAKEAWAY*
# 1) Prime comes w/ major sacrifices in sensitivity per gene
# 2) The sensitivity drop comes w/ a loss of dynamic range for genes w/in a cell
# 3) The specs are unlikely to move me off Xenium V1 for most of our work
# 4) Time to think hard about 5K vs HD for borader profiling

# I think we also expect to see this with CosMx 6k vs 1k

path <- "/oldvol/apattison/Data/Spatial/Comparing_technologies/xenium/"
# Load the Xenium data
# tar -xvf cell_feature_matrix.tar.gz 
xenium.obj <- ReadXenium(data.dir = path, type =c('centroids', 'segmentation'))

str(xenium.obj)

# Make a seurat object
xenium.obj_seu = CreateSeuratObject(counts = xenium.obj$matrix$`Gene Expression`, assay="Spatial")

xenium.obj$centroids[1:5,1:3]

coord.df = data.frame(x=xenium.obj$centroids$x, y=xenium.obj$centroids$y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = xenium.obj$centroids$cell

# Barcodes match up
min(as.numeric(rownames(coord.df) ))

xenium.obj_seu@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)

xenium.obj_seu@assays$Spatial$counts[1:5,1:5]

head(xenium.obj_seu@meta.data)

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- xenium.obj_seu@meta.data%>%
  mutate(Barcode = colnames(xenium.obj_seu))%>%
  mutate(m = median(nFeature_Spatial))%>%
  mutate(s = mad(nFeature_Spatial))%>%
  mutate(robzscore_nFeature_Spatial = abs((nFeature_Spatial - m) / (s)))%>%
  mutate(m = median(nCount_Spatial))%>%
  mutate(s = mad(nCount_Spatial))%>%
  mutate(robzscore_nCount_Spatial = abs((nCount_Spatial - m) / (s)))%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
xenium.obj_seu@meta.data <- md

min_QC_robz <- 3

VlnPlot(xenium.obj_seu, features = c("nCount_Spatial", "robzscore_nFeature_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)

dim(xenium.obj_seu)

# Subset down based on robust Z score cutoffs for each sample
xenium.obj_seu <- subset(xenium.obj_seu, subset = robzscore_nFeature_Spatial < min_QC_robz & robzscore_nCount_Spatial < min_QC_robz & nCount_Spatial > 50)

VlnPlot(xenium.obj_seu, features = c("nCount_Spatial", "robzscore_nFeature_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)

dim(xenium.obj_seu)

# Log normalise
xenium.obj_seu <- NormalizeData(xenium.obj_seu)
xenium.obj_seu <- FindVariableFeatures(xenium.obj_seu, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(xenium.obj_seu), 20) 

t20

# Scale data and run PCA
xenium.obj_seu <- ScaleData(xenium.obj_seu)
xenium.obj_seu <- RunPCA(xenium.obj_seu, features = VariableFeatures(object = xenium.obj_seu))

print(xenium.obj_seu[["pca"]], dims = 1:5, nfeatures = 5)

# Find neighbour cells (in PCA space, not real space)
xenium.obj_seu <- FindNeighbors(xenium.obj_seu, dims = 1:30)

# Run with default cluster params
xenium.obj_seu <- FindClusters(xenium.obj_seu)

# Run UMAP
xenium.obj_seu <- RunUMAP(xenium.obj_seu, dims = 1:30)

# Add on the annotations for the MMR atlas
predictions <- SingleR(test=xenium.obj_seu@assays$Spatial$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       aggr.ref = T, num.threads = 10)

xenium.obj_seu$SingleR_pred <- predictions$labels 

DimPlot(xenium.obj_seu, label = T)+ NoLegend()
DimPlot(xenium.obj_seu, group.by = "SingleR_pred", label = T, cols = cell_type_colors)

# NEAT1 is still present
FeaturePlot(xenium.obj_seu, features = c("nCount_Spatial"))

# 
clust2 <- FindMarkers(xenium.obj_seu, ident.1 = "2", ident.2 = "1")
# clust18 <- FindMarkers(xenium.obj_seu, ident.1 = "18")

pred_cluster <- DimPlot(xenium.obj_seu, group.by = "SingleR_pred", label = T, cols = cell_type_colors)

ggsave(plot = pred_cluster, filename = paste0(outdir, "/plots/Xenium_predictions_vs_clusters.pdf"), width = 10, height = 8)

fovs <- DimPlot(xenium.obj_seu, group.by = "fov", label = T)+ NoLegend()

ggsave(plot = fovs, filename = paste0(outdir, "/plots/Xenium_FOVs.pdf"), width = 10, height = 8)

FeaturePlot(xenium.obj_seu, features = "EPCAM")
FeaturePlot(xenium.obj_seu, features = "FN1")
FeaturePlot(xenium.obj_seu, features = "CD8A")
FeaturePlot(xenium.obj_seu, features = "CD4")
FeaturePlot(xenium.obj_seu, features = "ANXA1")
FeaturePlot(xenium.obj_seu, features = "AREG")

SpatialFeaturePlot(xenium.obj_seu, features = "nCount_Spatial", pt.size.factor = 1.2) +
  ggtitle("Spatial Distribution of Total Counts")

marks <- FindAllMarkers(xenium.obj_seu)

# Vizgen merscope
# https://info.vizgen.com/ffpe-showcase
# I took colon cancer 2 from here as it had higher qc scores:

# I had to download off google cloud. It was very annoying. 
# Read in the object using centroids only
vizgen.obj <- ReadVizgen(data.dir = "/oldvol/apattison/Data/Spatial/Comparing_technologies/visium/files/", 
                         type= "centroids", metadata = c("volume", "fov"))

str(vizgen.obj)

# Make a seurat object
vizgen.obj_seu = CreateSeuratObject(counts = vizgen.obj$transcripts, assay="Spatial")

vizgen.obj$centroids[1:5,1:3]

str(vizgen.obj_seu)

coord.df = data.frame(x=vizgen.obj$centroids$x, y=vizgen.obj$centroids$y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) = vizgen.obj$centroids$cell

# Barcodes match up
min(as.numeric(rownames(coord.df) ))

vizgen.obj_seu@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)

SpatialFeaturePlot(vizgen.obj_seu, features = "nCount_Spatial", pt.size.factor = 1.2) +
  ggtitle("Spatial Distribution of Total Counts")

vizgen.obj_seu@assays$Spatial$counts[1:5,1:5]

colnames(vizgen.obj_seu)[1:5]

head(vizgen.obj_seu@meta.data)

head(vizgen.obj$metadata)

md_df <- vizgen.obj$metadata%>%
  rownames_to_column("Barcode")

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- vizgen.obj_seu@meta.data%>%
  mutate(Barcode = colnames(vizgen.obj_seu))%>%
  # Join on the volume and FOV metadata
  left_join(md_df)%>%
  mutate(m = median(nFeature_Spatial))%>%
  mutate(s = mad(nFeature_Spatial))%>%
  mutate(robzscore_nFeature_Spatial = abs((nFeature_Spatial - m) / (s)))%>%
  mutate(m = median(nCount_Spatial))%>%
  mutate(s = mad(nCount_Spatial))%>%
  mutate(robzscore_nCount_Spatial = abs((nCount_Spatial - m) / (s)))%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
vizgen.obj_seu@meta.data <- md

min_QC_robz <- 3

VlnPlot(vizgen.obj_seu, features = c("nCount_Spatial", "robzscore_nFeature_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)

dim(vizgen.obj_seu)

# Subset down based on robust Z score cutoffs for each sample
vizgen.obj_seu <- subset(vizgen.obj_seu, subset = robzscore_nFeature_Spatial < min_QC_robz & robzscore_nCount_Spatial < min_QC_robz & nCount_Spatial > 50)

VlnPlot(vizgen.obj_seu, features = c("nCount_Spatial", "robzscore_nFeature_Spatial", "robzscore_nCount_Spatial"), pt.size = -1)

# Look at the first 1000 FOVs for now
vizgen.obj_seu_fov1 <- vizgen.obj_seu[,vizgen.obj_seu$fov %in% unique(vizgen.obj_seu$fov)[1:1000]]

dim(vizgen.obj_seu_fov1)

# Log normalise
vizgen.obj_seu_fov1 <- NormalizeData(vizgen.obj_seu_fov1)
vizgen.obj_seu_fov1 <- FindVariableFeatures(vizgen.obj_seu_fov1, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(vizgen.obj_seu_fov1), 20) 

t20

# Scale data and run PCA
vizgen.obj_seu_fov1 <- ScaleData(vizgen.obj_seu_fov1)
vizgen.obj_seu_fov1 <- RunPCA(vizgen.obj_seu_fov1, features = VariableFeatures(object = vizgen.obj_seu_fov1))

print(vizgen.obj_seu_fov1[["pca"]], dims = 1:5, nfeatures = 5)

# Find neighbour cells (in PCA space, not real space)
vizgen.obj_seu_fov1 <- FindNeighbors(vizgen.obj_seu_fov1, dims = 1:30)

# Run with default cluster params
vizgen.obj_seu_fov1 <- FindClusters(vizgen.obj_seu_fov1)

# Run UMAP
vizgen.obj_seu_fov1 <- RunUMAP(vizgen.obj_seu_fov1, dims = 1:30)

# Add on the annotations for the MMR atlas
predictions <- SingleR(test=vizgen.obj_seu_fov1@assays$Spatial$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       aggr.ref = T, num.threads = 10)

vizgen.obj_seu_fov1$SingleR_pred <- predictions$labels 

DimPlot(vizgen.obj_seu_fov1, label = T)+ NoLegend()
DimPlot(vizgen.obj_seu_fov1, group.by = "SingleR_pred", label = T, cols = cell_type_colors)

# NEAT1 is still present
FeaturePlot(vizgen.obj_seu_fov1, features = c("nCount_Spatial"))

# 
# clust3 <- FindMarkers(vizgen.obj_seu_fov1, ident.1 = "3")
# clust18 <- FindMarkers(vizgen.obj_seu_fov1, ident.1 = "18")

pred_cluster <- DimPlot(vizgen.obj_seu_fov1, group.by = "SingleR_pred", label = T, cols = cell_type_colors)

ggsave(plot = pred_cluster, filename = paste0(outdir, "/plots/Merscope_ct_predictions_vs_clusters.pdf"), width = 10, height = 8)

fovs <- DimPlot(vizgen.obj_seu_fov1, group.by = "fov", label = T)+ NoLegend()

ggsave(plot = fovs, filename = paste0(outdir, "/plots/Merscope_FOVs.pdf"), width = 10, height = 8)

FeaturePlot(vizgen.obj_seu_fov1, features = "EPCAM")
FeaturePlot(vizgen.obj_seu_fov1, features = "FN1")
FeaturePlot(vizgen.obj_seu_fov1, features = "CD8A")
FeaturePlot(vizgen.obj_seu_fov1, features = "CD4")

marks <- FindAllMarkers(vizgen.obj_seu_fov1)




