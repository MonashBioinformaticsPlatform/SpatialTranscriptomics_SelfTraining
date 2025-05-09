library(tidyverse)
library(Seurat)
library(future)
library(LoomExperiment)
library(sf)
library(readr) # For read_file, in case it's not a perfectly formatted GeoJSON
library(jsonlite) # For parsing JSON if needed

# For quicly loading datasets
plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

# Load up the MMR atlas as a singleR reference
# Rerence dataset is from this paper: https://doi.org/10.1016/j.cell.2021.08.003
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

# Have a look at the daat structures
test_trans <- read_csv("/oldvol/apattison/Data/Spatial/raw_to_publish/Run5654_Tumor_A/Run5654_Tumor_A_tx_file.csv", n_max = 100)
test_md <- read_csv("/oldvol/apattison/Data/Spatial/raw_to_publish/Run5654_Tumor_A/Run5654_Tumor_A_metadata_file.csv", n_max = 100)

head(test_md)

# Look at the data structure for Baysor
test_trans[1:5,1:9]

genes <- unique(test_trans$target)

genes

# Calculate radius for each cell
test_md_with_radius <- test_md %>%
  mutate(radius_px = sqrt(Area / pi))

# Calculate typical radius values
mean_radius <- mean(test_md_with_radius$radius_px, na.rm = TRUE)
median_radius <- median(test_md_with_radius$radius_px, na.rm = TRUE)

cat("Mean estimated cell radius (pixels):", mean_radius, "\n")
cat("Median estimated cell radius (pixels):", median_radius, "\n")

# Add vertical lines for mean and median
ggplot(test_md_with_radius, aes(x = radius_px)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  geom_vline(aes(xintercept = mean_radius), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median_radius), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Estimated Cell Radii",
       x = "Estimated Radius (pixels)",
       y = "Frequency") +
  annotate("text", x = mean_radius * 1.1, y = 5, label = paste("Mean =", round(mean_radius,1)), color = "red", angle=90, vjust=-0.5) + # Adjust y based on your histogram height
  annotate("text", x = median_radius * 0.9, y = 5, label = paste("Median =", round(median_radius,1)), color = "blue", angle=90, vjust=1.5) + # Adjust y
  theme_minimal()

# Start with FOVs 5/6/7/8

# --- Configuration ---
full_tx_file <- "/oldvol/apattison/Data/Spatial/raw_to_publish/Run5654_Tumor_A/Run5654_Tumor_A_tx_file.csv"
# 2. Path for the new smaller (test) transcript CSV file
test_tx_file <- "/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx_Abud_Baysor/Run5654_Tumor_A_tx_file_fovs_5-8_test.csv"
# 3. FOVs to include in the test dataset
selected_fovs <- c(5, 6, 7, 8)

# Maybe try 14-18 as this is also single sample and has less clear cell types

# --- Script ---
# Read the full transcript data
full_data <- read_csv(full_tx_file, show_col_types = FALSE)

# Filter for the selected FOVs
test_data <- full_data %>%
  filter(fov %in% selected_fovs)

# Write the filtered data to the new test file
write_csv(test_data, test_tx_file)

# Try and cluster the data before Baysor
slide <- LoadNanostring(data.dir = "/oldvol/apattison/Data/Spatial/raw_to_publish/Run5654_Tumor_A/", fov = "fov_all")

slide$Barcode <- colnames(slide)

head(slide$Barcode)

# Read the metadata and add a barcode
meta <- read_csv("/oldvol/apattison/Data/Spatial/raw_to_publish/Run5654_Tumor_A/Run5654_Tumor_A_metadata_file.csv")%>%
  mutate(Barcode = paste0(cell_ID, "_", fov))


slide_md <- slide@meta.data%>%
  left_join(meta)

rownames(slide_md) <- slide_md$Barcode


slide@meta.data <- slide_md

slide <- slide[,slide$fov %in% c(5, 6, 7, 8)]

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- slide@meta.data%>%
  mutate(m = median(Area))%>%
  mutate(s = mad(Area))%>%
  mutate(robzscore_Area = abs((Area - m) / (s)))%>%
  mutate(m = median(nFeature_Nanostring))%>%
  mutate(s = mad(nFeature_Nanostring))%>%
  mutate(robzscore_nFeature_Nanostring = abs((nFeature_Nanostring - m) / (s)))%>%
  mutate(m = median(nCount_Nanostring))%>%
  mutate(s = mad(nCount_Nanostring))%>%
  mutate(robzscore_nCount_Nanostring = abs((nCount_Nanostring - m) / (s)))%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
slide@meta.data <- md

VlnPlot(slide, features = c("nCount_Nanostring", "robzscore_nCount_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("nFeature_Nanostring", "robzscore_nFeature_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("nFeature_Nanostring", "robzscore_nFeature_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("Area", "robzscore_Area"), pt.size = -1)

dim(slide)

min_QC_robz <- 2

# Subset down based on robust Z score cutoffs for each sample
slide <- subset(slide, 
                        subset = robzscore_nFeature_Nanostring < min_QC_robz & 
                          robzscore_nCount_Nanostring < min_QC_robz & 
                          robzscore_Area < min_QC_robz & 
                          nCount_Nanostring > 20)

dim(slide)

VlnPlot(slide, features = c("nCount_Nanostring", "robzscore_nCount_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("nFeature_Nanostring", "robzscore_nFeature_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("nFeature_Nanostring", "robzscore_nFeature_Nanostring"), pt.size = -1)
VlnPlot(slide, features = c("Area", "robzscore_Area"), pt.size = -1)

# Run without cluster info
# scDblFinder is for scRNA-Seq but might remove the cells with 
# the worst mixing
dbs <- scDblFinder(slide@assays$Nanostring$counts,
                   BPPARAM=MulticoreParam(20))

# Look at the doublet rate
table(dbs$scDblFinder.class)

slide$scDblFinder.class <- dbs$scDblFinder.class

# Keep score for inspection
slide$scDblFinder.score <- dbs$scDblFinder.score

# Drop doublets
#slide <- slide[,slide$scDblFinder.class == "singlet"]

# Visualize score distribution and potential cutoff
#VlnPlot(slide, features = "scDblFinder.score", pt.size = 0.1) + NoLegend()

# Log normalsie
slide <- NormalizeData(slide)
slide <- FindVariableFeatures(slide, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(slide), 20) 

# Scale data and run PCA
slide <- ScaleData(slide)
slide <- RunPCA(slide, features = VariableFeatures(object = slide))

print(slide[["pca"]], dims = 1:5, nfeatures = 5)

# Find neighbour cells (in PCA space, not real space)
slide <- FindNeighbors(slide, dims = 1:30)

# Run with default cluster params
slide <- FindClusters(slide)

# Run UMAP
slide <- RunUMAP(slide, dims = 1:30)

# Add on the annotations for the MMR atlas
predictions <- SingleR(test=slide@assays$Nanostring$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       num.threads = 30, aggr.ref = T)

# delta.next, a numeric vector containing the difference between the best and next-best score.
slide$SingleR_pred <- predictions$labels 
slide$SingleR_delta.next <- predictions$delta.next 
slide$SingleR_pred_pruned <- predictions$pruned.labels 

DimPlot(slide, label = T)+ NoLegend()

DimPlot(slide, group.by = "SingleR_pred", label = T, cols = cell_type_colors)
DimPlot(slide, group.by = "scDblFinder.class", label = F)

clust0 <- FindMarkers(slide, ident.1 = "0", ident.2 = "8")

ImageDimPlot(slide, axes = TRUE, cols = cell_type_colors, group.by = "SingleR_pred", fov = "8")

# Try running the same pipeline again with the Baysor results ----
baysor_loom <- import("/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx_Abud_Baysor/baysor_out/segmentation_counts.loom")

baysor_counts <- assays(baysor_loom)$matrix%>%
  as.matrix()

baysor_counts[1:5, 1:5]

dim(baysor_counts)

cell_metadata <- colData(baysor_loom)%>%
  data.frame()%>%
  mutate(Barcode = Name)
head(cell_metadata)

gene_metadata <- rowData(baysor_loom)%>%
  data.frame()

rownames(baysor_counts) <- gene_metadata$Name
colnames(baysor_counts) <- cell_metadata$Name

# Try and cluster the data before Baysor
Baysor_seu <- CreateSeuratObject(baysor_counts, assay = "Baysor")

Baysor_seu$Barcode <- colnames(Baysor_seu)

Baysor_md <- Baysor_seu@meta.data%>%
  left_join(cell_metadata)

rownames(Baysor_md) <- Baysor_md$Barcode

Baysor_seu@meta.data <- Baysor_md

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- Baysor_seu@meta.data%>%
  dplyr::rename(Area= area)%>%
  mutate(Area = replace(Area, is.na(Area), 0))%>%
  mutate(m = median(Area))%>%
  mutate(s = mad(Area))%>%
  mutate(robzscore_Area = abs((Area - m) / (s)))%>%
  mutate(m = median(nFeature_Baysor))%>%
  mutate(s = mad(nFeature_Baysor))%>%
  mutate(robzscore_nFeature_Baysor = abs((nFeature_Baysor - m) / (s)))%>%
  mutate(m = median(nCount_Baysor))%>%
  mutate(s = mad(nCount_Baysor))%>%
  mutate(robzscore_nCount_Baysor = abs((nCount_Baysor - m) / (s)))%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
Baysor_seu@meta.data <- md

VlnPlot(Baysor_seu, features = c("nCount_Baysor", "robzscore_nCount_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("nFeature_Baysor", "robzscore_nFeature_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("nFeature_Baysor", "robzscore_nFeature_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("Area", "robzscore_Area"), pt.size = -1)

dim(Baysor_seu)

min_QC_robz <- 2

# Subset down based on robust Z score cutoffs for each sample
Baysor_seu <- subset(Baysor_seu, 
                subset = robzscore_nFeature_Baysor < min_QC_robz & 
                  robzscore_nCount_Baysor < min_QC_robz & 
                  robzscore_Area < min_QC_robz & 
                  nCount_Baysor > 20)

dim(Baysor_seu)

VlnPlot(Baysor_seu, features = c("nCount_Baysor", "robzscore_nCount_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("nFeature_Baysor", "robzscore_nFeature_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("nFeature_Baysor", "robzscore_nFeature_Baysor"), pt.size = -1)
VlnPlot(Baysor_seu, features = c("Area", "robzscore_Area"), pt.size = -1)

# Run without cluster info
# scDblFinder is for scRNA-Seq but might remove the cells with 
# the worst mixing
dbs <- scDblFinder(Baysor_seu@assays$Baysor$counts,
                   BPPARAM=MulticoreParam(20))

# Look at the doublet rate
table(dbs$scDblFinder.class)

Baysor_seu$scDblFinder.class <- dbs$scDblFinder.class

# Keep score for inspection
Baysor_seu$scDblFinder.score <- dbs$scDblFinder.score

# Drop doublets
#Baysor_seu <- Baysor_seu[,Baysor_seu$scDblFinder.class == "singlet"]

# Visualize score distribution and potential cutoff
#VlnPlot(Baysor_seu, features = "scDblFinder.score", pt.size = 0.1) + NoLegend()

# Log normalsie
Baysor_seu <- NormalizeData(Baysor_seu)
Baysor_seu <- FindVariableFeatures(Baysor_seu, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(Baysor_seu), 20) 

# Scale data and run PCA
Baysor_seu <- ScaleData(Baysor_seu)
Baysor_seu <- RunPCA(Baysor_seu, features = VariableFeatures(object = Baysor_seu))

print(Baysor_seu[["pca"]], dims = 1:5, nfeatures = 5)

# Find neighbour cells (in PCA space, not real space)
Baysor_seu <- FindNeighbors(Baysor_seu, dims = 1:30)

# Run with default cluster params
Baysor_seu <- FindClusters(Baysor_seu)

# Run UMAP
Baysor_seu <- RunUMAP(Baysor_seu, dims = 1:30)

# Add on the annotations for the MMR atlas
predictions <- SingleR(test=Baysor_seu@assays$Baysor$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       num.threads = 30, aggr.ref = T)

# delta.next, a numeric vector containing the difference between the best and next-best score.
Baysor_seu$SingleR_pred <- predictions$labels 
Baysor_seu$SingleR_delta.next <- predictions$delta.next 
Baysor_seu$SingleR_pred_pruned <- predictions$pruned.labels 

DimPlot(Baysor_seu, label = T)+ NoLegend()


DimPlot(Baysor_seu, group.by = "scDblFinder.class", label = F)

clust0 <- FindMarkers(Baysor_seu, ident.1 = "0", ident.2 = "8")

ImageDimPlot(Baysor_seu, axes = TRUE, cols = cell_type_colors, group.by = "SingleR_pred", fov = "8")

# Compare methods
DimPlot(Baysor_seu, group.by = "SingleR_pred", label = T, cols = cell_type_colors)
DimPlot(slide, group.by = "SingleR_pred", label = T, cols = cell_type_colors)

FeaturePlot(slide, features = "CD4")
FeaturePlot(Baysor_seu, features = "CD4")

FeaturePlot(slide, features = "CD3E")
FeaturePlot(Baysor_seu, features = "CD3E")

# Try and plot the baysor coords
# Define the path to your polygons file.
# Baysor might output this as 'segmentation_boundaries.geojson', 'polygons.geojson',
# or 'segmentation_polygons.csv' (if each row in CSV is a JSON string, which is less common for true GeoJSON).
# Let's assume it's a GeoJSON file for now.
polygon_file_path <- "/oldvol/apattison/Data/Spatial/Comparing_technologies/CosMx_Abud_Baysor/baysor_out/segmentation_polygons_2d.json" #
# **ACTION: Please verify the exact name of your polygon file from Baysor's output.**

polygons_sf <- st_read(polygon_file_path)
polygons_sf <- st_set_crs(polygons_sf, NA_crs_) 

print(head(polygons_sf))
print(colnames(polygons_sf))

polygons_sf$id <- as.character(polygons_sf$id)

common_cells <- intersect(colnames(Baysor_seu), polygons_sf$id)
head(common_cells)
dim(Baysor_seu)
# All cells are there as expected
length(common_cells)

Baysor_seu <- Baysor_seu[, common_cells]

# Original polygon processing
polygons_sf_ordered <- polygons_sf[match(common_cells, polygons_sf$id), ]
rownames(polygons_sf_ordered) <- polygons_sf_ordered$id

# Add the 'cell' column required by CreateFOV's 'molecules' argument for sf objects
polygons_sf_ordered$cell <- rownames(polygons_sf_ordered)

centroids_df <- data.frame(x = Baysor_seu$x,
                                                         y = Baysor_seu$y,
                                                         row.names = Baysor_seu$Barcode)

sum(rownames(centroids_df) == colnames(Baysor_seu)) == ncol(Baysor_seu)

current_cell_names <- colnames(Baysor_seu) # Get cell names from the current state of Baysor_seu

# Ensure centroids_df is correctly subsetted and has rownames matching these cell names
# Your current centroids_df creation seems to do this already if Baysor_seu$Barcode is correct
# but let's be hyper-explicit for centroids_df used in CreateFOV:
centroids_for_fov <- centroids_df[current_cell_names, c("x", "y")] # Ensure order and subset, and correct columns

centroids_for_fov$Cells <- current_cell_names

str(centroids_for_fov)

fov_obj <- CreateFOV(
  coords = centroids_for_fov,  
  type = "segmentation",         
  assay = DefaultAssay(Baysor_seu)
)

Baysor_seu@images$fov_all <- fov_obj

ImageDimPlot(Baysor_seu, fov = "fov_all", group.by = "seurat_clusters")
