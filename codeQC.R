##conda activate sc_rnaseq_202410
rm(list = ls())
source("helper_v2.R")
source("utilsQC.R")

##---------------- QC
sample_names = c("AVID17", "AVID18", "AVID19", "AVID20")
cellranger_data_dir = "/media/hdd1/06_BLS_RNA_Seq/02_cellranger"
outputDir = "./"

ref_list = lapply(sample_names, function(sample_name) {
	so <- Read10X(data.dir = file.path(cellranger_data_dir, sample_name, "outs/filtered_feature_bc_matrix"))
	so <- CreateSeuratObject(counts = so, project = sample_name)
	
	sce = as.SingleCellExperiment(so)
	sce = scDblFinder(sce, verbose = F)
	
	so[["doublet_status"]] <- as.character(sce$scDblFinder.class)
	so[["doublet_score"]] <- sce$scDblFinder.score
	so[["percent_mt"]] = PercentageFeatureSet(so, pattern = "^mt-")
	so[["percent_ribo"]] = PercentageFeatureSet(so, pattern = "^Rps|^Rpl")
	so[["barcodes"]] = rownames(so@meta.data)
	so[["sample"]] = so$orig.ident
	
	so = NormalizeData(so)
	so = CellCycleScoring(so, s.features = mm10_S_genes, g2m.features = mm10_G2M_genes, set.ident = T)
	so[["cell_cycle_phase"]] = so[["Phase"]]
	
	Idents(so) = so$orig.ident
	so[["Phase"]] <- NULL; so[["old.ident"]] <- NULL; so[["barcodes"]] = NULL
	
	so <- SCTransform(so, verbose = F, variable.features.n = 3000)
	so <- RunPCA(so, features = VariableFeatures(object = so), verbose = F)
	so <- FindNeighbors(so, dims = 1:20, verbose = F)
	so <- RunUMAP(so, dims = 1:20, verbose = F)
	so
})
gc()
names(ref_list) = sample_names

##------------- Plotting QC matrices
outputDir = "01_pre_filtering/"
merged_so <- merge(ref_list[[1]], y = ref_list[c(2:length(ref_list))], add.cell.ids = names(ref_list), project = "TME")
merged_so$sample = factor(merged_so$sample, levels = sample_names)

dir.create(outputDir, showWarnings = F, recursive = T)
customQCPlots(merged_so, prefix = outputDir, customWidth = 7, customHeight = 6)

## Optimal cut-off for nGeneMinCutOff = 800
plotQC(ref_list, prefix = outputDir, percentMtCutOff = 5, nGeneMinCutOff = 400, nGeneMaxCutOff = NULL, nUmiMinCutOff = 500)
saveRDS(ref_list, paste0("ref_list.rds"))

##----------- Filtering Cut-offs
gc()
merged_so_subset = subset(merged_so, subset = (nFeature_RNA > 400 & percent_mt < 15 & doublet_status == "singlet"))

DefaultAssay(merged_so_subset) = "RNA"
merged_so_subset[["SCT"]] = NULL
ref_list_flt = SplitObject(merged_so_subset, split.by = "sample")
saveRDS(ref_list_flt, paste0("ref_list_flt.rds"))

outputDir = "02_post_filtering/"
dir.create(outputDir, showWarnings = F, recursive = T)

#~ ##---------- Statistics Table and POST-FILTERED QC plot
#~ df_cellstats = data.frame(table(merged_so$sample))
#~ colnames(df_cellstats) = c("sample_name", "pre_filtered_counts")
#~ temp_df = data.frame(table(merged_so_subset$sample))
#~ df_cellstats$post_filtered_counts = temp_df$Freq
#~ temp_df = as.data.frame(table(merged_so$sample, merged_so$doublet_status), stringsAsFactors = F)
#~ df_cellstats = cbind(df_cellstats, temp_df[temp_df$Var2 == "doublet",])
#~ df_cellstats$Var1 = NULL; df_cellstats$Var1 = NULL; df_cellstats$Var2 = NULL; df_cellstats$Var2 = NULL
#~ colnames(df_cellstats) = c("sample_name", "pre_filtered_counts", "post_filtered_counts", "number_of_doublets") #, "number_of_singlet"
#~ df_cellstats$number_of_doublets = paste0(df_cellstats$number_of_doublets, " (", round(df_cellstats$number_of_doublets/df_cellstats$pre_filtered_counts*100,2), "%)")
#~ write.table(df_cellstats, file = paste0(outputDir, "sample_wise_cell_stats.tsv"), sep = "\t", quote = F, row.names = F)

#~ customQCPlots(merged_so_subset, prefix = outputDir, customWidth = 7, customHeight = 6)

#~ ref_list_flt = lapply(ref_list_flt, function(so) {
#~ 	so <- SCTransform(so, verbose = F, variable.features.n = 3000)
#~ 	so <- RunPCA(so, features = VariableFeatures(object = so), verbose = F)
#~ 	so <- FindNeighbors(so, dims = 1:20, verbose = F)
#~ 	so <- RunUMAP(so, dims = 1:20, verbose = F)
#~ 	so
#~ })
#~ names(ref_list_flt) = sample_names
#~ plotQC(ref_list_flt, prefix = outputDir, percentMtCutOff = NULL, nGeneMinCutOff = NULL, nGeneMaxCutOff = NULL, nUmiMinCutOff = NULL, status = "POST-FILTERED")

##--------------- Integration
#~ DefaultAssay(merged_so_subset) = "RNA"
#~ merged_so_subset[["SCT"]] = NULL
#~ ref_list_flt = SplitObject(merged_so_subset, split.by = "sample")
#~ saveRDS(ref_list_flt, paste0("ref_list_flt.rds"))
