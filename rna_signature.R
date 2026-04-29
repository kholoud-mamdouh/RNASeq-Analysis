# Load libraries
library(readxl)
library(dplyr)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(ggplot2)


#  Load expression data and MMC5 signatures

expr <- read.csv("log2.cpm.filtered.df.csv")
mmc5 <- read.csv("mmc5.csv")

# Make gene names unique and set rownames
expr$geneID <- make.unique(expr$geneID)
expr <- expr %>% column_to_rownames("geneID")


# Subset samples and define groups

samples <- c("Oli_Ctrl_1","Oli_Ctrl_2","Oli_Ctrl_3","Oli_CC_2",
             "Oli_CC_3","Oli_CC_4")
expr <- expr[, samples]

group <- data.frame(
     Sample = samples,
     Condition = factor(c("Control","Control","Control","CoCulture","CoCulture","CoCulture"),levels = c("Control","CoCulture"))
)
rownames(group) <- group$Sample


#  Prepare gene sets

gene_sets <- lapply(mmc5, function(x) unique(na.omit(x)))
names(gene_sets) <- colnames(mmc5)

# Filter expression to only genes in MMC5
common_genes <- intersect(rownames(expr), unique(unlist(gene_sets)))
expr_filt <- expr[common_genes, ]


# Compute signature scores

signature_scores <- sapply(gene_sets, function(genes){
     genes <- intersect(genes, rownames(expr_filt))
     if(length(genes) > 2){
          colMeans(expr_filt[genes, , drop=FALSE])
     } else {
          rep(NA, ncol(expr_filt))
     }
})
signature_scores <- t(signature_scores)
colnames(signature_scores) <- colnames(expr_filt)



sum(is.na(as.matrix(dist(signature_scores))))


# # Scale signature scores
# signature_scores_scaled <- t(scale(t(signature_scores)))
# 
# #-------------------------------
# # 5. Heatmap of MMC5 signatures
# #-------------------------------
# ann <- data.frame(Condition = group$Condition)
# rownames(ann) <- rownames(group)
# # --- 1. Replace NAs with row mean ---
# signature_scores_scaled <- t(apply(signature_scores_scaled, 1, function(x) {
#      x[is.na(x)] <- mean(x, na.rm = TRUE)
#      return(x)
# }))
# 
# # --- 2. Remove rows with zero variance ---
# row_var <- apply(signature_scores_scaled, 1, var)
# signature_scores_scaled <- signature_scores_scaled[row_var > 0, ]
# 
# # --- 3. Check remaining rows ---
# if(nrow(signature_scores_scaled) < 2) {
#      warning("Too few rows to cluster heatmap. Plotting as single row heatmap or barplot.")
#      
#      # Optional: just plot as a barplot
#      barplot(as.numeric(signature_scores_scaled[1,]),
#              names.arg = colnames(signature_scores_scaled),
#              las = 2,
#              main = "MMC5 Signature Activity (single row)",
#              col = "steelblue")
#      
# } else {
#      # --- 4. Plot heatmap safely ---
#      pheatmap(signature_scores_scaled,
#               annotation_col = ann,
#               clustering_distance_rows = "correlation",
#               clustering_distance_cols = "correlation",
#               fontsize = 9,
#               main = "MMC5 Signature Activity Heatmap")
# }


giveNAs = which(is.na(as.matrix(dist(signature_scores))),arr.ind=TRUE)
head(giveNAs)

signature_scores[c(1,26),]

tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
     sum(is.na(as.matrix(dist(signature_scores[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
signature_scores = signature_scores[-as.numeric(rmv),]
nrow(signature_scores)
ncol(signature_scores)
str(signature_scores)



p <- pheatmap(
     signature_scores)








library(pheatmap)
library(RColorBrewer)

signature_scores_scaled <- t(scale(t(signature_scores)))

# Replace any remaining NAs with row mean
signature_scores_scaled <- t(apply(signature_scores_scaled, 1, function(x) {
     x[is.na(x)] <- mean(x, na.rm = TRUE)
     return(x)
}))

row_var <- apply(signature_scores_scaled, 1, var)
signature_scores_scaled <- signature_scores_scaled[row_var > 0, ]


ann <- data.frame(Condition = group$Condition)
rownames(ann) <- rownames(group)


my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

p <- pheatmap(
     signature_scores_scaled,
     color = my_palette,
     cluster_rows = TRUE,
     cluster_cols = TRUE,
     clustering_distance_rows = "correlation",
     clustering_distance_cols = "correlation",
     scale = "none",               
     annotation_col = ann,         
     show_rownames = TRUE,
     show_colnames = TRUE,
     fontsize = 10,
     fontsize_row = 8,
     fontsize_col = 9,
     angle_col = 45,
     border_color = NA,
     main = "MMC5 Signature Activity Heatmap"
)


