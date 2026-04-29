## https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/DE_analysis.md
## https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq-count-input

rm(list = ls())
library(DESeq2)
library(ggplot2)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(reshape2)
library(ggrepel)
library(enrichR)
library(stringr)

setEnrichrSite('Enrichr')

dbs = c('MSigDB_Hallmark_2020', 'GO_Biological_Process_2023', 'Reactome_2022', 'BioPlanet_2019', 'WikiPathway_2023_Human', 'Panther_2016', 'KEGG_2021_Human')

colours = colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)

get_head_tail <- function(x, n=10){
	return(rbind(head(x,n = n), tail(x, n = n)))
}

customCol_v2 = RColorBrewer::brewer.pal(name = 'Dark2', n = 8)
#customCol_v2 = c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D')

#~ countTable = read.csv(file = '06_quantification/raw_counts.csv', header = T, stringsAsFactors=F, check.names=F)
#~ rownames(countTable) = countTable$ensembl_id
#~ countTable$ensembl_id = NULL

countTable = read.csv(file = '06_quantification/featureCounts.deseq.input.txt', sep = '\t', header = T, stringsAsFactors=F, check.names=F)
rownames(countTable) = countTable$Geneid
countTable$Geneid = NULL
metaData = read.table('sample_metadata.tsv', header = T, sep = '\t', stringsAsFactors=F)
rownames(metaData) = metaData$sample
metaData$group = as.factor(metaData$group)

## Removing outliers
countTable = countTable[,metaData$sample_id]
countTable = countTable[rowSums(countTable) > 0,]

customCol_v2 = customCol_v2[c(1:length(unique(metaData$group)))]
names(customCol_v2) = unique(metaData$group)
customCol_v2 = customCol_v2[!is.na(names(customCol_v2))]

######
#~ gtf.file = '/opt/db/mm10_index/mm10.gtf'
#~ gtf.gr = rtracklayer::import(gtf.file) # creates a GRanges object
#~ gtf.df = as.data.frame(gtf.gr)
#~ gtf.df = gtf.df[!(gtf.df$gene_biotype %in% c('IG_C_pseudogene', 'IG_D_pseudogene', 'IG_pseudogene', 'IG_V_pseudogene', 'processed_pseudogene', 'polymorphic_pseudogene', 'pseudogene', 'TR_J_pseudogene', 'TR_V_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'transcribed_unitary_pseudogene', 'unitary_pseudogene', 'translated_unprocessed_pseudogene', 'unprocessed_pseudogene')),]
#~ genes = unique(gtf.df[ ,c('gene_id','gene_name', 'width')])
#~ write.csv(genes, file='gene_id_gene_name.csv')
genes = read.csv('gene_id_gene_name.csv', row.names = 1)
genes = genes[!duplicated(genes$gene_id),]

countTable = countTable[rownames(countTable) %in% genes$gene_id,]

dds <- DESeqDataSetFromMatrix(countData = countTable, colData = metaData, design = ~group)

## One can also omit this step entirely and just rely on the independent filtering procedures available in results()
#dds <- dds[ rowSums(counts(dds)) > 5, ]

prefix = '07_DGE_wo_pseudogenes/'
dir.create(prefix, showWarnings = F, recursive = T)

# vsd <- vst(dds, blind=FALSE)
vsd <- rlog(dds, blind=FALSE)
normcounts <- as.data.frame(assay(vsd))
normcounts$gene_id = rownames(normcounts)
normcounts = merge(x = normcounts, y = genes, by = 'gene_id', all.x = TRUE)
rownames(normcounts) = normcounts$gene_id
normcounts$width=NULL; normcounts$gene_id=NULL
write.table(x = normcounts, file = paste0(prefix, '/normalized_count.tsv'), sep = '\t', quote = F)
normcounts$gene_name = NULL
normcounts$gene_id = rownames(normcounts)

pcaData <- plotPCA(vsd, intgroup = 'group', returnData = TRUE)
percentVar <- round(100 * attr(pcaData, 'percentVar'))

p1 <- ggplot(pcaData, aes(PC1, PC2, color=group, label=name)) +
	scale_color_manual(values = customCol_v2) +
	geom_point(size=3) + geom_text(hjust=0.5, vjust=-0.5, size = 2) + 
	xlab(paste0('PC1: ', percentVar[1], '% variance')) +
	ylab(paste0('PC2: ', percentVar[2], '% variance')) +
	theme_bw() + theme(panel.grid=element_blank())
ggsave(paste0(prefix, '/pca_all_samples.png'), p1, height = 5, width = 6)
ggsave(paste0(prefix, '/pca_all_samples.pdf'), p1, height = 5, width = 6)

## distance matrix
sampleDists <- dist(t(assay(vsd)))
metadata_subset = metaData[,c('sample_id', 'group')]
rownames(metadata_subset) = metadata_subset$sample_id
metadata_subset$sample_id = NULL

sampleDistMatrix <- as.matrix(sampleDists)

ha = HeatmapAnnotation(condition = metadata_subset$group,
	col = list(condition = customCol_v2), 
	simple_anno_size = unit(0.6, 'cm'),
	gap = unit(1, 'mm'),
	border = TRUE)
ra = rowAnnotation(condition = metadata_subset$group,
	col = list(condition = customCol_v2), 
	simple_anno_size = unit(0.6, 'cm'),
	gap = unit(1, 'mm'),
	border = TRUE,
	show_legend = FALSE)
	
heatmap_col = circlize::colorRamp2(c(0, 100), c('white', 'black'))

p1 <- Heatmap(sampleDistMatrix, rect_gp = gpar(col = "white", lwd = 1), col = heatmap_col, name = 'dist', border = TRUE, cluster_columns = T, cluster_rows=T, row_title = NULL, show_column_names = T, show_row_names = T, row_dend_width = unit(4, 'cm'), column_dend_height = unit(4, 'cm'), top_annotation = ha, left_annotation=ra)

pdf(paste0(prefix, '/Distance_heatmap.pdf'), width = 9, height = 8)
print(p1)
dev.off()

png(paste0(prefix, '/Distance_heatmap.png'), width = 900, height = 800)
print(p1)
dev.off()

### Comparisions
compList = list('EO771_vs_PY230' = c('EO771', 'PY230'), 'EO771_vs_PY8119' = c('EO771', 'PY8119'), 'PY230_vs_PY8119' = c('PY230', 'PY8119'))

log2fc_cutoff = 1
n_top = 10
p_val_cutoff = 0.01

for(comp in names(compList)) {
	output_dir = paste0(prefix, '/', comp, '/')
	dir.create(output_dir, showWarnings = F, recursive = T)
	
	temp_metaData = metaData[metaData$group %in% compList[[comp]],]
	rownames(temp_metaData) = temp_metaData$sample_id
	temp_countTable = countTable[,temp_metaData$sample_id]
	#all(rownames(temp_metaData) == colnames(temp_countTable))
	temp_metaData$condition = factor(temp_metaData$group, levels = compList[[comp]])
	
	dds <- DESeqDataSetFromMatrix(countData = temp_countTable, colData = temp_metaData, design = ~condition)
	##dds <- dds[ rowSums(counts(dds)) > 5, ]
	
	## PCA Plot
	vsd <- rlog(dds, blind=FALSE)
	pcaData <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
	percentVar <- round(100 * attr(pcaData, 'percentVar'))
	
	p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
		geom_point(size=3) + geom_text(hjust=0, vjust=0) + 
		scale_color_manual(values = customCol_v2) +
		xlab(paste0('PC1: ', percentVar[1], '% variance')) +
		ylab(paste0('PC2: ', percentVar[2], '% variance')) +
		theme_bw() + theme(panel.grid=element_blank())
	ggsave(paste0(output_dir, '/pca.png'), p1, height = 5, width = 6)
	ggsave(paste0(output_dir, '/pca.pdf'), p1, height = 5, width = 6)
	
	## DESeq and Export data
	dds <- DESeq(dds, quiet=T)
	
	res <- results(dds)
	dge <- as.data.frame(res)
	dge$gene_id = rownames(res)
	dge = merge(x = dge, y = genes, by = 'gene_id', all.x = TRUE)
	dge$width=NULL
	
	df_export = merge(x = dge, y = normcounts[,c(temp_metaData$sample_id, 'gene_id')], by = 'gene_id', all.x = TRUE)
	write.table(x = df_export, file = paste0(output_dir, '/deseq2_res_w_normalized_counts.tsv'), sep = '\t', quote = F, row.names = F)
	
	## Volcano Plot
	volcanoData = data.frame(log2FoldChange = dge$log2FoldChange, pValue = dge$pvalue, pAdj = dge$padj, GeneNames = dge$gene_name)
	volcanoData = volcanoData[!is.na(volcanoData$pAdj),]
	volcanoData$Regulation = ifelse((volcanoData$pAdj < p_val_cutoff & volcanoData$log2FoldChange > 1), 'Up Regulated', ifelse((volcanoData$pAdj < p_val_cutoff & volcanoData$log2FoldChange < -1), 'Down Regulated', 'Other'))
	volcanoData$log2FoldChange <- as.numeric(volcanoData$log2FoldChange) 
	volcanoData$pValue <- as.numeric(volcanoData$pValue) 
	volcanoData$pAdj <- as.numeric(volcanoData$pAdj) 
	volcanoData <- volcanoData[is.finite(volcanoData$log2FoldChange),]
	volcanoData$log10Pval = -log(volcanoData$pAdj,10)
	volcanoData$log10Pval[volcanoData$log10Pval > 15] = 15
	df_genenames = volcanoData[order(-volcanoData$log2FoldChange), ]
	df_genenames = df_genenames[(df_genenames$pAdj < p_val_cutoff) & (abs(df_genenames$log2FoldChange)>1),] %>% do(get_head_tail(.,n=10))
	
	p1 <- ggplot(volcanoData, aes(log2FoldChange, log10Pval)) + 
		geom_point(aes(color = Regulation), size = 1.0) + 
		geom_text_repel(data = df_genenames, mapping = aes(log2FoldChange, log10Pval, label = GeneNames), size = 4, nudge_y = .2, min.segment.length = unit(0, 'lines')) + 
		xlab(expression('log'[2]*'FC')) +
		ylab(expression('-log'[10]*'PAdj')) +
		scale_color_manual(values = c('blue', 'gray50', 'red')) + 
		guides(colour = guide_legend(override.aes = list(size=1.5))) + 
		theme_classic() + 
		theme(panel.border = element_rect(colour = 'black', fill=NA, size=0.5), axis.text=element_text(colour = 'black',size=12), axis.title=element_text(size=14, colour = 'black', face='bold'), plot.title = element_text(hjust = 0.5)) +
		ggtitle('')
	ggsave(paste0(output_dir, '/volcano_plot.png'), p1, height = 5, width = 6)
	ggsave(paste0(output_dir, '/volcano_plot.pdf'), p1, height = 5, width = 6)
	
	## Heatmap
	temp_df = na.omit(df_export[(df_export$padj < p_val_cutoff) & abs(df_export$log2FoldChange) > log2fc_cutoff,])
	temp_df = temp_df[order(-temp_df$log2FoldChange), ]
	heatmap_input = temp_df[,c(9:ncol(temp_df))]
	heatmap_input = heatmap_input[,temp_metaData$sample_id]
	#heatmap_input = log2(heatmap_input + 1)
	
	temp_gene_names <- data.frame('gene_names' = temp_df$gene_name)
	temp_gene_names = temp_gene_names %>% group_by(gene_names) %>% mutate(gene_names = paste(gene_names, seq_along(gene_names), sep = '#'))
	temp_gene_names = gsub(pattern = '#1', replacement = '', x = temp_gene_names$gene_names)
	temp_gene_names = gsub(pattern = '#(\\d+)', replacement = '\\.\\1', x = temp_gene_names, perl = T)
	
	## paste0(temp_df$gene_name, seq_along(temp_df$gene_name))
	rownames(heatmap_input) = temp_gene_names
	heatmap_input = t(scale(t(heatmap_input)))		
	
	rowAnno = temp_df %>% do(get_head_tail(.,n=20))
	ha_row = rowAnnotation(gene_names = anno_mark(at = which(rownames(heatmap_input) %in% rowAnno$gene_name), which = 'row', labels = rownames(heatmap_input)[rownames(heatmap_input) %in% rowAnno$gene_name], padding = unit(1, 'mm'), labels_gp = gpar(fontsize = 4), link_width = unit(10, 'mm')), width = unit(4, 'cm'))
	
	color_grads = customCol_v2[names(customCol_v2) %in% unique(temp_metaData$condition)]
	
	p1 <- Heatmap(heatmap_input, name = 'exp', clustering_distance_rows = 'pearson', clustering_method_rows = 'average', column_title = paste0(''), row_title = NULL, row_dend_width = unit(4, 'cm'), show_row_names = FALSE, top_annotation = HeatmapAnnotation('Condition' = temp_metaData$condition, show_annotation_name = FALSE, col = list('Condition' = color_grads), simple_anno_size = unit(0.4, 'cm')), right_annotation = ha_row)
	
	png(paste0(output_dir, '/expression_heatmap.png'), width = 800, height = 800)
	print(p1)
	dev.off()
	pdf(paste0(output_dir, '/expression_heatmap.pdf'), width = 8, height = 8)
	print(p1)
	dev.off()
	
	## Pathway Enrichment
	for (x in c('Up', 'Down')) {
		if(x == 'Up') {
			geneSymbols = df_export[(df_export$pval < p_val_cutoff) & (df_export$log2FoldChange > log2fc_cutoff),]
		} else {
			geneSymbols = df_export[(df_export$pval < p_val_cutoff) & (df_export$log2FoldChange < -log2fc_cutoff),]
		}
		geneSymbols = na.omit(geneSymbols)
		if(nrow(geneSymbols) > 500) {
			geneSymbols = geneSymbols[order(-geneSymbols$log2FoldChange),] %>% do(head(.,n=500))
		}
		geneSymbols = geneSymbols$gene_name
		
		output_dir2 = paste0(output_dir, '/Pathway_for_', x, '_Reg_genes/')
		dir.create(output_dir2, showWarnings = F, recursive = T)
		
		enriched <- enrichr(geneSymbols, dbs)
		
		for (i in names(enriched)) {
			enrich_df = enriched[[i]]
			temp = stringr::str_split_fixed(enrich_df$Overlap, '/', 2)
			enrich_df$NoOfGenes = as.numeric(temp[,1])
			enrich_df$totalGene = as.numeric(temp[,2])
			
			enrich_df$GeneRatio = (enrich_df$NoOfGenes/enrich_df$totalGene) #*100
			enrich_df$P.value2 = -log10(enrich_df$P.value)
			enrich_df = enrich_df[order(enrich_df$P.value),]
			enrich_df$Db = i
			
			write.csv(x = enrich_df, file = paste0(output_dir2, '/', i, '.csv'))
			
			df = enrich_df[order(-enrich_df$P.value2),]
			df = df[c(1:n_top),]
			df = df[df$P.value < 0.05,]
			
			if(nrow(df) > 5) {
				max_term_len = max(unlist(lapply(df$Term, function(x){length(unlist(strsplit(x, split = '')))})))
				custom_w = min(8, round(max_term_len/10) + 1)
				
				p1 <- ggplot(df, aes(x = Db, y = reorder(Term, -P.value), size = NoOfGenes, color = P.value2)) + geom_point(alpha = 0.8) + theme_bw() + scale_color_gradient(low = 'mediumblue',  high = 'red2', limits = c(-log10(0.05), max(df$P.value2))) + labs(size='No. of Genes',col='-log10(pvalue)') + theme(plot.title = element_text(hjust = 0.5)) + xlab('') + ylab('') + theme_bw() + theme(panel.grid = element_blank()) #+ scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 80))
				ggsave(filename = paste0(output_dir2, '/', i, '.pdf'), plot = p1, width = custom_w, height = 6)
				ggsave(filename = paste0(output_dir2, '/', i, '.png'), plot = p1, width = custom_w, height = 6)
			}
		}
	}
}
