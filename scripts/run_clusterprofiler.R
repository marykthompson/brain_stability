# Run ClusterProfiler enrichment on a cluster of genes
# https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#ora-algorithm
# The background distribution by default is all the genes that have annotation.

# to run from the commandline:
# Rscript run_clusterprofiler.R -sub ../clust_1.csv -bg ../clust_all.csv -outdir .
library(tidyverse)
library(AnnotationHub)
hub <- AnnotationHub()
library(clusterProfiler)
library(argparse)

parser <- ArgumentParser()

# if you write a double flag, i.e. --flag, then this will be the name in args
parser$add_argument("-sub", help="csv file with subset genes")
parser$add_argument("-bg", help="csv file with background genes")
parser$add_argument("-outdir", help="output directory")
parser$add_argument("-sim_co", type="double", help="similarity cutoff to use for GO trimming")

# To use for non-interactive use:
# setwd('/Users/mk/Desktop/Davislab/3.10_brain_stability/brain_stability/Figures/CTS/genesets/')
# sub <- 'CTS_10stable_genes.csv'
# bg <- 'bg_genes.csv'
# outdir <-'../Figures/CTS/'
# sim_co <-0.5
# args <- parser$parse_args(c("-sub", sub, "-bg", bg, "-outdir", outdir, "-sim_co", sim_co))

args <- parser$parse_args()
sub_genes <- read.table(args$sub, header=FALSE)$V1
bg_genes <- read.table(args$bg, header=FALSE)$V1
Dmel <- hub[["AH100407"]]
pval_co = 0.05
cats <- c("BP", "MF", "CC")

print(paste0('argssub', args$sub))
clust <- unlist(strsplit(basename(args$sub), split=".csv"))
print(paste0('clust name ', clust))

# Convert Flybase IDs to ENTREZ
gene.df_fb <- bitr(bg_genes, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = Dmel)
unmapped_genes_fb <- setdiff(bg_genes, gene.df_fb$FLYBASE)
# Convert Flybase IDs to ENSEMBL
gene.df_ens <- bitr(bg_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = Dmel)
unmapped_genes_ens <- setdiff(bg_genes, gene.df_ens$ENSEMBL)
# Get final conversion, with priority given to the Flybase mapping
df1 <- data.frame(bg_genes)
fb_mapped <- dplyr::left_join(df1, gene.df_fb, by=c("bg_genes"="FLYBASE"))
ens_mapped <- dplyr::left_join(df1, gene.df_ens, by=c("bg_genes"="ENSEMBL"))
final_mapped <- data.frame(fb_mapped)

# this will replace nans in the fb mapping with ones from the ensembl mapping
# https://stackoverflow.com/questions/34697032/fill-in-missing-values-nas-with-values-from-another-dataframe-in-r
final_mapped$ENTREZID[is.na(final_mapped$ENTREZID)] <- ens_mapped$ENTREZID[match(final_mapped$bg_genes, ens_mapped$bg_genes)][which(is.na(final_mapped$ENTREZID))]
# Create converted gene lists
sub_mapped <- final_mapped[final_mapped$bg_genes %in% sub_genes,'ENTREZID']
bg_mapped <- final_mapped$ENTREZID

# Print some stats about name conversion
# Num unmapped
print(paste('final unmapped:', sum(is.na(final_mapped$ENTREZID))))

# check if assignment done correctly
test <- left_join(fb_mapped, final_mapped, by = c("bg_genes"="bg_genes"))
test <- left_join(test, ens_mapped, by = c("bg_genes"="bg_genes"))
test['yes'] = test$ENTREZID.x == test$ENTREZID.y
# view the ones that were nan on fb_mapped
test[is.na(test['yes']),]

# make sure that there are no FALSE ones (mismatch between Fb and final)
not_na <- !is.na(test$yes)
print(paste('no mismatch:', sum(test[not_na,'yes']) == length(test[not_na, 'yes'])))

for (cat in cats) {
  res <- enrichGO(sub_mapped, Dmel, pvalueCutoff = pval_co, universe = bg_mapped, keyType = "ENTREZID", ont = cat)
  #simplify the results
  res_simp <- clusterProfiler::simplify(res, cutoff = args$sim_co, by = "p.adjust", select_fun = min)
  filtered_res <- filter(res_simp@result, p.adjust < pval_co)
  outfile <- file.path(args$outdir, paste0(clust, "_", cat, "_", args$sim_co, ".csv"))
  write.csv(filtered_res, outfile, row.names = FALSE)
}

#Run Kegg pathway enrichment analysis
#convert NCBI IDs to Kegg IDs
kegg_df <- bitr_kegg(final_mapped$ENTREZID, fromType = "ncbi-geneid", toType = "kegg", organism = "dme")
kegg_sub <- kegg_df[kegg_df$`ncbi-geneid` %in% sub_mapped, 'kegg']
kegg_bg <- kegg_df$kegg

res_kegg <- enrichKEGG(kegg_sub, organism = 'dme', pvalueCutoff = 0.05, universe = kegg_bg, keyType = "kegg")
filtered_res <- filter(res_kegg@result, p.adjust < pval_co)
outfile <- file.path(args$outdir, paste0(clust, "_kegg.csv"))
write.csv(filtered_res, outfile, row.names = FALSE)