### Run RBP enrichment tests ###
# Run RBP enrichment tests for the significant AUC genes
# Run the tests for the given transcript region, e.g. 5'UTR, CDS, 3'UTR

library(Biostrings)
library(transite)
library(tidyverse)
library(argparse)
library(dplyr)

parser <- ArgumentParser()
#if you write a double flag, i.e. --flag, then this will be the name in args
parser$add_argument("-sub", help="csv file with subset genes")
parser$add_argument("-bg", help="csv file with background genes")
parser$add_argument("-outdir", help="output directory")
parser$add_argument("-fasta", help="fasta file with the sequences of the corresponding region")
parser$add_argument("-name", help="analysis name, e.g. threeprime")
parser$add_argument("-runtype", help="tsma or spma")
args <- parser$parse_args()

seqs <- readDNAStringSet(args$fasta)

outdir <- args$outdir
dir.create(outdir, showWarnings = FALSE)

bg_file <- args$bg
bg_genes = read.csv(bg_file)$gene

seqs <- gsub("T", "U", seqs)
avail_seqs <- names(seqs)
bg_shared <- intersect(avail_seqs, bg_genes)
bg_seqs <- paste(seqs[bg_genes[bg_genes %in% bg_shared]])

if (args$runtype == 'tsma') {
    sub_file <- args$sub
    sub_genes <- read.csv(sub_file)$gene
    subset_shared <- intersect(avail_seqs, sub_genes)
    subset_seqs <- paste(seqs[sub_genes[sub_genes %in% subset_shared]])

    # Write the gene names used in the subset:
    subset_names_used <- sub_genes[sub_genes %in% subset_shared]
    subset_df <- data.frame(subset_names_used)
    colnames(subset_df) <- c("gene")
    write.csv(subset_df, file.path(outdir, paste0(args$name, '_subset.csv')), row.names=FALSE)

    #Run tsma analysis:
    #This function requires passing the significant genes in a list because you could test more than one gene set
    results <- run_matrix_tsma(list(subset_seqs), bg_seqs, cache=FALSE)
    #Save the motifs by gene
    motif_cols <- select(results$foreground_scores[[1]]$df, 'motif_id', 'motif_rbps')
    #Convert the hit list to dataframe so we can see how many hits individual transcripts had
    #https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
    foreground_l <- results$foreground_scores[[1]]$absolute_hits
    hits_df <- data.frame(matrix(unlist(foreground_l), nrow = length(foreground_l), byrow = TRUE))
    colnames(hits_df) <- sub_genes[sub_genes %in% subset_shared]
    rownames(hits_df) <- motif_cols$motif_id
    write.csv(hits_df, file.path(outdir, paste0(args$name, "_tsma", "_hits_by_gene.csv")), row.names = TRUE)
    #Save the enrichment df
    write.csv(results$enrichment_dfs[[1]], file.path(outdir, paste0(args$name, "_tsma", "_enrich_scores.csv")), row.names = FALSE)

} else {
    results <- run_matrix_spma(bg_seqs, cache=FALSE)
    # write the overall results
    write.csv(results$spectrum_info_df, file.path(outdir, paste0(args$name, "_spma", "_enrich_scores.csv")), row.names = FALSE)
    # write the bin-wise enrichment
    bins_df <- dplyr::bind_rows(results$enrichment_dfs, .id="column_label")
    write.csv(bins_df, file.path(outdir, paste0(args$name, "_spma", "_enrich_bins.csv")), row.names = FALSE)
}

#Save the RDS file
outfile <- file.path(outdir, paste0(args$name, "_", args$runtype, ".rds"))
saveRDS(results, file = outfile)

# Write the gene names used in the background:
bg_names_used <- bg_genes[bg_genes %in% bg_shared]
bg_df <- data.frame(bg_names_used)
colnames(bg_df) <- c("gene")
write.csv(bg_df, file.path(outdir, paste0(args$name, "_bg.csv")), row.names=FALSE)
#Can reload the data like this:
#d = readRDS(infile)
#The saved results RDS file does not seem to include the names of the sequences that were used, need to save these separately
