#Run notebooks and Rscripts to create the NMJ figures
#$ bash run_notebooks.sh TestFigs/
#Use cmd + / to uncomment and run individual blocks

# conda activate pretty
# ### Prep data & annotations ####
# jupyter nbconvert --to pdf --execute notebooks/Prep_filter_lowExp.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/Prep_rate_SDs.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/Prep_genelists.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/Prep_summarize_exps.ipynb --output-dir $1
# ##################

# # ### Figure 1 and S2: Data overview and localized RNAs ####
# Intron levels in 4sU and total libraries
# jupyter nbconvert --to pdf --execute notebooks/Overview_intronfrac.ipynb --output-dir $1
# # Histograms and scatter plots for synthesis and decay rates
# jupyter nbconvert --to pdf --execute notebooks/Overview_histograms.ipynb --output-dir $1
# # ################
# #
# # ### Figure S1: Incubation conditions ####
# jupyter nbconvert --to pdf --execute notebooks/Inc_control.ipynb --output-dir $1

# ### Figure 2 and S3: GO slim analysis of decay rates  ####
# jupyter nbconvert --to pdf --execute notebooks/GO_stability_and_function.ipynb --output-dir $1
# #################

# ### Figure 3 and S5: Stability of celltype specific RNAs  ####
# jupyter nbconvert --to pdf --execute notebooks/CTS_genes_corrales.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/CTS_genes_dillon.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/CTS_stability.ipynb --output-dir $1
# #################

### Figure 4 and S7: Stability of developmentally regulated RNAs, identified by xtin marks ###
# Get the genes surrounded by H3K27me3 marks
# jupyter nbconvert --to pdf --execute notebooks/Devreg_xtin_states.ipynb --output-dir $1
# # Plot stability of genes with H3K27me3 marks
# jupyter nbconvert --to pdf --execute notebooks/Devreg_xtin_stab.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute notebooks/Devreg_php.ipynb --output-dir $1

### Figure 5: Examples of developmentally regulated RNAs and their stability
# jupyter nbconvert --to pdf --execute notebooks/Examples_devTFs.ipynb --output-dir $1

### Figure S6: Relationships between RNA decay and gene attributes
# Show the overlap between H3K27me3, CTS, and TF
# jupyter nbconvert --to pdf --execute notebooks/Reg_me_vs_TF.ipynb --output-dir $1
# # Compare overlapping gene groups-- TF, CTS, and me3
# jupyter nbconvert --to pdf --execute notebooks/Reg_compare_groups.ipynb --output-dir $1

### Figure S4: GO Analysis of the most and least stable CTS RNAs ###
Rscript scripts/run_clusterprofiler.R -sub Figures/CTS/genesets/CTS_10stable_genes.csv -bg Figures/CTS/genesets/bg_genes.csv -sim_co 0.5 -outdir Figures/CTS/genesets/
Rscript scripts/run_clusterprofiler.R -sub Figures/CTS/genesets/CTS_10unstable_genes.csv -bg Figures/CTS/genesets/bg_genes.csv -sim_co 0.5 -outdir Figures/CTS/genesets/
jupyter nbconvert --to pdf --execute notebooks/CTS_goplot.ipynb --output-dir $1

### Figure S5: transiteMotif analysis:
# jupyter nbconvert --to pdf --execute notebooks/Motifs_genesets.ipynb --output-dir $1

# Rscript scripts/run_transite.R -sub Figures/Motifs/CTS_TF_genes.csv -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_ctstf_tmsa_3pUTR/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna_threeprime -runtype tsma
## Rscript scripts/run_transite.R -sub Figures/Motifs/unstable_5perc.csv -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_subset/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna
## Rscript scripts/run_transite.R -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_spectrum_3pUTR/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna_threeprime -runtype spma

### Scrape Interactive Fly to get the segment polarity genes ###
# jupyter nbconvert --to pdf --execute notebooks/Examples_scrape_segment_polarity.ipynb --output-dir $1