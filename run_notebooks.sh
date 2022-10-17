#Run notebooks and Rscripts to create the NMJ figures
#$ bash run_notebooks_pdf.sh TestFigs/
#Use cmd + / to uncomment and run individual blocks

#datadir=/Users/mkthompson/Desktop/Davislab/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_results/run_200915/
# This doesn't work anymore
# source ~/opt/miniconda3/etc/profile.d/conda.sh

# conda activate pretty
# ### Prep data & annotations ####
# jupyter nbconvert --to pdf --execute Drafting/Prep_filter_lowExp.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute Drafting/Prep_rate_SDs.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute Drafting/Prep_genelists.ipynb --output-dir $1
# ##################
# #
# # # ### Get gene attributes ####
# # # jupyter nbconvert --to pdf --execute Drafting/calc_CAI.ipynb --output-dir $1
# # # jupyter nbconvert --to pdf --execute Drafting/calc_attributes.ipynb --output-dir $1
# # # ##################
#
# # ### Figure 1 and S5: Data overview and localized RNAs ####
# # Intron levels in 4sU and total libraries
# jupyter nbconvert --to pdf --execute Drafting/Overview_intronfrac.ipynb --output-dir $1
# # Histograms and scatter plots for synthesis and decay rates
# jupyter nbconvert --to pdf --execute Drafting/Overview_histograms.ipynb --output-dir $1
# ################
#
# ### Figure S1: Incubation conditions ####
# jupyter nbconvert --to pdf --execute Drafting/Inc_control.ipynb --output-dir $1

# ### Figure 2 and S2: GO slim analysis of decay rates  ####
# jupyter nbconvert --to pdf --execute Drafting/GO_stability_and_function.ipynb --output-dir $1
# #################
#
# ### Figure 3: Stability of celltype specific RNAs  ####
# jupyter nbconvert --to pdf --execute Drafting/CTS_stability.ipynb --output-dir $1
# #################

### Figure 4 and S4: Stability of developmentally regulated RNAs, identified by xtin marks ###
# Get the genes surrounded by H3K27me3 marks
# jupyter nbconvert --to pdf --execute Drafting/Devreg_xtin_states.ipynb --output-dir $1
# # Plot stability of genes with H3K27me3 marks
# jupyter nbconvert --to pdf --execute Drafting/Devreg_xtin_stab.ipynb --output-dir $1
# jupyter nbconvert --to pdf --execute Drafting/Devreg_php.ipynb --output-dir $1

### Figure S3: GO Analysis of the most and least stable CTS RNAs ###
# Rscript notebooks/run_clusterprofiler.R -sub Figures/CTS/genesets/CTS_10stable_genes.csv -bg Figures/CTS/genesets/bg_genes.csv -sim_co 0.5 -outdir Figures/CTS/genesets/
# Rscript notebooks/run_clusterprofiler.R -sub Figures/CTS/genesets/CTS_10unstable_genes.csv -bg Figures/CTS/genesets/bg_genes.csv -sim_co 0.5 -outdir Figures/CTS/genesets/

### Figure S5: transiteMotif analysis:
# Rscript scripts/run_transite.R -sub Figures/Motifs/unstable_5perc.csv -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_subset/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna
# Rscript scripts/run_transite2.R -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_spectrum_3pUTR/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna_threeprime -runtype spma
# Rscript scripts/run_transite2.R -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_spectrum_5pUTR/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_fiveprime.fa -name mrna_fiveprime -runtype spma
# Rscript scripts/run_transite2.R -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_spectrum_cds/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_cds.fa -name mrna_cds -runtype spma
Rscript scripts/run_transite2.R -sub Figures/Motifs/CTS_TF_genes.csv -bg Figures/Motifs/bg_genes.csv -outdir Figures/Motifs/transite_cts_tf_tsma/ -fasta  ../../../Davislab_old/3.4_NMJ_4Tu_4sU/3.4e_pipeline_dev/nmj_figures/resources/region_fastas/longest_threeprime.fa -name mrna_threeprime -runtype tsma
