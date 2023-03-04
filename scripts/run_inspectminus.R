# Run INSPEcT on Kallisto transcript quantitations
# Calculate first guess rates
# This version will run INPSEcT- rate estimation for experiment with 
# no matched 4sU librares. This analysis was done to address a reviewer's 
# question about whether the RNA dynamics is changing before and after brain
# incubation

library('INSPEcT')

setwd('/Users/mk/Desktop/Davislab/3.10_brain_stability/brain_results/')

# input files
indir <- '/Users/mk/Desktop/Davislab/3.10_brain_stability/brain_results/301122_brain_inctest/inspect'
tot_exon_file <- file.path(indir, 'tot_exon_tpm.csv')
tot_intron_file <- file.path(indir, 'tot_intron_tpm.csv')
tot_exp_des_file <- file.path(indir, 'tot_expDes.csv')

# output files
outdir <- 'brain_inc_inspectMinus'
dir.create(outdir)
mature_file <- file.path(outdir, 'mature.csv')
premature_file <- file.path(outdir, 'premature.csv')
mature_var_file <- file.path(outdir, 'mature_var.csv')
premature_var_file <- file.path(outdir, 'premature_var.csv')
rdata_file <- file.path(outdir, 'brain_inc_inspectMinus.RData')

# Read in total exonic and intronic counts and convert to matrix
totexon_ma <- as.matrix(read.csv(tot_exon_file, row.names = 'gene'))
totintron_ma <- as.matrix(read.csv(tot_intron_file, row.names = 'gene'))

# Read in experimental design (tpts, then reps): e.g. t1_1 t2_1, t3_1, t1_2, t2_2, t3_2
# or e.g.: wt_1, syp_1, wt_2, syp_2,...
tot_exp_des <- read.csv(tot_exp_des_file)$conditions

totL <- list('exonsAbundances' = totexon_ma, 'intronsAbundances' = totintron_ma)

totExp_plgem<-quantifyExpressionsFromTrAbundance(trAbundaces = totL,
                                                 experimentalDesign = tot_exp_des)

# There is a bug in INSPEcT in which it reduces one condition matrix to vector.
# Transpose the row/columns in order to fix the matrix subsetting problem which
# becomes a problem with one condition.

if (ncol(totExp_plgem$exonsExpressions) != length(unique(tot_exp_des))) {
  tot_exonsExpressions2 = t(totExp_plgem$exonsExpressions)
  tot_exonsVariance2 = t(totExp_plgem$exonsVariance)
  tot_intronsExpressions2 = t(totExp_plgem$intronsExpressions)
  tot_intronsVariance2 = t(totExp_plgem$intronsVariance)
  
  rownames(tot_exonsVariance2) <- rownames(tot_exonsExpressions2)
  rownames(tot_intronsVariance2) <- rownames(tot_intronsExpressions2)

  totExp_plgem <- list(exonsExpressions = tot_exonsExpressions2, exonsVariance = tot_exonsVariance2,
                       intronsExpressions = tot_intronsExpressions2, intronsVariance = tot_intronsVariance2)
}

num_tpts <- length(unique(tot_exp_des))
tpts <- tot_exp_des[1: num_tpts]
# for steady-state data, timepoints will actually be conditions
matureInspObj<-newINSPEcT(tpts = tpts
                           ,labeling_time = NULL
                           ,nascentExpressions = NULL
                           ,matureExpressions = totExp_plgem)

# There doesn't seem to be any rates stored in this object, but the k3/k2 values
# should be able to be calculated by using the ratio premrna/mature, as showin
# in Fig 1B of Furlan et al., 2020

print('writing output files')

# write the rates first guess:
write.csv(matureInspObj@mature, file=mature_file)
write.csv(matureInspObj@premature, file=premature_file)
write.csv(matureInspObj@matureVar, file=mature_var_file)
write.csv(matureInspObj@prematureVar, file=premature_var_file)

saveRDS(matureInspObj, file = rdata_file)
