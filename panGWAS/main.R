rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ape)
library(scales)

source("ridge_regression.R")
source("tree_plot_functions.R")
source("utilities.R")

## Parameters --------------------------------------------------------------------

## Input files
snp.file <- "~/Desktop/gwas/data/Resultat_SNP_INDEL_filtered_ann_binaire.txt"
geno.file <- "~/Desktop/gwas/data/gene_presence_absence.Rtab"
tree.file <- "~/Desktop/gwas/data/RAxML_bestTree.inference_minidataset.nwk"
pheno.file <- "~/Desktop/gwas/data/phenotype.txt"

## Output files
gemma.cov      <- "gemma_covariates.txt"
gemma.kinship  <- "gemma_kinship.txt"
gemma.variants <- "gemma_variants.txt"
gemma.pheno    <- "gemma_pheno.txt"
gemma.results  <- "gemma_results"

## Other parameters
kinship.type <- "SNP" ## add options to compute kinship matrix from SNPs or Tree
gemma <- "gemma" ## Path to gemma executable file
outdir <- "./gemma"

## Import and format data -------------------------------------------------------
data <- read_data(snp.file, geno.file, tree.file, pheno.file,
                  kinship.type = kinship.type)

## Export data ------------------------------------------------------------------
export_gemma_data(data = data, outdir = outdir,
                  gemma.cov = gemma.cov, gemma.kinship = gemma.kinship,
                  gemma.variants = gemma.variants, gemma.pheno = gemma.pheno)

## Launch gemma using precomputed kinship matrix --------------------------------
system(paste0(gemma,
              " -g ", file.path(outdir, gemma.variants),   ## gemma.variants,
              " -p ", file.path(outdir, gemma.pheno),      ## phenotype
              " -c ", file.path(outdir, gemma.cov),        ## covariates
              " -k ", file.path(outdir, gemma.kinship),    ## gemma kinship matrix
              " -lmm 1",                                   ## Wald test for significance
              " -outdir ", outdir,                         ## Output directory
              " -o ", gemma.results                        ## Gemma output file
              ),
       ignore.stdout = TRUE,
       ignore.stderr = TRUE)

## Format GEMMA results and correct values for peculiar SNPs.
results <- format_gemma_results(data$variants.description,
                                data$pattern.description,
                                gemma.results = file.path(outdir, gemma.results))

## Compute principal components (PC) and test them using Wald's test -----------

## Wald test for PCs. Note that Wald's test is useful only for SNP-type kinship matrix
## Another test (in progress) should be used for tree-type kinship matrix
wald.values <- pc_wald(y = data$phenotype,                                  ## phenotype
                       x = scale_patterns(pm = data$pattern.matrix,         ## scaled centered
                                          pd = data$pattern.description),   ## pattern matrix
                       w = data$covariates,                                 ## covariates
                       K = data$kinship                                     ## kinship matrix
                       )

## select (at most 6) significant PCs and sort them in decreasing significance order
significant.pcs <- select_pc(wald.values = wald.values,
                             npcs = 6,
                             signif.cutoff = 0.05) ## 5 PCs turn out to be significant

## Graphical representations ---------------------------------------------------

## Patterns --------------------------

## Low level pattern plot
plot_pheno_pattern(tree = data$tree,
                   pheno = data$phenotype,
                   pattern = data$pattern.matrix[1, ])
## High level pattern plot
plot_pattern("pattern20", data = data, lmm.results = results)

## Select and plot all patterns with adjusted p-values < 0.05
significant.patterns <- results %>% filter(pval <= 0.05) %>% arrange(pval) %>%
  pull(var = "pattern") %>% unique() ## 8 patterns

par(mfrow = c(3, 3)) ## Split the plot to create 3 by 3 small graphics, numbered from left to right, top to bottom
for (pattern in significant.patterns) {
  par(mar = c(1, 1, 1, 1)) ## less white space around plots
  plot_pattern(pattern,
               data = data,
               lmm.results = results,
               add.legend = TRUE)
}

par(mfrow = c(1, 1))

## Principal components ---------------

serovar <- fread(pheno.file) ## recuperer d'autre donnee pour la legende couleur des graph acp
serovar <- setNames(object = serovar$Serovar, nm = serovar$Strain)

## Plot two most significants PCs
plot_pc(tree = data$tree,
        scores = significant.pcs,
        missing.color = "gray50",
        pcs = 1:2,
        fatten.by = "size",
        color.tip.by = serovar, ## data$phenotype à la place de serovar pour afficher en fonction du phenotype source
        add.legend = TRUE)

## Plot all significants PCs
plot_pc(tree = data$tree,
        scores = significant.pcs,
        missing.color = "gray50",
        fatten.by = "size",
        color.tip.by = data$phenotype,
        add.legend = TRUE)
