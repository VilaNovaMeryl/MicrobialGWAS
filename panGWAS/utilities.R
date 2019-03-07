#' Keep and order siginificant PCs
#'
#' @param wald.values Required. List as produced by pc_wald
#' @param npcs Optional. Maximum number of significant PCs to retain. Defaults to 6.
#' @param signif.cutoff Optional. Threshold used for significance. Defaults to 0.01
#'
#' @return Siginificant PCs (at mosts npcs) ordered by increasing p-values. P-values are Bonferroni corrected. P-values are also stored in the \code{pval} attribute of the result.
#' @export
#'
#' @examples
select_pc <- function(wald.values, npcs = 6, signif.cutoff = 0.01) {
  test.stat <- wald.values$Eg^2 / wald.values$Vg ## Wald-stats
  log10pval <- -log10(exp(1)) * pchisq(test.stat, df = 1, lower.tail = FALSE, log.p = TRUE)
  ## order PCs
  pc.rank <- order(-log10pval)
  ## significant PCs
  npcs <- min(npcs, length(log10pval))
  pval.thresh <- -log10(signif.cutoff / npcs)
  signif.pcs <- (log10pval > pval.thresh)
  ## Keep at most npcs significant pc in order
  selected.pcs <- pc.rank[signif.pcs]
  selected.pcs <- selected.pcs[1:min(npcs, length(selected.pcs))]
  ## return results
  pc.scores <- wald.values$PC.scores[ , selected.pcs, drop = FALSE]
  evidence <- log10pval[selected.pcs]
  attr(pc.scores, "evidence") <- evidence
  colnames(pc.scores) <- paste0(colnames(pc.scores), " (evid. = ", prettyNum(evidence, big.mark = ","), ")")
  return(pc.scores)
}

#' Read and format files produced by iVarCall2 for use with GEMMA
#'
#' @param snp.file Required. File path to binary snp file
#' @param geno.file Required. File path to genomic profile file
#' @param tree.file Required. File path to the tree file (in newick or nexus format)
#' @param pheno.file Required. File path to the phenotype file
#' @param kinship.type Optional. Kinship type. Defaults to \code{SNP} (kinship computed from binary snps)
#'                     but support for \code{tree} (kinship computed from phylogenetic tree) will be added later.
#' @param gemma Optional. File path to GEMMA executable. Defaults to "gemma" (assumes that gemma's path in the environment PATH)
#' @param process.variants Optional. Defaults to TRUE. Should variants be processed to identify and remove (i) invariants patterns and (ii) patterns perfectly correlated to the phenotype (can cause numerical errors in GEMMA)
#'
#' @return A list with several items \begin{itemize}
#' \item \code{kinship} kinship matrix
#' \item \code{phenotype} Phenotype named vector
#' \item \code{tree} Phylogenetic tree (object of class \code{\link{class-phylo}})
#' \item \code{variants.description}: a data.frame with one row per feature and several columns:
#'       - the original feature-level descriptors (Feature, Chr, Pos, ...)
#'       - one additional \code{pattern} column with the unique ID of the pattern corresponding to that feature
#' \item \code{pattern.description}: a data.frame with one row per pattern and five columns:
#'       \code{pattern} (unique pattern ID),
#'       \code{Feature} (first feature with that pattern),
#'       \code{SNP.count} (number of SNPs with that pattern),
#'       \code{Gene.count} (number of genes with that pattern) and
#'       \code{Comment} (comment for that pattern)
#' \item \code{pattern.matrix}: a binary matrix coding with one row per unique pattern. Rownames corresponds to
#'       pattern IDs and colnames to strain names.
#' \item \code{gemma.variants}: a gemma formatted pattern data.frame. Only patterns passing the filters are included
#' \end{itemize}
#'
#' @export
#'
#' @examples
read_data <- function(snp.file = NULL,
                      geno.file = NULL,
                      tree.file = NULL,
                      pheno.file = NULL,
                      kinship.type = "SNP", ## add support for tree later
                      ## gemma = "gemma", ## Path to gemma executable file
                      process.variants = TRUE
) {

  ## Read data
  snps  <- fread(snp.file)
  geno  <- fread(geno.file)
  tree  <- read.tree(tree.file) %>% phangorn::midpoint()
  pheno <- fread(pheno.file) %>% as.data.frame()

  ## Format for GEMMA ------------------------------------------------------------------

  ## Variants
  snps <- snps %>% rename(Feature = Variant) %>%
    separate(col = "Feature", into = c("Chr", "Pos", "Allele"), sep = "_", remove = FALSE) %>%
    separate(col = "Allele",  into = c("REF", "ALT"), sep = "/", remove = TRUE) %>%
    mutate(Type = "SNP")
  geno <- geno %>% rename(Feature = Gene) %>%
    mutate(Chr = NA, Pos = NA, REF = "-", ALT = "+", Type = "Gene")
  ## Concatenate SNP and gene presence/absence profile
  variants <- bind_rows(snps, geno)

  ## Phenotype (in same order as in variant file)
  strain.order <- variants %>% select(-one_of("Feature", "Chr", "Pos", "REF", "ALT", "Type")) %>% colnames()
  ii <- match(strain.order, pheno$Strain)
  pheno <- pheno[ii, ]
  ## TODO correct phenotype selection
  phenotype <- setNames(pheno[, 3], pheno$Strain)

  ## Covariates
  ## TODO account for more covariates
  covariates <- rep(1, length(strain.order))


  ## Compress patterns -----------------------------------------------------------------
  compressed.patterns <- compress_patterns(variants)
  variants.description <- compressed.patterns$variants.description
  pattern.description <- compressed.patterns$pattern.description
  pattern.matrix <- compressed.patterns$pattern.matrix

  ## Filter patterns -------------------------------------------------------------------
  if (process.variants) {
    message("Filtering invariant and perfect variants")
    pattern.description$Comment <- filter_patterns(pattern.matrix, phenotype)
  }

  ## Create a vcf-like table with only valid patterns
  filtered.variants <- inner_join(pattern.description, variants, by = "Feature") %>%
    filter(is.na(Comment)) %>%
    select(-SNP.count, -Geno.count, -Feature, -Comment, -Chr, -Pos, -Type)

  ## Compute kinship matrix ------------------------------------------------------------
  kinship <- compute_kinship(pattern.matrix, pattern.description, tree, kinship.type)

  ## Export ----------------------------------------------------------------------------

  ## Return data.frames as list to the user
  return(list("kinship"              = kinship,
              "phenotype"            = phenotype,
              "covariates"           = covariates,
              "tree"                 = tree,
              "variants.description" = variants.description,
              "pattern.description"  = pattern.description,
              "pattern.matrix"       = pattern.matrix,
              "gemma.variants"       = filtered.variants))
}


#' Export gemma-formatted tables to specified files
#'
#' @param data Required. A list containing the phenotype, covariates, variants and kinship matrix for matrix.
#'             As produced by \code{\link{read_data}}
#' @param outdir Optional. Folder where gemma-formatted files should be stored. Defaults to "./gemma/"
#' @param gemma.cov Required. File name for covariates
#' @param gemma.kinship Required. File name for kinship matrix
#' @param gemma.variants Required. File name for variants
#' @param gemma.pheno Required. File name for phenotype
#'
#' @return Nothing, the function is used for its side effect of writing several gemma formatted files:
#' \begin{itemize}
#' \item A tabulated file with covariates to pass on to GEMMA (intercept only for now)
#' \item A tabulated file with the kinship matrix to pass on to GEMMA
#' \item A tabulated file with the variants matrix to pass on to GEMMA. Variants are compressed, to reduce the amounts of test to perform
#' \item A tabulated file with the phenotype to pass on to GEMMA. Discrete phenotype (e.g. Bovine/Swine) are recoded as integers.
#' \end{itemize}
#'
#' @export
#'
#' @examples
export_gemma_data <- function(data,
                              outdir = "./gemma",
                              gemma.cov,
                              gemma.kinship,
                              gemma.variants,
                              gemma.pheno
                              ) {
  ## Create output dir (invisibly)
  dir.create(outdir, showWarnings = FALSE)

  ## Covariates
  write.table(data$covariates,
              file = file.path(outdir, gemma.cov),
              sep="\t", row=F, col=F, quote=F)

  ## Phenotype
  ## convert character to integer
  if (is.character(data$phenotype)) { data$phenotype <- factor(data$phenotype) }
  if (is.factor(data$phenotype)) { data$phenotype <- as.integer(data$phenotype)-1 }
  write.table(data$phenotype,
              file = file.path(outdir, gemma.pheno),
              sep="\t", row=F, col=F, quote=F)

  ## Variants
  data$gemma.variants %>%
    fwrite(file = file.path(outdir, gemma.variants),
           col.names = FALSE,
           sep = ",",
           quote = "auto"
    )

  ## kinship
  write.table(data$kinship,
              file = file.path(outdir, gemma.kinship),
              sep = "\t",
              col.names = FALSE, row.names = FALSE)

}


#' Format gemma results and correct values for peculiar SNPs
#'
#' @param vd Required. Variant description data.frame, as produced by \code{compress_pattern}
#' @param pd Required. Pattern description data.frame, as produced by \code{compress_patterns}
#' @param gemma.results Filename of results produced by GEMMA
#'
#' @return A data.frame with GEMMA results and additional feature level informations.
#'         The column \code{pval} contains p-values adjusted by BH correction.
#' @export
#'
#' @examples
format_gemma_results <- function(vd, pd, gemma.results) {
  ## read gemma results
  results <- fread(paste0(gemma.results, ".assoc.txt")) %>%
    select(-one_of(c("chr", "ps", "n_miss", "allele0", "allele1"))) %>%
    rename(pattern = rs)
  ## Merge with pattern description and correct values for peculiar patterns
  ## Correction for beta values
  format_beta <- function(x, Comment) {
    res <- x
    res[Comment == "Perf. correlated"] <- 1
    res[Comment == "Perf. anti-correlated"] <- -1
    return(res)
  }
  results <- left_join(pd, results, by = "pattern") %>%
    select(-Feature, -SNP.count, -Geno.count) %>%
    ## correction for correlated sites
    mutate(beta   = format_beta(beta, Comment),
           se     = ifelse(is.na(Comment), se, 0),
           p_wald = ifelse(is.na(Comment), p_wald, 0)
    ) %>%
    ## correction for invariant sites
    mutate_at(.vars = c("beta", "se", "p_wald"),
              .funs = function(x, condition) { if_else(condition, true = NA_real_, false = x, missing = x) },
              condition = (quo(Comment == "Invariant")))
  ## Merge with variant description
  results <- vd %>% inner_join(results, by = "pattern") %>% mutate(pval = p.adjust(p_wald))

  return(results)
}


#' Compress a binary pattern data.frame for efficient downstream computation.
#'
#' @title compress_patterns
#' @param variants Required. A binary feature data.frame, as produced by extracting the GT component
#'                 from a vcf file with two additional columns \code{Feature} (unique ID for each feature) and
#'                 \code{Type} (feature type, typically SNP or Gene).
#'
#' @return A list with three components: \begin{itemize}
#' \item pattern.description: a data.frame with one row per pattern and five columns:
#'       \code{pattern} (unique pattern ID),
#'       \code{Feature} (first feature with that pattern),
#'       \code{SNP.count} (number of SNPs with that pattern),
#'       \code{Gene.count} (number of genes with that pattern) and
#'       \code{Comment} (comment for that pattern)
#' \item variants.description: a data.frame with one row per feature and several columns:
#'       - the original feature-level descriptors (Feature, Chr, Pos, ...)
#'       - one additional \code{pattern} column with the unique ID of the pattern corresponding to that feature
#' \item pattern.matrix: a binary matrix coding with one row per unique pattern. Rownames corresponds to
#'       pattern IDs and colnames to strain names.
#' \end{itemize}
#' @export
#'
#' @examples
compress_patterns <- function(variants) {

  ## variants description
  variants.description <- variants %>%
    unite(col = "pattern", -one_of(c("Feature", "Chr", "Pos", "REF", "ALT", "Type")), sep = "") %>%
    mutate(pattern = paste0("pattern", as.integer(factor(pattern, levels = unique(pattern)))))

  ## pattern description
  pattern.description <- variants.description %>% group_by(pattern) %>%
    summarize(SNP.count = sum(Type == "SNP"),   ## Number of SNP with that pattern
              Geno.count = sum(Type == "Gene"), ## Number of profile with that pattern
              Feature = first(Feature),         ## First feature with that pattern
              Comment = NA                      ## Is pattern special (typically used for invariant
              ##    and "perfect" patterns)
    )

  ## pattern matrix
  pm <- inner_join(pattern.description, variants, by = "Feature") %>%
    select(-Feature, -SNP.count, -Geno.count, -Chr, -Pos, -Type, -Comment) %>%
    select(-REF, -ALT)
  pattern.matrix <- pm %>% select(-pattern) %>% data.matrix()
  rownames(pattern.matrix) <- pm$pattern

  ## return results
  return(list("variants.description" = variants.description,
              "pattern.description"  = pattern.description,
              "pattern.matrix"       = pattern.matrix))

}


#' Filter invariant and highly correlated to the phenotype patterns
#'
#' @param pattern.matrix Required. Binary pattern matrix
#' @param phenotype Optional. Phenotype used for testing high phenotype/pattern correlations.
#'                  If NULL (default), no test is performed.
#'
#' @return A named character with levels \code{NA} (pattern passed the filter),
#'         \code{"Invariant"} (pattern is invariant),
#'         \code{"Perf. correlated"} (pattern is almost perfectly correlated to (rho > 0.999) to the phenotype)
#'         \code{"Perf. anti-correlated"} (pattern is almost perfectly anti-correlated to (rho < -0.999) to the phenotype)
#' @export
#'
#' @examples
filter_patterns <- function(pattern.matrix, phenotype = NULL)
{
  pattern.status <- setNames(object = rep(NA, nrow(pattern.matrix)),
                             nm = rownames(pattern.matrix))

  ## Invariant patterns
  inv.pattern <- apply(pattern.matrix, 1, sd) == 0
  if (any(inv.pattern)) {
    message(paste0("Pattern ", rownames(pattern.matrix)[inv.pattern], " is invariant."))
    pattern.status[inv.pattern] <- "Invariant"
  }

  ## Correlated patterns
  if (!is.null(phenotype)) {
    ## Transform character to factor
    if (is.character(phenotype)) { phenotype <- as.factor(phenotype) }
    ## transform factor to numeric
    if (is.factor(phenotype)) { phenotype <- as.numeric(phenotype) - 1 }
    ## correlations
    cor.values <- apply(pattern.matrix, 1, function(x) { suppressWarnings(cor(phenotype, x)) } )
    cor.values[is.na(cor.values)] <- 0
    perfect.patterns <- (abs(cor.values) > 1 - 1e-07)
    if (any(perfect.patterns)) {
    message(paste0("Pattern(s) ",
                   paste(rownames(pattern.matrix)[perfect.patterns],collapse = " "),
                   " is (are) perfectly correlated to phenotype."))
    pattern.status[perfect.patterns] <- paste0("Perf. ", ifelse(cor.values[perfect.patterns] > 0,
                                                                "correlated",
                                                                "anti-correlated"))
    }
  }

  return(pattern.status)
}


#' Scale and center binary pattern matrix. Patterns are scaled by the square of their count
#'
#' @param pm Required. Binary pattern matrix, as produced by \code{\link{compress_patterns}}
#' @param pd Required. Pattern description data.frame, as produced by \code{\link{compress_patterns}}
#' @param SNP.only Optional. Should the scaling counts be used by considering only SNPs (TRUE, default) or
#'                 both SNPs and Genes (FALSE).
#'
#' @return A numeric matrix. Pattern with 0 count are removed from the matrix.
#' @export
#'
#' @examples
scale_patterns <- function(pm, pd, SNP.only = TRUE) {
  ## If SNP.only, use only patterns found in SNPs, else aggregate Genes and SNPs.
  if (!SNP.only) {
     pd <- pd %>% mutate(SNP.count = SNP.count + Gene.count)
  }
  ## Use only patterns found at least once in SNPs
  snp.patterns <- pd %>% filter(SNP.count > 0)
  pattern.counts <- snp.patterns %>% pull(var = SNP.count)
  pattern.names <- snp.patterns %>% pull(var = pattern)
  ## Extract SNP features
  snp.pm <- pm[pattern.names, ] %>% t()
  ## scale (with SNP counts) centered pattern matrix
  return(scale(snp.pm, center = TRUE, scale = 1/sqrt(pattern.counts)))
}

#' Compute the kinship matrix from either binary SNPs or a phylogenetic tree.
#'
#' @param pattern.matrix Required. Binary pattern matrix, as produced by \code{\link{compress_patterns}}
#' @param pattern.description Required. Pattern description data.frame, as produced by \code{\link{compress_pattern}}
#' @param tree Optional. A tree of class \code{\link{phylo-class}}
#' @param kinship.type One of \code{SNP} (default) or \code{tree}. Method used to compute the kinship matrix.
#'
#' @return A square numeric kinship matrix.
#' @export
#'
#' @examples
compute_kinship <- function(pattern.matrix = NULL,
                            pattern.description = NULL,
                            tree = NULL,
                            kinship.type = c("SNP", "tree")) {
  kinship.type <- match.arg(kinship.type)
  ## specific functions to compute kinship
  kinship_tree <- function(pm, tree) {
    K <- vcv(tree)
    ii <- match(colnames(pm), colnames(K))
    return(K[ii, ii])
  }
  kinship_snps <- function(pm, pd) {
    ## scaled centered pattern matrix
    scpm <- scale_patterns(pm, pd, SNP.only = TRUE)
    pattern.counts <- apply(scpm, 2, function(x) { diff(range(x))^2 })
    ## Compute kinship matrix
    return(tcrossprod(scpm) / sum(pattern.counts))
  }

  if (kinship.type == "tree") { kinship <- kinship_tree(pattern.matrix, tree) }
  if (kinship.type == "SNP")  { kinship <- kinship_snps(pattern.matrix, pattern.description) }
  return(kinship)
}

#' Helper function to extract the binary pattern of a feature / pattern given by its name
#'
#' @param string Required. Feature or pattern name
#' @param vd Required. Variant description data.frame, as produced by \code{\link{compress_patterns}}
#' @param pm Required. Pattern matrix, as produced by \code{\link{compress_patterns}}
#'
#' @return A binary named pattern vector
#' @export
#'
#' @examples
get_pattern <- function(string, vd, pm) {
  if (!pmatch("pattern", string, nomatch = FALSE)) { ## string is a feature name
    string <- vd %>% filter(Feature == string) %>% pull("pattern")
    if (length(string) == 0) {
      stop(paste0(string, " is neither a valid feature nor pattern name."))
    }
  }
  pattern <- pm[string, ]
  return(pattern)
}

































################################################################################################
## Get Bayesian Wald Test inputs.
## @fit.lmm: ridge regression results
## @pca: principal component analysis
## @svd.XX: single value decomposition of biallelic SNPs
##
## Outputs:
## pca.Ebeta
## pca.Vbeta
################################################################################################

get_wald_input <- function(fit.lmm = NULL,
                           pca = NULL,
                           svd.XX = NULL,
                           y = NULL,
                           npcs = NULL,
                           XX = NULL){

  # Get full posterior covariance matrix for Bayesian Wald Test
  # Need the full posterior covariance matrix for the Bayesian Wald test,
  # to get the posterior uncertainty for each point estimate
  lambda = fit.lmm$lambda_MLE
  Cstar = diag(lambda * svd.XX$d^2 / (lambda * svd.XX$d^2 + 1))
  # For the null model, the posterior mean and variance (may be slow!)
  astar = t(y)%*%y - t(y)%*%XX%*%fit.lmm$Ebeta
  dstar = length(y)
  tau = as.numeric((dstar-2)/astar)

  rotation = t(pca$rotation[,1:npcs])
  # Based on the PCA rotations of raw genetic diversity
  pca.Ebeta = rotation %*% fit.lmm$Ebeta

  rtr = tcrossprod(rotation,rotation); # Should be n (sample size) by n
  rv = rotation %*% svd.XX$v; # Should be n by n
  pca.Vbeta = rv %*% Cstar; # Should be n by n
  pca.Vbeta = tcrossprod(pca.Vbeta,rv); # Should be n by n
  pca.Vbeta = lambda/tau * (rtr - pca.Vbeta); # Should be n by n

  return(list("Ebeta" = pca.Ebeta, "Vbeta" = pca.Vbeta))


}

################################################################################################
## Get order of principal components by Bayesian Wald Test results.
## @fit.lmm: ridge regression results
## @pca: principal component analysis
## @svd.XX: single value decomposition of biallelic SNPs
##
## Outputs:
## pca.Ebeta
## pca.Vbeta
################################################################################################

get_pc_order <- function(p.pca.bwt = NULL,
                         signif_cutoff = NULL){
  o <- order(p.pca.bwt,decreasing=T)
  # How many are above a significance cut-off
  if(length(which(p.pca.bwt>= signif_cutoff))>20){
    pc.lim <- 1:20
  } else {
    if(length(which(p.pca.bwt>= signif_cutoff))>0){
      pc.lim <- 1:length(which(p.pca.bwt>= signif_cutoff))
    } else {
      pc.lim <- NULL
    }
  }

  return(list("pc_order" = o, "pc.lim" = pc.lim))
}

################################################################################################
## Do PCA
################################################################################################

do_pca <- function(pcs = NULL, XX = NULL, XX.ID = NULL){
  if(is.null(pcs)){
    pca <- prcomp(XX)
  } else {
    ## Read in PCs
    pca <- read.table(pcs, header = T, as.is = T)
    if(any(XX.ID != rownames(pca))){
      m <- match(XX.ID, rownames(pca))
      pca <- pca[m, ]
    }
    pca <- list("x" = pca, "rotation" = NULL)
  }
  return(list("pca" = pca))
}




################################################################################################
## Do Wald test
################################################################################################

wald_test <- function(y = NULL,
                      XX = NULL,
                      svd.XX = NULL,
                      lambda = NULL,
                      XX.all = NULL,
                      prefix = NULL,
                      npcs = NULL,
                      pca = NULL){


  fit.lmm <- ridge_regression(y, XX, svdX=svd.XX,
                              lambda_init=as.numeric(lambda)/sum(XX.all$bippat),
                              maximize=FALSE, skip.var=TRUE)

  # Fit the grand null model

  fit.0 <- lm(y~1)

  # LRT for the LMM null vs grand null
  LRTnullVgrand <- -log10(pchisq(2*(fit.lmm$ML - as.numeric(logLik(fit.0))), 1, low=F)/2)
  cat(paste0("## LRT for the LMM null vs grand null = ", LRTnullVgrand),
      file = paste0(prefix, "_logfile.txt"), sep="\n", append = TRUE)

  # Heritability
  fit.lmm.ypred <- XX %*% fit.lmm$Ebeta
  cat(paste0("## Heritability (R^2) = ", cor(fit.lmm.ypred,y)^2),
      file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)

  # Get full posterior covariance matrix for Bayesian Wald Test
  # Need the full posterior covariance matrix for the Bayesian Wald test,
  # to get the posterior uncertainty for each point estimate

  wald_input <- get_wald_input(fit.lmm = fit.lmm, pca = pca, svd.XX = svd.XX,
                               y = y, npcs = npcs, XX = XX)

  # Bayesian Wald Test
  pca.bwt <- wald_input$Ebeta^2/diag(wald_input$Vbeta)

  p.pca.bwt <- -log10(exp(1))*pchisq(pca.bwt, 1, low=F, log=T)

  cat(paste0("## Bayesian Wald Test for PCs range = ", paste(range(p.pca.bwt), collapse=" ")),
      file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)

  write.table(p.pca.bwt, file = paste0(prefix, "_Bayesian_Wald_Test_negLog10.txt"),
              sep="\t", row=T, col = F, quote=F)

  # Get order of PCs by Wald test results
  signif_cutoff <- -log10(0.05/npcs)

  pc_order <- get_pc_order(p.pca.bwt = p.pca.bwt, signif_cutoff = signif_cutoff)
  # Predict phenotype using effect sizes
  effect <- t(t(XX) * as.vector(fit.lmm$Ebeta))
  pred <- rowSums(effect)


  return(list("pc_order" = pc_order, "p.pca.bwt" = p.pca.bwt, "pred" = pred,
              "signif_cutoff" = signif_cutoff))

}


ridge_regression = function(y,x,w=matrix(1,length(y),1),K_type="uncentred",svdX=NULL,lambda_init=200,prefer.REML=TRUE,skip.var=FALSE,maximize=TRUE) {
  if(lambda_init<=0) stop("ridge_regression: lambda_init must be positive")
  n = length(y)
  # Process x as uncentred, centred or standardized, and label the result X
  # (not yet implemented: take x at face value)
  if(K_type!="uncentred") stop("Currently the data x must be pre- centred or standardized if desired")
  X = x
  #	Note that after centring/standardizing X, K = X t(X) (i.e. K_ij = sum_l X_il X_jl)
  #	and that by singular value decomposition (SVD), X = U D V`, with U and V orthogonal,
  #	i.e. V^(-1) = V`. So one can write
  #	K = X t(X) = U D V` V` D U` = U D^2 U`.
  #	If you define R = U D, then K = X t(X) = R t(R)
  #	Calculating K via SVD will be faster in general because R is only an n x n matrix
  # Do the singular value decomposition if necessary
  if(is.null(svdX)) {
    svdX = svd(X)
  }
  rk = length(svdX$d)
  if(rk<n) {
    svdX$d = c(svdX$d,rep(0,n-rk))
    svdX$u = cbind(svdX$u,matrix(0,n,n-rk))
    svdX$v = cbind(svdX$v,matrix(0,nrow(svdX$v),n-rk))
  }
  # Calculate the relatedness matrix
  svdR = svdX$u %*% diag(svdX$d)
  #	svdR = matrix(0,n,n)
  #	svdR[1:n,1:length(svdX$d)] = svdX$u %*% diag(svdX$d)
  K = (svdR %*% t(svdR)) # Unlike GEMMA, not dividing by number of loci

  # Precompute other quantities
  Wx = as.matrix(w)
  c = ncol(w)

  # Estimate lambda using maximum likelihood
  lambda_loglik = function(log_lambda) {
    lambda = exp(log_lambda)
    H = lambda * K + diag(n)
    Hinv = solve(H)
    Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
    ret = 0.5*n*log(0.5*n/pi) - 0.5*n - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*n*log(t(y) %*% Px %*% y)
    if(!is.finite(ret)) return(-.Machine$double.xmax)
    return(ret)
  }
  #(lambda_MLE = optimize(lambda_loglik,c(1e-5,1e5),maximum=TRUE)$maximum)
  #	nlm_MLE = suppressWarnings(nlm(function(x) -lambda_loglik(x),log(lambda_init)))
  if(maximize) {
    nlm_MLE = suppressWarnings(nlm(function(x) tryCatch(-lambda_loglik(x),error=function(e) -Inf),log(lambda_init)))
  } else {
    nlm_MLE = list("estimate"=log(lambda_init),"minimum"=-lambda_loglik(log(lambda_init)))
  }
  lambda_MLE = exp(nlm_MLE$estimate)
  ML = -nlm_MLE$minimum
  # Estimate lambda using REML
  lambda_logpartiallik = function(log_lambda) {
    lambda = exp(log_lambda)
    H = lambda * K + diag(n)
    Hinv = solve(H)
    Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
    # Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
    ret = 0.5*(n-c)*log(0.5*(n-c)/pi) - 0.5*(n-c) + 0.5*as.numeric(determinant(t(Wx)%*%Wx,log=TRUE)$modulus) - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*as.numeric(determinant(t(Wx)%*%Hinv%*%Wx,log=TRUE)$modulus) - 0.5*(n-c)*log(t(y) %*% Px %*% y)
    if(!is.finite(ret)) return(-.Machine$double.xmax)
    return(ret)
  }
  #(lambda_REML = optimize(lambda_logpartiallik,c(1e-5,1e5),maximum=TRUE)$maximum)
  #	nlm_REML = suppressWarnings(nlm(function(x) -lambda_logpartiallik(x),log(lambda_init)))
  if(maximize) {
    nlm_REML = suppressWarnings(nlm(function(x) tryCatch(-lambda_logpartiallik(x),error=function(e) -Inf),log(lambda_init)))
  } else {
    nlm_REML = list("estimate"=log(lambda_init),"minimum"=-lambda_logpartiallik(log(lambda_init)))
  }
  lambda_REML = exp(nlm_REML$estimate)
  REML = -nlm_REML$minimum

  # Estimate a and tau given lambda using maximum likelihood
  lambda = lambda_MLE
  H = lambda * K + diag(n)
  Hinv = solve(H)
  a_MLE = solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y
  Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
  tau_MLE = n / (t(y) %*% Px %*% y)
  # Estimate tau given lambda using REML
  lambda = lambda_REML
  H = lambda * K + diag(n)
  Hinv = solve(H)
  a_MLE = solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y
  Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
  # Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
  tau_REML = (n - c) / (t(y) %*% Px %*% y)

  # Calculate the posterior mean and variance in effect size for each SNP
  lambda = lambda_REML
  tau = tau_REML
  if(!prefer.REML) {
    lambda = lambda_MLE
    tau = tau_MLE
  }
  # Note that with prior mean m=0, the posterior mean is the ridge regression estimate (Eq 11.57 of O`Hagan & Forster with tau=sigma^(-2))
  #	mstar	= solve(t(X) X + 1/lambda I_L) t(X) y
  # where L is the number of loci. This would involve inverting a very large L x L matrix. It is more efficient to use SVD:
  #	mstar	= solve(X` X + 1/lambda I_L) X` y
  # next introduce identity matrices
  #			= (V solve(V)) solve(X` X + 1/lambda I_L) (solve(V`) V`) X` y
  # then use the rule solve(A) solve(B) solve(C) = solve(C B A)
  #			= V solve(V` (X` X + 1/lambda I_L) V) V` X` y
  # next expand the brackets
  #			= V solve(V` X` X V + 1/lambda V` I_L V) V` X` y
  # and note that V` I_L V = I_n
  #			= V solve(V` X` X V + 1/lambda I_n) V` X` y
  # using SVD, write X = U D V` and X` = V D U`
  #			= V solve(V` V D U` U D V` V + 1/lambda I_n) V` V D U` y
  # and use the property that V` V = I_n
  #			= V solve(D U` U D + 1/lambda I_n) D U` y
  # then introduce R = U D and R` = D U`
  #			= V solve(R` R + 1/lambda I_n) R` y
  # finally define C = solve(R` R + 1/lambda I_n) to give
  #			= V C R` y
  # Since R`R = D U` U D, using the property U`U = I_n, we get R`R = D^2
  C = solve(t(svdR) %*% svdR + (1/lambda * diag(n)))
  post.mean = (svdX$v %*% (C %*% t(svdR) %*% y))
  # The posterior variance, Wstar = solve(solve(W) + tau t(X) X) (Eq 11.61 of O`Hagan & Forster, p327 with tau=sigma^(-2))
  # Assuming a prior variance W = lambda/tau I_L,
  #	Wstar	= lambda/tau solve(I_L + lambda t(X) X)
  # Applying the SVD,
  #			= lambda/tau solve(I_L + lambda V D U` U D V`)
  #			= lambda/tau solve(I_L + lambda V D D V`)
  # substituting E = lambda D^2 (a diagonal matrix)
  #			= lambda/tau solve(I_L + lambda V E V`)
  # Following section 1.5 of Henderson and Searle 1981 SIAM Review 23:53-60, Eq. 17, with A = I_n, B = E and U = V. In their notation,
  #			solve(A + U B U`) = solve(A) - solve(A) U solve(solve(B) + U` solve(A) U) U` solve(A)
  # gives in our notation
  #			solve(I_L + V E V`) = I_L - I_L V solve(solve(E) + V` I_L V) V` I_L = I_L - V solve(solve(E) + I_n) V`
  # so
  #	Wstar	= lambda/tau (I_L - V Cstar V`)
  # where
  #	Cstar	= lambda D^2/(lambda D^2 + 1) I_n
  Cstar = diag(lambda * svdX$d^2 / (lambda * svdX$d^2 + 1))
  if(skip.var) {
    post.var = NA
  } else {
    #		post.var = lambda/tau * apply(svdX$v,1,function(V) 1 - V %*% Cstar %*% V)
    #		post.var = lambda/tau * sapply(1:nrow(svdX$v),function(i) as.vector(1 - svdX$v[i,] %*% Cstar %*% svdX$v[i,]))
    #		post.var = lambda/tau * sapply(1:nrow(svdX$v),function(i) 1 - sum(svdX$v[i,]^2*diag(Cstar)))
    post.var = lambda/tau * (1 - colSums(t(svdX$v^2)*diag(Cstar)))
  }

  ret = list("lambda_MLE"=lambda_MLE,"lambda_REML"=lambda_REML,"a_MLE"=a_MLE,"tau_MLE"=tau_MLE[1,1],"tau_REML"=tau_REML[1,1],"Ebeta"=post.mean,"Vbeta"=post.var,"ML"=ML,"REML"=REML,"func"=function(x){x})
  return(ret)
}

