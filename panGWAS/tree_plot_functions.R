#' Simple wrapper around \code{\link{read.tree}} to renumber tips of a tree. Ensure that tips are numbered
#' from 1 to ntaxa(tree) (from bottom to top). Most useful after rerooting a tree.
#'
#' @title reorder_tips
#' @param phy Required. A object of class \code{\link{phylo-class}}.
#'
#' @return A tree with correctly numbered tips
#' @export
#'
#' @examples
reorder_tips <- function(phy) {
  tree.string <- write.tree(tree, file = "")
  return(read.tree(text = tree.string))
}

#' Compute incidence matrix of a tree
#'
#' @title incidence_matrix
#' @param phy Required. A \code{phylo} class object
#' @return An incidence matrix M of size nedges(tree) x ntaxa(tree) where
#'         M[i,j] is set to 1 if taxa derives j from edge i and 0 otherwise.
#' @note incidence_matrix assumes that the tree is rooted. If not, it roots it
#'       arbitrarily (at taxa 1).
incidence_matrix <- function(phy) {
  if (!is.rooted(phy)) {
    warning("Tree is not rooted, incidence matrix may be meaningless")
  }
  ## Construct incidence matrix of the tree (taxa x edge matrix)
  ## All taxa descending from an edge are set to 1, all others to -1
  ntaxa <- length(phy$tip.label)
  nedges <- nrow(phy$edge)
  incidence <- matrix (0,
                       nrow = ntaxa,
                       ncol = nedges,
                       dimnames = list(phy$tip.label, ## taxa names
                                       phy$edge[, 2]  ## edge number (child node)
                       ))
  ## Incidence of internal edges
  phy.part <- prop.part(phy) ## clade composition indexed by (shifted) node number
  for (i in 2:length(phy.part)) { ## first clade corresponds to root node
    edge <- which(phy$edge[, 2] == i + ntaxa) ## node numbers are shifted by ntaxa
    incidence[phy.part[[i]] , edge] <- 1
  }
  ## Incidence of pendant edges
  ## pendant.edges[i] is the edge leading to tip i.
  pendant.edges <- match(seq_len(ntaxa), phy$edge[ , 2])
  for (i in seq_len(ntaxa)) {
    incidence[i, pendant.edges[i]] <- 1
  }
  attr(incidence, "pendant.edges") <- pendant.edges
  return(incidence)
}


#' Function for plotting a tree with edge fattened and colored according to
#' their correlation with a PC scores
#'
#' @title plot_pc
#' @param tree Required. A \code{\link{phylo-class}} class tree.
#' @param scores Required. PC scores of the taxa
#' @param pcs Optional. PC to be plotted. Defaults to all PC in scores.
#' @param width.lim Optional. Numeric. Minimum and maximal edge after fattening.
#'                  Defaults to c(0.5, 4). Set to c(0, x) for true linear scaling.
#' @param fatten.by Required. Aesthetics used for fattening.
#'                        Subset of 'size' and 'alpha'.
#' @param fatten.tips Optional. Logical. Should tips be fattened like pendant edges.
#'                    Defaults to FALSE.
#' @param color Optional. Color vector of length 2 used to distinguish positive
#'              loading edges (color[1]) and negative loading edges (color[2]).
#'              Defaults to c("#EF8A62", "#67A9CF")
#' @param color.tip.by Optional. A factor with the same length as the number of taxa in the tree.
#'              Defaults to NULL. If NULL, nothing happens.
#' @param missing.color Optional. Color used for edges with low correlation with the PC.
#'                      Defaults to "gray50". If NULL, nothing happens. Use "white", "transparent",
#'                      or par("bg") to remove them from plot.
#' @param cor.cutoff Optional. Edge with absolute correlation lower than cor.cutoff (default 0.2)
#'                   are considered missing.
#' @param add.legend Optional. Shoud a legend be added to the plot (TRUE, default) or not (FALSE)
#' @param ... Additional arguments passed on to plot.phylo
#'
#' @return Nothing, this function is used for its side effect of plotting a tree
#' @seealso \code{\link{plot.pc}}
plot_pc <- function(tree,
                    scores,
                    pcs = 1:ncol(scores),
                    fatten.by = c("size"), fatten.tips = FALSE,
                    width.lim = c(0.5, 4), color = c("#EF8A62", "#67A9CF"),
                    color.tip.by = NULL, missing.color = "gray50", cor.cutoff = 0.2,
                    add.legend = TRUE,
                    ...) {
  ## Check fatten.by
  if (any( !fatten.by %in% c("size", "alpha"))) {
    stop("fatten.by must be a subset of c(\"size\", \"alpha\")")
  }
  ## Extract scores
  scores <- scores[ , pcs, drop = FALSE]
  ## Compute correlation between scores and tree edges
  edge.patterns <- incidence_matrix(tree)
  edge.patterns <- edge.patterns[rownames(scores), ]
  x <- cor(edge.patterns, scores)
  ## Extract tree
  pendant.edges <- match(1:length(tree$tip.label) , tree$edge[ , 2])
  ## Get edge.color
  edge.color <- ifelse(x > 0, color[1], color[2])
  ## Get fattening factor
  fattening.factor <- rescale(abs(x), from = c(0, 1), to = width.lim)
  ## Fatten edges by size
  if ("size" %in% fatten.by) {
    edge.width <- fattening.factor
  } else {
    edge.width <- rep(1, length(x))
  }
  dim(edge.width) <- dim(x)
  ## If alpha, fatten edges by alpha
  if ("alpha" %in% fatten.by) {
    scaled.alpha <- fattening.factor / width.lim[2]
    edge.color <- alpha(edge.color, scaled.alpha)
    dim(edge.color) <- dim(x)
  }
  ## Color tips
  if (is.null(color.tip.by)) {
    ## if no color.tip.by, color tip like pendant edges
    tip.color <- edge.color[pendant.edges, , drop = FALSE]
  } else {
    if (is.null(names(color.tip.by))) {
      warning("No names provided for color.tip.by, assuming that colors are given in tree order.")
    } else {
      color.tip.by <- color.tip.by[tree$tip.label]
    }
    tip.levels <- unique(color.tip.by)
    legend.palette <- hue_pal()(length(tip.levels))
    tip.color <- col_factor(legend.palette, domain = tip.levels)(color.tip.by)
    tip.color <- matrix(rep(tip.color, ncol(x)), ncol = ncol(x))
  }
  ## Fatten tips
  if (fatten.tips) {
    tip.width <- edge.width[pendant.edges, , drop = FALSE]
    if ("alpha" %in% fatten.by) {
      dim(scaled.alpha) <- dim(x)
      scaled.alpha <- scaled.alpha[pendant.edges, , drop = FALSE]
      tip.color <- alpha(tip.color, scaled.alpha)
      dim(tip.color) <- dim(scaled.alpha)
    }
  }
  ## Set missing edges to missing.color
  if (!is.null(missing.color)) {
    edge.color[abs(x) < cor.cutoff] <- missing.color
  }
  ## Plot tree

  ## Prepare arguments for plot phylo
  args <- list(x = tree, edge.width = 1, edge.color = "black")
  args <- c(args, list(...))

  ## Prepare layout if there are multiple trees
  oldmar <- par("mar")
  n.axis <- ncol(x)
  layout.nrow <- floor(sqrt(n.axis))
  layout.ncol <- ceiling(n.axis / layout.nrow)
  par(mfrow = c(layout.nrow, layout.ncol), mar = c(1, 0, 1, 1))
  for (i in 1:n.axis) {
    args[["edge.width"]] <- edge.width[ , i]
    args[["edge.color"]] <- edge.color[ , i]
    args[["tip.color"]] <- tip.color[ , i]
    args[["traits"]] <- setNames(scores[, i], nm = rownames(scores))
    if (fatten.tips) {
      args[["cex"]] <- tip.width[ , i]
    }
    do.call("plot_tree_barplot", args)
    plot.title <- if (is.null(colnames(scores))) { paste("PC", pcs[i]) } else { colnames(scores)[i] }
    title(plot.title, line = -1)
    if (add.legend) {
      legend(x = "bottomleft",
             inset = 0.05,
             legend = c("0", "0.5", "1",
                        paste0(">", cor.cutoff), paste0(cor.cutoff, ">, >", -cor.cutoff), paste0(-cor.cutoff, ">")),
             fill = c(rep("transparent", 3), color[1], missing.color, color[2]),
             border = c(rep("transparent", 3), rep("black", 3)),
             ncol = 2,
             xjust = 0,
             title = "PC-edge correlation\n(abs. value / sign)",
             box.lwd = 0,
             lwd = c(rescale(c(0, 0.5, 1), from = c(0, 1), to = width.lim), 0, 0, 0),
             merge = TRUE,
             seg.len = 1)
    }
  }
  par(mar = oldmar)
}

#' Title plot_tree_barplot
#'
#' @param x Required. An object of class \code{\link{phylo-class}}
#' @param traits Required. A named vector of trait values.
#' @param ... Optional. Arguments passed on to \code{\link{plot.phylo}}. Argument \code{xlim} and \code{ylim}
#'             are replaced internally
#'
#' @return Plots a tree with an associated bar plot for a continuously valued character at the tips.
#' @export
#'
#' @comment Code heavily inspired from Paul's Bastide plot function in PhylogeneticEM.
#' @examples
plot_tree_barplot <- function(x, traits, ...) {
  ## give more meaningful names to variables
  tree <- x

  ## Sanity check for trait values
  if (!all(tree$tip.label %in% names(traits))) {
    stop("Trait values are missing for some taxa.")
  }
  x <- traits[tree$tip.label]
  ## Preliminary computation, devote 3/4 of the space to the tree and 1/4 to the barplot
  ntaxa <- length(x)

  ## Capture args
  args <- c(list(x = tree), list(...))

  ## Fake plot to compute required coordinates
  args[["x.lim"]] <- NULL
  args[["y.lim"]] <- NULL
  args[["plot"]]  <- FALSE
  do.call("plot.phylo", args)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ## Padding: add extra space on the right to add barplot and on top and bottom to add scale
  tree_height <- lastPP$x.lim[2]
  barplot_width <- tree_height / 3
  x.lim.max <- tree_height + barplot_width
  y.lim.min <- -ntaxa/10
  y.lim.max <- ntaxa + ntaxa/10

  ## Plot tree with extra space on the right
  args[["x.lim"]] <- c(0, x.lim.max)
  args[["y.lim"]] <- c(y.lim.min, y.lim.max)
  args[["plot"]] <- TRUE

  par(new = TRUE) ## necessary to ensure that the new call to plot.phylo overrides the previous
  ## call used for getting coordinates
  do.call("plot.phylo", args)

  ## Barplot at tip
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  barplot.x.min <- max(lastPP$xx) + max(strwidth(tree$tip.label, cex = lastPP$cex))
  barplot.x.max <- lastPP$x.lim[2]
  offset <- (barplot.x.max - barplot.x.min)/10
  barplot.range <- c(barplot.x.min + offset, barplot.x.max)

  ## Plots characters
  trait.scaled <- rescale(x,
                          to = barplot.range,
                          from = range(x))
  trait.start <- rescale(0, to = barplot.range, from = range(x))
  min.trait <- min(x)
  max.trait <- max(x)

  ## 0 bar
  segments(x0 = trait.start, y0 = y.lim.min + ntaxa/15,
           x1 = trait.start, y1 = y.lim.max - ntaxa/10,
           lty = 3)
  ## axis
  axis(1, at = barplot.range,
       labels = round(c(min.trait, max.trait), digits = 2),
       pos = y.lim.min + ntaxa/15 #,
       # las = axis_las
  )
  ## trait values
  segments(x0  = trait.start,
           x1  = trait.scaled,
           y0  = lastPP$yy[1:ntaxa],
           col = if (is.null(args$tip.color)) { "black" } else { args$tip.color })
}

#' Plot binary phenotype and binary pattern side by side at the tips of the tree.
#'
#' @param tree Required. A object of class \code{\link{phylo-class}}
#' @param pheno Required. A named binary phenotype to plot at the tree tips.
#' @param pattern  Required. A named binary pattern to plot at the tree tips.
#' @param add.legend Optional. Shoud a legend be added to the plot (TRUE, default) or not (FALSE)
#' @param ... Optional. Additional arguments passed on to plot.phylo
#'
#' @return Nothing. Used for the plot side effect
#' @export
#'
#' @examples
plot_pheno_pattern <- function(tree, pheno, pattern, add.legend = TRUE, ...) {
  ## Only works for binary phenotypes
  if (length(unique(pheno)) != 2) {
    stop("Only binary phenotypes are supported for now. Consider using plot_tree_bar for quantitative phenotypes")
  } else {
    pheno <- factor(pheno)
  }
  ## HACK: expand tip labels by three spaces to have space to add symbols
  tmp.tree <- tree; tmp.tree$tip.label <- paste0("   ", tree$tip.label)
  plot(tmp.tree, ...)
  pal.pheno <- setNames(c("#EF8A62", "#67A9CF"), nm = levels(pheno))
  pal.pattern <- setNames(c("white", "black"), c("REF/-", "ALT/+"))
  tiplabels(pch = 21, bg = pal.pheno[pheno], tip = match(names(pheno), tree$tip.label),
            adj = c(0.5, 0.5), offset = -1.5*strwidth(" "))
  tiplabels(pch = 22, bg = pal.pattern[1+pattern], tip = match(names(pattern), tree$tip.label),
            adj = c(0.5, 0.5), offset = 1.5*strwidth(" "))
  if (add.legend) {
    legend(x = "bottomleft",
           pch = c(21, 21, 22, 22),
           pt.bg = c(pal.pheno, pal.pattern),
           legend = c(names(pal.pheno), names(pal.pattern)),
           ncol = 2, title = c("Pheno. / Pattern"))
  }
}

#' High level function to plot a pattern side-by-side with a phenotype at the tips of the tree.
#' Patterns/Features are called by name.
#'
#' @param pattern.name Required. Feature or pattern name
#' @param data Required. A list containing the phenotype, covariates, variants and kinship matrix for matrix.
#'             As produced by \code{\link{read_data}}
#' @param lmm.results Required. A list containing the results of Gemma,
#'                    as produced by \code{\link{format_gemma_results}}
#' @param add.legend Optional. Shoud a legend be added to the plot (TRUE, default) or not (FALSE)
#' @param ... Additional arguments passed on to plot.phylo
#'
#' @return Nothing. The function is used for its plot side effect
#' @export
#'
#' @examples
plot_pattern <- function(pattern.name, data, lmm.results, add.legend = TRUE, ...) {
  bip <- get_pattern(pattern.name, vd = data$variants.description, pm = data$pattern.matrix)
  bip.info <- results %>% filter(pattern == !!pattern.name | Feature == !!pattern.name) %>% slice(1)
  ## construct title
  plot.title <- paste0("Feature/Pattern ", pattern.name,
                       " (p = ", prettyNum(bip.info["pval"], digits = 2), ")")
  plot_pheno_pattern(tree = data$tree, pheno = data$phenotype, pattern = bip, add.legend = add.legend, ...)
  title(plot.title,
        line = -1, ## make sure title is in plot region
        adj = 0    ## left align title
        )
}
