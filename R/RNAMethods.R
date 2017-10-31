#' @include RNAClass.R
#' @title sepSpike
#' @name sepSpike
#' @rdname sepSpike-methods
#' @description
#' Separate out spikein from genes based upon rownames of tpm.value slot of tpm object.
#' @param obj A \code{tpm} object.
#' @param invert A \code{logical} object, indicating whether to keep spike-in only (TRUE) or genes only (FALSE).
#' @export sepSpike
setGeneric(name = 'sepSpike',
  def = function(obj, invert = FALSE) {
    standardGeneric("sepSpike")
  }
)

#' @rdname sepSpike-methods
setMethod(f = "sepSpike",
  signature = c("tpm", "logical"),
  def = function(obj, invert) {
    tpm.value <- obj@tpm.value[grep("spike", rownames(obj@tpm.value), ignore.case = TRUE, invert = invert),]
    return(tpm.value)
  }
)

#' @title rmLow
#' @name rmLow
#' @rdname rmLow-methods
#' @description
#' Filter of genes or transcripts with expression levels above given threshold in one or more samples.
#' @param obj A \code{tpm} object.
#' @param thresh A TPM threshold.
#' @export rmLow
setGeneric(name = "rmLow",
  def = function(obj, thresh) {
    standardGeneric("rmLow")
  }
)

#' @rdname rmLow-methods
setMethod(f = "rmLow",
  signature = c("tpm", "numeric"),
  definition = function(obj, thresh) {
    kpt.idx <- apply(obj@tpm.value, 1,
      function(vec)any(vec>thresh)
    )
    stopifnot(sum(kpt.idx) > 0)
    obj@tpm.value <- obj@tpm.value[kpt.idx, ]
    return(obj)
  }
)

#' @title rmNonVar
#' @name rmNonVar
#' @rdname rmNonVar-methods
#' @description
#' Filter genes or transcripts with expression variabilities above a given quantile.
#' @param obj A \code{tpm} object.
#' @param probs A probability threshold for removing invariably expressed genes.
#' @export rmNonVar
setGeneric(name = "rmNonVar",
  def = function(obj, probs = 0.1) {
    standardGeneric("rmNonVar")
  }
)

#' @rdname rmNonVar-methods
setMethod(f = "rmNonVar",
  signature = c("tpm", "numeric"),
  definition = function(obj, probs) {
    var <- apply(obj@tpm.value, 1, var)
    thresh <- quantile(var, probs = probs, names = FALSE)
    kpt.idx <- var >= thresh
    obj@tpm.value <- obj@tpm.value[kpt.idx, ]
    return(obj)
  }
)

#' @title lmAdjCovar
#' @name lmAdjCovar
#' @rdname lmAdjCovar-methods
#' @description
#' Linear regression to adjust for batch effects.
#' @param x A matrix of gene expression data with genes in the rows and samples in the cols.
#' @param covar A data.frame covariates.
#' @param add.mean A logical vector indicating whether to add group mean to the residuals
#' @export lmAdjCovar
setGeneric(name = "lmAdjCovar",
  def = function(x, covar, add.mean = TRUE) {
    standardGeneric("lmAdjCovar")
  }
)

#' @rdname lmAdjCovar-methods
setMethod(f = "lmAdjCovar",
  signature = c("matrix", "data.frame", "logical"),
  definition = function(x, covar, add.mean = TRUE){
    stopifnot(ncol(x) == nrow(covar))
    colnames(covar) <- paste0('V', 1:ncol(covar))
    x <- t(x)
    t(lm(x~., data = covar)$residuals) + if(add.mean){colMeans(x)}else{0}}
)

#' @title loessnorm
#' @name loessnorm
#' @rdname loessnorm-methods
#' @description
#' Loess normalization of gene or transcript expressions to expressions of spike-in.
#' @param obj A \code{tpm} object.
#' @param small A numeric vector indicating the adjustment to the TPM values before log2 transformation.
#' @export loessnorm
setGeneric(name="loessnorm",
  def = function(obj, small = 0.05) {
    standardGeneric("loessnorm")
  }
)

#' @rdname loessnorm-methods
setMethod(f = "loessnorm",
  signature = c("tpm", "ANY"),
  definition = function(obj, small = 0.05){
    tpm.norm <- normalize.loess(obj@tpm.value + small,
      subset = grep("spikein", row.names(tpm)))
    return(tpm.norm)
  }
)

#' @title rdcntnorm
#' @name rdcntnorm
#' @rdname rdcntnorm-methods
#' @description
#' Normalization of gene or transcript expressions to expressions of spike-in based upon read counts in a table.
#' @param obj A \code{tpm} object.
#' @param stats A data.frame with named columns of 'spikein' and 'mm10' read counts.
#' @export rdcntnorm
setGeneric(name="rdcntnorm",
  def = function(obj, stats) {
    standardGeneric("rdcntnorm")
  }
)

#' @rdname rdcntnorm-methods
setMethod(f = "rdcntnorm",
  signature = c("tpm", "data.frame"),
  definition = function(obj, stats) {
    frac <- stats$spikein/(stats$spikein + stats$mm10)
    normfac <- mean(frac)/frac
    tpm.norm <- data.frame(t(t(obj@tpm.value)*normfac))
    return(tpm.norm)
  }
)

#' @title distplot
#' @name distplot
#' @rdname distplot-methods
#' @description
#' tpm.value either normalised to spikein or not
#' tpm.value either spikein, non-spikein, or all
#' idx for selecting differnetially expressed genes.
#' @param obj A \code{tpm} object.
#' @param ylab A character string specify y-labl on the plot.
#' @param pdffout A character string specify pdf output file.
#' @param probs A probability specifying the top quantile to use a ceilling in plotting.
#' @export distplot
setGeneric(name="distplot",
  def = function(obj, ylab, pdffout, probs = 0.85) {
    standardGeneric("distplot")
  }
)

#' @rdname distplot-methods
setMethod(f = "distplot",
  signature = c("tpm", "character", "character", "numeric"),
  definition = function(obj, ylab, pdffout, probs) {
    dat <- data.table(obj@tpm.value)
    ldat <- melt(dat)
    max.y <- quantile(ldat$value, probs = probs)
    p1 <- ggplot(ldat, aes_(x = ~variable, y = ~value)) +
      geom_boxplot(aes_(fill = factor(~variable))) +
      scale_y_continuous(labels = comma) +
      labs(x = "", y = ylab) +
      theme(legend.title = element_blank(), legend.position="top") +
      coord_cartesian(ylim=c(0, max.y))
      pdf(pdffout)
      theme_set(theme_grey(base_size=15))
    multiplot(p1,cols=1)
    dev.off()
  }
)

#' @title corplot
#' @name corplot
#' @rdname corplot-methods
#' @description
#' Scatter plot of pairwise sample comparisons
#' @param obj A \code{tpm} object.
#' @param pdffout A character string specify pdf output file.
#' @export corplot
setGeneric(name = "corplot",
  def = function(obj, pdffout) {
    standardGeneric("corplot")
  }
)

#' @rdname corplot-methods
setMethod(f = "corplot",
  signature = c("tpm", "character"),
  definition = function(obj, pdffout) {
    tpm.value <- obj@tpm.value
    thresh <- floor(quantile(as.matrix(tpm.value), probs = 0.999))
    pm <- ggpairs(data.table(tpm.value))
    pm2 <- pm
    for(i in 2:pm$nrow) {
      for(j in 1:(i-1)) {
        pm2[i,j] <- pm[i,j] + coord_cartesian(xlim = c(0,thresh),ylim = c(0,thresh)) +
          scale_x_continuous(breaks = c(0, as.numeric(thresh))) +
          scale_y_continuous(breaks = c(0, as.numeric(thresh)))
      }
    }
    pdf(pdffout, height = 7 * (ncol(obj@tpm.value)/10), width = 7 *(ncol(obj@tpm.value)/10))
    print(pm2)
    dev.off()
  }
)

#' @title hireplot
#' @name hireplot
#' @rdname hireplot-methods
#' @description
#' Hierarchical clustering of sample correlation coefficients
#' @param obj A \code{tpm} object.
#' @param pdffout A character string specify pdf output file.
#' @export hireplot
setGeneric(name = "hireplot",
  def = function(obj, pdffout) {
    standardGeneric("hireplot")
  }
)

#' @rdname hireplot-methods
setMethod(f = "hireplot",
  signature = c("tpm"),
  definition = function(obj, pdffout) {
    corstats <- cor(obj@tpm.value,method="spearman")
    grps <- factor(obj@grps,levels=unique(obj@grps),ordered=T)
    tip.col <- if(length(levels(grps)) < 3) {
      brewer.pal(3, "Dark2")[1:2]
    } else if (length(levels(grps)) <= 8) {
      brewer.pal(length(levels(grps)), "Dark2")
    } else if (length(levels(grps)) <= 12) {
      brewer.pal(length(levels(grps)), "Paired")
    } else if (length(levels(grps)) > 12) {
      colorRampPalette(brewer.pal(12, "Paired"))(length(levels(grps)))
    }
    tip.col <- tip.col[as.numeric(grps)]
    pdf(pdffout, height = 7 * (ncol(obj@tpm.value)/10), width = 7 *(ncol(obj@tpm.value)/10))
    plot(as.phylo(hclust(as.dist(1-corstats), method = "average")),
      cex = 2, label.offset = 0, tip.color = tip.col)
    dev.off()
  }
)

#' @title heatcorplot
#' @name heatcorplot
#' @rdname heatcorplot-methods
#' @description
#' Correlation coefficient heatmap of samples
#' @param obj A \code{tpm} object.
#' @param pdffout A character string specify pdf output file.
#' @export heatcorplot
# alternative code:
# corstats=cor(tpm.value[idx,],method="spearman")
# ana.col=data.frame(Group=grps)
# rownames(ana.col)=colnames(tpm.value)
# png(pngfout,width=3000,height=3000,res=300)
# pheatmap(corstats,colorRampPalette(rev(brewer.pal(n=7, name ="RdYlBu")))(100),breaks=seq(-1,1,length.out=101),annotation_col=ana.col)
# dev.off()
setGeneric(name = "heatcorplot",
  def = function(obj, pdffout) {
    standardGeneric("heatcorplot")
  }
)

#' @rdname heatcorplot-methods
setMethod(f = "heatcorplot",
  signature = "tpm",
  definition = function(obj, pdffout) {
    pm <- ggcorr(obj@tpm.value, method = c("pairwise","spearman"),
      label = TRUE, label_alpha = TRUE)
    pdf(pdffout)
    print(pm)
    dev.off()
  }
)

#' @title bplot
#' @name bplot
#' @rdname bplot-methods
#' @description
#' Boxplot of tpm.value slot of tpm object.
#' @param obj A \code{tpm} object.
#' @param title A title for boxplot.
#' @param pdffout A character string specify pdf output file.
#' @param probs A probability specifying the top quantile to use a ceilling in plotting.
#' @param ylab A string specifying the ylab in the plot.
#' @param isLog A logical value indicating whether TPM values are already log2 transformed.
#' @param small A numeric vector indicating the adjustment to the TPM values before log2 transformation.
#' @export bplot
setGeneric(name = "bplot",
  def = function(obj, title, pdffout, probs = 0.80,
    ylab = expression(paste(log[2], "(TPM)")),
    isLog = FALSE, small = 0.05) {
    standardGeneric("bplot")
  }
)

#' @rdname bplot-methods
setMethod(f = "bplot",
  signature = "tpm",
  definition = function(obj, title, pdffout, probs, ylab, isLog, small) {
    tpm.value <- obj@tpm.value
    if(!isLog) {
      tpm.value <- log2(tpm.value + small)
    }
    map.it <- obj@grps
    names(map.it) <- colnames(tpm.value)
    ldat <- melt(data.table(tpm.value))
    ldat$grps <- factor(ldat$variable, levels = colnames(tpm.value))
    levels(ldat$grps) <- obj@grps
    min.y <- min(ldat$value)
    max.y <- as.numeric(quantile(ldat$value, probs = probs))
    p1 <- ggplot(ldat, aes_(x = ~variable, y = ~value, fill = ~grps))+
      geom_boxplot()+
      coord_cartesian(ylim = c(min.y, max.y))+
      theme(legend.title = element_blank(), legend.position = "top") +
      labs(x = "", y = ylab) +
      ggtitle(title)
    pdf(pdffout, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)

#' @title PCAplot
#' @name PCAplot
#' @rdname PCAplot-methods
#' @description
#' PCA plot (PC1, PC2 and PC3) of samples
#' @param obj A \code{tpm} object.
#' @param pdffout A character string specify pdf output file.
#' @param fout A text output file of top N most variably expressed genes.
#' @param excl.col A numeric vector indicating columns to be excluded from the \code{tpm} slot of a \code{tpm} object.
#' @param ntop A numeric value indicating selection of the top N most variably expressed genes, or Inf (default).
#' @param isLog A logical value of whether the \code{tpm} data is already log2 transformed.
#' @param small A numeric vector indicating the adjustment to the TPM values before log2 transformation.
#' @export PCAplot
setGeneric(name = "PCAplot",
  def = function(obj, pdffout, fout = NULL, excl.col = NULL,
    ntop = Inf, isLog = FALSE, small = 0.05) {
    standardGeneric("PCAplot")
  }
)

#' @rdname PCAplot-methods
setMethod(f = "PCAplot",
  signature = "tpm",
  definition = function(obj, pdffout, fout, excl.col, ntop, isLog, small) {
    out <- function(dat, fout) {
      tmp <- data.frame(gene = rownames(dat), dat)
      write.table(tmp, file = fout, row.names = F,
        col.names = T, quote = F, sep = "\t")
    }
    if (is.null(excl.col)) {
      dat <- as.matrix(obj@tpm.value)
    } else {
      dat <- as.matrix(obj@tpm.value[, -excl.col])
    }
    if(ncol(dat) <= 3) {
      stop('The number of samples must be larger than 3\n')
    }
    ridx <- apply(dat, 1, function(vec)any(vec > 0))
    stopifnot(any(ridx))
    dat <- dat[ridx,]
    if(!isLog) {
      dat <- log2(dat + small)
    }
    dat <- dat[apply(dat, 1, var) > 0, ]
    if (ntop < nrow(dat)){
       vars <- apply(dat, 1, var)
       dat <- dat[order(vars, decreasing = TRUE), ]
       dat <- dat[1:ntop, ]
    }
    pca <- prcomp(t(dat), center=TRUE, scale=TRUE, retx=TRUE)
    ve <- pca$sdev^2/sum(pca$sdev^2) #variance explained by each PC
    ve13 <- sprintf("%.1f",ve[1:3]*100)
    x <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3])
    rownames(x) <- colnames(dat)
    if (!is.infinite(ntop) && !is.null(fout)) {
      out(dat, fout)
    }
    uniq.cols <- rainbow(length(unique(obj@grps)))
    cols <- uniq.cols[as.numeric(factor(obj@grps, levels = unique(obj@grps), ordered = TRUE))]
    x$colors <- cols
    pch <- as.numeric(factor(rownames(x)), ordered = TRUE)
    pdf(pdffout, pointsize = 14, height = max(7, 7 * (ncol(obj@tpm.value)/30)), width = max(7, 7 *(ncol(obj@tpm.value)/30)))
    par(mar = c(1, 1, 1, 1))
    s3d <- scatterplot3d(x[ , 1:3], grid = FALSE, box = FALSE, mar = c(3, 3, 2, 2), pch = "")
    addgrids3d(x[, 1:3], grid = c("xy", "xz", "yz"))
    #s3d$points3d(x[, 1:3], pch=16, col=x$colors)
    text(s3d$xyz.convert(x[, 1:3] + 1), labels=1:nrow(x), col=x$colors)
    legend(par('usr')[1] - 0.5,par('usr')[4] + 0.3,
      legend = paste0(1:nrow(x), ": ", rownames(x)),
      pch = 20, col = x$colors, xpd = TRUE,
      ncol = 3, cex = 0.7)
    dev.off()
  }
)

#' @title MAchart
#' @name MAchart
#' @rdname MAchart-methods
#' @description
#' plot of log2FC over pvalue during differential analysis.
#' @param dd A \code{data.frame} object, with columns of 'P.Value', 'logFC', 'nlogpval' and 'DEG'.
#' @param pdffout A character string specify pdf output file.
#' @export MAchart
setGeneric(name = "MAchart",
  def = function(dd, pdffout) {
    standardGeneric("MAchart")
  }
)

#' @rdname MAchart-methods
setMethod(f = "MAchart",
  signature = c("data.frame", "character"),
  definition = function(dd, pdffout) {
    nd <- dd
    nd$nlogpval <- -log10(nd$P.Value)
    ylab <- "P"
    p1 <- ggplot(nd, aes_(x = ~logFC, y = ~nlogpval, color = ~DEG)) + geom_point() +
      labs(y = bquote(-log[10](.(ylab)))) +
      theme(legend.title = element_blank(), legend.position = "top")
    pdf(pdffout, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)

#' @title BICplot
#' @name BICplot
#' @rdname BICplot-methods
#' @description
#' Plot BIC values over cluster size.
#' @param g A numeric vector of cluster sizes.
#' @param BIC A numeric vector of BIC values.
#' @param pdffout A character string specify pdf output file.
#' @export BICplot
setGeneric(name = "BICplot",
  def = function(g, BIC, pdffout) {
    standardGeneric("BICplot")
  }
)

#' @rdname BICplot-methods
setMethod(f = "BICplot",
  signature = c("numeric", "numeric", "character"),
  definition = function(g, BIC, pdffout) {
    dd <- data.frame(g = g ,BIC = BIC)
    p1 <- ggplot(dd, aes_(x = g, y = BIC)) + geom_point(colour = "#FF9999") +
      xlab("Cluster Size")+ylab("BIC")
    pdf(pdffout, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)

#' @title diffHeatmap
#' @name diffHeatmap
#' @rdname diffHeatmap-methods
#' @description
#' heatmap during or after differential analysis
#' @param tpm.value A \code{tpm} matrix.
#' @param col.idx Numeric or logical indices specifying columns to keep of the \code{tpm} matrix.
#' @param row.idx Numeric or logical indices specifying rows to keep of the code{tpm} matrix.
#' @param pdffout A character string specify pdf output file.
#' @param cutreek A logical value indicating whether to perform clustering.
#' @param cut.alg A string value of selecting the clustering algorithm "pam","hclust" or "emmix".
#' @param rank.man A logical value indicating whether to perform manual ranking on the row.
#' @param log.it.already A logical value indicating whether the \code{tpm} matrix is already log2 transformed.
#' @param scale.it A logical value indicating whether to row standardize the \code{tpm} matrix.
#' @param cluster_columns_par A logical value indicating whether to cluster the columns of the \code{tpm} matrix.
#' @param cluster_rows_par A logical value indicating whether to cluster the rows of the \code{tpm} matrix.
#' @param show_column_dend_par A logical value indicating whether to show the column dendrogram.
#' @param show_row_dend_par A logical value indicating whether to show the row dendrogram.
#' @param small A numeric value indicating the adjustment to the TPM values before log2 transformation.
#' @param ... Additional arguments to be passed to methods.
#' @export diffHeatmap
setGeneric(name = "diffHeatmap",
  def = function(tpm.value, col.idx, row.idx, pdffout,
    cutreek = NULL, cut.alg = NULL, rank.man = FALSE, log.it.already = FALSE,
    scale.it = TRUE, cluster_columns_par = TRUE, cluster_rows_par = TRUE,
    show_column_dend_par = TRUE, show_row_dend_par = FALSE, small = 0.05, ...) {
    standardGeneric("diffHeatmap")
  }
)

#' @rdname diffHeatmap-methods
setMethod(f = "diffHeatmap",
  signature = c("matrix"),
  definition = function(tpm.value, col.idx, row.idx, pdffout, cutreek,
    cut.alg, rank.man, log.it.already, scale.it, cluster_columns_par,
    cluster_rows_par, show_column_dend_par, show_row_dend_par, small, ...) {
    if (!is.null(cut.alg)) {
      cut.alg <- match.arg(cut.alg, c("pam","hclust","emmix"))
    }
    tpm.value <- tpm.value[row.idx, col.idx]
    if (!log.it.already) {
      tpm.value <- log2(tpm.value+small)
    }
    if (scale.it) {
      norm <- t(scale(t(tpm.value)))
    } else {
      norm <- tpm.value
    }
    if (rank.man) {
      norm <- norm[order(norm[,1], -norm[,2], decreasing=T),];
      cluster_rows <- FALSE
    }
    if (!is.null(cutreek)) {
      if (cut.alg == "hclust") {
        corstats <- cor(t(norm),method="spearman");
        clusters <- cutree(hclust(as.dist(1-corstats), method = 'average'), k = cutreek)
      } else if (cut.alg == "emmix") {
        if (!is.null(clusters)) {
          clusters <- clusters
        }
      } else if (cut.alg == "pam") {
        clusters <- pam(norm,cutreek)$clustering
      }
      pm <- Heatmap(norm, cluster_columns = cluster_columns_par,
        cluster_rows = cluster_rows_par, show_row_names = FALSE,
        show_row_dend = show_row_dend_par, show_column_dend = show_column_dend_par,
        heatmap_legend_param = list(title = "", color_bar = "continuous"),
        clustering_distance_rows = "spearman", clustering_method_rows = "average",
        clustering_distance_columns = "spearman", clustering_method_columns = "average",
        split = paste0("Cluster", clusters))
    } else {
      pm <- Heatmap(norm, cluster_columns = cluster_columns_par,
        cluster_rows = cluster_rows_par, show_row_names = FALSE,
        show_row_dend = show_row_dend_par, show_column_dend = show_column_dend_par,
        heatmap_legend_param = list(title="", color_bar="continuous"),
        clustering_distance_rows = "spearman", clustering_method_rows = "average",
        clustering_distance_columns = "spearman", clustering_method_columns = "average")
    }
    pdf(pdffout, pointsize = 14)
    draw(pm)
    dev.off()
    if (exists("clusters")) {
      invisible(clusters)
    } else {
      invisible(NULL)
    }
  }
)

#' @title clusing
#' @name clusing
#' @rdname clusing-methods
#' @description
#' Determine the optimum number of k-mean clusters, and return cluster membership.
#' @param dat A \code{data.frame} of \code{tpm} values.
#' @param pdffout A character string specify the diagnosis pdf file of the optimum number of clusters.
setGeneric(name="clusing",
  def=function(dat, pdffout) {
    standardGeneric("clusing")
  }
)

#' @rdname clusing-methods
setMethod(f = "clusing",
  signature=c(dat = "data.frame", pdffout="character"),
  definition=function(dat, pdffout){
    if(!is.null(dirname(pdffout)) || dirname(pdffout)==".") dir.create(dirname(pdffout), showWarnings = FALSE)
    pdf(pdffout)
    opt <- Optimal_Clusters_KM(dat, max_clusters = min(10,ncol(dat)), plot_clusters=TRUE, criterion = 'distortion_fK', fK_threshold = 0.85, initializer= 'optimal_init', tol_optimal_init = 0.2)
    dev.off()
    km_mb <- MiniBatchKmeans(dat, clusters = opt, batch_size = 20, num_init = 5, max_iters = 100,
    init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,verbose = F)
    pr_mb <- predict_MBatchKMeans(dat, km_mb$centroids)
    return(pr_mb)
  }
)


#' @title kHeat
#' @name kHeat
#' @rdname kHeat-methods
#' @description
#' K-mean clustering and heatmap
#' @param obj A \code{tpm} object.
#' @param pdffout A character string specify pdf output file.
#' @param k A numeric value specifying the number of clusters.
#' @param log2.it A logical value specifying whether to log2 transform the tpm value data.
#' @param scale.it A logical value specifying whether to row standardize the tpm value data.
#' @param small A numeric value indicating the adjustment to the TPM values before log2 transformation.
#' @export kHeat
setGeneric(name = "kHeat",
  def = function(obj, pdffout, k = NULL, log2.it = TRUE, scale.it = TRUE, small = 0.05) {
    standardGeneric("kHeat")
  }
)

#' @rdname kHeat-methods
setMethod(f = "kHeat",
  signature = c("tpm", "ANY"),
  definition = function(obj, pdffout, k, log2.it, scale.it, small) {
    mat <- obj@tpm.value
    if (log2.it) {
      mat <- log2(mat + small)
    }
    if (scale.it) {
      mat <- t(apply(mat, 1, scale))
    }
    colnames(mat) <- colnames(obj@tpm.value)
    if (is.null(k)) {
        message("clusing function when k is null.")
        cl <- clusing(data.frame(mat), pdffout = gsub(".pdf","_optK.pdf", pdffout))
        message("Heatmap function when k is null.")
    } else {
        message("heatmap function when k is ", k)
        set.seed(888)
        cl <- kmeans(mat, k)$cluster
    }
    ht_list <- Heatmap(mat, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE,
      show_row_dend = FALSE,  cluster_columns = FALSE, show_column_dend = FALSE,
      heatmap_legend_param = list(title = "", color_bar = "continuous"),
      clustering_distance_rows = "spearman", clustering_method_rows = "average",
      clustering_distance_columns = "spearman", clustering_method_columns = "average",
      split = factor(cl), gap = unit(3, "mm"))
    pdf(pdffout)
    draw(ht_list)
    dev.off()
    message("running prepare function...")
    ht.obj <- prepare(ht_list)
    message("done with prepare function...")
    row_list <- row_order(ht_list)
    dd <- do.call("rbind", lapply(seq_along(row_list),
      function(i){data.frame(cluster = i, ori.idx = row_list[[i]], idx = seq_along(row_list[[i]]))}))
    dat <- data.frame(gene = rownames(mat)[dd$ori.idx], dd, mat = mat[dd$ori.idx,])
    dat <- dat[order(dat$ori.idx), grep("ori.idx", colnames(dat), invert = TRUE, fixed = TRUE)]
    colnames(dat)[-c(1:3)] <- colnames(mat)
    write.table(dat,
      file = gsub("pdf", "txt", pdffout),
      row.names = FALSE, col.names = TRUE,
      sep = "\t", quote = FALSE)
  }
)

#' @title limmaDiff
#' @name limmaDiff
#' @rdname limmaDiff-methods
#' @description
#' tpm.value either normalised to spikein or not
#' filtering of genes must be done before hand
#' @param obj A \code{tpm} object.
#' @param dout A character string specifying the output directory.
#' @param pat A character string specifying the prefix of file(s) without the directory name.
#' @param MA.it A logical value indicating whether to generate an MA plot.
#' @param HEAT.it A logical value specifying whether to draw a heatmap for each pairwise differnetially analysis.
#' @param GO.it A logical value specifying whether to do gene ontology analysis.
#' @param DiffOut.it A logical value specifying whether to write out differential analysis data to file.
#' @param logFCthresh A numeric value specifying the log2 fold change threshold for data to be called differential.
#' @param PValthresh A numeric value specifying the pvalue threshold for data to be called differential.
#' @param log2.it A logical value specifying whether to perform log2 transformation.
#' @param small A numeric value indicating the adjustment to the TPM values before log2 transformation.
#' @export limmaDiff
setGeneric(name = "limmaDiff",
  def = function(obj, dout, pat,
    MA.it = TRUE, HEAT.it = TRUE, GO.it = TRUE, DiffOut.it = TRUE,
    logFCthresh = 1, PValthresh = 0.05, log2.it = TRUE, small = 0.05) {
    standardGeneric("limmaDiff")
  }
)

#' @rdname limmaDiff-methods
setMethod(f = "limmaDiff",
  signature = c("tpm"),
  definition = function(obj, dout, pat, MA.it, HEAT.it, GO.it, DiffOut.it, logFCthresh, PValthresh, log2.it, small) {
    stopifnot(all(obj@tpm.value>=0))
    grps <- factor(obj@grps)
    if (log2.it) {
      tpm.value <- log2(obj@tpm.value + small)
    } else {
      tpm.value <- obj@tpm.value
    }
    contrast <- apply(combn(levels(grps), 2), 2,paste0,collapse="-")
    design <- model.matrix(~ 0 + grps)
    colnames(design) <- sub('^grps', '', colnames(design))
    fit <- lmFit(tpm.value, design)
    contrast.matrix <- makeContrasts(contrasts = contrast, levels = design)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)
    dat.list <- lapply(contrast,function(coef){
      dd <- topTable(fit, adjust.method = "BH", number = Inf, coef = coef, sort.by = 'none')
      dd <- dd[,c("logFC","P.Value","adj.P.Val")]
      query.population <- gsub(
        "[^\\|]+\\|([^\\|]+)\\|[^\\|]+(\\|[^\\|]+){0,1}",
        "\\1", rownames(dd))
      up.idx <- which(dd[,"logFC"] >= logFCthresh & dd[,"P.Value"] <= PValthresh)
      down.idx <- which(dd[,"logFC"] <= (-logFCthresh) & dd[,"P.Value"] <= PValthresh)
      if (length(up.idx) > 0 || length(down.idx) > 0) {
        dd$DEG <- "NDiff"
        dd[up.idx, "DEG"] <- "Up"
        dd[down.idx, "DEG"] <- "Down"
        if (MA.it) {
          MAchart(dd, pdffout = file.path(dout, paste0(pat, "_", coef, "_MA.pdf")))
        }
        if (HEAT.it) {
          sm1 <- gsub("([^\\-]+)\\-([^\\-]+)", "\\1", coef)
          sm2 <- gsub("([^\\-]+)\\-([^\\-]+)", "\\2", coef)
          diffHeatmap(tpm.value, col.idx = c(grep(sm1, as.character(grps)), grep(sm2, as.character(grps))),
            row.idx = which(dd$DEG != "NDiff"), pdffout = file.path(dout, paste0(pat, "_", coef, "_diffHeatmap.pdf")),
            cutreek = NULL, log.it.already = log2.it)
        }
        if (GO.it) {
          DEG <- data.frame(gene = query.population[c(up.idx, down.idx)],
            group = dd$DEG[c(up.idx, down.idx)], stringsAsFactors = FALSE)
          Res <- msigdb.gsea(DEG, query.population = query.population, background = 'query',
            genesets = c('C2.CP','C5.BP','C5.CC','C5.MF'), name.x='DEGs', name.go='MSigDB', species='mouse')
          write.table(Res, file = file.path(dout, paste0(pat, "_", coef, "_GO.xls")), row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")
        }
      };
      colnames(dd) <- paste(colnames(dd), coef, sep="_")
      return(dd)
    })
    dat <- do.call("cbind",dat.list)
    rownames(dat) <- rownames(tpm.value)
    gid <- gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}",
      "\\1", rownames(dat))
    gname <- gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}",
      "\\2", rownames(dat))
    gtype <- gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}",
      "\\3", rownames(dat))
    dat <- data.frame(gid = gid, gname = gname, gtype = gtype, dat, obj@tpm.value)
    if (DiffOut.it) {
      write.table(dat, file = file.path(dout, paste0(pat, "_DiffAna.xls")),
        row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
    return(dat)
  }
)

#' @title lineGraph
#' @name lineGraph
#' @rdname lineGraph-methods
#' @description
#' Line graphs of clustering
#' @param obj A \code{tpm} object.
#' @param col.idx A numeric or logical vector specifying columns to keep.
#' @param row.idx A numeric or logical vector specifying rows to keep.
#' @param clusters A numeric vector indicating cluster membership.
#' @param pdffout A character string specify pdf output file.
#' @param log.it.already A logical value specifying whether tpm values were already log2 transformed.
#' @param mean.it A logical value indicating whether to generate group values from individual samples.
#' @param small A numeric value indicating the adjustment to the TPM values before log2 transformation.
#' @export lineGraph
setGeneric(name = "lineGraph",
  def = function(obj, col.idx, row.idx, clusters, pdffout, log.it.already = FALSE, mean.it = FALSE, small = 0.05) {
    standardGeneric("lineGraph")
  }
)

#' @rdname lineGraph-methods
setMethod(f = "lineGraph",
  signature = c("tpm"),
  definition = function(obj, col.idx, row.idx, clusters, pdffout, log.it.already, mean.it, small) {
    tpm.value <- obj@tpm.value[row.idx, col.idx]
    if (!log.it.already) {
      tpm.value <- log2(tpm.value + small)
    }
    norm <- t(scale(t(tpm.value)))
    grps <- factor(obj@grps, levels = unique(obj@grps))
    norm.grp <- t(apply(norm, 1, function(vec)tapply(vec, grps, mean)))
    if (mean.it) {
      dat <- data.table(norm.grp, clusters = clusters)
      dat <- data.table(dat[, lapply(.SD, mean, na.rm = TRUE), by = clusters, .SDcols = levels(grps)], id="mean")
      setcolorder(dat, c(2:(ncol(dat)-1), 1, ncol(dat)))
    } else {
      dat <- data.table(norm.grp, clusters = clusters, id = rownames(norm.grp))
    }
    ldat <- melt(dat, id.vars = c(1:ncol(dat))[!1:ncol(dat) %in% 1:ncol(norm.grp)])
    ylab <- "TPM"
    p1 <- ggplot(data = ldat, aes_(x = ~variable, y = ~value, group = ~id, colour = "#FF9999" )) + geom_line() +
      geom_point() + facet_wrap(~ clusters) + labs(y = bquote(paste("Standardized ", log[2](.(ylab))))) + theme(legend.position="none")
    pdf(pdffout, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)
