#' @include RNAClass.R
#' @title rmSpike
#' @name rmSpike
#' @rdname rmSpike-methods
#' @export rmSpike
setGeneric(name = 'sepSpike',
  def = function(obj, invert = FALSE) {
    standardGeneric("sepSpike")
  }
)

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
#' @export rmLow
setGeneric(name = "rmLow",
  def = function(obj, thresh) {
    standardGeneric("rmLow")
  }
)

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
#' @export rmNonVar
setGeneric(name = "rmNonVar",
  def = function(obj, probs = 0.1) {
    standardGeneric("rmNonVar")
  }
)

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
#' @export lmAdjCovar
setGeneric(name = "lmAdjCovar",
  def = function(x, covar, add.mean = TRUE) {
    standardGeneric("lmAdjCovar")
  }
)

setMethod(f = "lmAdjCovar",
  signature = c("numeric", "data.frame", "logical"),
  definition = function(x, covar, add.mean){
    #function to correct for covariates using linear regression
    #x, a matrix of gene expression data with genes in the rows and samples in the cols
    #covar, a matrix/data.frame or vector of covariates
  	stopifnot(ncol(x) == nrow(covar))
  	colnames(covar) <- paste0('V', 1:ncol(covar))
  	x <- t(x)
  	t(lm(x~., data = covar)$residuals) +
      if(add.mean){colMeans(x)} else {0}
  }
)

#' @title loessnorm
#' @name loessnorm
#' @rdname loessnorm-methods
#' @export loessnorm
setGeneric(name="loessnorm",
  def = function(obj, small = 0.05) {
    standardGeneric("loessnorm")
  }
)

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
#' @export rdcntnorm
setGeneric(name="rdcntnorm",
  def = function(obj, stats) {
    standardGeneric("rdcntnorm")
  }
)

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
#' idx for selecting differnetially expressed genes
#' @export distplot
setGeneric(name="distplot",
  def = function(obj, ylab, pdffout, probs = 0.85) {
    standardGeneric("distplot")
  }
)

setMethod(f = "distplot",
  signature = c("tpm", "character", "character"),
  definition = function(obj, ylab, pdffout, probs) {
    dat <- data.table(obj@tpm.value)
    ldat=melt(dat)
    max.y=quantile(ldat$value, probs = probs)
    p1=ggplot(ldat, aes(x = variable, y = value)) +
      geom_boxplot(aes(fill = factor(variable))) +
      scale_y_continuous(labels = comma) +
      labs(x = "", y = ylab) +
      theme(legend.title = element_blank(), legend.position="top") +
      coord_cartesian(ylim=c(0,max.y))
      pdf(pdffout)
      theme_set(theme_grey(base_size=15))
    multiplot(p1,cols=1)
    dev.off()
  }
)

#' @title corplot
#' @name corplot
#' @rdname corplot-methods
#' @export corplot
setGeneric(name = "corplot",
  def = function(obj, pdffout) {
    standardGeneric("corplot")
  }
)

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
    pdf(pdffout)
    print(pm2)
    dev.off()
  }
)

#' @title hireplot
#' @name hireplot
#' @rdname hireplot-methods
#' @export hireplot
setGeneric(name = "hireplot",
  def = function(obj, pdffout) {
    standardGeneric("hireplot")
  }
)

setMethod(f = "hireplot",
  signature = c("tpm"),
  definition = function(obj, pdffout) {
    corstats <- cor(obj@tpm.value,method="spearman")
    grps <- factor(obj@grps,levels=unique(obj@grps),ordered=T)
    tip.col <- brewer.pal(length(levels(grps)), "Dark2")[as.numeric(grps)]
    pdf(pdffout)
    plot(as.phylo(hclust(as.dist(1-corstats), method = 'average')),
      cex = 2, label.offset = 0, tip.color = tip.col)
    dev.off()
  }
)

#' @title heatcorplot
#' @name heatcorplot
#' @rdname heatcorplot-methods
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
#' @export bplot
setGeneric(name = "bplot",
  def = function(obj, title, pdffout, maxPcnt = 0.80,
    ylab = expression(paste(log[2], "(TPM)")),
    isLog = FALSE, small = 0.05) {
    standardGeneric("bplot")
  }
)

setMethod(f = "bplot",
  signature = "tpm",
  definition = function(obj, title, pdffout, maxPcnt, ylab, isLog, small) {
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
    max.y <- as.numeric(quantile(ldat$value, probs = maxPcnt))
    p1 <- ggplot(ldat, aes(x = variable, y = value, fill = ldat$grps))+
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
#' @export PCAplot
setGeneric(name = "PCAplot",
  def = function(obj, pdffout, fout = NULL, excl.col = NULL,
    ntop = Inf, isLog = FALSE, small = 0.05) {
    standardGeneric("PCAplot")
  }
)

setMethod(f = "PCAplot",
  signature = "tpm",
  definition = function(obj, pdffout, fout, excl.col, ntop, isLog, small) {
    out <- function(dat,fout) {
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
    if (!is.infinite(ntop)) {
      out(dat, fout)
    }  
    uniq.cols <- rainbow(length(unique(obj@grps)))
    cols <- uniq.cols[as.numeric(factor(obj@grps, levels = unique(grps), ordered = TRUE))]
    x$colors <- cols
    pch <- as.numeric(factor(rownames(x)), ordered = TRUE)
    pdf(pdffout, pointsize = 14)
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

#' @title MAplot
#' @name MAplot
#' @rdname MAplot-methods
#' @description
#' plot of log2FC over pvalue during differential analysis
#' @export MAplot
setGeneric(name = "MAplot",
  def = function(dd, pdffout) {
    standardGeneric("MAplot")
  }
)

setMethod(f = "MAplot",
  signature = c("data.frame", "character"),
  definition = function(dd, pdffout) {
    nd <- dd
    nd$nlogpval <- -log10(nd$P.Value)
    ylab <- "P"
    p1 <- ggplot(nd, aes(x = logFC, y = nlogpval, color = DEG)) + geom_point() +
      labs(y = bquote(-log[10](.(ylab)))) +
      theme(legend.title = element_blank(), legend.position = "top")
    png(pngfout, width = 3000, height = 3000, res = 300, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)

#' @title BICplot
#' @name BICplot
#' @rdname BICplot-methods
#' @export BICplot
setGeneric(name = "BICplot",
  def = function(g, BIC, pdffout) {
    standardGeneric("BICplot")
  }
)

setMethod(f = "BICplot",
  signature = c("numeric", "numeric", "character"),
  definition = function(g, BIC, pdffout) {
    dd <- data.frame(g = g ,BIC = BIC)
    p1 <- ggplot(dd, aes(x = g, y = BIC)) + geom_point(colour = "#FF9999") +
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
#' @export diffHeatmap
setGeneric(name = "diffHeatmap",
  def = function(tpm.value, col.idx, row.idx, pdffout,
    cutreek = NULL, cut.alg, rank.man = FALSE, log.it.already = FALSE,
    scale.it = TRUE, cluster_columns_par = TRUE, cluster_rows_par = TRUE,
    show_row_dend_par = FALSE, small = 0.05, ...) {
    standardGeneric("diffHeatmap")
  }
)

setMethod(f = "diffHeatmap",
  signature = c("matrix"),
  definition = function(tpm.value, col.idx, row.idx, pdffout, cutreek,
    cut.alg, rank.man, log.it.already, scale.it, cluster_columns_par,
    cluster_rows_par, show_row_dend_par, small, ...) {
    cut.alg <- match.arg(cut.alg, c("pam","hclust","emmix"))
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
      norm <- norm[order(norm[,1],-norm[,2],decreasing=T),];
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

#' @title kHeat
setGeneric(name = "kHeat",
  def = function(obj, pdffout, log2.it = TRUE, scale.it = TRUE, small = 0.05) {
    standardGeneric("kHeat")
  }
)

setMethod(f = "kHeat",
  signature = c("tpm", "ANY"),
  definition = function(obj, pdffout, log2.it, scale.it, small) {
    mat <- obj@tpm.value
    if (log2.it) {
      mat <- log2(mat + small)
    }
    if (scale.it) {
      mat <- t(apply(mat, 1, scale))
    }
    pr_mb <- clusing(mat, pdffout = gsub(".pdf","_optK.pdf", pdffout))
    ht_list <- Heatmap(mat, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE,
      show_row_dend = FALSE,  cluster_columns = FALSE, show_column_dend = FALSE,
      heatmap_legend_param = list(title = "", color_bar = "continuous"),
      clustering_distance_rows = "spearman", clustering_method_rows = "average",
      clustering_distance_columns = "spearman", clustering_method_columns = "average",
      split=factor(pr_mb), gap = unit(3, "mm"))
    # png(pngfout,width=2500*2,height=2500,res=300)
    pdf(pdffout)
    draw(ht_list)
    dev.off()
    write.table(data.frame(mat, cluster = pr_mb),
      file = gsub("pdf", "txt", pdffout),
      row.names = TRUE, col.names = TRUE,
      sep = "\t", quote = FALSE)
  }
)

#' @title limmaDiff
#' @name limmaDiff
#' @rdname limmaDiff-methods
#' @description
#' tpm.value either normalised to spikein or not
#' filtering of genes must be done before hand
#' @export limmaDiff
setGeneric(name = "limmaDiff",
  def = function(obj, dout, pat,
    MA.it = TRUE, HEAT.it = TRUE, GO.it = TRUE, DiffOut.it = TRUE,
    logFCthresh = 1, PValtrhesh = 0.05, log2.it = TRUE, small = 0.05) {
    standardGeneric("limmaDiff")
  }
)

setMethod(f = "limmaDiff",
  signature = c("tpm"),
  definition = function(obj, dout, pat, MA.it, HEAT.it, GO.it, DiffOut.it, logFCthresh, PValtrhesh, log2.it, small) {
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
      dd <- topTable(fit, adjust = "BH", number = Inf, coef = coef, sort = 'none')
      dd <- dd[,c("logFC","P.Value","adj.P.Val")]
      query.population <- gsub(
        "[^\\|]+\\|([^\\|]+)\\|[^\\|]+(\\|[^\\|]+){0,1}",
        "\\1", rownames(dd))
      up.idx <- which(dd[,"logFC"] >= logFCthresh & dd[,"P.Value"] <= PValtrhesh)
      down.idx <- which(dd[,"logFC"] <= (-logFCthresh) & dd[,"P.Value"] <= PValtrhesh)
      if (length(up.idx) > 0 || length(down.idx) > 0) {
        dd$DEG <- "NDiff"
        dd[up.idx, "DEG"] <- "Up"
        dd[down.idx, "DEG"] <- "Down"
        if (MA.it) {
          MAplot(dd, pdffout = file.path(dout, paste0(pat, "_", coef, "_MA.pdf")))
        }
        if (HEAT.it) {
          sm1 <- gsub("([^\\-]+)\\-([^\\-]+)", "\\1", coef)
          sm2 <- gsub("([^\\-]+)\\-([^\\-]+)", "\\2", coef)
          diffHeatmap(tpm.value, col.idx = c(grep(sm1, as.character(grps)), grep(sm2, as.character(grps))),
            row.idx = which(dd$DEG != "NDiff"), pngfout = file.path(dout, paste0(pat, "_", coef, "_diffHeatmap.png")),
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
#' line graphs of clustering
#' @export lineGraph
setGeneric(name = "lineGraph",
  def = function(obj, col.idx, row.idx, clusters, pdffout, log.it.already = FALSE, mean.it = FALSE, small = 0.05) {
    standardGeneric("lineGraph")
  }
)

setMethod(f = "lineGraph",
  signature = c("obj"),
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
      dat <- data.table(dat[, lapply(.SD, mean, na.rm = TRUE), by = clusters, .SDcols = levels(grps), id="mean")
      setcolorder(dat, c(2:(ncol(dat)-1), 1, ncol(dat)))
    } else {
      dat <- data.table(norm.grp, clusters=clusters, id=rownames(norm.grp))
    }
    ldat <- melt(dat, id.vars = c(1:ncol(dat))[!1:ncol(dat) %in% 1:ncol(norm.grp)])
    ylab <- "TPM"
    p1 <- ggplot(data = ldat, aes(x = variable, y = value, group = id, colour = "#FF9999" )) + geom_line() +
      geom_point() + facet_wrap(~ clusters) + labs(y = bquote(paste("Standardized ", log[2](.(ylab))))) + theme(legend.position="none")
    pdf(pdffout, pointsize = 14)
    theme_set(theme_grey(base_size = 15))
    multiplot(p1, cols = 1)
    dev.off()
  }
)
