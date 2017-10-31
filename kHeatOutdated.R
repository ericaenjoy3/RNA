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
    if (is.null(k)) {
        message("clusing function when k is null.")
        pr_mb <- clusing(data.frame(mat), pdffout = gsub(".pdf","_optK.pdf", pdffout))
        message("Heatmap function when k is null.")
        ht_list <- Heatmap(mat, show_row_names = TRUE, show_column_names = TRUE, cluster_rows = TRUE,
          show_row_dend = FALSE,  cluster_columns = FALSE, show_column_dend = FALSE,
          heatmap_legend_param = list(title = "", color_bar = "continuous"),
          clustering_distance_rows = "spearman", clustering_method_rows = "average",
          clustering_distance_columns = "spearman", clustering_method_columns = "average",
          split = factor(pr_mb), gap = unit(3, "mm"))
    } else {
        message("heatmap function when k is ", k)
        set.seed(888)
        ht_list <- Heatmap(mat, show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE,
          show_row_dend = FALSE,  cluster_columns = FALSE, show_column_dend = FALSE,
          heatmap_legend_param = list(title = "", color_bar = "continuous"),
          clustering_distance_rows = "spearman", clustering_method_rows = "average",
          clustering_distance_columns = "spearman", clustering_method_columns = "average",
          km = k, gap = unit(3, "mm"))
    }
    pdf(pdffout)
    draw(ht_list)
    dev.off()
    message("running prepare function...")
    ht.obj <- prepare(ht_list)
    message("done with prepare function...")
    idx.dat <- do.call("rbind",
      lapply(seq_along(ht.obj@row_order_list),function(i){
        data.frame(cluster = i, idx = ht.obj@row_order_list[[i]])
      })
    )
    idx.dat <- idx.dat[order(idx.dat[, "idx"]), ]
    dat <- data.frame(gene = rownames(mat), idx.dat, mat = mat)
    write.table(dat,
      file = gsub("pdf", "txt", pdffout),
      row.names = FALSE, col.names = TRUE,
      sep = "\t", quote = FALSE)
  }
)
