#!/usr/bin/env Rscript

# RNA Library
library(RNA)

# Other library
libs <- c("data.table", "limma", "affy", "ggplot2", "scatterplot3d",
  "RColorBrewer", "grDevices", "scales", "GGally", "GOtest", "ape",
  "pheatmap", "ComplexHeatmap", "circlize", "multiplot", "addgrids3d")
sapply(libs, library, character.only = TRUE, quietly = TRUE)

fin <- file.path(
  "/home/liuyiyua/athena/RNA/RNA_seq/DF5154_2017_08_25/hera",
  "TPM_rd_merge.txt"
)
dout <- dirname(fin)
tpm.value <- SepTPMCnt(fin)$tpm

tpm.obj <- new("tpm",
  tpm.value = tpm.value,
  grps = gsub("_\\d", "", colnames(tpm.value)))
gene.obj <- new("tpm",
  tpm.value = sepSpike(tpm.obj, invert = TRUE),
  grps = tpm.obj@grps)
spike.obj <- new("tpm",
  tpm.value = sepSpike(tpm.obj, invert = FALSE),
  grps = tpm.obj@grps)

distplot(spike.obj, ylab="TPM",
  pdffout = file.path(dout, "Daf_spikein_distplot.pdf"))
distplot(gene.obj, ylab="TPM",
  pdffout = file.path(dout, "Daf_genes_distplot.pdf"))
bplot(spike.obj, title="Spike-in",
  pdffout = file.path(dout, "Daf_spikein_bplot.pdf"),
  probs = 0.80, ylab=expression(paste(log[2], "(TPM)")),
  isLog = FALSE, small = 0.05)
bplot(gene.obj, title = "Genes",
  pdffout = file.path(dout,"Daf_genes_bplot.pdf"), probs = 0.80,
  ylab = expression(paste(log[2], "(TPM)")),
  isLog = FALSE, small = 0.05)

distplot(gene.obj, ylab = "TPM",
  pdffout=file.path(dout, "Daf_distplot.pdf"))
corplot(gene.obj,
  pdffout = file.path(dout, "Daf_corplot.pdf"))
hireplot(gene.obj,
  pdffout = file.path(dout, "Daf_hireplot.pdf"))
heatcorplot(gene.obj,
  pdffout = file.path(dout, "Daf_heatcorplot.pdf"))
PCAplot(gene.obj,
  pdffout = file.path(dout,"Daf_PCAplot.png"),
  fout = NULL, excl.col = NULL, ntop = Inf, isLog = FALSE, small = 0.05)

genefilter.obj <- rmLow(gene.obj, thresh = 1)
genevar.obj <- rmNonVar(genefilter.obj, probs = 0.1)
dat <- limmaDiff(genefilter.obj, dout, pat = "Daf",
  MA.it = TRUE, HEAT.it = TRUE, GO.it = TRUE, DiffOut.it = TRUE,
  logFCthresh = 1, PValthresh = 10^(-5), log2.it = TRUE, small = 0.05)

# differential expressed genes and also in genevar.obj -------------------------
ridx <- apply(dat[,grep("DEG", colnames(dat))], 1 ,
  function(vec)any(vec != 'NDiff'))
com.obj <- new("tpm",
  tpm.value = genefilter.obj@tpm.value[rownames(dat)[ridx] %in% rownames(genevar.obj@tpm.value), ],
  grps = genefilter.obj@grps)
ridx <- dat[, "DEG_ESC.MEF"] != 'NDiff' & apply(dat[, grep("^(ESC|MEF)", colnames(dat))], 1, function(vec)any(vec>10))
esc.obj <- new("tpm",
  tpm.value = genefilter.obj@tpm.value[ridx, ],
  grps = genefilter.obj@grps)


# differential expressed genes between ESC and MEF------------------------------
kHeat(com.obj, pdffout = file.path(dout, "Daf_comDiff_kHeat.pdf"), k = 9)
kHeat(esc.obj, pdffout = file.path(dout, "Daf_escDiff_k4Heat.pdf"), k = 4)
