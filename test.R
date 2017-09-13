#!/usr/bin/env Rscript

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
  grps=gsub("_\\d","",colnames(tpm)))
gene.obj <- new("tpm",
  tpm.value = sepSpike(tpm.obj, invert = FALSE),
  grps = tpm.obj@grps)
spike.obj <- new("tpm",
  tpm.value = sepSpike(tpm.obj, invert = TRUE),
  grps = tpm.obj@grps)

distplot(spike.obj, grps, ylab="TPM",
  pdffout = file.path(dout, "Daf_spikein_distplot.pdf"))
distplot(gene.obj, grps, ylab="TPM",
  pdffout = file.path(dout, "Daf_genes_distplot.pdf"))
bplot(spike.obj, grps, title="Spike-in",
  pdffout = file.path(dout, "Daf_spikein_bplot.pdf"),
  maxPcnt = 0.80, ylab=expression(paste(log[2], "(TPM)")),
  isLog = FALSE, small = 0.05)
bplot(gene.obj, grps, title = "Genes",
  pdffout = file.path(dout,"Daf_genes_bplot.pdf"), maxPcnt = 0.80,
  ylab = expression(paste(log[2], "(TPM)")), isLog = FALSE, small = 0.05)

distplot(gene.obj, grps, ylab = "TPM",
  pdffout=file.path(dout, "Daf_distplot.pdf"))
corplot(gene.obj, grps,
  pdffout = file.path(dout, "Daf_corplot.pdf"))
hireplot(gene.obj, grps,
  pdffout = file.path(dout, "Daf_hireplot.pdf"))
heatcorplot(gene.obj, grps,
  pdffout = file.path(dout, "Daf_heatcorplot.pdf"))
PCAplot(gene.obj, grps,
  pdffout = file.path(dout,"Daf_PCAplot.png"),
  fout = NULL, excl.col = NULL, ntop = Inf, isLog = FALSE, small = 0.05)

genefilter.obj <- rmLow(gene.obj, thresh = 1)
genevar.obj <- rmNonVar(genefilter.obj, probs = 0.1)
dat <- limmaDiff(genefilter.obj, dout, pat = "Daf",
  MA.it = TRUE, HEAT.it = TRUE, GO.it = TRUE, DiffOut.it = TRUE,
  logFCthresh = 1, PValtrhesh = 10^(-5), log2.it = TRUE, small = 0.05)

# differential expressed genes and also in genevar.obj -------------------------

# differential expressed genes between ESC and MEF------------------------------
mat = structure(list(dat[ridx, cidx]), row.names=paste(dat$gid, dat$gname, dat$gtype, sep="|")[ridx], class="data.frame")
kHeat(mat, pdffout="Daf_kHeat.pdf")
