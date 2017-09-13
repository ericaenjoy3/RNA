#' @include RNAClass.R
#' @title read in Salmon differential expression file
#' @name SepTPMCnt
#' @rdname SepTPMCnt-methods
#' @export SepTPMCnt
SepTPMCnt <- function(fin) {
  dat <- fread(fin)
  tpm <- apply(dat[,-1,with <- F], 2, function(x)as.numeric(gsub("(.+);(.+)", "\\1", x)))
  cnt <- apply(dat[,-1,with <- F], 2, function(x)round(as.numeric(gsub("(.+);(.+)", "\\2", x)), digits <- 0))
  rpm <- t(t(cnt)/colSums(cnt) * 1e6)
  rownames(tpm) <- dat$gene
  rownames(cnt) <- dat$gene
  rownames(rpm) <- dat$gene
  grps <- gsub("(.+?)[_-]*[^-_]+$", "\\1", colnames(tpm))
  tpm.grp <- t(apply(tpm, 1,
    function(vec)
    tapply(vec, factor(grps, levels <- unique(grps)), mean)
  ))
  return(list(tpm <- tpm, cnt <- cnt, rpm <- rpm, tpm.grp <- tpm.grp))
}
