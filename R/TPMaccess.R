#!/usr/bin/env Rscript

#' Read in TPM from Salmon output and differential analysis
#'
#' @author Erica Liu \email{ericaenjoy3@@gmail.com}
#' import data.table
#' import limma
#' import affy
#' import ggplot2
#' import scatterplot3d
#' import RColorBrewer
#' import grDevices
#' import scales
#' import GGally
#' import GOtest
#' import ape
#' import pheatmap
#' import ComplexHeatmap
#' import circlize
#' import multiplot
#' import addgrids3d
###

#' @title read in Salmon differential expression file
#' @name SepTPMCnt
#' @rdname SepTPMCnt-methods
#' @export SepTPMCnt
SepTPMCnt=function(fin) {
  dat=fread(fin)
  tpm=apply(dat[,-1,with=F],2,function(x)as.numeric(gsub("(.+);(.+)","\\1",x)))
  cnt=apply(dat[,-1,with=F],2,function(x)round(as.numeric(gsub("(.+);(.+)","\\2",x)),digits=0))
  rpm=t(t(cnt)/colSums(cnt)*1e6)
  rownames(tpm)=dat$gene
  rownames(cnt)=dat$gene
  rownames(rpm)=dat$gene
  grps=gsub("(.+?)[_-]*[^-_]+$","\\1",colnames(tpm))
  tpm.grp=t(apply(tpm,1,function(vec)tapply(vec,factor(grps,levels=unique(grps)),mean)))
  return(list(tpm=tpm,cnt=cnt,rpm=rpm,tpm.grp=tpm.grp))
}

#' @title lmAdjCovar
#' @name lmAdjCovar
#' @rdname lmAdjCovar-methods
#' @export lmAdjCovar
lmAdjCovar=function(x,covar,add.mean=TRUE){
#function to correct for covariates using linear regression
#x, a matrix of gene expression data with genes in the rows and samples in the cols
#covar, a matrix/data.frame or vector of covariates
	covar=as.data.frame(covar)
	stopifnot(ncol(x)==nrow(covar))
	colnames(covar)=paste0('V',1:ncol(covar))
	x=t(x)
	t(lm(x~.,data=covar)$residuals)+if(add.mean){colMeans(x)}else{0}
}

#' @title loessnorm
#' @name loessnorm
#' @rdname loessnorm-methods
#' @export loessnorm
loessnorm=function(tpm,small=0.05){
  tpm.norm=normalize.loess(tpm+small,subset=grep("spikein", row.names(tpm)))
  return(tpm.norm)
}

#' @title rdcntnorm
#' @name rdcntnorm
#' @rdname rdcntnorm-methods
#' @export rdcntnorm
rdcntnorm=function(tpm,stats) {
  frac=stats$spikein/(stats$spikein+stats$mm10)
  normfac=mean(frac)/frac
  tpm.norm=data.frame(t(t(tpm)*normfac))
  return(tpm.norm)
}

#' @title distplot
#' @name distplot
#' @rdname distplot-methods
#' @description
#' tpm.value either normalised to spikein or not
#' tpm.value either spikein, non-spikein, or all
#' idx for selecting differnetially expressed genes
#' @export distplot
distplot=function(tpm.value,grps,ylab,pngfout) {
  dat=data.table(tpm.value)
  ldat=melt(dat)
  max.y=quantile(ldat$value, probs=0.85)
  p1=ggplot(ldat, aes(x=variable, y=value))+geom_boxplot(aes(fill = factor(variable)))+ scale_y_continuous(labels=comma)+labs(x="",y=ylab)+theme(legend.title=element_blank(), legend.position="top")+coord_cartesian(ylim=c(0,max.y))
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  theme_set(theme_grey(base_size=15))
  multiplot(p1,cols=1)
  dev.off()
}

#' @title corplot
#' @name corplot
#' @rdname corplot-methods
#' @export corplot
corplot=function(tpm.value,grps,pngfout) {
  tpm.value=tpm.value
  thresh=floor(quantile(as.matrix(tpm.value),probs=0.999))
  pm=ggpairs(tpm.value)
  pm2=pm
  for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- pm[i,j]+coord_cartesian(xlim=c(0,thresh),ylim=c(0,thresh))+scale_x_continuous(breaks=c(0,as.numeric(thresh)))+scale_y_continuous(breaks=c(0,as.numeric(thresh)))
    }
  }
  png(pngfout,width=3000,height=3000,res=300)
  print(pm2)
  dev.off()
}

#' @title hireplot
#' @name hireplot
#' @rdname hireplot-methods
#' @export hireplot
hireplot=function(tpm.value,grps,pngfout) {
  corstats=cor(tpm.value,method="spearman")
  if(!is.factor(grps)) {grps=factor(grps,levels=unique(grps),ordered=T)}
  tip.col=brewer.pal(length(levels(grps)),"Dark2")[as.numeric(grps)]
  png(pngfout,width=3000,height=3000,res=300)
  plot(as.phylo(hclust(as.dist(1-corstats),method='average')),cex=2,label.offset=0,tip.color=tip.col)
  dev.off()
}

#' @title heatcorplot
#' @name heatcorplot
#' @rdname heatcorplot-methods
#' @export heatcorplot
heatcorplot=function(tpm.value,grps,pngfout) {
  # corstats=cor(tpm.value[idx,],method="spearman")
  # ana.col=data.frame(Group=grps)
  # rownames(ana.col)=colnames(tpm.value)
  # png(pngfout,width=3000,height=3000,res=300)
  # pheatmap(corstats,colorRampPalette(rev(brewer.pal(n=7, name ="RdYlBu")))(100),breaks=seq(-1,1,length.out=101),annotation_col=ana.col)
  # dev.off()
  pm=ggcorr(tpm.value,method=c("pairwise","spearman"),label=TRUE,label_alpha=TRUE)
  png(pngfout,width=3000,height=3000,res=300)
  print(pm)
  dev.off()
}

#' @title bplot
#' @name bplot
#' @rdname bplot-methods
#' @export bplot
bplot=function(tpm.value,grps,title,pngfout,maxPcnt=0.80,ylab=expression(paste(log[2], "(TPM)")),isLog=FALSE,small=0.05) {
  if(!isLog) dat=log2(tpm.value+small)
  map.it=grps
  names(map.it)=colnames(tpm.value)
  ldat=melt(data.table(tpm.value))
  ldat$grps=factor(ldat$variable,levels=colnames(tpm.value))
  levels(ldat$grps)=grps
  min.y=min(ldat$value)
  max.y=as.numeric(quantile(ldat$value,probs=maxPcnt))
  p1=ggplot(ldat, aes(x=variable, y=value, fill=grps))+
  geom_boxplot()+
  coord_cartesian(ylim=c(min.y, max.y))+
  theme(legend.title=element_blank(), legend.position="top")+
  labs(x="",y=ylab)+
  ggtitle(title)
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  theme_set(theme_grey(base_size=15))
  multiplot(p1,cols=1)
  dev.off()
}

#' @title PCAplot
#' @name PCAplot
#' @rdname PCAplot-methods
#' @export PCAplot
PCAplot=function(tpm.value,grps,pngfout,fout=NULL,excl.col=NULL,ntop=Inf,isLog=FALSE,small=0.05) {
  tpm.value=tpm.value
  out=function(dat,fout) {
    tmp=data.frame(gene=rownames(dat),dat)
    write.table(tmp,file=fout,row.names=F,col.names=T,quote=F,sep="\t")
  }
  if (is.null(excl.col)) {
    dat=as.matrix(tpm.value)
  } else {
    dat=as.matrix(tpm.value[,-excl.col])
  }
  if(ncol(dat)<=3) stop('The number of samples must be larger than 3\n')
  ridx=apply(dat,1,function(vec) any(vec>0))
  stopifnot(any(ridx))
  dat=dat[ridx,]
  if(!isLog) dat=log2(dat+small)
  dat=dat[apply(dat,1,var)>0,]
  if(ntop<nrow(dat)){
     vars=apply(dat,1,var)
     dat=dat[order(vars,decreasing=T),]
     dat=dat[1:ntop,]
  }
  pca=prcomp(t(dat), center=TRUE, scale=TRUE, retx=TRUE)
  ve=pca$sdev^2/sum(pca$sdev^2) #variance explained by each PC
  ve13=sprintf("%.1f",ve[1:3]*100)
  x=data.frame(PC1=pca$x[,1],PC2=pca$x[,2],PC3=pca$x[,3])
  rownames(x)=colnames(dat)
  if (!is.infinite(ntop)) {out(dat,fout)}
  uniq.cols=rainbow(length(unique(grps)))
  cols=uniq.cols[as.numeric(factor(grps,levels=unique(grps),ordered=TRUE))]
  x$colors=cols
  pch=as.numeric(factor(rownames(x)),ordered=TRUE)
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  par(mar=c(1,1,1,1))
  s3d=scatterplot3d(x[,1:3], grid=FALSE, box=FALSE, mar=c(3,3,2,2),pch="")
  addgrids3d::addgrids3d(x[, 1:3], grid = c("xy", "xz", "yz"))
  #s3d$points3d(x[, 1:3], pch=16, col=x$colors)
  text(s3d$xyz.convert(x[, 1:3]+1),labels=1:nrow(x),col=x$colors)
  legend(par('usr')[1]-0.5,par('usr')[4]+0.3, legend=paste0(1:nrow(x),": ",rownames(x)), pch=20,col=x$colors,xpd=TRUE, ncol=3,cex=0.7)
  dev.off()
}

#' @title MAplot
#' @name MAplot
#' @rdname MAplot-methods
#' @description
#' plot of log2FC over pvalue during differential analysis
#' @export MAplot
MAplot=function(dd, pngfout) {
  nd=dd
  nd$nlogpval=-log10(nd$P.Value)
  ylab="P"
  p1=ggplot(nd, aes(x=logFC, y=nlogpval, color=DEG))+geom_point()+labs(y=bquote(-log[10](.(ylab))))+theme(legend.title=element_blank(), legend.position="top")
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  theme_set(theme_grey(base_size=15))
  multiplot(p1,cols=1)
  dev.off()
}

#' @title BICplot
#' @name BICplot
#' @rdname BICplot-methods
#' @export BICplot
BICplot=function(g,BIC,pngfout) {
  dd=data.frame(g=g,BIC=BIC)
  p1=ggplot(dd, aes(x=g, y=BIC)) + geom_point(colour="#FF9999")+xlab("Cluster Size")+ylab("BIC")
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  theme_set(theme_grey(base_size=15))
  multiplot(p1,cols=1)
  dev.off()
}

#' @title diffHeatmap
#' @name diffHeatmap
#' @rdname diffHeatmap-methods
#' @description
#' heatmap during or after differential analysis
#' @export diffHeatmap
diffHeatmap=function(tpm.value,col.idx,row.idx,pngfout,cutreek=NULL,cut.alg=c("pam","hclust","emmix"),rank.man=FALSE,log.it.already=FALSE,scale.it=TRUE,cluster_columns_par=TRUE,cluster_rows_par=TRUE,show_row_dend_par=FALSE,small=0.05,...) {
  tpm.value=tpm.value[row.idx,col.idx]
  if (!log.it.already) {tpm.value=log2(tpm.value+small)}
  if (scale.it) {norm=t(scale(t(tpm.value)))} else {norm=tpm.value}
  if (rank.man) {
    norm=norm[order(norm[,1],-norm[,2],decreasing=T),];
    cluster_rows=FALSE
  }
  if (!is.null(cutreek)) {
    if (cut.alg=="hclust") {
      corstats=cor(t(norm),method="spearman");
      clusters=cutree(hclust(as.dist(1-corstats),method='average'),k=cutreek)
    } else if (cut.alg=="emmix") {
      if (!is.null(clusters)) {clusters=clusters}
    } else if (cut.alg=="pam") {
      clusters=pam(norm,cutreek)$clustering
    }
    pm=ComplexHeatmap::Heatmap(norm,cluster_columns=cluster_columns_par,cluster_rows=cluster_rows_par,show_row_names=FALSE,show_row_dend=show_row_dend_par,show_column_dend=show_column_dend_par,heatmap_legend_param=list(title="",color_bar="continuous"),clustering_distance_rows="spearman",clustering_method_rows="average",clustering_distance_columns="spearman",clustering_method_columns="average",split=paste0("Cluster",clusters))
  } else {
    pm=ComplexHeatmap::Heatmap(norm,cluster_columns=cluster_columns_par,cluster_rows=cluster_rows_par,show_row_names=FALSE,show_row_dend=show_row_dend_par,show_column_dend=show_column_dend_par,heatmap_legend_param=list(title="",color_bar="continuous"),clustering_distance_rows="spearman",clustering_method_rows="average",clustering_distance_columns="spearman",clustering_method_columns="average")
  }
  png(pngfout,width=3000,height=4500,res=300,pointsize=14)
  draw(pm)
  dev.off()
  if (exists("clusters")) {
    return(clusters)
  } else {
    return(NULL)
  }
}

#' @title limmaDiff
#' @name limmaDiff
#' @rdname limmaDiff-methods
#' @description
#' tpm.value either normalised to spikein or not
#' filtering of genes must be done before hand
#' @export limmaDiff
limmaDiff=function(tpm.value.ori,grps,dout,pat,MA.it=TRUE,HEAT.it=TRUE,GO.it=TRUE,DiffOut.it=TRUE,logFCthresh=1,PValtrhesh=0.05,log2.it=T,small=0.05) {
  stopifnot(all(tpm.value.ori>=0))
  if (!is.factor(grps)) {grps=factor(grps)}
  if (log2.it) {
    tpm.value=log2(tpm.value.ori+small)
  } else {
    tpm.value=tpm.value.ori
  }
  contrast=apply(combn(levels(grps),2),2,paste0,collapse="-")
  design <- model.matrix(~0+grps)
  colnames(design) <- sub('^grps','',colnames(design))
  fit <- lmFit(tpm.value,design)
  contrast.matrix <- makeContrasts(contrasts=contrast,levels=design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  dat.list=lapply(contrast,function(coef){dd=topTable(fit,adjust="BH",number=Inf,coef=coef,sort='none');
  dd=dd[,c("logFC","P.Value","adj.P.Val")]
  query.population=gsub("[^\\|]+\\|([^\\|]+)\\|[^\\|]+(\\|[^\\|]+){0,1}","\\1",rownames(dd));
  up.idx=which(dd[,"logFC"]>=logFCthresh & dd[,"P.Value"]<=PValtrhesh);
  down.idx=which(dd[,"logFC"]<=(-logFCthresh) & dd[,"P.Value"]<=PValtrhesh);
  if (length(up.idx)>0 || length(down.idx)>0) {
    dd$DEG="NDiff"
    dd[up.idx,"DEG"]="Up";
    dd[down.idx,"DEG"]="Down";
    if (MA.it) {MAplot(dd,pngfout=file.path(dout,paste0(pat,"_",coef,"_MA.png")))}
    if (HEAT.it) {
      sm1=gsub("([^\\-]+)\\-([^\\-]+)","\\1",coef)
      sm2=gsub("([^\\-]+)\\-([^\\-]+)","\\2",coef)
      diffHeatmap(tpm.value,col.idx=c(grep(sm1,as.character(grps)),grep(sm2,as.character(grps))),row.idx=which(dd$DEG!="NDiff"),pngfout=file.path(dout,paste0(pat,"_",coef,"_diffHeatmap.png")),cutreek=NULL,log.it.already=log2.it)
    }
    if (GO.it) {
      DEG=data.frame(gene=query.population[c(up.idx,down.idx)],group=dd$DEG[c(up.idx,down.idx)],stringsAsFactors=FALSE)
      Res=msigdb.gsea(DEG,query.population=query.population, background='query',genesets=c('C2.CP','C5.BP','C5.CC','C5.MF'),name.x='DEGs',name.go='MSigDB',species='mouse')
      write.table(Res,file=file.path(dout,paste0(pat,"_",coef,"_GO.xls")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
    }
  };
  colnames(dd)=paste(colnames(dd),coef,sep="_");
  return(dd)})
  dat=do.call("cbind",dat.list)
  rownames(dat)=rownames(tpm.value)
  gid=gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}","\\1",rownames(dat));
  gname=gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}","\\2",rownames(dat));
  gtype=gsub("([^\\|]+)\\|([^\\|]+)\\|([^\\|]+)(\\|[^\\|]+){0,1}","\\3",rownames(dat));
  dat=data.frame(gid=gid,gname=gname,gtype=gtype,dat,tpm.value.ori)
  if (DiffOut.it) {
    write.table(dat,file=file.path(dout,paste0(pat,"_DiffAna.xls")),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  }
  return(dat)
}

#' @title lineGraph
#' @name lineGraph
#' @rdname lineGraph-methods
#' @description
#' line graphs of clustering
#' @export lineGraph
lineGraph=function(tpm.value,col.idx,row.idx,clusters,pngfout,log.it.already=FALSE,mean.it=FALSE,small=0.05) {
  tpm.value=tpm.value[row.idx,col.idx]
  if (!log.it.already) {tpm.value=log2(tpm.value+small)}
  norm=t(scale(t(tpm.value)))
  grps=gsub("_*\\d*","",colnames(norm))
  norm.grp=t(apply(norm,1,function(vec)tapply(vec,factor(grps,levels=c("M","EG","LG")),mean)))
  if (mean.it) {
    dat=data.table(norm.grp,clusters=clusters)
    dat=data.table(dat[,lapply(.SD,mean,na.rm=TRUE),by=clusters,.SDcols=c("M", "EG","LG")],id="mean")
    setcolorder(dat,c(2:(ncol(dat)-1),1,ncol(dat)))
  } else {
    dat=data.table(norm.grp,clusters=clusters,id=rownames(norm.grp))
  }
  ldat=melt(dat,id.vars=c(1:ncol(dat))[!1:ncol(dat) %in% 1:ncol(norm.grp)])
  ylab="TPM"
  p1=ggplot(data=ldat, aes(x=variable, y=value, group=id,colour="#FF9999"))+geom_line()+geom_point()+facet_wrap(~ clusters)+labs(y=bquote(paste("Standardized ",log[2](.(ylab)))))+theme(legend.position="none")
  png(pngfout,width=3000,height=3000,res=300,pointsize=14)
  theme_set(theme_grey(base_size=15))
  multiplot(p1,cols=1)
  dev.off()
}
