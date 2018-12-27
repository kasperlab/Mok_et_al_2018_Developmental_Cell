# Load source script (RaceID3_StemID2_class.R), available on https://github.com/dgrun/RaceID3_StemID2
source("~/RaceID3_StemID2_class_GO.R")

# Load required packages
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
install.packages("XML")
install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit"))
install.packages("gProfileR")
install.packages("randomForest")
install.packages (c("vegan", "Rtsne", "scran", "DESeq2", "Go.db"))
install.packages("caTools")
library("biomaRt")

require(colorspace)

# Load .csv file (available on GEO (GSE122026))
prdata <- read.csv("~/GSE122026_E15.0_Mesenchymal_cells_scRNAseq_Final.csv", sep = "\t")

# Remove spike-in genes and mitochondrial genes
prdata <- prdata[grep("ERCC-",rownames(prdata),invert=TRUE),]
prdata <- prdata[grep("_chrM",rownames(prdata),invert=TRUE),]

# RaceID3 parameters
params <- TRUE

if ( params ){
  init          <- TRUE
  celltypeid    <- TRUE
  s             <- "Run"
  do.downsample <- F
  dsn           <- 1
  mintotal      <- 6000
  minexpr       <- 4
  minnumber     <- 2
  maxexpr       <- Inf
  metric        <- "manhattan" 
  cln           <- 0
  do.gap        <- FALSE
  sat           <- TRUE
  clustnr       <- 30
  bootnr        <- 50
  SE.method     <- "Tibs2001SEmax"
  SE.factor     <- .25
  B.gap         <- 50
  rseed         <- 17000
  rseed.tsne    <- 1337
  outminc       <- 25
  probthr       <- 1e-6
  outlg         <- 3
  outdistquant  <- .95
  qmark         <- .95
  genegroups    <- list(APO="APO",Collagen="COL",Keratin="KRT",Ribo="RPL|RPS",Mito="MTRNR",C1Q="C1Q(A|B)")
  genes         <- c("Gene name")
  FUNcluster    <- "hclust"
  old           <- FALSE
  comp          <- TRUE
  
  lineagetree   <- TRUE
  bootcells     <- FALSE
  cthr          <- 10
  
  cellproj      <- TRUE
  pdiback       <- TRUE
  pdishuf       <- 100
  pthr          <- .01
  
  monocle       <- FALSE
  
  lineagegenes  <- TRUE
  ksom          <- 50
  locreg        <- FALSE
  alpha         <- .25
  nbsom         <- 25
  minexprbr     <- 5
  minnumberbr   <- 2
  cellorder     <- "dist" # "som"  "tspc" "entr" "tspg" "dist"
  logscale      <- FALSE
  cty           <- list(ENT=c(3,2,1))
  
  lineagenetwork <- TRUE
  lgenes        <- c("Hmgn2","Atf4","Sox9","Fos","Jun")
  netk          <- 1
  netgenes      <- "TF"
  netcthr       <- .5
  corsig        <- FALSE
  
  markerid       <- FALSE
  stringent      <- FALSE
  mlist          <- list()
  nset           <- 2
  markermincount <- 10
  clmedian       <- TRUE
  minsize        <- 2
  wd 			 <- getwd()
  norm			 <- "n"
  if (do.downsample) {
    norm 		 <- "d"
  }	
}

cdiff <- list()
sco   <- list()
entr  <- list()

## RaceID3
## Create a directory where all the plots etc. will be saved
dir.create(paste(wd,"/RaceID3",s,"_",norm,mintotal,"_out",outlg,"_cln",cln,"minexp",minexpr,sep=""),showWarnings= FALSE)
setwd(paste(wd,"/RaceID3",s,"_",norm,mintotal,"_out",outlg,"_cln",cln,"minexp",minexpr,sep=""))

## Initialize SCseq object with transcript counts
sc <- SCseq(prdata)

## Filter expression data
sc <- filterdata(sc, mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=do.downsample, sfn=FALSE, hkn=FALSE, dsn=dsn, rseed=rseed, CGenes=NULL, FGenes=NULL, ccor=.4)

# Correct for cell cycle, proliferation
require(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
g   <- sub("__chr.+","",rownames(sc@fdata));
k   <- getBM(attributes = c("external_gene_name", "go_id","name_1006"),filters="external_gene_name",values=g,mart=mart)
gCC <- name2id( k$external_gene_name[k$name_1006 == "cell cycle"],rownames(sc@fdata)) 
gCP <- name2id( k$external_gene_name[k$name_1006 == "cell proliferation"],rownames(sc@fdata))
vset <- list(gCC,gCP)
x <- CCcorrect(sc@fdata,vset=vset,CGenes=NULL,ccor=.4,nComp=NULL,pvalue=.05,quant=.01,mode="pca")
# Number of principal components that have been removed
x$n
# Loadings of the first principal component that has been removed
y <- x$pca$rotation[,x$n[1]]
# Genes from vset are either enriched in the head or the tail of this list
tail(y[order(y,decreasing=TRUE)],10)
# Reassign the corrected expression matrix to sc@fdata
sc@fdata <- x$xcor

## K-means clustering
sc <- clustexp(sc,clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,sat=sat,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed,FUNcluster=FUNcluster,FSelect=FALSE)

## Compute t-SNE map
fast <- FALSE
sc <- comptsne(sc,rseed=15555,sammonmap=FALSE,initial_cmd=FALSE,fast=fast,perplexity=30)

# Detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=3,probthr=1e-7,thr=2**-(1:40),outdistquant=.95)

# Reassign clusters based on random forest
sc <- rfcorrect(sc,rfseed=12345,final=TRUE,nbfactor=5)

# Cell breakdown for final clusters
dir.create("clust_txt",showWarnings=F) 
x <- data.frame(CELLID=names(sc@cpart),cluster=sc@cpart)
write.table(x[order(x$cluster,decreasing=FALSE),],"clust_txt/cell_clust.xls",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

## Diagnostic plots

# Gap statistics
if ( do.gap ){
  fname <- paste("gap_clust_cells",s,sep="_")
  plot.2.file(fname)
  plotgap(sc)
  dev.off()
}

if ( length(unique(sc@cluster$kpart)) > 1 ){
  # Silhouette of k-means clusters
  fname <- paste("silhouette",s,sep="_")
  plot.2.file(fname)
  plotsilhouette(sc)
  dev.off()
  
  # Jaccard's similarity of k-means clusters
  fname <- paste("jaccard",s,sep="_")
  plot.2.file(fname)
  plotjaccard(sc)
  dev.off()
  
  if ( sat ){
    fname <- paste("saturation",s,sep="_")
    plot.2.file(fname)
    plotsaturation(sc)
    dev.off()
    
    fname <- paste("dispersion",s,sep="_")
    plot.2.file(fname)
    plotsaturation(sc,disp=TRUE)
    dev.off()
  }
}

# Barchart of outlier probabilities
fname <- paste("outlierprob",s,sep="_")
plot.2.file(fname)
plotoutlierprobs(sc)
dev.off()

# Regression of background model
fname <- paste("var_mean_cell_orig",s,sep="_")
plot.2.file(fname)
plotbackground(sc)
dev.off()

# Dependence of outlier number on probability threshold (probthr)
fname <- paste("par_sensitivity",s,sep="_")
plot.2.file(fname)
plotsensitivity(sc)
dev.off()

# Highlight final clusters in t-SNE map
fname <- paste("tsne_map_cells_c",s,"final_v2",sep="_")
plot.2.file(fname)
plottsne3(sc,final=TRUE,ax=F,text=T)
dev.off()

# Highlight k-means clusters in t-SNE map
fname <- paste("tsne_map_cells_c",s,"orig",sep="_")
plot.2.file(fname)
plottsne3(sc,final=FALSE,ax=TRUE)
dev.off()

# Highlight cell groups (experiments) in t-SNE map
fname <- paste("tsne_map_cells_types_symbol",s,sep="_")
plot.2.file(fname)
types <- sapply(names(sc@ndata),function(x) strsplit(x,"_")[[1]][2])
plotsymbolstsne3(sc,types=types,ax=F,seed=6) ## seed = 8 for dev_order figure ##
dev.off()

## Extensions
cdata   <- sc@ndata[apply(sc@ndata >= minexpr,1,sum) >= minnumber,]
part <- sc@cpart

set.seed(rseed)
minn <- min(apply(sc@expdata[,names(sc@cpart)],2,function(x) sum(round(x,0))))
dsd <- if ( do.downsample ) sc@ndata else downsample(sc@expdata[,names(sc@cpart)],n=minn,dsn=dsn)
probs <- t(t(dsd)/apply(dsd,2,sum))
entr[[s]] <- -apply(probs*log(probs),2,sum)

probs <- list()
for ( i in 1:max(part) ) probs[[i]] <- if ( sum( part == i ) == 1 ) dsd[,part == i]/sum(dsd[,part == i]) else  apply(dsd[,part == i],1,mean)/sum(apply(dsd[,part == i],1,mean))
ent <- c()
for ( p in probs ) ent <- append(ent,-sum(p[p>0]*log(p[p>0])))

fname <- paste("barplot_entropy",s,sep="_")
plot.2.file(fname)
plot(1:length(ent),ent,axes=FALSE,ylim=c(min(ent,na.rm=TRUE) - .1, max(ent,na.rm=TRUE) + .1),xlab="Cluster",ylab="Entropy",cex=0)
rect(1:length(ent) - .4,min(ent) - .1,1:length(ent) + .4,ent,col="grey")
axis(2)
axis(1,at=1:length(ent),cex.axis=.75)
dev.off()


tsd <- sc@tsne
co  <- sc@fcol

l <- list(entr=entr[[s]],number_t=apply(sc@expdata[,names(sc@cpart)],2,sum),number_g=apply(sc@expdata[,names(sc@cpart)],2,function(x) sum(x>1)),number_t_log=log2(apply(sc@expdata[,names(sc@cpart)],2,sum)))
nl <- list(entr="Entropy",number_t="Number of transcripts",number_g="Number of genes",number_t_log="log2 Number of transcripts")
for ( i in 1:length(l) ){
  mi <- min(l[[i]],na.rm=TRUE)
  ma <- max(l[[i]],na.rm=TRUE)
  if ( mi < ma ){
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l[[i]] - mi)/(ma - mi)*99 + 1,0)
    fname <- paste("tsne_map",names(l)[i],s,sep="_")
    plot.2.file(fname)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(3,5,2.5,2))
    plot(tsd,xlab="Dim 1",ylab="Dim 2",main=paste("tsne: ",s," (",nl[[i]],")",sep=""),pch=20,cex=0,col="grey")
    for ( k in 1:length(v) ){
      points(tsd[k,1],tsd[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
    }
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
    dev.off()
    
    ll <- list()
    pv <- c()
    qv <- c()
    ll[["all"]] <- l[[i]]
    for ( k in 1:max(part) ){
      ll[[paste("cl",k,sep=".")]] <- l[[i]][part == k]
      wi <- if (min(length(ll[["all"]]),length(ll[[paste("cl",k,sep=".")]])) > 5 ) wilcox.test(ll[["all"]],ll[[paste("cl",k,sep=".")]])$p.value else NA
      pv <- append(pv,wi)
      qv <- append(qv, quantile(ll[[paste("cl",k,sep=".")]],.95))
    }
    fname <- paste("boxplot_cell_clusters",names(l)[i],s,sep="_")
    plot.2.file(fname)
    boxplot(ll,xlab="Cluster",main=s,ylab=nl[[i]],ylim=c(min(l[[i]]),max(l[[i]]) + ( max(l[[i]]) - min(l[[i]]) )*.1),names=sub("cl\\.","",names(ll)),cex.axis=.5,cex=.5,pch=20)
    f <- !is.na(pv) & pv < .05 & pv >= 0.001
    if ( sum(f) > 0 ) text((2:(max(part) + 1))[f],max(l[[i]]) + ( max(l[[i]]) - min(l[[i]]) )*.05,rep("*",max(part))[f])
    f <- !is.na(pv) & pv < .001
    if ( sum(f) > 0 ) text((2:(max(part) + 1))[f],max(l[[i]]) + ( max(l[[i]]) - min(l[[i]]) )*.05,rep("**",max(part))[f])
    abline(a=median(l[[i]]),b=0,col="red",lty=2)
    dev.off()
  }
}

for ( i in 1:max(part) ){
  x  <- if ( sum(part != i) == 1 ) cdata[,part != i] else apply(cdata[,part != i],1,quantile,probs=qmark)
  #x  <- if ( sum(part != i) == 1 ) cdata[,part != i] else apply(cdata[,part != i],1,max)
  y  <- if ( sum(part == i) == 1 ) cdata[,part == i] else apply(cdata[,part == i],1,median)
  sd <- if ( sum(part != i) == 1 ) x else apply(cdata[,part != i],1,mad)
  z  <- (y - x)/sd
  n <- head(z[order(z,decreasing=TRUE)],10)
  n[n>100] <- 100
  names(n) <- sub("\\_\\_chr\\w+","",names(n))
  fname <- paste("cluster_marker_genes",s,"clust",i,sep="_")
  dir.create("marker",showWarnings= FALSE)
  plot.2.file(paste("marker",fname,sep="/"))
  b <- barplot(n,cex.names=.5,main=paste(s,": Cluster ",i,sep=""),ylab="z-score",ylim=c(min(0,n)-( max(n,0) - min(0,n) )/2,max(0,n)),names.arg=FALSE)
  text(b,min(0,n) - ( max(0,n) - min(0,n) )/3,names(n),srt=90)
  dev.off()
}

# Create tSNEs in separate folder
dir.create("exp", showWarnings = FALSE)
for ( g in genes ) {
  facs <- F
  ax <- FALSE
  # fname <- paste("tsne_map",s,g,"log","noaxis",sep="_")
  # plot.2.file(paste("exp",fname,sep="/"))
  # tryCatch(plotexptsne2(sc,g,logsc=TRUE,n=paste(g,"log",sep="_"),ax=ax,facs=facs),error=function(err) 0 )
  # dev.off()
  fname <- paste("tsne_map",s,g,"noaxis",sep="_")
  plot.2.file(paste("exp",fname,sep="/"))
  tryCatch(plotexptsne2(sc,g,ax=ax,facs=facs),error=function(err) 0)
  dev.off()
}

## StemID2

# Initialization
ltr <- Ltree(sc)

# Computation of entropy
ltr <- compentropy(ltr)

# Computation of the projections for all cells
ltr <- projcells(ltr,cthr=2,nmode=FALSE)

# Computation of the projections for all cells after randomization
ltr <- projback(ltr,pdishuf=2000,nmode=FALSE,fast=FALSE,rseed=17000)

# Assembly of the lineage tree
ltr <- lineagetree(ltr,pthr=0.05,nmode=FALSE,fast=FALSE)

# Determination of significant differentiation trajectories
ltr <- comppvalue(ltr,pethr=0.05,nmode=FALSE,fast=FALSE)

## Diagnostic plots
# Histogram of ratio between cell-to-cell distances in the embedded and the input space
fname <- paste("stemid_distanceratio",s,sep="_")
plot.2.file(fname)
plotdistanceratio(ltr)
dev.off()

# t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
fname <- paste("stemid_map",s,sep="_")
plot.2.file(fname)
plotmap(ltr)
dev.off()

# Visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
fname <- paste("stemid_mapprojections",s,sep="_")
plot.2.file(fname)
plotmapprojections(ltr)
dev.off()

# Lineage tree showing the projections of all cells in t-SNE space
fname <- paste("stemid_tree_projections",s,sep="_")
plot.2.file(fname)
plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=.3)
dev.off()

# lineage tree without showing the projections of all cells
fname <- paste("stemid_tree",s,sep="_")
plot.2.file(fname)
plottree(ltr,showCells=FALSE,nmode=FALSE,scthr=.3)
dev.off()

# Heatmap of the enrichment p-values for all inter-cluster links
fname <- paste("stemid_link_enrich",s,sep="_")
plot.2.file(fname)
plotlinkpv(ltr)
dev.off()

# Heatmap of the link score for all inter-cluster links
fname <- paste("stemid_link_score",s,sep="_")
plot.2.file(fname)
plotlinkscore(ltr)
dev.off()

# Heatmap showing the fold enrichment (or depletion) for significantly enriched or depleted links
fname <- paste("stemid_projenrichment",s,sep="_")
plot.2.file(fname)
projenrichment(ltr)
dev.off()

## Computing the StemID2 score
x <- compscore(ltr,nn=1,scthr=0)

# Plotting the StemID2 score
fname <- paste("stemid_score",s,sep="_")
plot.2.file(fname)
plotscore(ltr,nn=1,scthr=0)
dev.off()

# Retrieve cells from branch in pseudo-temporal order as inferred by the projection coordinates
n <- cellsfromtree(ltr,c(4,5,6,7,8))

# Filter out lowly expressed genes
fs  <- filterset(ltr@sc@ndata,n=n$f,minexpr=3,minnumber=1)
#Ccompute self organizing map (SOM) of co-expressed genes
s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
ps  <- procsom(s1d,corthr=.85,minsom=3)

# Coloring scheme for clusters (vector with colors)
fcol <- ltr@sc@fcol
y    <- ltr@sc@cpart[n$f]

# Plot average z-score for all modules derived from the SOM
fname <- paste("stemid_SOM_zscore_modules",s,sep="_")
plot.2.file(fname)
plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

# Plot z-score profile of each module
fname <- paste("stemid_SOM_zscore_genes",s,sep="_")
plot.2.file(fname)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

# Plot normalized expression profile of each module
fname <- paste("stemid_SOM_norm_genes",s,sep="_")
plot.2.file(fname)
plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

# Plot binarized expression profile of each module (z-score < -1, -1 < z-score < 1, z-score > 1)
fname <- paste("stemid_SOM_binary_genes",s,sep="_")
plot.2.file(fname)
plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

# Extract all genes from one module (2) of the SOM
stemnode<-2
g <- names(ps$nodes)[ps$nodes == stemnode]
# plot average expression profile of these genes along the trajectory
plotexpression(fs,y,g,n$f,k=25,col=fcol,name=paste("node",stemnode,sep=" "),cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)

# Plot expression of a single gene
plotexpression(fs,y,"Gene",n$f,k=25,col=fcol,cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)

# Plot average expression profile of these genes along the trajectory, highlighting batch origin
stemnode<-37
g <- names(ps$nodes)[ps$nodes == stemnode]
plotexpression(fs,y,g,n$f,k=25,col=fcol,name=paste("node",stemnode,sep=" "),cluster=FALSE,locreg=TRUE,alpha=.5,types=sub("\\_\\d+","",n$f))


