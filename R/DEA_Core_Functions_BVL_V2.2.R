#' Core Differential Expression analysis Function
#'
#'@param Design.file File with design matrix SampleID / Condition / Batch
#'@param counts.dir Directory with counts per gene for each sample(one file per sample)
#'@param DEComparison.file File with pairwise comparisons to be performed
#'@param name.prefix a name to be appended to generated files
#'@param OutputDir output directory 
#'@param anot.file annotation file full path 
#'@param source.format annotation format G
#'@param SelFeature
#'@param fields
#'@param adjpval
#'@param cleandescription T/F (default= T) change ascii code in decrition to char (comon in gff/gtf from http://tritrypdb.org/)
#'@param correct.batch
#'@param correct.batch.inmodel
#'@param filter.lcg filter low count genes
#'@return PCA/ HeatMap /Non Zero GEnes
#'@return annotated DE genes results table
#'@return BED Track with DE genes
#'@export

DEAnalysis=function(Design.file,counts.dir,DEComparison.file,name.prefix,OutputDir,anot.file,source.format,SelFeature,fields,adjpval,cleandescription,correct.batch,correct.batch.inmodel,filter.lcg)
{
  Design=read.table(Design.file,header=T)
  DEComparisons=read.table(DEComparison.file)
  nprefix=basename(Design.file)
  conds=Design$C
  Batch=Design$B
  print("Read the first count file to fetch genes names") 
  setwd(counts.dir)
  tmp=read.delim(as.character(paste(Design[1,1],".count",sep="")),header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE,comment.char = "#")
  
  
  print("consolidate CountsTable from individual count files")
  CountsTable=matrix(nrow=dim(tmp)[1],ncol=dim(Design)[1])
  rownames(CountsTable)=tmp[,1]
  colnames(CountsTable)=Design[,1]
  CountsTable[,1]=tmp[,2]
  for(i in 2:dim(Design)[1])
  {
    #print(as.character(Design[i,1]))
    CountsTable[,i]=read.delim(as.character(paste(Design[i,1],".count",sep="")),header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE,comment.char = "#",)[,2]   
  }
  #print(head(CountsTable))
  # write count table into a file
  dir.create(paste("../",OutputDir,sep=""), showWarnings = FALSE)
  setwd(paste("../",OutputDir,sep=""))
  
  ## add option writeCountTable = T/F
  ## if (writeCountTable = T) {}
  print("# write count table into a file ")
  write.table(cbind(GeneID=rownames(CountsTable),CountsTable),paste(name.prefix,"_CountsTable.xls",sep=""),sep="\t",row.names = F)
  
  print("# find number of expressed gene / total reads per condition")

  NZG_TR=data.frame(SampleID=colnames(CountsTable),NonZeroGenes=colSums(CountsTable>1),Count=colSums(CountsTable)*1e-6,Condition=Design$C,Batch=Design$B)
  write.table(NZG_TR,file=paste(name.prefix,"_","NZG_CPM_data.xls",sep=""),sep="\t",row.names=F)
  
  library(ggplot2)
  mycpmp=ggplot(NZG_TR, aes(Count,NonZeroGenes, color=Condition, shape=Batch))+ geom_point(stat="identity",size=5) + theme(axis.text.x=element_text(angle=-90)) + geom_text(aes(label=SampleID), size=3,vjust=-1)
  
  ggsave(mycpmp, file=paste(name.prefix,"_","NZG_CPM_plot.png",sep=""), dpi=600)
  
  
  if (filter.lcg) {
    print("Filter low count genes")
    ntgenes=dim(CountsTable)[1]
    CountsTable=CountsTable[rowSums(CountsTable) > dim(CountsTable)[2],]
    nfgenes=dim(CountsTable)[1]
    print(paste("Total Number Of Genes: ",ntgenes," / Number Of Genes After Filtering: ",nfgenes, " / Percentage Of Genes after filtering: ", round(100*nfgenes/ntgenes,digits=2),"%"))
    ### add number of gene  before and after filtering and percentage 
    
  }
  print("Plot Library size (millions)")
  library(ggplot2)
  tmp=cbind(Design[,-3],Count=colSums(CountsTable)*1e-6)
  mycpmp=ggplot(tmp, aes(HPGL_ID,Count, fill=Condition))+ geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=-90))
  ggsave(mycpmp, file=paste(name.prefix,"_","Count_plot.png",sep=""), dpi=600)
  
  print("# generate HeatMap & PCA raw data")
  require.auto("cbcbSEQ")
  plot.cluster2(CountsTable,Design,name.prefix) 
  print("generate table with correlation of (Condition,Batch) to Principal component")
 
  
  print("# generate HeatMap & PCA normalized log2cpm data") 
  qcounts = qNorm(CountsTable)
  # convert to log(counts per million)
  res = log2CPM(qcounts)
  y = res$y
  plot.cluster2(y,Design,paste(name.prefix,"_Norm_log2CPM_Heatmap.png",sep=""))
  
  print("# generate PCA to Condition/Batch Correlation table row data ") 
  cSVD=makeSVD(CountsTable)
  cpcRes=pcRes(cSVD$v,cSVD$d,Design$Condition,Design$Batch)
  write.table(cpcRes,paste(name.prefix,"_PCcorCondBatch_raw_Table.xls",sep=""),sep="\t",row.names = F)
  print("# generate PCA to Condition/Batch Correlation table normalized log2CPM data ") 
  cSVD=makeSVD(y)
  cpcRes=pcRes(cSVD$v,cSVD$d,Design$Condition,Design$Batch)
  write.table(cpcRes,paste(name.prefix,"_PCcorCondBatch_log2CPM_Table.xls",sep=""),sep="\t",row.names = F)
  if (correct.batch)
  {
    print("Perform Batch Correction")  
    CountsTable2=batchSEQ(CountsTable,model.matrix(~0+Design$C),Design$B,Design$C)
    CountsTable2 = CountsTable2$elist
    colnames(CountsTable2$design)=sapply(strsplit(colnames(CountsTable2$design),split="\\$"), "[[", 2)   
    colnames(CountsTable2$design)=substr(colnames(CountsTable2$design),2,nchar(colnames(CountsTable2$design)))
    plot.cluster2(CountsTable2$E,Design,paste(name.prefix,"Batch_corrected"))
    cSVD=makeSVD(CountsTable2)
    cpcRes=pcRes(cSVD$v,cSVD$d,Design$Condition,Design$Batch)
    write.table(cpcRes,paste(name.prefix,"_PCcorCondBatch_BatchCorrected_Table.xls",sep=""),sep="\t",row.names = F)
  } else {
    CountsTable2= voom(CountsTable, model.matrix(~0+Design$C), plot=TRUE)
    colnames(CountsTable2$design)=sapply(strsplit(colnames(CountsTable2$design),split="\\$"), "[[", 2)   
    colnames(CountsTable2$design)=substr(colnames(CountsTable2$design),2,nchar(colnames(CountsTable2$design)))
  }
  
  
  if (correct.batch.inmodel)
  {
  print("correct for batch in limma model")  
  Design$Batch<-as.character(Design$Batch)
  mod=model.matrix(~0+Design$C+Design$B)  #limma
  colnames(mod)=substr(colnames(mod),9,nchar(colnames(mod)))
  if (correct.batch){CountsTable2$design= mod}else{CountsTable2= voom(CountsTable, mod, plot=TRUE)}
  print("create Contrast matrix from DEComparison file")
  kont=paste("makeContrasts(",paste("ctr",1,"=", DEComparisons[1,1], "-", DEComparisons[1,2],sep=""))
  for (i in 2:nrow(DEComparisons))
  {
    kont=paste(kont,paste("ctr",i,"=", DEComparisons[i,1], "-", DEComparisons[i,2],sep=""),sep=",")
  }
  kont=paste(kont,", levels=CountsTable2$design)")
  kont
  assign("matC",eval(parse(text=kont)))
  print("Contrast matrix created")
  print(matC)
  
  fit0 = lmFit(CountsTable2)
  fit = contrasts.fit(fit0, matC)  
  eb = eBayes(fit)
  
  for (i in 1:dim(DEComparisons)[1])
  {
    cond1=as.character(DEComparisons[i,1])
    cond2=as.character(DEComparisons[i,2])
    print(paste("performing DE analysis :",cond1,"Vs.",cond2,sep=" "))
    top1 = topTable(eb, coef=paste("ctr",i,sep=""), n=nrow(CountsTable2$E),sort.by="M")
    
    print("MAPlot generation")
    madata=data.frame(AvgExp=rowMeans(CountsTable2$E[top1$ID,]),LogFC=top1$logFC,AdjPVal=top1$adj.P.Val)
    mymap=ggplot(madata, aes(AvgExp,LogFC, color=(AdjPVal<adjpval))) +
                 geom_hline(yintercept=c(-1,1),color="Red",size=2) +
                 geom_point(stat="identity",size=3) +
                 theme(axis.text.x=element_text(angle=-90)) +
                 xlab("Average Count(Million Reads)") +ylab("log FC")
                 
                 ggsave(mymap, file=paste(name.prefix,"_",cond1,"_", cond2,"_MA_plot.png",sep=""), dpi=600)
                 
    print("Plot P-value Distribution")
    png(paste(name.prefix,"_",cond1,"_", cond2,"_pval_plot.png",sep=""))
    hist(top1$P)
    dev.off()
    FC=2^top1$logFC
    FC_hr=FC
    FC_hr[which(FC_hr<1)]=-1/FC_hr[which(FC_hr<1)]  
    E1=2 * top1$AveExpr / (1+FC)
    E2=FC * E1
    
    # cont1= A - B <=> B ir ref  
    res=cbind(id=top1$ID,E1,E2,AveExpr=top1$AveExpr,FC,FC_hr,top1[,c("logFC","P.Value","adj.P.Val")])
    names(res)[match("E1", names(res))] = cond2
    names(res)[match("E2", names(res))] = cond1
    print("Annotate full table")
    resAnot=annotate.data(res,anot.file,source.format,SelFeature,fields,cleandescription)
    write.table(resAnot,paste(name.prefix,"_",cond1,"_", cond2,"_res_anot.xls",sep=""),sep="\t",row.names = F)
    print("Annotate resSig table")
    resSig=res[ res$adj.P.Val < adjpval, ]
    if (dim(resSig)[1]>0)
    {
      print(dim(resSig))
      resSigAanot=annotate.data(resSig,anot.file,source.format,SelFeature,fields,cleandescription)
      write.table(resSigAanot,paste(name.prefix,"_",cond1,"_", cond2,"_res_sig_anot.xls",sep=""),sep="\t",row.names = F)
      print("generate Bed Track")
      sv2bed(resSigAanot,cond1,cond2,nprefix)  
    }
    
  } 
}
  else
  {
    print("NO correction for batch in limma model") 
    #limma
    print("create Contrast matrix from DEComparison file")
    kont=paste("makeContrasts(",paste("ctr",1,"=", DEComparisons[1,1], "-", DEComparisons[1,2],sep=""))
    for (i in 2:nrow(DEComparisons))
    {
      kont=paste(kont,paste("ctr",i,"=", DEComparisons[i,1], "-", DEComparisons[i,2],sep=""),sep=",")
    }
    kont=paste(kont,", levels=unique(unlist(Design$C)))")
    kont
    assign("matC",eval(parse(text=kont)))
    print("Contrast matrix created")
    print(matC)
    
    fit0 = lmFit(CountsTable2)
    fit = contrasts.fit(fit0, matC)  
    eb = eBayes(fit)
    
    for (i in 1:dim(DEComparisons)[1])
    {
      cond1=as.character(DEComparisons[i,1])
      cond2=as.character(DEComparisons[i,2])
      print(paste("performing DE analysis :",cond1,"Vs.",cond2,sep=" "))
      top1 = topTable(eb, coef=paste("ctr",i,sep=""), n=nrow(CountsTable2$E),sort.by="M")
      print("MAPlot generation")
      madata=data.frame(AvgExp=rowMeans(CountsTable2$E[top1$ID,]),LogFC=top1$logFC,AdjPVal=top1$adj.P.Val)
      mymap=ggplot(madata, aes(AvgExp,LogFC, color=(AdjPVal<adjpval))) +
        geom_point(stat="identity",size=3) +
        geom_hline(yintercept=c(-1,1),color="Red",size=2) +
        theme(axis.text.x=element_text(angle=-90)) +
        xlab("Average Count(Million Reads)") +ylab("log FC")
      print(mymap)
      ggsave(mymap, file=paste(name.prefix,"_",cond1,"_", cond2,"_MA_plot.png",sep=""), dpi=600)
      
      print("Plot P-value Distribution")
      png(paste(name.prefix,"_",cond1,"_", cond2,"_pval_plot.png",sep=""))
      hist(top1$P)
      dev.off()
      FC=2^top1$logFC
      FC_hr=FC
      FC_hr[which(FC_hr<1)]=-1/FC_hr[which(FC_hr<1)]  
      E1=2 * top1$AveExpr / (1+FC)
      E2=FC * E1
      
      # cont1= A - B <=> B ir ref  
      res=cbind(id=top1$ID,E1,E2,AveExpr=top1$AveExpr,FC,FC_hr,top1[,c("logFC","P.Value","adj.P.Val")])
      names(res)[match("E1", names(res))] = cond2
      names(res)[match("E2", names(res))] = cond1
      print("Annotate full table")
      resAnot=annotate.data(res,anot.file,source.format,SelFeature,fields,cleandescription)
      write.table(resAnot,paste(name.prefix,"_",cond1,"_", cond2,"_res_anot.xls",sep=""),sep="\t",row.names = F)
      print("Annotate resSig table")
      resSig=res[ res$adj.P.Val < adjpval, ]
      if (dim(resSig)[1]>0)
      {
        print(dim(resSig))
        resSigAanot=annotate.data(resSig,anot.file,source.format,SelFeature,fields,cleandescription)
        write.table(resSigAanot,paste(name.prefix,"_",cond1,"_", cond2,"_res_sig_anot.xls",sep=""),sep="\t",row.names = F)
        print("generate Bed Track")
        sv2bed(resSigAanot,cond1,cond2,nprefix)  
      }
      
    } 
  }
 
  
  
}



#' Annotate Differentially expressed gene table 
#'
#'@param datTable DE results table 
#'@param anot.file
#'@param source.format
#'@param SelFeature 
#'@param fields
#'@param cleandescription  
#'@return annotated DE genes results table
#'@export

annotate.data = function (datTable, anot.file,source.format="gff",SelFeature="gene",fields=c("ID","Name","description"),cleandescription) 
{
  if (source.format=="gff")
  {
    require.auto("genomeIntervals")  
    gff <- readGff3(anot.file,isRightOpen=FALSE)
    anotattribs=unique(as.character(annotation(gff)$type))
    if ((!SelFeature%in%anotattribs) & length(intersect(fields,anotattribs) > 0 )) { stop("SelFeature not found or none of the element of fields found in gff column 9 !") } 
    else {
      tmp=gff[annotation(gff)$type==SelFeature]
      myanot=as.data.frame(cbind(as.character(annotation(tmp)$seq_name),tmp[,1:2],getGffAttribute(tmp,fields)))
      colnames(myanot)[1:3]=c("chr","start","end")
      if (cleandescription) {
        for (i in fields) {myanot=ascii2char(myanot,i)}
      }
      #apply( myanot$description, 1, function(x) gsub("%2C",",",x) )
      #myanot=getGffAttribute(gff[annotation(gff)$type==SelFeature],fields) 
      return(resanot=cbind(datTable,myanot[match(tolower(datTable$id),tolower(myanot$ID)),]))
    }  
  }  
  # anotate from txt tab delimited txt file 
  # key column name in annotation file must be "ID"
  if (source.format=="txt")
  {
    myanot=read.delim(anot.file,header=T)
    return(resanot=cbind(datTable,myanot[match(tolower(datTable$id),tolower(myanot$ID)),])) 
  }
  
  
}   





#' convert results from DE table to bed track for visualization in IGV  
#'
#'@param res DE results table 
#'@param Cond1
#'@param Cond2
#'@param Design.file name of the design file 
#'@param lfct log fold-change threshold 
#'@return Bed file with DE genes results
#'@export

sv2bed=function(res,cond1,cond2,Design.file,lfct=1)
{
  nprefix=basename(Design.file)
  bedcolor=rep("0,0,0",dim(res)[1])
  bedcolor[which(res[,6]>lfct)]="102,204,0" #UP
  bedcolor[which(res[,6]>0 & res[,6]<lfct)]="204,255,153" #UP
  bedcolor[which(res[,6]>-lfct & res[,6]<=0)]="255,153,153" #DOWN
  bedcolor[which(res[,6]<=-lfct)]="255,0,0" #DOWN
  
  score=rep(0,dim(res)[1])
  blank=rep("",dim(res)[1])
  svplotd=cbind( res[,c(10,11,12)],
                 info=paste(as.character(cond1),round(res[,3], digits = 0), as.character(cond2),round(res[,4], digits = 0),"LFC",round(res[,6], digits = 2),"_AdjP",round(res[,8], digits = 2),sep="_"),
                 score,blank,res[,c(11,12)],
                 bedcolor)
  
  
  
  svplotd <- data.frame(lapply(svplotd, as.character), stringsAsFactors=FALSE)
  svplotd=svplotd[order(svplotd[,1],svplotd[,2],svplotd[,3]),]
  write.table(svplotd,paste(name.prefix,"_",as.character(cond1),"_", as.character(cond2),".bed",sep=""),col.names=F,row.names=F,quote = F,sep="\t")
}


#' Plot heatmap and 2D scatter plot of principal components 
#'
#'@param CountsTable count table 
#'@param Design design ( sample id / Condition / Batch )
#'@param name.prefix a name for the generated plot files 
#'@return Heatmap plot file adn PCA plot file in png format
#'@export


######
## Plot2 on count + design
######
plot.cluster2 = function(CountsTable,Design,name.prefix="sample"){
  
  require.auto("gplots") 
  require.auto("ggplot2") 
  require.auto("RColorBrewer")  
  
  # HeatMap data
  print("# HeatMap Plot")
  if (length(as.integer(as.factor(Design$Batch))) >= 2) 
    rc=brewer.pal(12,"Set3")[as.integer(as.factor(Design$Batch))] else rc = rep("green",length(Design$Batch))
  
  if (length(as.integer( as.factor(Design$Condition))) >= 2) 
    cc=brewer.pal(9,"Set1")[as.integer( as.factor(Design$Condition))] else cc =  rep("red",length(Design$Batch))
  
  
  dists = dist( t( CountsTable ) )
  mat = as.matrix( dists )
  #hc = hclust(dist(t(CountsTable)))
  
  png(file=paste(name.prefix,"_Heatmap.png",sep=""))
  par(mar = c(1,1,1,1))
  hv <- heatmap.2(mat, 
                  RowSideColors=rc, ColSideColors=cc, margin=c(6, 6), 
                  #main=name.prefix,
                  trace="none",key=F,
                  #="Conditions",labRow="Batch"
  )
  legend(x="topleft", legend=unique(Design$Condition), col=unique(cc), pch=15)
  legend(x="topright", legend=unique(Design$Batch), col=unique(rc), pch=15)
  dev.off()
  
  print("# PCA plot")
  pca=makeSVD(CountsTable)
  pcVar <- round((pca$d^2)/sum(pca$d^2)*100,2)
  xl <- sprintf("PC1: %.2f%% variance", pcVar[1])  
  yl <- sprintf("PC2: %.2f%% variance", pcVar[2]) 
  pcaData=data.frame(SampleID=colnames(CountsTable),PC1=pca$v[,1],PC2=pca$v[,2],Condition=Design$C,Batch=Design$B)
  mypcap=ggplot(pcaData, aes(PC1,PC2, color=Condition, shape= Batch)) + geom_point(stat="identity",size=5) + geom_text(aes(label=SampleID),angle=45, size=4,vjust=2) + xlab(xl) +ylab(yl) +theme(axis.ticks = element_blank(),axis.text.x=element_text(angle=-90))
  print(mypcap)
  ggsave(mypcap, file=paste(name.prefix,"_pca_plot.png",sep=""), dpi=600)
  
}




#'clean up description by converting ascii code to char in the annotation
#'
#'@param myanot a dataframe containing results
#'@param field name of the field that require conversion
#'@return a dataframe with cleaned field
#'@export

ascii2char=function(myanot,field)
{
  print("Clean-up description")
  
  myanot[,which(colnames(myanot)==field)]=gsub("%21","!",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%23","#",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%24","$",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%25","%",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%26","&",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%28","(",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%29",")",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%2A","*",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%2C",",",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%2D","-",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%2F","/",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3A",":",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3B",";",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3C","<",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3D","=",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3E",">",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("%3F","?",myanot[,which(colnames(myanot)==field)],fixed = TRUE)
  myanot[,which(colnames(myanot)==field)]=gsub("+"," ",myanot[,which(colnames(myanot)==field)],fixed = TRUE)  

  return(myanot)
}



#' auto install of packages from Cran / Bioconductor and their dependecies (if they are not installed)
#' Adapted from http://sbamin.com/2012/11/05/tips-for-working-in-r-automatically-install-missing-package/
#' 
#' @param x package to be installed 
#' 
require.auto=function(x)
{  
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
  eval(parse(text=paste("require(", x, ")", sep="")))
} else {
  update.packages(ask=F) # update dependencies, if any.
  eval(parse(text=paste("install.packages('", x, "')", sep="")))
}
if(isTRUE(x %in% .packages(all.available=TRUE))) {
  eval(parse(text=paste("require(", x, ")", sep="")))
} else {
  source("http://bioconductor.org/biocLite.R")
  biocLite(character(), ask=FALSE) # update dependencies, if any.
  eval(parse(text=paste("biocLite('", x, "')", sep="")))
  eval(parse(text=paste("require(", x, ")", sep="")))
}
}