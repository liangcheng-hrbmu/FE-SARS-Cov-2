# The function HieClust can re-encode the nucleic acid sequence of the SARS-CoV-2 and perform hierarchical clustering.
# Parameter fasta.file: Fasta file path.
# Parameter sample.class: We built sample.class in folder TestData where the sample information of the continent classification is stored.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
# Parameter plot: It control whether the output of HieClust is a hierarchical clustering result or drawing a hierarchical clustering diagram.
# Dependent package: dendextend, factoextra, flexclust, RColorBrewer.
HieClust<-function(fasta.file="",sample.class,refseqID="NC_045512.2",plot=F){
  genes<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  
  test<-readLines(fasta.file)
  wz<-which(grepl(">",test)==T)
  samples<-test[wz]
  test<-tolower(test)
  
  xl<-NULL
  for(i in 1:length(wz)){
    if(i<length(wz)){
      a<-wz[i]
      b<-wz[i+1]
      xl<-c(xl,paste(test[(a+1):(b-1)],collapse=""))
    }else{
      a<-wz[i]
      b<-length(test)
      xl<-c(xl,paste(test[(a+1):b],collapse=""))
    }
  }
  
  mRNAxl<-list()
  for(i in 1:length(xl)){
    a<-unlist(strsplit(xl[i],split=""))
    mRNAxl[[i]]<-a
  }
  
  samples<-gsub(">","",samples)
  samples<-gsub(" ","",samples)
  pp<-which(samples==refseqID)
  refseq<-mRNAxl[[pp]] 
  mRNAxl<-mRNAxl[-pp]
  samples<-samples[-pp]
  refseqindex<-NULL 
  j=0
  for (i in 1:length(refseq)) {
    if(refseq[i]=="-"){
      refseqindex<-c(refseqindex,NA)
    }else{
      j<-j+1
      refseqindex<-c(refseqindex,j)
    }
  }
  mRNAxl<-do.call("rbind",mRNAxl)
  
  mRNAxl[mRNAxl=="a"]<-1
  mRNAxl[mRNAxl=="t"]<-2
  mRNAxl[mRNAxl=="c"]<-3
  mRNAxl[mRNAxl=="g"]<-4
  mRNAxl<-apply(mRNAxl, 2, as.numeric)
  
  codeindex<-sapply(c(1:nrow(codematrix)), function(i){
    return(c(codematrix[i,1]:codematrix[i,2]))
  })
  codeindex<-unlist(codeindex)

  pp<-match(codeindex,refseqindex)
  mRNAxl_cds<-mRNAxl[,pp]
  row.names(mRNAxl_cds)<-samples
  
  result <- dist(mRNAxl_cds, method = "euclidean")
  result_hc <- hclust(d = result, method = "ward.D2")
  
  if(plot==T){
    dend <- as.dendrogram(result_hc)
    col_car_type<-rep(1,nrow(mRNAxl_cds))
    for(i in 1:length(sample.class)){
      col_car_type[na.omit(match(sample.class[[i]],row.names(mRNAxl_cds)))]<-c("yellow","red","blue")[i]
    }
  
    par(mar = c(12,4,1,1))
    labels(dend) <- ""
    plot(dend)
    colored_bars(col_car_type,dend,rowLabels="")
  }else{
    return(result_hc)
  }
}