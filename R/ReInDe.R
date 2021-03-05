# The function ReInDe is used to remove strains that contain insertions, deletions, and degenerate bases in their sequences.
# Parameter fasta.file: Fasta file path.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
ReInDe<-function(fasta.file="",refseqID="NC_045512.2"){
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
  
  mRNAxl[which(mRNAxl=="r")]<-"-"
  mRNAxl[which(mRNAxl=="y")]<-"-"
  mRNAxl[which(mRNAxl=="m")]<-"-"
  mRNAxl[which(mRNAxl=="k")]<-"-"
  mRNAxl[which(mRNAxl=="s")]<-"-"
  mRNAxl[which(mRNAxl=="w")]<-"-"
  mRNAxl[which(mRNAxl=="h")]<-"-"
  mRNAxl[which(mRNAxl=="b")]<-"-"
  mRNAxl[which(mRNAxl=="v")]<-"-"
  mRNAxl[which(mRNAxl=="d")]<-"-"
  mRNAxl[which(mRNAxl=="n")]<-"-"
  
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  refcdsindex<-NULL
  for(i in 1:nrow(codematrix)){
    refcdsindex<-c(refcdsindex,c(codematrix[i,1]:codematrix[i,2]))
  }
  refcdsindex<-unique(refcdsindex)
  
  qs<-which(refseqindex==1)
  zz<-which(refseqindex==max(na.omit(refseqindex)))
  blindex<-c(qs:zz)
  blindex<-match(refcdsindex,refseqindex)
  qcr<-NULL
  qcc<-NULL
  for (i in 1:length(blindex)) {
    if(refseq[blindex[i]]=="-"){
      qcc<-c(qcc,blindex[i])
      qcr<-c(qcr,which(mRNAxl[,blindex[i]]!="-"))
    }else{
      qcr<-c(qcr,which(mRNAxl[,blindex[i]]=="-"))
    }
  }
  
  if(length(qcr)>0){
    mRNAxl<-mRNAxl[-qcr,]
    samples<-samples[-qcr]
  }
  
  jg<-NULL
  for(i in 1:length(samples)){
    a<-paste(">",samples[i],collapse = "")
    b<-paste(mRNAxl[i,],collapse = "")
    jg<-c(jg,a,b)
  }
 
  return(jg)
}