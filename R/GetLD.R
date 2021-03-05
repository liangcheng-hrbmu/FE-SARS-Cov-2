# The function GetLD can obtain the input required by software Haploview according to the base sequence position.
# Parameter fasta.file: Fasta file path.
# Parameter seq.index: A numeric vector is the position of the mutation site in the reference sequence where you want to calculate the linkage disequilibrium.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
GetLD<-function(fasta.file="",seq.index,refseqID="NC_045512.2"){
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
  
  
  refseqindex<-NULL
  j=0
  for (i in 1:length(refseq)){
    if(refseq[i]=="-"){
      refseqindex<-c(refseqindex,NA)
    }else{
      j<-j+1
      refseqindex<-c(refseqindex,j)
    }
  }
  seq.index<-match(seq.index,refseqindex)
  mRNAxl<-do.call("rbind",mRNAxl)
  mRNAxl<-mRNAxl[,seq.index]
  
  mRNAxl[which(is.element(mRNAxl,c("a","t","c","g"))==F)]<-0
  mRNAxl[which(mRNAxl=="a")]<-1
  mRNAxl[which(mRNAxl=="c")]<-2
  mRNAxl[which(mRNAxl=="g")]<-3
  mRNAxl[which(mRNAxl=="t")]<-4
  
  for(j in 1:ncol(mRNAxl)){
    jj_table<-table(mRNAxl[,j])
    if(length(jj_table)>2){
      jj_table_px<-sort(jj_table,decreasing = T)
      jj_table_px<-jj_table_px[-c(1:2)]
      for(y in 1:length(jj_table_px)){
        mRNAxl[which(mRNAxl[,j]==names(jj_table_px)[y]),j]<-0
      }
    }
  }
  
  xl<-NULL
  for(i in 1:ncol(mRNAxl)){
    a<-cbind(mRNAxl[,i],mRNAxl[,i])
    xl<-cbind(xl,a)
  }
  seq1<-cbind(1:nrow(mRNAxl),1:nrow(mRNAxl),rep(0,nrow(mRNAxl)),rep(0,nrow(mRNAxl)),rep(1,nrow(mRNAxl)),
              rep(2,nrow(mRNAxl)),xl)
  
  pp<-refseqindex[seq.index]
  SNP_index<-cbind(pp,pp)
  jg<-list(samples=seq1,SNP_sites=SNP_index)
  return(jg)
}