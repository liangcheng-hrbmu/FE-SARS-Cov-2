# The function SeqDivision can classify strains according to the time of diagnosis and the continent of patients infected with SARS-CoV-2.
# Parameter fasta.file: Fasta file path.
# Parameter divide: Sort sequence data according to time or continent.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
SeqDivision<-function(fasta.file="",divide="time",refseqID="NC_045512.2"){
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
  
  
  if(divide=="time"){
    samplecf<-NULL
    for(i in 1:length(samples)){
      samplecf[i]<-unlist(strsplit(samples[i],"-"))[2]
    }
    sj<-NULL
    for(i in 1:length(samplecf)){
      sj<-c(sj,substr(samplecf[i],1,6))
    }
    
    hf<-names(table(sj))
    
    jg<-list()
    for(y in 1:length(hf)){
      yy<-NULL
      pp<-which(sj==hf[y])
      all_data<-mRNAxl[pp,]
      all_sample<-samples[pp]
      for(i in 1:nrow(all_data)){
          a<-paste(">",all_sample[i],sep = "")
          b<-paste(all_data[i,],collapse = "")
          yy<-c(yy,a,b)
      }
      yy<-c(paste(">",refseqID,sep = ""),paste(refseq,collapse = ""),yy)
      jg[[y]]<-yy
    }
    names(jg)<-hf
    
  }
  
  if(divide=="continent"){
    countrys<-list(Europe=c("SLO","DEN","NET","NOR","RUS","DEN","GER","CZE","SWE","FIN","ITA","SPA",
              "TUR","BEL","LAT","LUX","AT","WAL","SWI","FRA","ENG","SCO"),Asia=c("CHI","BJ","ZJ","YN","WH","SZ","SD","SC","JX","JS","HZ","GZ","GD","FS","CQ",
            "FY","JZ","NC","XY","SR","TM","JJ","SIN","SK","JAP","TW","HK","NON","NEP","CAM","VIE",
            "PAK","JOR","SRI","MAL","IND","QAT","KUW","GEO","SAU","ISR","THA"),America=c("USA","MEX","CAN","BRA"))
    
    jg<-list()
    for(c in 1:3){
      pp<-NULL
      for(i in 1:length(countrys[[c]])){
          pp<-c(pp,grep(countrys[[c]][i],samples))
      }
      pp<-unique(pp)
      
      gj_f<-NULL
      gj_data<-mRNAxl[pp,]
      ls_sample<-samples[pp]
      for(i in 1:nrow(gj_data)){
          a<-paste(">",ls_sample[i],sep = "")
          b<-paste(gj_data[i,],collapse = "")
          gj_f<-c(gj_f,a,b)
      }
      gj_f<-c(paste(">",refseqID,sep = ""),paste(refseq,collapse = ""),gj_f)
      jg[[c]]<-gj_f
    }
    names(jg)<-names(countrys)
  }
  
  return(jg)
}