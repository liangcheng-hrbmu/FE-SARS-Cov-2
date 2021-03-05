# The GetContinentS function can get the strain ID of each continent grouped in the nucleic acid sequence.
# Parameter fasta.file: Fasta file path.
GetContinentS<-function(fasta.file="",refseqID="NC_045512.2"){
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
      
      ls_sample<-samples[pp]
      
      jg[[c]]<-ls_sample
    }
    names(jg)<-names(countrys)
  return(jg)
}