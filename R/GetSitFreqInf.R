# The function GetSitFreqInf can get the change of mutation site frequency over time.
# Parameter SeqDivision.data: The result of the function SeqDivision.
# Parameter seq.index: A numeric vector is the position of the mutation site in the reference sequence where you want to calculate the linkage disequilibrium.
GetSitFreqInf<-function(SeqDivision.data,seq.index){
  jg<-NULL
  for(i in 1:length(seq.index)){
    wdxx<-NULL
    wd<-seq.index[i]
    for(j in 1:length(SeqDivision.data)){
      tb_table<-SeqDivision.data[[j]][[1]]
      pp<-match(wd,tb_table[,2])
      if(is.na(pp)==T){
        wdxx<-c(wdxx,0)
      }else{
        ck<-tb_table[pp,10]
        jjxl<-tb_table[pp,11]
        jjxl<-unlist(strsplit(jjxl,"///"))
        ck_index<-unlist(grep(ck,jjxl))
        jjxl1<-toupper(c(jjxl[ck_index],jjxl[-ck_index]))
        jjxl1<-paste(jjxl1,collapse = "/")
        wdxx<-c(wdxx,jjxl1)
      }
    }
    jg<-rbind(jg,wdxx)
  }
  row.names(jg)<-seq.index
  colnames(jg)<-names(SeqDivision.data)
  return(jg)

}