# The GetSpSite function can obtain specific mutation sites of different strains.
# Parameter con.mut.list: It is the tabular data composed of grouped nucleic acid sequences of different strains.
# Parameter RaTG13.seq.file: It is the fasta data after alignment with the reference sequence NC_045512.2.
# Parameter refseqID: The ID of the reference nucleic acid sequence.

#获取特异性位点信息和序列
GetSpSite<-function(con.mut.list,RaTG13.seq.file="",refseqID="NC_045512.2"){
  cds.list<-list()
  for(i in 1:length(con.mut.list)){
    cds.list[[i]]<-con.mut.list[[i]][which(con.mut.list[[i]][,3]!="Noncoding_region"),2]
  }
  
  s_s<-list()
  s_wd<-list()
  for (j in 1:length(cds.list)) {
    t<-setdiff(cds.list[[j]],Reduce(union, cds.list[-j]))
    s_s[[j]]<-con.mut.list[[j]][match(t,con.mut.list[[j]][,2]),]
    s_wd[[j]]<-t
  }
  names(s_s)<-names(con.mut.list)
  
  test<-readLines(RaTG13.seq.file)
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
  pp<-which(samples==refseqID)
  RaTG13refseq<-mRNAxl[[pp]] 
  mRNAxl<-mRNAxl[-pp]
  samples<-samples[-pp]
  RaTG13refseqindex<-NULL 
  j=0
  for (i in 1:length(RaTG13refseq)) {
    if(RaTG13refseq[i]=="-"){
      RaTG13refseqindex<-c(RaTG13refseqindex,NA)
    }else{
      j<-j+1
      RaTG13refseqindex<-c(RaTG13refseqindex,j)
    }
  }
  RaTG13seq<-do.call("rbind",mRNAxl)
  
  
  jg_seq<-NULL
  for(s in 1:length(s_wd)){
    ty_wd_seq1<-NULL
    for(i in 1:length(s_wd)){
      swd<-s_wd[[i]]
      if(i == s){
        s_dq<-s_s[[i]]
          for(j in 1:length(swd)){
            pp<-which(s_dq[,2]==swd[j])
            ck<-s_dq[pp,10]
            tb<-s_dq[pp,11]
            tb<-unlist(strsplit(tb,"///"))
            tb<-tb[-grep(ck,tb)]
            if(length(tb)>=2){
              tb_jj1<-unlist(strsplit(tb,"_"))
            blindex<-which.max(as.numeric(tb_jj1[c(1:length(tb_jj1))%%2==0]))
              tb<-tb[blindex]
              tb<-unlist(strsplit(tb,"_"))[1]
            }else{
              tb<-unlist(strsplit(tb,"_"))[1]
            }
            ty_wd_seq1<-c(ty_wd_seq1,tb)
          }

      }else{
        ckindex<-match(swd,RaTG13refseqindex)
        ty_wd_seq1<-c(ty_wd_seq1,RaTG13refseq[ckindex])
      }
    }
    jg_seq<-rbind(jg_seq,ty_wd_seq1)
  }
  colnames(jg_seq)<-unlist(s_wd)
  
  pp<-match(unlist(s_wd),RaTG13refseqindex)
  NC_seq<-RaTG13refseq[pp]
  RaTG13_seq<-RaTG13seq[pp]
  jg_seq<-rbind(jg_seq,NC_seq)
  jg_seq<-rbind(jg_seq,RaTG13_seq)
  
  
  row.names(jg_seq)<-c(names(s_s),"NC_045512","Bat_RaTG13")
  
  jg<-list(Specific_site_information=s_s,Specific_site_sequence=jg_seq)
  return(jg)
}


