# The function GetSNP can identify the mutation site according to the nucleic acid sequence of the SARS-CoV-2 and count the mutation site information.
# Parameter fasta.file: Fasta file path.
# Parameter mutation.ft: The mutation frequency threshold of the mutation site. The mutation site can be selected. The default value: 0.01.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
# Parameter parallel.sz: Number of processors to use when doing the calculations in parallel (default value: 1). If parallel.sz=0, then it will use all available core processors unless we set this argument with a smaller number.
# Dependent package: parallel.
GetSNP<-function(fasta.file="",mutation.ft=0.01,refseqID="NC_045512.2",parallel.sz=0){
  genes<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  mimz<-data.frame(code=c("ttt","ttc","tta","ttg","tct","tcc","tca","tcg","tat","tac","taa","tag","tgt","tgc","tga","tgg","ctt","ctc","cta","ctg","cct","ccc","cca","ccg","cat","cac","caa","cag","cgt","cgc","cga","cgg","att","atc","ata","atg","act","acc","aca","acg","aat","aac","aaa","aag","agt","agc","aga","agg","gtt","gtc","gta","gtg","gct","gcc","gca","gcg","gat","gac","gaa","gag","ggt","ggc","gga","ggg"),
                   amino=c("F","F","L","L","S","S","S","S","Y","Y","stopcode","stopcode","C","C","stopcode","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"),stringsAsFactors =F)
  
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
  
  if(parallel.sz==0){
    nCores <- parallel::detectCores()
  }else{
    nCores <- parallel.sz
  }
  cl <- makeCluster(nCores)
  Mutation_v<-parLapply(cl,1:ncol(mRNAxl),function(n,mutation.ft,mRNAxl,refseq,refseqindex,genes,codematrix,mimz){
    jj_t<-table(mRNAxl[,n])
    qcindex<-which(is.element(names(jj_t),c("a","t","c","g"))==F)
    if(length(jj_t)>1){
      ref_jj<-refseq[n]
      jj_pl<-vector("numeric",length(jj_t) )
       for(i in 1:length(jj_t)){
        jj_pl[i]<-signif(jj_t[i]/sum(jj_t),digits = 2) 
       }
      if(length(qcindex)==0){
        pd<-min(jj_pl)>mutation.ft|names(jj_t)[which.min(jj_t)]==ref_jj
      }else{
        if(length(jj_pl)==2){
          pd<-min(jj_pl)>mutation.ft|names(jj_t)[which.min(jj_t)]==ref_jj
        }else{
          pd<-min(jj_pl[-qcindex])>mutation.ft|names(jj_t[-qcindex])[which.min(jj_t[-qcindex])]==ref_jj
        }
      }
      
      if(pd){ 
        ref_wz<-refseqindex[n] 
        if(is.na(ref_wz)==T){
         wz_pd<-0 
        }else{
          wz_pd<-sapply(1:nrow(codematrix), function(o){
                            x<-codematrix[o,]
                            if(ref_wz>=x[1]&ref_wz<=x[2]){
                                  dj<-(ref_wz-x[1])%%3
                                  if(dj==0){
                                      xq<-c(ref_wz,ref_wz+1,ref_wz+2)
                                  }
                                  if(dj==1){
                                      xq<-c(ref_wz-1,ref_wz,ref_wz+1)
                                  }
                                  if(dj==2){
                                      xq<-c(ref_wz-2,ref_wz-1,ref_wz)
                                  }
                                  aj_wd<-(xq[3]-x[1]+1)/3 
                                  ppindex<-match(xq,refseqindex)
                                  sjj<-apply(mRNAxl[,ppindex], 1, function(aj){paste(aj,collapse = "")})
                                  sjjpp<-match(sjj,mimz$code)
                                  sjjpp<-mimz$amino[sjjpp]
                                  sjjpp<-table(sjjpp)
                                  refajs<-mimz$amino[match(paste(refseq[ppindex],collapse = ""),mimz$code)]
                                  if(length(sjjpp)>1){
                                    tbhy<-("nonsynonymy")
                                  }else{
                                    tbhy<-("synonymy")
                                  }
                                  zh<-c(tbhy,paste(refseq[ppindex],collapse = ""),paste(names(table(sjj)),collapse = ", "),refajs,paste(names(sjjpp),collapse = ", "),aj_wd)
                              }else{
                                  return(0)
                              }})
        }
        if(all(wz_pd=="0")){
          qy<-"Noncoding_region"
          tbzt<-c(0,0,0,0,0,0)
        }else{
          qy<-genes[which(wz_pd!="0")] 
          tbzt<-unlist(wz_pd[which(wz_pd!="0")])
        }
        Mutation<-paste(names(jj_t),"_",jj_pl,sep = "",collapse = "///")
        if(length(qcindex)==0){
          jg<-c(n,ref_wz,qy,tbzt,refseq[n],Mutation,min(jj_pl))
        }else{
          jg<-c(n,ref_wz,qy,tbzt,refseq[n],Mutation,min(jj_pl[-qcindex]))
        }
        
        return(jg)
        
        }
      
      
    }
    
  },mutation.ft,mRNAxl,refseq,refseqindex,genes,codematrix,mimz)
  stopCluster(cl)
  cd<-NULL
  for(i in 1:length(Mutation_v)){
    cd[i]<-length(Mutation_v[[i]])
  }
  qc<-which(cd!=12)
  Mutation_v<-Mutation_v[-qc]
  Mutation_v<-do.call("rbind",Mutation_v)
  colnames(Mutation_v)<-c("Site","Reference sequence site","ORF","Mutation status","Reference sequence codon","Codon","Reference sequence amino acid","Mutant amino acid",
                          "Mutation site in ORF","Reference sequence nucleotide","Mutation frequency","Minimum mutation frequency")
  jg<-list(Mutation_v=Mutation_v,refseq=refseq,refseqindex=refseqindex)
  return(jg)
}