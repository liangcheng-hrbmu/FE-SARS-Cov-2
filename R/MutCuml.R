# The function MutCuml is used to count the cumulative mutation frequency of mutation sites over time.
# Parameter SeqDivision.times: The result of the function SeqDivision obtained by dividing=time.
# Parameter refseqID: The ID of the reference nucleic acid sequence.
# Parameter parallel.sz: Number of processors to use when doing the calculations in parallel (default value: 1). If parallel.sz=0, then it will use all available core processors unless we set this argument with a smaller number.
# Dependent package: parallel.
MutCuml<-function(SeqDivision.times,refseqID="NC_045512.2",parallel.sz=0){
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  mimz<-data.frame(code=c("ttt","ttc","tta","ttg","tct","tcc","tca","tcg","tat","tac","taa","tag","tgt","tgc","tga","tgg","ctt","ctc","cta","ctg","cct","ccc","cca","ccg","cat","cac","caa","cag","cgt","cgc","cga","cgg","att","atc","ata","atg","act","acc","aca","acg","aat","aac","aaa","aag","agt","agc","aga","agg","gtt","gtc","gta","gtg","gct","gcc","gca","gcg","gat","gac","gaa","gag","ggt","ggc","gga","ggg"),
                   amino=c("F","F","L","L","S","S","S","S","Y","Y","stopcode","stopcode","C","C","stopcode","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"),stringsAsFactors =F)
  codeindex<-sapply(c(1:nrow(codematrix)), function(i){
    return(c(codematrix[i,1]:codematrix[i,2]))
  })
  codeindex<-unlist(codeindex)
  
  jg1<-vector("numeric",length = 3)
  for(i in 1:length(SeqDivision.times)){
    a<-SeqDivision.times[[i]]
    samples<-a[c(1:length(a)) %% 2!=0]
    samples<-gsub(">","",samples)
    
    seqs<-a[c(1:length(a)) %% 2==0]
    mRNAxl<-list()
    for(s in 1:length(seqs)){
      mRNAxl[s]<-strsplit(seqs[s],"")
    }
    mRNAxl<-do.call("rbind",mRNAxl)
    
    pp<-which(samples==refseqID)
    refseq<-mRNAxl[pp,]
    mRNAxl<-mRNAxl[-pp,]
    samples<-samples[-pp]
    
    refseqindex<-NULL 
    j=0
    for (c in 1:length(refseq)) {
      if(refseq[c]=="-"){
        refseqindex<-c(refseqindex,NA)
      }else{
        j<-j+1
        refseqindex<-c(refseqindex,j)
      }
    }
    
    if(parallel.sz==0){
    nCores <- parallel::detectCores()
  }else{
    nCores <- parallel.sz
  }
  cl <- makeCluster(nCores)
  
  Mutation_v<-parLapply(cl,1:nrow(mRNAxl),function(m,codeindex,mRNAxl,refseq,refseqindex,codematrix,mimz){
    sampleseq<-mRNAxl[m,]
    jg<-c(0,0,0)
    for(n in 1:length(codeindex)){
    
    wd<-codeindex[n]
    sjwd<-which(refseqindex==wd)
    ajj<-sampleseq[sjwd]
    ckjj<-refseq[sjwd]
    if(ajj!=ckjj){
      jg[1]<-jg[1]+1
      
      ajzb<-sapply(1:nrow(codematrix), function(o){
        x<-codematrix[o,]
        if(wd>=x[1]&wd<=x[2]){
          dj<-(wd-x[1])%%3
          if(dj==0){
                xq<-c(wd,wd+1,wd+2)
              }
          if(dj==1){
                xq<-c(wd-1,wd,wd+1)
              }
          if(dj==2){
                xq<-c(wd-2,wd-1,wd)
              }
        }else{
          xq<-0
        }
        return(xq)
      })
      ajzb<-unlist(ajzb)
      ajzb<-ajzb[which(ajzb!=0)]
      
      dywz<-which(ajzb==wd)
      ppindex<-match(ajzb,refseqindex)
      ckajjj<-refseq[ppindex]
      ckaj<-mimz[which(mimz[,1]==paste(ckajjj,collapse  = "")),2]
      
        ckajjj[dywz]<-ajj
        tbaj<-mimz[which(mimz[,1]==paste(ckajjj,collapse  = "")),2]
        if(length(ckaj)!=0&length(tbaj)!=0){
          if(ckaj==tbaj){
            jg[2]<-jg[2]+1
          }else{
            jg[3]<-jg[3]+1
          }
        }

    }
    
    }
    return(jg)
  },codeindex,mRNAxl,refseq,refseqindex,codematrix,mimz)
  Mutation_v<-do.call("rbind",Mutation_v)
  stopCluster(cl)
    
   jg1<-rbind(jg1,c(sum(Mutation_v[,1])/nrow(Mutation_v),sum(Mutation_v[,2])/nrow(Mutation_v),sum(Mutation_v[,3])/nrow(Mutation_v))) 
  }
  jg1<-jg1[-1,]
  row.names(jg1)<-names(SeqDivision.times)
  colnames(jg1)<-c("Mutation","SynonymousMutation","Non-synonymousMutation")
  return(jg1)
}