# The function GetqzM can get potential mutations. It is an internal program used by function GetDistMutORF.
GetqzM<-function(wdseq,refseq,refseqindex,parallel.sz=0){
  genes<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  mimz<-data.frame(code=c("ttt","ttc","tta","ttg","tct","tcc","tca","tcg","tat","tac","taa","tag","tgt","tgc","tga","tgg","ctt","ctc","cta","ctg","cct","ccc","cca","ccg","cat","cac","caa","cag","cgt","cgc","cga","cgg","att","atc","ata","atg","act","acc","aca","acg","aat","aac","aaa","aag","agt","agc","aga","agg","gtt","gtc","gta","gtg","gct","gcc","gca","gcg","gat","gac","gaa","gag","ggt","ggc","gga","ggg"),
                   amino=c("F","F","L","L","S","S","S","S","Y","Y","stopcode","stopcode","C","C","stopcode","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G"),stringsAsFactors =F)
  
  jjxl<-c("a","t","c","g")
  
  if(parallel.sz==0){
    nCores <- parallel::detectCores()
  }else{
    nCores <- parallel.sz
  }
  cl <- makeCluster(nCores)
  Mutation_v<-parLapply(cl,1:length(wdseq),function(n,wdseq,refseq,refseqindex,codematrix,mimz,jjxl){
     wd<-wdseq[n]
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
      
      if(length(ajzb)==3){
      dywz<-which(ajzb==wd)
      
      ppindex<-match(ajzb,refseqindex)
      ckajjj<-refseq[ppindex]
      ckaj<-mimz[which(mimz[,1]==paste(ckajjj,collapse  = "")),2]
      ckajjj1<-ckajjj
      jjxl1<-setdiff(jjxl,ckajjj[dywz])
      
      jg<-c(0,0)
      for(i in 1:length(jjxl1)){
        ckajjj1[dywz]<-jjxl1[i]
        tbaj<-mimz[which(mimz[,1]==paste(ckajjj1,collapse  = "")),2]
        
        if(ckaj==tbaj){
          jg[1]<-jg[1]+1
        }else{
          jg[2]<-jg[2]+1
        }
      }
      }else{
        jg<-c(0,0)
        for(j in 1:(length(ajzb)/3)){
          ajzb1<-ajzb[c((1+3*(j-1)):(j*3))]
          dywz<-which(ajzb1==wd)
          
           ppindex<-match(ajzb1,refseqindex)
           ckajjj<-refseq[ppindex]
           ckaj<-mimz[which(mimz[,1]==paste(ckajjj,collapse  = "")),2]
           ckajjj1<-ckajjj
           jjxl1<-setdiff(jjxl,ckajjj[dywz])
           
           for(i in 1:length(jjxl1)){
            ckajjj1[dywz]<-jjxl1[i]
            tbaj<-mimz[which(mimz[,1]==paste(ckajjj1,collapse  = "")),2]
        
            if(ckaj==tbaj){
              jg[1]<-jg[1]+1
            }else{
              jg[2]<-jg[2]+1
            }
         }
        }
      }
      return(jg)
  },wdseq,refseq,refseqindex,codematrix,mimz,jjxl)
  Mutation_v<-do.call("rbind",Mutation_v)
  colnames(Mutation_v)<-c("ty","fty")
  stopCluster(cl)
  Mutation_v<-Mutation_v/3
  return(Mutation_v)
}
