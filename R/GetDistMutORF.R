# The function GetDistMutORF can get the distribution and significance of mutation sites in each ORF region. 
# Parameter GetSNP.data: The result of the function GetSNP.
# Parameter parallel.sz: Number of processors to use when doing the calculations in parallel (default value: 1). If parallel.sz=0, then it will use all available core processors unless we set this argument with a smaller number.
# Dependent package: parallel.
GetDistMutORF<-function(GetSNP.data,parallel.sz=0){
  genes<-c("ORF1a","ORF1b","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
  codematrix<-matrix(c(266,13483,13468,21555,21563,25384,25393,26220,26245,26472,26523,27191,27202,27387,27394,27759,27756,27887,27894,28259,28274,29533,29558,29674),
                     ncol = 2,byrow = T)
  Allmutation<-GetSNP.data[[1]]
  refseq<-GetSNP.data[[2]]
  refseqindex<-GetSNP.data[[3]]
  
  qztb<-list()
  for(i in 1:length(genes)){
    cs<-GetqzM(c(codematrix[i,1]:codematrix[i,2]),refseq,refseqindex,parallel.sz)
    qztb[[i]]<-cs
  }
  all<-do.call("rbind",qztb)
  
  ty_Pvalue<-NULL
  fty_Pvalue<-NULL
  all_Pvalue<-NULL
  for(i in 1:length(genes)){
    cs<-qztb[[i]]
    ty_Pvalue<-c(ty_Pvalue,fisher.test(matrix(c(length(which(Allmutation[,3]==genes[i]&Allmutation[,4]=="synonymy")),sum(cs[,1]),length(which(Allmutation[,4]=="synonymy")),sum(all[,1])),nrow=2))$p.value)
    fty_Pvalue<-c(fty_Pvalue,fisher.test(matrix(c(length(which(Allmutation[,3]==genes[i]&Allmutation[,4]=="nonsynonymy")),sum(cs[,2]),length(which(Allmutation[,4]=="nonsynonymy")),sum(all[,2])),nrow=2))$p.value)
    all_Pvalue<-c(all_Pvalue,fisher.test(matrix(c(length(which(Allmutation[,3]==genes[i])),length(c(codematrix[i,1]:codematrix[i,2])),length(which(Allmutation[,3]!="Noncoding_region")),29280),nrow=2))$p.value)

  }
  ORF_trend1<-data.frame(ORF=genes,ty_Pvalue,fty_Pvalue,all_Pvalue,log10ty_pvalue=(-log10(ty_Pvalue)),log10fty_pvalue=(-log10(fty_Pvalue)),log10all_pvalue=(-log10(all_Pvalue)),stringsAsFactors = F)
  colnames(ORF_trend1)<-c("ORFs","Synonymous_Pvalues","NonSynonymous_Pvalues","Mutation_Pvalues","(-log10(Synonymous_Pvalues))","(-log10(NonSynonymous_Pvalues))","(-log10(Mutation_Pvalues))")
  
  ty_qz<-NULL
  fty_qz<-NULL
  bftyftytb_qz<-NULL
  for(i in 1:length(genes)){
    cs<-qztb[[i]]
    ty_qz<-c(ty_qz,length(which(Allmutation[,3]==genes[i]&Allmutation[,4]=="synonymy"))/sum(cs[,1]))
    fty_qz<-c(fty_qz,length(which(Allmutation[,3]==genes[i]&Allmutation[,4]=="nonsynonymy"))/sum(cs[,2]))
    bftyftytb_qz<-c(bftyftytb_qz,length(which(Allmutation[,3]==genes[i]))/length(c(codematrix[i,1]:codematrix[i,2])))
  }
  bl<-data.frame(ORF=genes,ty_qz,fty_qz,bftyftytb_qz,stringsAsFactors = F)
  colnames(bl)<-c("ORFs","Synonymous_substitution_rate","Nonsynonymous_substitution_rate",
                "Mutation_rate")
  jg<-list(Significance=ORF_trend1,Percentage=bl)
  return(jg)
}