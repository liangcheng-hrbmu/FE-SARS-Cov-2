# The function plotMtutSamp can display the number of derived strains with individual mutations.
# Parameter GetSNP.data: The result of function GetSNP.
# Parameter sample.number: The number of strain samples is used to obtain GetSNP.data.
# Dependent package: reshape2, ggplot2.
plotMtutSamp<-function(GetSNP.data,sample.number){
  GetSNP.data<-GetSNP.data[[1]]
  GetSNP.data<-GetSNP.data[which(GetSNP.data[,3]!="Noncoding_region"),]
  
  gs<-NULL
  for(i in 1:nrow(GetSNP.data)){
    sefjj<-GetSNP.data[i,10]
    jj<-unlist(strsplit(GetSNP.data[i,11],"///"))
    jj<-jj[-grep(sefjj,jj)]
    sn<-as.numeric(unlist(strsplit(jj,"_"))[2])*sample.number
    gs[i]<-sn
  }
  gs<-ceiling(gs)
  
  gs1<-sort(unique(gs))
  sy<-NULL
  nsy<-NULL
  for(j in 1:length(gs1)){
    pp<-which(gs==gs1[j])
    
    sy[j]<-length(which(GetSNP.data[pp,4]=="synonymy"))
    nsy[j]<-length(which(GetSNP.data[pp,4]=="nonsynonymy"))
    
  }
  data<-data.frame(Var1=gs1,nonsynonymy=nsy,synonymy=sy,stringsAsFactors = F)
  
  data$Var1 <- as.factor(data$Var1)
  data <- melt(data)
  
  par <- ggplot(data,aes(x=Var1,y=value,fill=variable))+geom_bar(stat = "identity",position = "dodge")+
    xlab("Frequency of derived mutation")+ylab("Number of mutations")+
    theme_bw()+theme(panel.grid =element_blank())+
    theme(text=element_text(size=16,  family="serif"))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
    geom_text(aes(label=value, y=value+5), position=position_dodge(0.8), vjust=0,size = 3.5)
  par
  
}