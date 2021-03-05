# The function plotLD can display the linkage disequilibrium between the mutation sites as a scatter plot.
# Parameter LD.file: The output data of Haploview.
# Dependent package: ggplot2.
plotLD<-function(LD.file=""){
  LD.table<-read.table(LD.file,stringsAsFactors = F,header = T)
  
  LD.table$LOD <- log10(LD.table$LOD)

  par7 <- ggplot(data = LD.table, mapping = aes(x = r.2, y = LOD)) + geom_point(color= "blue",size = 1.75,alpha=0.5)+
    xlab("r^2")+ylab("log10(LOD)")+
    theme_bw()+theme(panel.grid =element_blank())+
    theme(text=element_text(size=16,  family="serif"))+
    theme(legend.position = "none")+
    geom_hline(aes(yintercept=2), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=0.95), colour="red", linetype="dashed")+
    geom_rug()


  par7
  
}
