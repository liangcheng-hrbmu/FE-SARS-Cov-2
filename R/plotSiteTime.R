# The function plotSiteTime visualizes the results of GetSitFreqInf.
# Parameter sitfreqinf.data: The results of GetSitFreqInf.
# Dependent package: ggplot2, reshape2.
plotSiteTime<-function(sitfreqinf.data){
  sitfreqinf.data<-as.data.frame(sitfreqinf.data,optional = T,stringsAsFactors = F)
  for(i in 1:ncol(sitfreqinf.data)){
    for(j in 1:nrow(sitfreqinf.data))
      if(sitfreqinf.data[j,i]!="0"){
        sitfreqinf.data[j,i]<-unlist(strsplit(unlist(strsplit(sitfreqinf.data[j,i],"/"))[2],"_"))[2]
      }
    sitfreqinf.data[,i]<-as.numeric(sitfreqinf.data[,i])
    
  }
  
  sitfreqinf.data$Location<-row.names(sitfreqinf.data)
  
  sitfreqinf.data$col<-rep("black",nrow(sitfreqinf.data))
  
  
  data<- melt(sitfreqinf.data,id.vars = c("Location","col"))
  
  par8<- ggplot(data,aes(x=Location,y=variable,color=value))+geom_point(aes(size=value))+

  scale_color_gradientn(values = seq(0,1,0.2),colours = c('cyan','blue','green','orange','red'))+
  xlab("Location")+ylab("Month")+
  theme_bw()+theme(panel.grid =element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  theme(text=element_text(size=16,  family="serif"))

  par8
}