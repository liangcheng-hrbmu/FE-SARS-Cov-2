# The Function plotAcMu can visualize the results of function MutCuml.
# Parameter MutCuml.data: The result of the function MutCuml.
# Parameter axis.text.size: The size of axis.text.
# Parameter axis.title.size: The size of axis.title.
# Parameter legend.text.size: The size of legend.text.
# Dependent package: ggplot2.
plotAcMu<-function(MutCuml.data,axis.text.size=15,axis.title.size=20,legend.text.size=10){
  
  MutCuml.data1<-data.frame(times=rep(row.names(MutCuml.data),3),values=c(MutCuml.data[,1],MutCuml.data[,2],MutCuml.data[,3]),
                            class=c(rep("Mutation",nrow(MutCuml.data)),rep("SynMutation",nrow(MutCuml.data)),rep("NonSynMutation",nrow(MutCuml.data))))
  
  ggplot(data=MutCuml.data1, aes(x=times, y=values,colour=class,group=class))+geom_line(size=1.5)+geom_point(size=4)+
  labs(x = "Month", y = "Average")+theme(axis.text=element_text(size = axis.text.size),
                                         axis.title=element_text(size = axis.title.size),
                                         legend.text=element_text(size = legend.text.size),
                                         legend.title=element_blank()
                                         )
  
}