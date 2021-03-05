# The Function plotAcMu can visualize the results of function MutCuml.
# Parameter GetDistMutORF.data: The result of the function GetDistMutORF.
# Dependent package: ggplot2, cowplot, showtext.
plotDistMutORFHist<-function(GetDistMutORF.data){
  win.graph()
  plot1<-ggplot(GetDistMutORF.data[["Significance"]],aes(x=ORFs, y=(-log10(Mutation_Pvalues)), fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round((-log10(Mutation_Pvalues)),3)),vjust=-0.2)+
    labs(x = "All mutation", y = "-log10(P value)")
  
  plot2<-ggplot(GetDistMutORF.data[["Significance"]],aes(x=ORFs, y=(-log10(Synonymous_Pvalues)), fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round((-log10(Synonymous_Pvalues)),3)),vjust=-0.2)+
    labs(x = "Synonymous mutation", y = "-log10(P value)")
  
  plot3<-ggplot(GetDistMutORF.data[["Significance"]],aes(x=ORFs, y=(-log10(NonSynonymous_Pvalues)), fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round((-log10(Synonymous_Pvalues)),3)),vjust=-0.2)+
    labs(x = "Non-synonymous mutation", y = "-log10(P value)")
  
  plot4<-ggplot(GetDistMutORF.data[["Percentage"]],aes(x=ORFs, y=Mutation_rate, fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round(Mutation_rate,3)),vjust=-0.2)+
    labs(x = "All mutation", y = "Percentage(%)")
  
  plot5<-ggplot(GetDistMutORF.data[["Percentage"]],aes(x=ORFs, y=Synonymous_substitution_rate, fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round(Synonymous_substitution_rate,3)),vjust=-0.2)+
    labs(x = "Synonymous mutation", y = "Percentage(%)")
  
  plot6<-ggplot(GetDistMutORF.data[["Percentage"]],aes(x=ORFs, y=Nonsynonymous_substitution_rate, fill=ORFs))+geom_bar(stat='identity')+
    geom_text(aes(label=round(Nonsynonymous_substitution_rate,3)),vjust=-0.2)+
    labs(x = "Non-synonymous mutation", y = "Percentage(%)")
  
  gg <- ggdraw() + 
    draw_plot(plot1, 0, 1, 1/3, 0.5) + draw_plot(plot2, 1/3, 1, 1/3, 0.5) +
    draw_plot(plot3, 1, 1, 1/3, 0.5)+draw_plot(plot4, 0, 0, 1/3, 0.5)+
    draw_plot(plot5, 1/3, 0, 1/3, 0.5)+draw_plot(plot6, 1, 0, 1/3, 0.5)
  p5 <- cowplot::plot_grid(plot1, plot2, plot3, plot4,plot5,plot6, nrow = 2)
  p5
}