#Virus sequence from May to December

AllSeqData<-ReInDe(fasta.file,refseqID="NC_045512.2")
# Output as fasta file.
# SeqData.fasta has been stored in the Nucleic Acid Sequence Data_May-October folder.
write.table(AllSeqData,file = ".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",sep ="\n",row.names = F,col.names = F,quote = F)

Time.class<-SeqDivision(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",divide="time",refseqID="NC_045512.2")
save(Time.class,file = ".../test_data/Time.class")
Continent.class<-SeqDivision(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",divide="continent",refseqID="NC_045512.2")
for(i in 1:length(Continent.class)){
  write.table(Continent.class[[i]],file = paste(names(Continent.class)[i],".fasta",sep=""),sep ="\n",row.names = F,col.names = F,quote = F)
}

library(parallel)
AllMutation<-GetSNP(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",mutation.ft=0,parallel.sz=0)
save(AllMutation,file = ".../test_data/AllMutation.Rdata")
AllMutation_0.01<-GetSNP(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",mutation.ft=0.01,parallel.sz=0)
save(AllMutation_0.01,file = ".../test_data/AllMutation_0.01.Rdata")

library(reshape2)
library(ggplot2)
plotMtutSamp(AllMutation,sample.number=2872)

library(parallel)
library(ggplot2)
MutCuml.data<-MutCuml(Time.class,refseqID="NC_045512.2",parallel.sz=0)
save(MutCuml.data,file = ".../test_data/MutCuml.data.Rdata")
plotAcMu(MutCuml.data)

library(parallel)
GetDistMutORF.data<-GetDistMutORF(AllMutation)
save(GetDistMutORF.data,file = ".../test_data/GetDistMutORF.data.Rdata")
library(ggplot2)
library(cowplot)
library(showtext)
plotDistMutORFHist(GetDistMutORF.data)

temp<-AllMutation_0.01[["Mutation_v"]][which(AllMutation_0.01[["Mutation_v"]][,3]!="Noncoding_region"),]
seq.index<-as.numeric(temp[,2])
LD<-GetLD(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",seq.index)
write.table(LD$samples,file = "LD.ped",quote = F,row.names = F,col.names = F)
write.table(LD$SNP_sites,file = "LD.info",quote = F,row.names = F,col.names = F)
ld_data<-read.table(".../test_data/ld1",stringsAsFactors=F,header=T)
ld_data<-ld_data[,c(4,5)]
write.table(ld_data,file = ".../test_data/ld_data.txt",row.names = F,quote=F)
library(ggplot2)
plotLD(".../test_data/ld_data.txt")

Time_seq_path<-".../test_data/Time_seq/"
Time_seq_files<-list.files(Time_seq_path)
Time_seq_list<-list()
for(i in 1:length(Time_seq_files)){
  Time_seq_list[[i]]<-GetSNP(fasta.file=paste(Time_seq_path,Time_seq_files[i],sep = ""),mutation.ft=0.01,parallel.sz=0)
}
names(Time_seq_list)<-Time_seq_files
sitfreqinf.data<-GetSitFreqInf(Time_seq_list,seq.index)
library(ggplot2)
library(reshape2)
plotSiteTime(sitfreqinf.data)

library(dendextend)
library(factoextra)
library(flexclust)
library(RColorBrewer)
sample.class<-GetContinentS(".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta")
HieClust.data<-HieClust(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",sample.class,refseqID="NC_045512.2",plot=F)
HieClust(fasta.file="AllSeqData.fasta",sample.class,refseqID="NC_045512.2",plot=T)