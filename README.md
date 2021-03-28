<h1 align="center">Functional alterations caused by mutations reflect evolutionary trends of SARS-Cov-2</h1>

## Introduction
FE-SARS-Cov-2 (Functional alterations caused by mutations reflect evolutionary trends of SARS-Cov-2) is a pipline based on the R language for the functional and evolutionary analysis of the nucleic acid sequence mutations of the SARS-Cov-2 strain. The main process is as follows:

![wf](https://github.com/chl198478/FE-SARS-Cov-2/blob/main/Figures/wf.png?raw=true)

&ensp; The procedure will be introduced below.

## Folder description
* Folder *R*: all R program files for  SARS-Cov-2 nucleic acid sequence analysis. The parameters contained in each function are described in the file.
* Folder *Software*: dependent software.
* Folder *test_data*: for theSARS-CoV-2 nucleic acid sequence data containing May to December, our analysis process.

## Sequence data preprocessing of SARS-CoV-2 strains
With [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) as the reference sequence, all SARS-CoV-2 nucleic acid sequence data are aligned by **MAFFT**. The function **ReInDe** is used to remove strains that contain insertions, deletions, and degenerate bases in their sequences. The code are as follows:
```
AllSeqData<-ReInDe(fasta.file,refseqID="NC_045512.2")
# Output as fasta file.
write.table(AllSeqData,file = ".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",sep ="\n",row.names = F,col.names = F,quote = F)
```
&ensp; AllSeqData has been saved in the folder *Nucleic Acid Sequence Data_May-October*.

&ensp; In addition, function **SeqDivision** can classify strains according to the time of diagnosis and the continent of patients infected with SARS-CoV-2. The code are as follows:
```
# The SARS-CoV-2 strains are classified according to the time the patient is diagnosed.
Time.class<-SeqDivision(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",divide="time",refseqID="NC_045512.2")
# The SARS-CoV-2 strains are classified according to the continent where the patient is located.
# The program only supports recognition in Asia, America and Europe.
Continent.class<-SeqDivision(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",divide="continent",refseqID="NC_045512.2")
# Output as fasta file.
for(i in 1:length(Continent.class)){
  write.table(Continent.class[[i]],file = paste(names(Continent.class)[i],".fasta",sep=""),sep ="\n",row.names = F,col.names = F,quote = F)
}
```

## The distribution of mutations in SARS-CoV-2 strains
The function **GetSNP** can identify the mutation site according to the nucleic acid sequence of the SARS-CoV-2 and count the mutation site information. The code are as follows:
```
library(parallel)
AllMutation<-GetSNP(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",mutation.ft=0,parallel.sz=0)
AllMutation_0.01<-GetSNP(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",mutation.ft=0.01,parallel.sz=0)
```
&ensp; The function **plotMtutSamp** can display the number of derived strains with individual mutations. The code are as follows:
```
library(reshape2)
library(ggplot2)
plotMtutSamp(AllMutation,sample.number=2872)
```

&ensp; The result is shown in Figure 1A.

![Figure 1](https://github.com/chl198478/FE-SARS-Cov-2/blob/main/Figures/fig1.tif?raw=true)
<p align="center"> Figure 1. A. The number of derived strains with individual mutations in SARS-CoV-2 virus; B. The average of accumulative mutations in each month.<p align="center"></p><br/>

&ensp; The function **MutCuml** is used to count the cumulative mutation frequency of mutation sites over time and **plotAcMu** can visualize the results. The code are as follows:
```
library(parallel)
library(ggplot2)
MutCuml.data<-MutCuml(Time.class,refseqID="NC_045512.2",parallel.sz=0)
plotAcMu(MutCuml.data)
```

&ensp; The result is shown in Figure 1B.

## The distribution of mutations in each of ORFs
The function **GetDistMutORF** can get the distribution and significance of mutation sites in each ORF region. The function **plotDistMutORFHist** will draw a histogram to display the results. The code are as follows:
```
library(parallel)
GetDistMutORF.data<-GetDistMutORF(AllMutation)
library(ggplot2)
library(cowplot)
library(showtext)
plotDistMutORFHist(GetDistMutORF.data)
```

&ensp; The result is shown in Figure 2.

![Figure 2](https://github.com/chl198478/FE-SARS-Cov-2/blob/main/Figures/fig2.tiff?raw=true)
<p align="center">igure 2. A. Significant score of the number of mutation locations in each of ORFs; B. Significant score of the number of synonymous mutation locations in each of ORFs; C. Significant score of the number of non-synonymous mutation locations in each of ORFs; D. Mutation rate in each of ORFs; E. Synonymous substitution rate in each of ORFs; F. Non-synonymous substitution rate in each of ORFs.<p align="center"></p><br/>

## Linkage and tendency of mutations occurred in 1% and more virus strains
The function **GetLD** can obtain the input required by **Haploview** according to the base sequence position. Then **Haploview** is able to identify the LD association between the mutation sites. The code are as follows:
```
temp<-AllMutation_0.01[["Mutation_v"]][which(AllMutation_0.01[["Mutation_v"]][,3]!="Noncoding_region"),]
seq.index<-as.numeric(temp[,2])
LD<-GetLD(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",seq.index)
write.table(LD$samples,file = "LD.ped",quote = F,row.names = F,col.names = F)
write.table(LD$SNP_sites,file = "LD.info",quote = F,row.names = F,col.names = F)
```

&ensp; Files LD.ped and LD.info are input into **Haploview** to get the linkage disequilibrium information of the mutation site.

&ensp; The function **plotLD** can display the linkage disequilibrium between the mutation sites as a scatter plot. The code are as follows:
```
ld_data<-read.table(".../test_data/ld1",stringsAsFactors=F,header=T)
ld_data<-ld_data[,c(4,5)]
library(ggplot2)
#LD.file is the output data of Haploview.
plotLD(".../test_data/ld_data.txt")
```

&ensp; The result is shown in Figure 3A.

![Figure 3](https://github.com/chl198478/FE-SARS-Cov-2/blob/main/Figures/fig3.tif?raw=true)
<p align="center">Figure 3. A. Scatter diagram of linkage disequilibrium between 36 SNPs. Horizontal axis and vertical axis represent r2 and LOD of pair-wise SNPs, respectively; B. The ratio of 36 mutations occurs in each month.<p align="center"></p><br/>

&ensp; The function **GetSitFreqInf** can get the change of mutation site frequency over time. The function **plotSiteTime** visualizes the results. The code are as follows:
```
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
```
##  Associations between SARS-CoV-2 strains among different continents
The function **HieClust** can re-encode the nucleic acid sequence of the SARS-CoV-2 and perform hierarchical clustering. The code are as follows:
```
library(dendextend)
library(factoextra)
library(flexclust)
library(RColorBrewer)
# Parameter plot can control whether the output of HieClust is a hierarchical clustering result or drawing a hierarchical clustering diagram.
# We built sample.class in folder TestData where the sample information of the continent classification is stored.
# Get hierarchical clustering result data.
sample.class<-GetContinentS(".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta")
HieClust.data<-HieClust(fasta.file=".../Nucleic Acid Sequence Data_May-October/AllSeqData.fasta",sample.class,refseqID="NC_045512.2",plot=F)
# Draw a hierarchical cluster diagram
HieClust(fasta.file="AllSeqData.fasta",sample.class,refseqID="NC_045512.2",plot=T)
```

&ensp; The hierarchical clustering diagram is shown in Figure 4.

![Figure 4.](https://github.com/chl198478/FE-SARS-Cov-2/blob/main/Figures/fig4.tiff?raw=true)
<p align="center">Figure 4. A. Scatter diagram of linkage disequilibrium. Horizontal axis and vertical axis represent r2 and LOD of pair-wise SNPs, respectively; B. The ratio of mutations occurs in each month.<p align="center"></p><br/>

&ensp; The **GetSpSite** function can obtain specific mutation sites of different strains.
```
Continent_seq_path<-".../Continent_seq/"
Continent_seq_files<-list.files(Continent_seq_path)
con.mut.list<-list()
for(i in 1:length(Continent_seq_files)){
  con.mut.list[[i]]<-GetSNP(fasta.file=paste(Continent_seq_path,Continent_seq_files[i],sep = ""),mutation.ft=0.01,parallel.sz=0)
}
names(con.mut.list)<-Continent_seq_files
# RaTG13.seq is the fasta data after alignment with the reference sequence NC_045512.2.
SpSite.data<-GetSpSite(con.mut.list,RaTG13.seq.file,refseqID="NC_045512.2")
tem<-NULL
for(i in 1:nrow(SpSite.data[["Specific_site_sequence"]])){
  a<-paste(SpSite.data[["Specific_site_sequence"]][i,],collapse  = "")
  b<-paste(">",row.names(SpSite.data[["Specific_site_sequence"]])[i],sep = "")
  tem<-c(tem,c(b,a))
}
write.table(tem,file = "SpSite.data.fasta",sep ="\n",row.names = F,col.names = F,quote = F)
# Draw VN diagram.
library(VennDiagram)
venn.diagram(SpSite.data[[1]],filename = "VN.tif",cat.cex=0.5,cex=2)
# Convert to phy format
library(phylotools)
dat <- read.fasta(".../SpSite.data.fasta")
dat2phylip(dat, outfile = "con_synonymy_seq.phy")
```

&ensp; Finally, the **GetSpSite** function can identify specific mutation sites in SARS-CoV-2 nucleic acid sequence groups and can integrate RaTG13 and NC_045512.2 sequences (The alignment sequence file of RaTG13 in NC_045512.2 has been stored in the folder test_data.) to output the nucleotide sequence composed of these sites. Then you can use the **PhyML** software to draw the evolutionary tree.

### Note:
Results have been updated to February 2021. The nucleotide sequence data is saved in the Nucleic Acid Sequence Data folder.

## References:
Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability, Mol Biol Evol 2013;30:772-780.

Guindon S, Dufayard J, Lefort V et al. New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0, Systematic Biology 2010;59:307-321.
