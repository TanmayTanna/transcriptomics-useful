library(ggplot2)
library(readr)
library(readxl)
library(reshape2)

RNAseqCountsPath = "NC_012947.1_RNAseq.txt"
RNAseqCounts = as.data.frame(read.table(RNAseqCountsPath, header = TRUE))
rownames(RNAseqCounts)=RNAseqCounts[,1]
genelengths<-RNAseqCounts$Length
RNAseqCounts=RNAseqCounts[,-c(1:6)]
for(i in 1:dim(RNAseqCounts)[2]){
  RNAseqCounts[,i]<-RNAseqCounts[,i]/genelengths
  RNAseqCounts[,i]<-RNAseqCounts[,i]*1000000/sum(RNAseqCounts[,i])
}

o=order(rowMeans(RNAseqCounts), decreasing = TRUE)
RNAseqCounts=RNAseqCounts[o,]
plasmid_countspaths=list.files(path = ".", pattern = "*Counts.txt")
Counts<-list()
for(i in 1:length(plasmid_countspaths)){
  plasmid_name = strsplit(plasmid_countspaths[i],".txt")[[1]][1]
  Counts[[plasmid_name]] = as.data.frame(read.table(plasmid_countspaths[i], header = TRUE))
  rownames(Counts[[plasmid_name]])=Counts[[plasmid_name]][,1]
  Counts[[plasmid_name]]=Counts[[plasmid_name]][,-c(1:6)]
}

for(i in 1:length(Counts)){
  Counts[[i]]<-Counts[[i]][match(rownames(RNAseqCounts), rownames(Counts[[i]])),]
  for(j in 1:length(colnames(Counts[[i]]))){
    if(grepl("_MC_",colnames(Counts[[i]])[j])){
    colnames(Counts[[i]])[j]= strsplit(colnames(Counts[[i]])[j], "2018")[[1]][length(strsplit(colnames(Counts[[i]])[j], "2018")[[1]])]
    colnames(Counts[[i]])[j]= paste0("Sample_2018",strsplit(colnames(Counts[[i]])[j], "_BL21_")[[1]][1])
    }
    print(colnames(Counts[[i]][j]))
    k=sum(Counts[[i]][,j])
    Counts[[i]][,j]=Counts[[i]][,j]/k
    Counts[[i]][,j]=cumsum(Counts[[i]][,j])
    }
  
  len=length(colnames(Counts[[i]]))+1
  Counts[[i]][,len]<-rowMeans(Counts[[i]])
  colnames(Counts[[i]])[len] = "Mean_Cumulative_transcript_Counts"
  len2=length(rownames(Counts[[i]]))
  write.csv(Counts[[i]], file = paste0("TranscriptionalBiasCounts_",as.character(names(Counts)[i]),".csv"))
  
  ## Uncomment this section if you want to plot only mean of replicates ##  
  # ggplot(Counts[[i]], aes(x=1:len2, y=Mean_Cumulative_transcript_Counts))+
  #   geom_smooth(col="DARKBLUE")+xlab("Genes sorted by RNA-seq counts")+ylab("Cumulative transcript count fraction")+geom_line(aes(x=1:len2, y=seq(1/len2,1,1/len2)), col="RED")+
  #   theme_classic()
  # ggsave(paste0("TranscriptionalBiasPlot_",as.character(names(Counts)[i]), ".svg"), height = 8.5, width = 10)
  
  ## Uncomment this section if you want to plot all replicates on the same plot ## 
  plot_counts<-data.frame(GeneNo=1:dim(Counts[[i]])[1], Counts[[i]])
  plot_counts<-melt(plot_counts, id.vars = "GeneNo")
  colnames(plot_counts)[2]<-"Replicate" 
  ggplot(plot_counts, aes(GeneNo, value, color=Replicate))+
    geom_smooth()+theme_classic()+xlab("Genes sorted by RNA-seq counts")+ylab("Cumulative transcript count fraction by replicate")
  ggsave(paste0("TranscriptionalBiasPlot_",as.character(names(Counts)[i]), ".svg"), height = 8.5, width = 10)
}

## Uncomment this section if you want to compare mean of all replicates across samples on the same plot ## 

Comparisons<-data.frame(row.names = rownames(Counts[[1]]))
Comparisons$GeneNo=1:dim(Comparisons)[1]
for(i in 1:length(Counts)){
  Comparisons[names(Counts)[i]]=Counts[[i]]$Mean_Cumulative_transcript_Counts
}

Comparisons<-melt(Comparisons, id.vars = "GeneNo")
colnames(Comparisons)[2]<-"Sample" 
  
ggplot(Comparisons, aes(GeneNo, value, color=Sample))+
  geom_smooth()+theme_classic()+xlab("Genes sorted by RNA-seq counts")+ylab("Cumulative transcript count fraction")
ggsave(paste0("TranscriptionalBiasPlot_Comparison.svg"), height = 8.5, width = 10)

  