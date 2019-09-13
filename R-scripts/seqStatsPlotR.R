library(ggplot2)
library(readr)
library(seqinr)
seqStatsPlotR<-function(
  input_path,
  output_path="."){
  files=list.files(path = input_path, pattern = "*\\.fasta$", full.names = TRUE)
  files1=list.files(path = input_path, pattern = "*\\.fa$", full.names = TRUE)
  files<-c(files,files1)
  GC_content_global<-c()
  Sequence_Length_global<-c()
  for(i in 1: length(files)){
    file_fasta=strsplit(files[i], '/')[[1]][length(strsplit(files[i], '/')[[1]])]
    filename=strsplit(file_fasta, '.f')[[1]][1]
    file_current<-read.fasta(files[1])
    GC_content<-c()
    Sequence_Length<-c()
    for(j in 1:length(file_current)){
      GC_content<-c(GC_content, round(GC(file_current[[j]])*100))
      Sequence_Length<-c(Sequence_Length, length(file_current[[j]]))
    }
    GC_content_global<-c(GC_content_global, GC_content)
    Sequence_Length_global<-c(Sequence_Length_global, Sequence_Length)
    GC_content<-as.data.frame(GC_content)
    Sequence_Length<-as.data.frame(Sequence_Length)
    Sequence_Length_freq<-data.frame(table(Sequence_Length))
    k<-sum(Sequence_Length_freq$Freq)
    Sequence_Length_freq$Freq<-Sequence_Length_freq$Freq*100/k
    GC_content_freq<-data.frame(GC_content=seq(0,99,by=5), frequency=0)
    for(l in 1:20){
      GC_content_freq$frequency[l]=as.numeric(length(which(GC_content$GC_content>GC_content_freq$GC_content[l]&GC_content$GC_content<(GC_content_freq$GC_content[l]+5))))
    }
    k = sum(GC_content_freq$frequency)
    GC_content_freq$frequency<-GC_content_freq$frequency*100/k
    ggplot(GC_content_freq, aes(x=GC_content, y=frequency))+
      geom_bar(stat="identity", fill="deepskyblue3")+
      theme_classic()+
      scale_x_continuous(breaks = seq(0, 100, by = 10))+
      xlab("GC content (%)")+
      ylab("Frequency (%)")+
      geom_vline(aes(xintercept=50), linetype=3)
    ggsave(paste0(output_path,"/", filename, "_GC_content.pdf"),  height = 8.5, width = 10)
    
    
    ggplot(Sequence_Length_freq, aes(x=Sequence_Length, y=Freq))+
      geom_bar(stat="identity", fill="deepskyblue3")+
      theme_classic()+
      xlab("sequence Length")+
      ylab("Frequency (%)")
    ggsave(paste0(output_path,"/", filename, "_sequence_Length.pdf"),  height = 8.5, width = 10)
  }


GC_content_global<-as.data.frame(GC_content_global)
Sequence_Length_global<-as.data.frame(Sequence_Length_global)
Sequence_Length_global_freq<-data.frame(table(Sequence_Length_global))
k<-sum(Sequence_Length_global_freq$Freq)
Sequence_Length_global_freq$Freq<-Sequence_Length_global_freq$Freq*100/k
GC_content_global_freq<-data.frame(GC_content_global=seq(0,99,by=5), frequency=0)
for(l in 1:20){
  GC_content_global_freq$frequency[l]=as.numeric(length(which(GC_content_global$GC_content_global>GC_content_global_freq$GC_content_global[l]&GC_content_global$GC_content_global<(GC_content_global_freq$GC_content_global[l]+5))))
}
k = sum(GC_content_global_freq$frequency)
GC_content_global_freq$frequency<-GC_content_global_freq$frequency*100/k
ggplot(GC_content_global_freq, aes(x=GC_content_global, y=frequency))+
  geom_bar(stat="identity", fill="deepskyblue3")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 100, by = 10))+
  xlab("GC content (%)")+
  ylab("Frequency (%)")+
  geom_vline(aes(xintercept=50), linetype=3)
ggsave(paste0(output_path,"/","Global_GC_content.pdf"),  height = 8.5, width = 10)


ggplot(Sequence_Length_global_freq, aes(x=Sequence_Length_global, y=Freq))+
  geom_bar(stat="identity", fill="deepskyblue3")+
  theme_classic()+
  xlab("Sequence Length")+
  ylab("Frequency (%)")
ggsave(paste0(output_path,"/", "Global_sequenceLength.pdf"),  height = 8.5, width = 10)

}
