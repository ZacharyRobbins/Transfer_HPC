
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

wd<-"C:/Users/345578/Documents/GitHub/FATES_BCI_transfer/extract_6_7/"
h1s=list.files(wd,pattern = "\\h1.extract.Rdata$")
h2s=list.files(wd,pattern = "\\h2.extract.Rdata$")
years<-5
###Need to make a dir named Figures.

#####################Opertations
for(cycle in 1:length(h1s)){
  print(h1s[cycle])
  load(paste0(wd,h1s[cycle]))
  LU<-as.data.frame(t(as.data.frame(var.res.arr)))
  colnames(LU)<-paste0("Replicate",seq(1,years))
  sub<-LU# [,sample(seq(1,100),50)]
  png(file=paste0(wd,"Figures/",h1s[cycle],"_6_7.png"),
      width = 10, height = 8, units = "in",res=100)
  plot(sub[,1],main=h1s[cycle],type="l")
  for(i in 1:years){
    lines(sub[,i],col=sample(col_vector, 1))
  }
  dev.off()
}
for(cycle in 1:length(h2s)){
  print(h1s[cycle])
  load(paste0(wd,h2s[cycle]))
  LU<-as.data.frame(t(as.data.frame(var.res.arr)))
  colnames(LU)<-paste0("Replicate",seq(1,years))
  sub<-LU# [,sample(seq(1,100),50)]
  png(file=paste0(wd,"Figures/",h2s[cycle],"_6_7.png"),
      width = 10, height = 8, units = "in",res=100)
  plot(sub[,1],main=h2s[cycle],type="l")
  for(i in 1:years){
    lines(sub[,i],col=sample(col_vector, 1))
  }
  dev.off()
}
