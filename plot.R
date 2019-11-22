#pdf(file = paste("PanCan_nooutlier","pdf",sep = "."),height = 6,width = 8)
pdf(file = paste("Firehose_nooutlier","pdf",sep = "."),height = 6,width = 8)
#setwd("/home/yehui/project/venn_gen/shi_genelist/tpye_tumor/TMB_coding/PanCan")
setwd("/home/yehui/project/venn_gen/shi_genelist/tpye_tumor/TMB_coding/Firehose")
myfiles <- Sys.glob("*.csv")
for(i in myfiles){
  #print(i)
  library(ggplot2)
  TMB_coding <- read.table(i,header = T,sep = ",")
  #print(TMB_coding)
  #TMB <- read.table(paste("/home/yehui/project/venn_gen/shi_genelist/tpye_tumor/TMB/PanCan",i,sep = "/"),header = T,sep = ",")
  TMB <- read.table(paste("/home/yehui/project/venn_gen/shi_genelist/tpye_tumor/TMB/Firehose",i,sep = "/"),header = T,sep = ",")
  #print(TMB)
  ff <- merge(TMB_coding,TMB,by="patient")
  q1 <- quantile(ff$TMB.x, 0.95)[[1]]
  q2 <- quantile(ff$TMB.y, 0.95)[[1]]
  ff$TMB.x <- ifelse(ff$TMB.x >q1,NA,ff$TMB.x)
  ff$TMB.y <- ifelse(ff$TMB.y >q2,NA,ff$TMB.y)
  #print(head(ff))
  n <- sub(pattern = ".csv",replacement = "",i)
  #dat.lm <- lm(ff[,3] ~ ff[,2], data = ff)
  #r2 <- sprintf("italic(R^2) == %.2f",summary(dat.lm)$r.squared)
  gg <- na.omit(ff)
  r2 <- sprintf("R^2 = %.4f",cor(gg[,3],gg[,2]))
  #print(r2)
  p <- ggplot(ff,aes(x=ff[,3],y=ff[,2])) + geom_point(shape=19,color="blue") + xlab("WES_TMB") + ylab("Coding_TMB") + geom_smooth(method = lm) + 
    #geom_abline(intercept = coef(dat.lm)[1],slope = coef(dat.lm)[2]) + 
    geom_text(aes(x = max(gg[,3])/2,y=max(gg[,2])/2,label=r2),size = 6) + labs(title = n) + theme(plot.title = element_text(hjust = 0.5))
  #pdf(file = paste(n,"pdf",sep = "."),height = 6,width = 8)
  plot(p)
  #dev.off()
}
dev.off()