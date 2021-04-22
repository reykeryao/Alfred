args = commandArgs(trailingOnly=TRUE)
dat<-read.table(args[1],col.names=c("Seq","Times"))
mat_seq<-"GCCGCTTCAGAGAGAAATCGC"
core_seq<-"CTTCAGAGAGAAA"
dat$Seq<-as.character(dat$Seq)
dat$ins_all<-0
dat$del_all<-0
dat$sub_all<-0
for (i in 1:dim(dat)[1]){
  dist<-drop(attr(adist(mat_seq,dat$Seq[i], count = TRUE), "counts"))
  dat$ins_all[i]<-dist[1]
  dat$del_all[i]<-dist[2]
  dat$sub_all[i]<-dist[3]
  }
dat$Len_all<-nchar(dat$Seq)
dat$Seq_core<-sapply(dat$Seq,FUN=function(x){
  unlist(regmatches(x,aregexec(core_seq,x,max=list(all = 4 ,del = 2, ins = 2, sub = 3))))})
dat$Seq_core<-as.character(dat$Seq_core)
dat$ins_core<-0
dat$del_core<-0
dat$sub_core<-0
for (i in 1:dim(dat)[1]){
  dist<-drop(attr(adist(core_seq,dat$Seq_core[i], count = TRUE), "counts"))
  dat$ins_core[i]<-dist[1]
  dat$del_core[i]<-dist[2]
  dat$sub_core[i]<-dist[3]
}
dat$Len_core<-nchar(dat$Seq_core)
write.table(dat,args[2],quote=F,sep="\t",row.names=F)
