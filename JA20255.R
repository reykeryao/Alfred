args = commandArgs(trailingOnly=TRUE)
con_r <- gzfile(args[1],open="r")
con_w <- gzfile(args[2],open="w")
# 3' NNN is not included, but will be analized in the following procedure

mat_seq<-"GCCGCTTCAGAGAGAAATCGC"

# process reads
# output is reads ID, Sequence (form previous processing if multiple TS), matched Templates, 5'sequences, 3' sequences, and FLAG
# multiple TS will out put multiple lines with the same ID
# maximum allowed ins/del is 2, substitution is 3, all errors is <=4
# FLAG is grouped as "TS" followed by a number, 
# indicates the round of TS for that read (1 is the rightmost one from the 5' end, larger number towards the 3'end)
# "TSL" means the length is long enough (>21 nt), but doesn't contain any template sequences
# "TSS" means the length is too short (< 22 nt) to be able to get any useful information, after remove template at least once
# "TS0" means the length is too short (< 22 nt), even before the TS, it's from primer-dimer
# in TSL and TSS condition, only ID and Sequene and FLAG is ouput, other field is empty

while (length(line <- readLines(con_r, n = 1, warn = FALSE)) > 0) {
  ID<-unlist(strsplit(line,split="\t"))[1]
  Seq<-unlist(strsplit(line,split="\t"))[2]
  n_TS<-0
  FLAG=""
  while ((nchar(Seq)>21) & FLAG!="TSL") {
    Pos<-aregexec(mat_seq,substr(Seq,1,30),max=list(all = 4 ,del = 2, ins = 2, sub = 3))
    if (unlist(Pos)[1]>=0){
      n_TS <- n_TS+1
      TSM<-unlist(regmatches(Seq, Pos))
      Seq_split<-unlist(strsplit(substr(Seq,1,30),TSM))
      if (length(Seq_split)==1){
        end5=""
        end3=paste0(Seq_split[1],substr(Seq,31,150))
      } else {
        end5=Seq_split[1]
        end3=paste0(Seq_split[2],substr(Seq,31,150))
      }
      FLAG<-paste0("TS",n_TS)
      out_line <- paste(ID,Seq,TSM,end5,end3,FLAG,collapse ="\t")
      writeLines(out_line,con_w)
      Seq<-end3
    } else {
      FLAG<-"TSL"
      out_line <- paste(ID,Seq,"","","",FLAG,collapse="\t")
      writeLines(out_line,con_w)
    }
  }
  if ((nchar(Seq)<=21) & n_TS==0 & FLAG=="") {
    FLAG<-"TS0"
    out_line <- paste(ID,Seq,"","","",FLAG,collapse="\t")
    writeLines(out_line,con_w)
  } else if ((nchar(Seq)<=21) & FLAG!="TSL" & n_TS!=0) {
    FLAG<-"TSS"
    out_line <- paste(ID,Seq,"","","",FLAG,collapse="\t")
    writeLines(out_line,con_w)
  }
} 

close(con_r)
close(con_w)

#calculate error rate
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

cat(sum(dat$Times*21),sum(dat$Times*dat$sub_all),sum(dat$Times*dat$ins_all),sum(dat$Times*dat$del_all),sep="\n")
cat(sum(dat$Times*13),sum(dat$Times*dat$sub_core),sum(dat$Times*dat$ins_core),sum(dat$Times*dat$del_core),sep="\n")

table(rep(dat$del_all,dat$Times))
table(rep(dat$sub_all,dat$Times))
table(rep(dat$ins_all,dat$Times))
table(rep(dat$del_core,dat$Times))
table(rep(dat$sub_core,dat$Times))
table(rep(dat$ins_core,dat$Times))