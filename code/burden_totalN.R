#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--casefile", action="store")
parser$add_argument("--casesize", action="store", type="integer")
parser$add_argument("--controlfile", action="store")
parser$add_argument("--controlsize", action="store", type="integer")
parser$add_argument("--outfile", action="store")


args <- parser$parse_args()

case.dat<-read.delim(args$casefile, header=T, stringsAsFactors=F, sep="\t")
names(case.dat)[1]<-"GENE"
control.dat<-read.delim(args$controlfile, header=T, stringsAsFactors=F, sep="\t")
names(control.dat)[1]<-"GENE"

dat<-merge(case.dat, control.dat, by="GENE", all.x=T, all.y=T)
dat[is.na(dat)]<-0

dat$P_DOM<-0
dat$P_REC<-0

dat$TOTAL_CASE<-0
dat$TOTAL_CONTROL <-0

for(i in 1:nrow(dat)){
  
  #Dominant model
  case_count<-dat[i,]$CASE_COUNT_HET+dat[i,]$CASE_COUNT_HOM
  control_count<-dat[i,]$CONTROL_COUNT_HET+dat[i,]$CONTROL_COUNT_HOM
  args_casesize <- dat[i,]$CASE_COUNT_HOMW+case_count
  args_controlsize<-round(dat[i,]$CONTROL_TOTAL_ANmean/2)

  if(args_casesize>args$casesize){
     args_casesize <- args$casesize
     dat[i,]$CASE_COUNT_HOMW <- args_casesize-case_count
  }
 

 
  if(case_count>args_casesize){
    case_count<-args_casesize
  }else if(case_count<0){
    case_count<-0
   }
  if(control_count>args_controlsize){
    control_count<-args_controlsize
  }else if(control_count<0){
    control_count<-0
   }

  if(args_casesize<=0){
    args_casesize <- args$casesize
  }
  if(args_controlsize<=0){
    args_controlsize <- args$controlsize
  }
  
  mat<-cbind(c(case_count, (args_casesize-case_count)), c(control_count, (args_controlsize-control_count)))
  dat[i,]$P_DOM<-fisher.test(mat, alternative="greater")$p.value
  
  
  #Recessive model
  case_count_rec<-dat[i,]$CASE_COUNT_CH+dat[i,]$CASE_COUNT_HOM
  control_count_rec<-dat[i,]$CONTROL_COUNT_HOM+(args_controlsize)*((dat[i,]$CONTROL_COUNT_HET)/(args_controlsize))^2
  
  if(control_count_rec<0){ control_count_rec<-0}
  if(case_count_rec>args_casesize){case_count_rec<-args_casesize}
  if(control_count_rec>args_controlsize){control_count_rec<-args_controlsize}
  control_count_rec<-round(control_count_rec, digits=0)
  
  mat_rec<-cbind(c(case_count_rec, (args_casesize-case_count_rec)), c(control_count_rec, (args_controlsize-control_count_rec)))
  dat[i,]$P_REC<-fisher.test(mat_rec, alternative="greater")$p.value
  dat[i,]$TOTAL_CASE <- args_casesize
  dat[i,]$TOTAL_CONTROL <- args_controlsize
}

write.table(x=dat[-9],file=args$outfile, sep="\t", quote=F, row.names=F, col.names=T)
