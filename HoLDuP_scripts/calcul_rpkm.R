#! /usr/bin/Rscript

args=commandArgs(trailingOnly=TRUE)

n_cell_lines=(length(args)-1)/2

count_file_paths=args[1:n_cell_lines]
rpkm_file_paths=args[(1:n_cell_lines)+n_cell_lines]

list_count_tabs=lapply(count_file_paths,read.table,header=FALSE,as.is=TRUE,sep="\t")
n_reads_total=sapply(1:n_cell_lines,function(i) {sum(list_count_tabs[[i]][,2])})

transcript_lengths=read.table(args[length(args)],header=FALSE,as.is=TRUE,sep="\t")
transcript_lengths=transcript_lengths[order(transcript_lengths$V1),]

list_rpkm_tabs=lapply(1:n_cell_lines,function(i) {
  tmp=list_count_tabs[[i]]
  tmp=tmp[1:(dim(tmp)[1]-5),]
  tmp=tmp[order(tmp$V1),]
  tmp$V2=((10^9)*as.double(tmp$V2))/(as.double(n_reads_total[i])*as.double(transcript_lengths$V2))
  return(tmp)
})

ctl=lapply(1:n_cell_lines,function(i) {
  write.table(list_rpkm_tabs[[i]],file=rpkm_file_paths[i],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
})
