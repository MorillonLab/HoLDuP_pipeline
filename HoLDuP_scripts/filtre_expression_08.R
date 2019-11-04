#! /usr/bin/Rscript

args=commandArgs(trailingOnly=TRUE)

n_cell_lines=length(args)/3

count_file_paths=args[1:n_cell_lines]
group_cell_lines=args[(1:n_cell_lines)+n_cell_lines]
rpkm_file_paths=args[(1:n_cell_lines)+2*n_cell_lines]

list_count_tabs=lapply(count_file_paths,read.table,header=FALSE,as.is=TRUE,sep="\t")
n_reads_total=sapply(1:n_cell_lines,function(i) {sum(list_count_tabs[[i]][,2])})
list_rpkm_tabs=lapply(rpkm_file_paths,read.table,header=FALSE,as.is=TRUE,sep="\t")

is_above_10_reads_at_least_once=rep(FALSE,dim(list_count_tabs[[1]])[1])
for(i in 1:n_cell_lines) {
  is_above_10_reads_at_least_once=is_above_10_reads_at_least_once|(list_count_tabs[[i]][,2]>=10)
}
is_above_10_reads_at_least_once=is_above_10_reads_at_least_once[1:(length(is_above_10_reads_at_least_once)-5)]

unique_group_cell_lines=unique(group_cell_lines)

is_above_08_rpkm_in_two_replicates_at_least_once=rep(FALSE,dim(list_rpkm_tabs[[1]])[1])
for(i in 1:length(unique_group_cell_lines)) {
  count_above_08_rpkm_within_group=rep(0,dim(list_rpkm_tabs[[1]])[1])
  i_cell_lines_of_group=which(group_cell_lines==unique_group_cell_lines[i])
  for(j in i_cell_lines_of_group) {
    count_above_08_rpkm_within_group=count_above_08_rpkm_within_group+(list_rpkm_tabs[[i]][,2]>=0.8)
  }
  is_above_08_rpkm_in_two_replicates_at_least_once=is_above_08_rpkm_in_two_replicates_at_least_once|(count_above_08_rpkm_within_group>=2)
}

cat(paste(list_rpkm_tabs[[1]][is_above_10_reads_at_least_once&is_above_08_rpkm_in_two_replicates_at_least_once,1],collapse="\n"))
cat("\n")

