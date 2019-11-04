#! /usr/bin/perl

$prev_gene_id="toto";

while(<>) {
  if(/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/) {
    $chr_=$1;
    $start_=$4;
    $stop_=$5;
    $strand_=$7;
    
    $gene_id_="NA";
    if($_ =~ /gene_id "([^"]+)";/) {
      $gene_id_=$1;
    }
    $gene_name_="NA";
    if($_ =~ /gene_name "([^"]+)";/) {
      $gene_name_=$1;
    }
    
#     print $chr_."\t".$start_."\t".$stop_."\t".$strand_."\t".$gene_id_."\n";
    if($prev_gene_id ne "toto") {
      if($gene_id_ ne $prev_gene_id) {
	print $chr."\t".$start."\t".$stop."\t".$gene_id."\t0\t".$strand."\t".$gene_name."\n";
	$chr=$chr_;
	$start=$start_;
	$stop=$stop_;
	$strand=$strand_;
	$gene_id=$gene_id_;
	$gene_name=$gene_name_;
      }
    }
    else {
      $chr=$chr_;
      $start=$start_;
      $stop=$stop_;
      $strand=$strand_;
      $gene_id=$gene_id_;
      $gene_name=$gene_name_;
    }
    
    if($start_<$start) {
      $start=$start_;
    }
    if($stop_>$stop) {
      $stop=$stop_;
    }
    
    $prev_gene_id=$gene_id_;
  }
}
