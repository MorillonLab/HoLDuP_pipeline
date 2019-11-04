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
    $gene_type_="NA";
    if($_ =~ /gene_type "([^"]+)";/) {
      $gene_type_=$1;
    }
    
    if($prev_gene_id ne "toto") {
      if($gene_id_ ne $prev_gene_id) {
	print $gene_id."\t".$gene_type."\n";
	$chr=$chr_;
	$start=$start_;
	$stop=$stop_;
	$strand=$strand_;
	$gene_id=$gene_id_;
	$gene_type=$gene_type_;
      }
    }
    else {
      $chr=$chr_;
      $start=$start_;
      $stop=$stop_;
      $strand=$strand_;
      $gene_id=$gene_id_;
      $gene_type=$gene_type_;
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
