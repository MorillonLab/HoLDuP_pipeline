#! /usr/bin/perl

$prev_gene_id="toto";
$n_exon=0;
my @exon_start;
my @exon_stop;

while(<>) {
  if(/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/) {
    $chr=$1;
    $source=$2;
    $type=$3;
    $start=$4;
    $stop=$5;
    $V6=$6;
    $strand=$7;
    $V8=$8;
    $V9=$9;
    if(($type eq "exon")) {
      if($V9 =~ /gene_id "([^"]+)";/) {
	$ID=$1;
	if(($ID ne $prev_gene_id)&&($prev_gene_id ne "toto")) {
	  my %nucl;
	  for($i=0;$i<$n_exon;$i++) {
	    for($j=$exon_start[$i];$j<=$exon_stop[$i];$j++) {
	      $nucl{$j}=1;
	    }
	  }
	  $exon_lengths= keys %nucl;
	  print($prev_gene_id."\t".$exon_lengths."\n");
	  undef %nucl;
	  @exon_start=();
	  @exon_stop=();
	  $n_exon=0;
	}
	$exon_start[$n_exon]=$start;
	$exon_stop[$n_exon]=$stop;
	$n_exon++;
	$prev_gene_id=$ID;
      }
    }
  }
}

my %nucl;
for($i=0;$i<$n_exon;$i++) {
  for($j=$exon_start[$i];$j<=$exon_stop[$i];$j++) {
    $nucl{$j}=1;
  }
}
$exon_lengths= keys %nucl;
print($prev_gene_id."\t".$exon_lengths."\n");
@exon_start=();
@exon_stop=();
$n_exon=0;

