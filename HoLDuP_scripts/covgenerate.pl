#!/usr/bin/perl
use POSIX
# ARGV[0]: input file
# ARGV[1]: read length
# ARGV[2]: strand +/-

open(FILE, $ARGV[0]);
my $start = -1;
my $end = 0;
my $len;
my $count;
my $num = 0;
my $id = 0;
my $maxnum = 0;
while (<FILE>) {
	chomp;
	$count = $num;
	$maxnum = $num;
	($chr, $pos, $num) = split("\t");
	if ($num > $maxnum) {
		$maxnum = $num;
	}
	if ($num == 0) {
		if ($count!=0) {
			$id = $id + 1;
			print "$chr\ttoto\t$start\t$start\t.\t$ARGV[2]\t.\tID=$id;";
        		$read = &max($count/$ARGV[1],$maxnum);
		        print "reads=$read;";
		        $cov = $count;
		        print "cov=$cov;\n";
		}
		$start = -1;
		next;
	}
	if ($start == -1) {
		$start = $pos;
		$len = 1;
	}
	else {
		$len = 2;
	}
	
	$end = $pos;
	print "$chr\t";
	print "contig\t";
	#print "unite_transcript\t";
	print "$start\t";
	$count = $count + $num;
	$id = $id + 1;
	while (<FILE>) {
		chomp; 	
		($chr, $pos, $num) = split("\t");
		last if ($num==0 || $pos!=$end+1);
		$len = $len + 1;
		$end = $pos;
		$count = $count + $num;		
		if ($num > $maxnum) {
			$maxnum = $num;
		}
	}
	if ($num==0) {
		$start = -1;
	}
	else {
		$start = $pos;
	}
	print "$end\t";
	print ".\t";
	print "$ARGV[2]\t";
	print ".\t";
	print "ID=$id;";
	$low = &max(ceil($count/($ARGV[1]+3)), $maxnum+floor($len/($ARGV[1]+3)-1));
	$high = floor($count/($ARGV[1]-3));

	$read = &max($count/$ARGV[1], $low);
	print "reads=$read;";
	print "$low $high;";
	$cov = $count/$len;
	print "cov=$cov;\n";
}
close(FILE);
exit;
 
sub max {
	if ($_[0] > $_[1]) {
		$_[0];
	}
	else {
		$_[1];
	}
}
