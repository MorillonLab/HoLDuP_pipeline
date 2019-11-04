#! /usr/bin/perl

while(<>) {chomp($_);

	@tab=split('\t', $_);
	
print "$tab[0]\t$tab[2]\t$tab[3]\t$tab[1]\t0\t$tab[5]\n";
}
