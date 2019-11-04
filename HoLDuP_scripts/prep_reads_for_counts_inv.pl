#! /usr/bin/perl

while(<>) {
  if(/XS:A:\+/) {
    chomp($_);
    @line=split('\t', $_);
    print($line[0].":P\t".'16'."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".'*'."\t".'0'."\t".'0'."\t".$line[9]."\t".$line[10]."\n");
  }
  if(/XS:A:-/) {
    chomp($_);
    @line=split('\t', $_);
    print($line[0].":M\t".'0'."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".'*'."\t".'0'."\t".'0'."\t".$line[9]."\t".$line[10]."\n");
  }
}
