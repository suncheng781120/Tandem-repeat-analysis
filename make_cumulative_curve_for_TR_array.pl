#!/usr/local/bin/perl-w
use strict;
#This script was used to make cumulative curve for the length of TR array.
my $usage = "perl  this_script   TRFlengthDistribution";
my $input = shift or die "$usage\n";
open IN, "<$input" or die;
open OUT, ">$input.CumulativeData.txt" or die;
print OUT "LociLength\tNumber\tCumulative\n";
my $totalNumber = 0;
my %length = ();
while(<IN>){
      chomp;
      if(/(\d+)\s+(\d+)/){
	  $totalNumber += $2;
	  $length{$1} = $2;
      }
}
close IN;
my $cumulative = 0;
foreach (sort {$a <=> $b} keys %length){
         $cumulative += $length{$_};
	 my $rate = $cumulative/$totalNumber;
	 print OUT "$_\t$length{$_}\t$rate\n";
}
close OUT;
exit;
      
