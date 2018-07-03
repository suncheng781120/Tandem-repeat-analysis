#!/usr/local/bin/perl-w
use strict;
#This script was used to summarize the profile, like copy number, repeat period length distribution of TR loci.
my $usage = "perl  this_script  Formated_TRFloci";
my $input = shift or die "$usage\n";
open IN, "<$input" or die;
#Chromosome
#NC_015762.1     1       31      2       15      CT
#NC_015762.1	150951	150984	9	3.777778	CTTTTCCCT
my %RepeatUnit = ();#This hash was used to keep the repeat unit, along with its frequency.
my %RepeatUnitLength = ();
while(<IN>){
      chomp;
      if(!/Chromosome/){
         my @array = split /\s+/;
         $RepeatUnit{$array[5]}++;
         $RepeatUnitLength{$array[3]}++;
      }
}
close IN;
open OUT1, ">$input.RepeatUnitDis.txt" or die;
print OUT1 "RepeatUnit\tNumber\n";
foreach(keys %RepeatUnit){
	print OUT1 "$_\t$RepeatUnit{$_}\n";
}
close OUT1;
open OUT2, ">$input.RepeatUnitLengthDis.txt" or die;
print OUT2 "RepeatUnitLength\tNumber\n";
foreach(keys  %RepeatUnitLength){
	print OUT2 "$_\t$RepeatUnitLength{$_}\n";
}
close OUT2;


