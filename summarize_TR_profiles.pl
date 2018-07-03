#!/usr/local/bin/perl-w
use strict;
#This script was used to summarize the profile, like copy number, length repeat period distribution of TR loci.
my $usage = "perl  this_script  Formated_TRFloci";
my $input = shift or die "$usage\n";
open IN, "<$input" or die;
#Chromosome
#NC_015762.1     1       31      2       15      CT
my %RepeatUnit = ();#This hash was used to keep the repeat unit, along with its frequency.
my %RepeatUnitLength = ();
my %Length = ();#This hash was used to keep the length distribution.
my %copyNumber = ();#This hash was used to keep the copy number distribution.
while(<IN>){
      chomp;
      if(!/Chromosome/){
         my @array = split /\s+/;
         $RepeatUnit{$array[5]}++;
         $RepeatUnitLength{$array[3]}++;
         my $length = $array[2]-$array[1]+1;
         $Length{$length}++;
         $copyNumber{$array[4]}++;
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
open OUT3, ">$input.LengthDis.txt" or die;
print OUT3 "ArrayLength\tNumber\n";
foreach (keys %Length){
	print OUT3 "$_\t$Length{$_}\n";
}
close OUT3;
open OUT4, ">$input.CopyNumber.txt" or die;
 print OUT4 "InternalCopyNumber\tNumber\n";
foreach (keys  %copyNumber){
	print OUT4 "$_\t$copyNumber{$_}\n";
}
close OUT4;


