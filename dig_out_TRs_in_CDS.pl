#!/usr/local/bin/perl-w
use strict;
#This script was used to calculation the number of TRs located in CDSs.
my $usage = "perl  this_script.pl  Gene_gff_file  TR_coordinates";
my $gff = shift or die "$usage\n";
my $SSR = shift or die "$usage\n";
my %gene = ();#This hash was used to keep gene coordinates.
open GFF, "<$gff" or die;
while(<GFF>){
      chomp;
      #NC_015762.1     Gnomon  gene    2279    18548 ID=cds0;Parent=rna0;Dbxref=GeneID:100649911,Genbank:XP_012171909.1;Name
      if(/^(\S+)\s+\S+\s+CDS\s+(\d+)\s+(\d+).+Genbank:(\w+\.\d+);Name/){
	   my $tem_string = "$2".'COOR'."$3".'COOR'."$4";
	   $gene{$1} .= "$tem_string".'SEPRATE';
      }
}
close GFF;
open SSR, "<$SSR" or die;
open OUT, ">$SSR.inCDS.txt" or die;
print OUT "Chr\tStart\tEnd\tUnitLength\tCopyNumber\tRepeatUnit\tSizeDiff\tProteinID\n";
my %distance = ();
LINE:while(<SSR>){
      chomp;
      #Chromosome Start	End PeriodLength CopyNumber
      if(!/Chromosome/){
	 my $last_distance = 100000000;#A number that longer than any bumblebee chromosome.
	 my @tem = split /\s+/;
	 my $middle = sprintf "%d", ($tem[1] + $tem[2])/2;#Here we use the middle position of TR locus as its location to determine its distance to genes.
	 foreach my $chr (keys %gene){
		 if($tem[0] eq $chr){
	            my @coordinates = split /SEPRATE/, $gene{$chr};
		    foreach my $coor (@coordinates){
			    my @position = split /COOR/,$coor;
                            if(($middle >= $position[0]) and ($middle <= $position[1])){
				 $last_distance = 0;
                                 print OUT "$_\t$position[2]\n";
				 $distance{$last_distance}++;
				 next LINE;
		            }
	            }
		 }
        }
     }
}
close SSR;
print OUT "Distance2CDS\tTRlocusNumber\n";
foreach (sort {$a <=> $b} keys %distance){
	print OUT "$_\t$distance{$_}\n";
}
close OUT;
print "This script was done successfully!\n";



