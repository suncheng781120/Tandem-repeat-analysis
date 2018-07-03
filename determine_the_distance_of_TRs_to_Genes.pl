#!/usr/local/bin/perl-w
use strict;
#This script was used to determine the distance betwween TRs and genes.
my $usage = "perl  this_script.pl  Gene_gff_file  TR_coordinates";
my $gff = shift or die "$usage\n";
my $SSR = shift or die "$usage\n";
my %gene = ();#This hash was used to keep gene coordinates.
open GFF, "<$gff" or die;
while(<GFF>){
      chomp;
      #NC_015762.1     Gnomon  gene    2279    18548
      if(/^(\S+)\s+\S+\s+gene\s+(\d+)\s+(\d+)/){
	   my $tem_string = "$2".'COOR'."$3";
	   $gene{$1} .= "$tem_string".'SEPRATE';
      }
}
close GFF;
open SSR, "<$SSR" or die;
open OUT, ">$SSR.GeneDistance.txt" or die;
print OUT "Distance2gene\tTR_locusNumber\n";
my %distance = ();
my $totalInputSSR = 0;
my $SSRonNW = 0;
LINE:while(<SSR>){
      chomp;
      #Chromosome Start	End PeriodLength CopyNumber
      if(!/Chromosome/){
	 $totalInputSSR++;
	 if(/^NW/){
	    $SSRonNW++;
	 }
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
				 $distance{$last_distance}++;
				 next LINE;
		            }
			    if($middle < $position[0]){
			       my $tem_distance = $position[0]-$middle;
			       if($tem_distance < $last_distance){
				  $last_distance = $tem_distance;
			       }
			    }
			    if($middle > $position[1]){
	  		       my $tem_distance = $middle - $position[1];
			       if($tem_distance < $last_distance){
				  $last_distance = $tem_distance;
		               }
			    }
	            }
		    $distance{$last_distance}++;
		    next LINE;
		 }
        }
     }
}
close SSR;
my %lengthDis = ();
foreach (sort {$a <=> $b} keys %distance){
	if($_ == 0){
           $lengthDis{genic_region} = $distance{$_};
	}
	if(($_ > 0) and ($_ <= 500)){
	    $lengthDis{500} += $distance{$_};
	}
	if(($_ > 500) and ($_ <= 1000)){
             $lengthDis{1000} += $distance{$_};
        }
        if(($_ > 1000) and ($_ <= 2000)){
	    $lengthDis{2000} += $distance{$_};
        }
        if($_ > 2000){
           $lengthDis{2001} += $distance{$_};
        }
}
print OUT "Total_TR $totalInputSSR\tTRonNW $SSRonNW\n";
foreach (sort {$a <=> $b} keys %lengthDis){
	print OUT "$_\t$lengthDis{$_}\n";
}
close OUT;
print "This script was done successfully!\n";

