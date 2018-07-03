#!/usr/local/bin/perl-w
use strict;
use warnings;
use Bio::DB::Fasta;
my $usage = "perl Genome-wide_Development_of_SSR_Molecular_Markers_v1.0.pl input_file\ninput_file, a Fasta file that will be mined for TR loci\nPlease run this program under the primer3-x.x.x/src folder\n";
my $input = shift or die "$usage\n";
###First, run Tandem Repeat Finder on target sequences.
######################################################
#Implementation Details:
#File = sequences input file
#Match  = matching weight
#Mismatch  = mismatching penalty
#Delta = indel penalty
#PM = match probability (whole number)
#PI = indel probability (whole number)
#Minscore = minimum alignment score to report
#MaxPeriod = maximum period size to report
#[options] = one or more of the following :
#          -m    masked sequence file
#          -f    flanking sequence
#          -d    data file
#          -h    suppress html output
#          -r    no redundancy elimination
# bless( {
#          'matchWeight' => 2,
#          'mismatchPenalty' => 7,
#          'delta' => 7,
#          'pm' => 80,
#          'pi' => 10,
#          'minScore' => 50,
#          'maxPeriod' => 500,
#          'maskedOutputFile' => '/users/bob/seq.fasta.masked',
#          'pathToEngine' => '/usr/local/bin/trf',
#          'workDir' => '/tmp',
#          'sequenceFile' => '/users/bob/seq.fasta'
#        },
#      'TRF' );
#####################################################
system("/mnt/cheng/software/trf $input 2 7 7 80 10 50 500 -h");
###Second, extract the coordinates of eligible TRs.
###Here we only keep TR loci that have a period between 2 and 12, with a total length of 30 or more.
my %SSRs =();#This hash was used to hold all the eligible TRs.
open SSR, "<$input.2.7.7.80.10.50.500.dat" or die;
open POSITION, ">coordinates.txt" or die;
my $key = '';
while(<SSR>){
      chomp;      
      if(/^Sequence:\s+(\S+)/){
         $key = $1;
      }
      elsif(/^(\d+)\s+(\d+)\s+(\d+)/){
            if(($3 > 1) and ($3 <= 12) and ($2 - $1 >= 30)){###Refine repeat unit between 2 and 12 and the loci length >= 30 bp.
               print POSITION "$key\t$1\t$2\n";
            }
      }
}
close SSR;
close POSITION;
###Third, extract 200 bp of flanking sequences of TRs, and design primers for those TR loci. 
my $db = Bio::DB::Fasta->new("$input");
open IN, "<coordinates.txt" or die;
open OUT, ">$input.SSRprimers.txt" or die;
no strict 'refs'; 
my $i = 1;
while(<IN>){
       chomp;
       my @array = split /\s+/;
       my $seq = $db->get_Seq_by_id("$array[0]");
       my $length = $seq->length;
       my $start = $array[1] - 200;
       my $stop = $array[2] + 200;
       my $ssr_length = $array[2] - $array[1];
       if(($start > 0) and ($stop < $length)){
           my $subseq = $seq->subseq($start=>$stop);
           ###Prepare primer3 software input file. Format is taken from default parameters (see Manual of Primer3).
           open Prime3input, ">primer3input.txt" or die;
           print Prime3input "SEQUENCE_ID=$array[0].$array[1].$array[2]\n";
           print Prime3input "SEQUENCE_TEMPLATE=$subseq\n";
           print Prime3input "PRIMER_TASK=pick_detection_primers\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\n";
           print Prime3input "SEQUENCE_TARGET=201,$ssr_length\n";#Begin from the start of SSR loci, and do not lie within this loci.
           print Prime3input "PRIMER_NUM_RETURN=5\nPRIMER_MAX_NS_ACCEPTED=1\n";
           print Prime3input "PRIMER_MIN_SIZE=18\nPRIMER_OPT_SIZE=20\nPRIMER_MAX_SIZE=27\n";
           print Prime3input "PRIMER_MIN_GC=20.0\nPRIMER_MAX_GC=80.0\nPRIMER_MIN_TM=55.0\nPRIMER_OPT_TM=60.0\nRIMER_MAX_TM=63.0\n";
           print Prime3input "=\n";
           close Prime3input;
           ###Run primer3.
           system("./primer3_core -output=primers_output primer3input.txt");
           ###Extract the best primer pair from primer3 output.
           my $handlename = 'primerout'."$i";
           open $handlename, "<primers_output" or die;
           print OUT ">$array[0]:$array[1]-$array[2]\n";
           while(<$handlename>){
                 chomp;
                 if((/PRIMER_LEFT_0_SEQUENCE/) or (/PRIMER_RIGHT_0_SEQUENCE/)){
                     print OUT "$_\n";
                 }
           }
           $i++;
           close $handlename;
       }
}
close IN;
close OUT;
system("rm ./coordinates.txt ./primers_output ./primer3input.txt $input.index");
print "This script was done successfully!\n";
exit;
