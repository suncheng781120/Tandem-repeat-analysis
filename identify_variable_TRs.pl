#!/usr/local/bin/perl-w
use strict;
use Bio::DB::Fasta;
use Bio::SearchIO;
use Bio::AlignIO;
#This script was used to identify varible and conserved Tandem repeats between two species based on pairwise alignments done by BLASTn.
my $usage = "perl  this_script  Formated_TRFinder_outfile   Bterrestris_softMasked_GenomeSeq";
my $input = shift or die "$usage\n";
my $genome = shift or die "$usage\n";
open IN, "<$input" or die "$usage\n";
open OUT, ">$input.variableTRs.txt" or die;
open OUT2, ">$input.poor-variableTRs.txt" or die;
open OUT3, ">$input.conservedTRs.txt" or die;
my $total_TRs = 0;
my $total_eligible_flank = 0;
my $total_good_hits = 0;
my $total_variable_TRs = 0;
my $total_conserved_TRs = 0;
my $poor_variable_TRs = 0;
###For each TR, extract 100bp of flank; check if there is any TR and N and soft_masked letters in its flanking;
LOOP:while(<IN>){
     chomp;
     if(!/Chromosome/){
        $total_TRs++;
        my @array = split /\s+/;
        my $array_length = $array[2]-$array[1]+1;
        my $db = Bio::DB::Fasta->new("$genome");
        my $seq = $db->get_Seq_by_id("$array[0]");
        my $length = $seq->length;
        my $start = $array[1] - 100;
        my $stop = $array[2] + 100;
        if(($start > 0) and ($stop < $length)){#Make sure coordinates are within sequence range.
            my $BLAST_seq = $seq->subseq($start=>$stop);
            my $up = $array[1]-1;
            my $flank_up = $seq->subseq($start=>$up);
            my $down = $array[2] + 1 ;
            my $flank_down = $seq->subseq($down=>$stop);
            my $flank_seq = "$flank_up"."$flank_down";##Concaternate up- and down-stream flanking sequences.
            if(($flank_seq !~ /N+/) and ($flank_seq !~ /[atcg]{10,}/)){###A total of 10 contignous low-complexity nucleotides at most will be allowed. 
                $total_eligible_flank++;
#Whenever there is not much soft-masked letters and gaps, the TRs, along with their flank, will be used as query to BLASTn against another species.
                open QUERY, ">BLAST_seq.txt" or die;
                print QUERY ">$array[0]:$start-$stop\n$BLAST_seq\n";
                close QUERY;
                open BLASTOUTPUTFILE, ">blastOutputFile" or die "$!";
                my $blastSysString = "/mnt/cheng/software/ncbi-blast-2.4.0+/bin/blastn -task blastn -dust no -db /mnt/cheng/blast_db/GCF_000188095.1_BIMP_2.0_genomic.fna -query BLAST_seq.txt -out blastOutputFile -evalue 0.00001";
                system($blastSysString);
                my $blastIn = new Bio::SearchIO(-format => 'blast', -file => 'blastOutputFile');
                while(my $result = $blastIn->next_result ){
                      while (my $hit = $result->next_hit ){
                             while (my $hsp = $hit->next_hsp){
                                    my $qStart = $hsp->start('query');
                                    my $qEnd = $hsp->end('query');
                                    my $qLength = $hsp->length('query');
                                    my $sStart = $hsp->start('hit');
                                    my $sEnd = $hsp->end('hit');
                                    my $qString = $hsp->query_string;
                                    my $sString = $hsp->hit_string;
####Check significant hits, and only hits with longer enough flank sequences will be considered for downstream analysis.
                                    if(($qStart <= 5) and ($qEnd >= $array_length + 200 - 5)){#Alignment can only miss 5 bp at both ends.
                                        $total_good_hits++;
                                        #if(($qEnd-$qStart-$sEnd+$sStart >= $array[3]) or ($sEnd-$sStart-$qEnd+$qStart >= $array[3])){
                                        #    $total_variable_TRs++;
                                        #    my $size_difference = $qEnd-$qStart-$sEnd+$sStart;
                                        #    print OUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$size_difference\n";
                                        #    my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Varied.hsp.msf"); 
                                        #    my $aln = $hsp->get_aln; 
                                        #    $alnIO->write_aln($aln);
                                        #    next LOOP;
                                        #}
                                        #else {
###Check the HSP and explore the TRs region to identify the presence of variable TRs.
                                        my @query = split //, $qString;
                                        my @hit = split //, $sString;
                                        my $position = $qStart;
                                        my $position_new = 0;
                                        my $query_TR = '';
                                        my $hit_TR = '';
                                        LINE:for(my $i = $qStart; $i < $qEnd; $i++){
                                            my $tem_query = shift @query;
                                            my $tem_hit = shift @hit;
                                            if($position <= 100){
                                               if($tem_query ne '-'){                                                      
                                                  $position++;
                                                  next LINE;
                                               }
                                            }
                                            if(($position > 100) and ($position_new <= $array_length)){
                                                   $query_TR .= $tem_query;
                                                   $hit_TR .= $tem_hit;
                                                   if($tem_query ne '-'){
                                                      $position_new++;
                                                   }
                                            }
                                        }
                                        my @query_TR = split //,$query_TR;
                                        my @hit_TR = split //, $hit_TR;
                                        my $query_ins = 0;
                                        my $hit_ins = 0;
                                        foreach (@query_TR){
                                                 if($_ eq '-'){
                                                    $query_ins++;
                                                 }
                                        }
                                        foreach(@hit_TR){
                                                if($_ eq '-'){
                                                   $hit_ins++;
                                                }
                                        }
                                        if(($query_ins >= 2) and ($hit_ins == 0)){
                                            my $size_diff = $hit_ins-$query_ins;
                                            print OUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$size_diff\n"; 
                                            my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Varied.hsp.msf");
                                            my $aln = $hsp->get_aln;
                                            $alnIO->write_aln($aln);
                                            $total_variable_TRs++;
                                            next LOOP;
                                        }
                                        elsif(($query_ins == 0) and ($hit_ins >= 2)){
                                                my $size_diff = $hit_ins-$query_ins;
                                                print OUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$size_diff\n"; 
                                                my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Varied.hsp.msf");
                                                my $aln = $hsp->get_aln;
                                                $alnIO->write_aln($aln);
                                                $total_variable_TRs++;
                                                next LOOP;
                                        }
                                        elsif(($query_ins >= 2) and ($hit_ins >= 2)){
                                                my $size_diff = $hit_ins-$query_ins;
                                                print OUT2 "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$size_diff\n";
                                                my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Varied.poorhsp.msf");
                                                my $aln = $hsp->get_aln;
                                                $alnIO->write_aln($aln);
                                                $poor_variable_TRs++;
                                                next LOOP;
                                        }
                                        elsif(($query_ins == 0) and ($hit_ins == 0)){
                                                my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Conserved.hsp.msf");
                                                my $aln = $hsp->get_aln;
                                                print OUT3 "$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\n"; 
                                                $alnIO->write_aln($aln);
                                                $total_conserved_TRs++;
                                                next LOOP;
                                        }
                                    }
                             }
                       }
                }
             }
          }
       }
}
print OUT "TotalTRs\t$total_TRs\nTotalTRsWithGoodFlank\t$total_eligible_flank\nTotalGoodHits\t$total_good_hits\nTotalVariableTRs\t$total_variable_TRs\nPoorVariable\t$poor_variable_TRs\nTotalConservedTRs\t$total_conserved_TRs\n";
close IN;
close OUT;
#close OUT2;
close OUT3;
system("rm BLAST_seq.txt blastOutputFile");
print "This script was done successfully!\n";
exit;            
