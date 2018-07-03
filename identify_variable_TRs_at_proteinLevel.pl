#!/usr/local/bin/perl-w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
#This script was used to identify varible and conserved Tandem repeats between two species based on pairwise alignment done by BLASTn.
my $usage = "perl  this_script  Bterrestris_genesHarboringVTRs";
my $input = shift or die "$usage\n";
open IN, "<$input" or die "$usage\n";
open OUT, ">$input.variableTRsAA.txt" or die;
open OUTprotein, ">$input.VTRHarboringProtein.fa";
my $total_TRs = 0;
my $total_good_hits = 0;
my $total_variable_TRs = 0;
###For each VTR-contained protein sequence, use it as query to do BLASTp against Bimpatiens protein database.
my $inseq = Bio::SeqIO -> new(-file => "<$input", -format => "fasta");
LOOP:while(my $seq = $inseq->next_seq){
      $total_TRs++;
      my $seq_length = $seq->length;
      my $name = $seq->id;
      my $sequences = $seq->seq;
      open QUERY, ">BLAST_seq.txt" or die;
      print QUERY ">$name\n$sequences\n";
      close QUERY;
      open BLASTOUTPUTFILE, ">blastOutputFile" or die "$!";
      my $blastSysString = "/mnt/cheng/software/ncbi-blast-2.4.0+/bin/blastp -db /mnt/cheng/blast_db/GCF_000188095.1_BIMP_2.0_protein -query BLAST_seq.txt -out blastOutputFile -evalue 0.0000000001";
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
                                    if(($qStart <= 5) and ($qEnd >= $seq_length - 5)){#Alignment can only miss 5 bp at both ends.
                                        $total_good_hits++;
                                        #if(($qEnd-$qStart-$sEnd+$sStart >= 1) or ($sEnd-$sStart-$qEnd+$qStart >= 1)){
                                            if(($qString =~ /-/) or ($sString =~ /-/)){
                                                $total_variable_TRs++;
                                                my $size_difference = $qEnd-$qStart-$sEnd+$sStart;
                                                print OUT "$name\t$size_difference\n";
                                                print OUTprotein ">$name\n$sequences\n";
                                                my $alnIO = Bio::AlignIO->new(-format => "msf", -file => ">>Varied.hsp_aa.msf"); 
                                                my $aln = $hsp->get_aln; 
                                                $alnIO->write_aln($aln);
                                                next LOOP;
                                            }
                                        #}
                                    }
                             }
                       }
                }
}
print OUT "TotalTRs\t$total_TRs\nTotalGoodHits\t$total_good_hits\nTotalVariableTRs\t$total_variable_TRs\n";
close IN;
close OUT;
system("rm BLAST_seq.txt blastOutputFile");
print "This script was done successfully!\n";
exit;            
