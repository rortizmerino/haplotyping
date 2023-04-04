#!/usr/bin/perl

use strict;
use warnings;

#use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::DB::SeqI;
use Bio::SeqIO;

my $usage = "\n Usage: compareAssemblies_ANI.pl <reference> <query>\n";
$usage .= " reference: fasta file; plotted along y axis\n";
$usage .= " query:     fasta file; plotted along x axis\n\n";

die $usage unless @ARGV > 0;

my $genome1 = $ARGV[0];
my $genome2 = $ARGV[1];

(my $ref_name = $genome1) =~ s/.fasta//g;
(my $qry_name = $genome2) =~ s/.fasta//g;
print "  $ref_name vs $qry_name\n";

@_ = split(/\//, $ref_name);
$ref_name = $_[$#_];

my @genomes = ($genome1, $genome2);
my ($inseq, $seq, $i, $out, $in);
$i = 0;

my $cmd;
$cmd = "nucmer --maxmatch $genome1 $genome2";
print "$cmd\n";
system($cmd);
$cmd = "mv out.delta $ref_name.comparison.out.delta";
print "$cmd\n";
system($cmd);
$cmd = "show-coords -c -L 1500 -r $ref_name.comparison.out.delta > $ref_name.comparison.coords";
print "$cmd\n";
system($cmd);

my %order;
foreach my $genome (@genomes){
 my $cnt=1;
 open $in, '<', $genome; 
 while (<$in>){
  if ($_ =~ /^>/){
   @_ = split(/\s/, $_);
   $_[0] =~ s/>//g;
   if ($cnt < 10) {$order{$_[0]} = "00"."$cnt";}
   elsif ($cnt >= 10 && $cnt < 100) {$order{$_[0]} = "0"."$cnt";}
   else {$order{$_[0]} = $cnt;}
   #print "$_[0] $cnt\n";
   $cnt++;
  }
 }
 close $in;
}

open $out, '>', "$ref_name.comparison.coords.txt";
print $out "S1 E1 S2 E2 PI CQ G1 G2 O1 O2\n";
open my $file1, '<', "$ref_name.comparison.coords";
while (<$file1>){
 chomp;
 if ($_ =~ /^\s*\d/){
  @_ = split(/\s+/, $_);
  #if ($_[13] > 50) { 
   print $out "$_[1] $_[2] $_[4] $_[5] $_[10] $_[13] $_[15] $_[16] $order{$_[15]} $order{$_[16]}\n";
  #}
 }
}
close $out;

foreach (@genomes){

 (my $out_file = $genomes[$i]) =~ s/.fasta/.coords/g;
 @_ = split(/\//, $out_file);
 $out_file = $_[$#_];
 open $out, '>', $out_file;
 #print $out "ID Len Ncount Nblocks Nblock_lens Nblock_coords\n";
 print $out "Order ID Len Ncount Nblocks Nblock_lens Nblock_coords\n";
 
 my $file = "$genomes[$i]";
 my $in_seq_db = Bio::DB::Fasta->new($file) ;#|| die "Could not index sequence database $!\n";
 $inseq = Bio::SeqIO->new(-file => "<$file", -format => "Fasta" );
 while ($seq = $inseq->next_seq) {
  
  my $ID = $seq->display_id;
  my $Len = $seq->length();
  my $seqstr = $seq->seq();
  my $desc = $seq->desc;

  print $out "$order{$ID} $ID $Len ";
  if ($seqstr =~ /N/){
   my $Ncount = ($seqstr =~ tr/N//);
   my $regex = 'N+';
   my @Nblocks_pos = all_match_positions($regex,$seqstr);
   my $Nblocks = @Nblocks_pos;
   print $out "$Ncount $Nblocks ";
   my (@Nblock_lens, @Nblock_coords);
   foreach my $pos(@Nblocks_pos){my $strt=@$pos[0]+1; my $stop=@$pos[1]+1;
    my $len;
    if ($stop==$strt) {$len = 1;}
    else {$len=$stop-$strt+1;}
    push @Nblock_lens, $len;
    push @Nblock_coords, "$stop-$strt"
   }
   print $out join(",", @Nblock_lens), " "; print $out join(",", @Nblock_coords), "\n";
# print gaps just in case
=pod  
   foreach my $pos(@Nblocks_pos){
    my $strt=@$pos[0]+1; my $stop=@$pos[1]+1; my $tmp_seq=$seq->subseq($strt,$stop); print  "$tmp_seq\n";
   }
=cut
  }
  else {print $out "0 0 0 0\n";}
 }
 $i++;
 close $out;
}

$cmd = "Rscript --vanilla ../../../scripts/DotPlot_ANI.R ";
$cmd.= "$qry_name ";
$cmd.= "$ref_name ";
$cmd.= "$qry_name"."_vs_"."$ref_name.ANI.jpeg";
print "$cmd\n";
system($cmd);

sub all_match_positions {
 my ($regex, $string) = @_;
 my @ret;
 while ($string =~ /($regex)/g) {
  push @ret, [(pos($string)-length $1),pos($string)-1];
 }
 return @ret
}


