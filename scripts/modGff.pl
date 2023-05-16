#!/usr/bin/perl

use strict;
use warnings;

my $usage = "\nUsage: perl modGff.pl <tag> <infile>\n\n";
$usage .= "\ttag: submission tag\n";
$usage .= "\tinfile: full path to gff3 file location; output will be written there\n";
$usage .= "\tdSfile: full path to file with dS values to label genes\n";
$usage .= "\te.g. perl modGff.pl LM07full_ annotation/FA_07/LM07LM07full/annotate_results/Saccharomyces_pastorianus_LM07.gff3 Scer_Seho_LM07full_dS.txt\n\n";

die $usage unless @ARGV > 1;

my $tag = $ARGV[0];
my $infile = $ARGV[1];
my $dSfile = $ARGV[2];
(my $inScffFile = $infile) =~ s/.gff3/.scaffolds.mod.fa/;
(my $inAnnFile = $infile) =~ s/.gff3/.annotations.txt/;
(my $outAnnFile = $infile) =~ s/.gff3/.annotations.mod.txt/;
(my $outGffFile = $infile) =~ s/.gff3/.mod.gff3/; 

my (%scffs_n_ids, %genes_n_IDs);
open my $in, '<', $inScffFile or die "$!";
while (<$in>){
 chomp;
 if ($_ =~ />/) {
  $_ =~ s/>//g;
  @_ = split(/\s/,$_);
  $scffs_n_ids{$_[0]} = $_[1];
 }
}

my $lab;
my %genes_n_lab;
open $in, '<', $dSfile or die "$!";
while (<$in>){
 if ($_ !~ /^Seho_/){
  @_ = split(/\s/,$_);
  $lab = "N";
  if ($_[2] < 0.25)                   {$lab = "A";} # ScerA
  if ($_[2] >= 0.25 && $_[2] <= 2.75) {$lab = "B";} # ScerB
  if ($_[2] > 2.75)                   {$lab = "C";} # SeubC
  $genes_n_lab{$_[1]} = $lab;
#  print "$_[1] $genes_n_lab{$_[1]} $lab\n";
 }
}

my $pre_ctg = "";
my $gene_cnt = 0;
my %genes_n_ids;
open my $out, '>', $outAnnFile or die "$!";
open $in, '<', $inAnnFile or die "$!";
while (<$in>){
 if ($_ =~ /^GeneID/) { print $out "$_"; }
 else {
  my @fields = split(/\t/,$_);
  if ($pre_ctg eq "")    {
   $pre_ctg = $fields[3];
  }
  if ($fields[3] eq $pre_ctg) {
   $pre_ctg = $fields[3]; my $tmp = $gene_cnt; $gene_cnt = $tmp + 10; 
  }
  if ($fields[3] ne $pre_ctg) {
   $pre_ctg = $fields[3]; $gene_cnt = 10;
  }

  if ($gene_cnt <  100) {my $tmp = $gene_cnt; $gene_cnt = "0000".$tmp;}
  if ($gene_cnt >= 100 && $gene_cnt < 1000) {
   my $tmp = $gene_cnt; $gene_cnt = "000".$tmp;
  }
  if ($gene_cnt >= 1000 && $gene_cnt < 10000) {
   my $tmp = $gene_cnt; $gene_cnt = "00".$tmp;
  }
  if ($gene_cnt >= 10000 && $gene_cnt < 100000) {
   my $tmp = $gene_cnt; $gene_cnt = "0".$tmp;
  }

  $genes_n_ids{$fields[0]} = $scffs_n_ids{$fields[3]}.$gene_cnt;
  my $tmp = $scffs_n_ids{$fields[3]}.$gene_cnt;
  if (exists $genes_n_lab{$fields[0]}) {
   $genes_n_ids{$fields[0]} .= $genes_n_lab{$fields[0]};
   $tmp .= $genes_n_lab{$fields[0]};
  }
  else {
   $genes_n_ids{$fields[0]} .= "N";
   $tmp .= "N";
  }
  $fields[0] = $tmp;
  $fields[1] = $tmp."-T1";

#  $genes_n_ids{$fields[1]} = $scffs_n_ids{$fields[3]}.$gene_cnt;
#  $fields[1] = $scffs_n_ids{$fields[3]}.$gene_cnt;

  $fields[3] = "$fields[3];$scffs_n_ids{$fields[3]}";
  
  print $out join("\t", @fields);

 }
}
close $in;
close $out;

open my $out2, '>', $outGffFile or die "$!";
open my $in2, '<', $infile or die "$!";
while (<$in2>){
 if ($_ =~ /^#/) { print $out2 $_; }
 else {
  my @fields = split(/\t/,$_);

  $fields[0] = $scffs_n_ids{$fields[0]};

  my @IDs = split (/;/,$fields[8]);
  (my $ID = $IDs[0]) =~ m/ID=\w+;/;
  my @tmp = split (/-/,$ID);
  (my $gene = $tmp[0]) =~ s/ID=//;
  $_ =~ s/$gene/$genes_n_ids{$gene}/g; 
  print $out2 $_;
 }
}
close $in2; close $out2;

