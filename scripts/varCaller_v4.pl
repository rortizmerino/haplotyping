#!/usr/bin/perl

use strict;

my $usage = "\nUsage: varCaller_v4.pl <species> <reference> <path lib> <path run> <threads> <run>\n";
$usage .= " species:   identifier for logfile, readgroups, and final joint vcf (e.g. Kmar) \n";
$usage .= " reference: fasta file to create indexes and use for mapping\n";
$usage .= " path lib:  path to fastq files\n";
$usage .= " path run:  path to run directory\n";
$usage .= " threads:   to be used for multithreaded steps\n";
$usage .= " run:       0, print commands to log and exit; 1, run commands and print to log\n\n";

die $usage unless @ARGV > 0;

my $cmd;

my $sp        = $ARGV[0];
my $ref_fa    = $ARGV[1];
my $path_lib  = $ARGV[2];
my $path_run  = $ARGV[3];
my $threads   = $ARGV[4];
my $clean_run = $ARGV[5];

(my $ref_name = $ref_fa) =~ s/\.fa//g;

open (LOG, '>', "varCaller_v4.".$sp.".log");
print LOG "cmd: varCaller_v4.pl $sp $ref_fa $path_lib $path_run $threads $clean_run\n";

print LOG "\n\n\n------>\n";
# initial index files
$cmd = "samtools faidx $ref_fa";
if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
$cmd = "bwa index $ref_fa";
if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
$cmd = "gatk ";
$cmd.= "CreateSequenceDictionary --REFERENCE $ref_fa --OUTPUT $ref_name.dict";
if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n\n";

my %files;
opendir(my $dh, $path_lib) || die "Can't open $path_lib: $!";
while (readdir $dh) {
 if ($_ =~ /_1/){
  my $name1 = $_; 
  (my $name2 = $name1) =~ s/_1/_2/g; 
  $files{$name1} = $name2;
 }
 if ($_ =~ /#/){
  my $name1 = $_;
  $files{$name1} = "";
 }
}

my $lib_count = 0;
my ($file1, $file2, $filename1, $filename2, @strains);
foreach $file1 (sort keys %files){
$lib_count++;
print LOG "\n\n\n-->$lib_count\n";
 @_ = split(/\//, $file1);
 # remove path here and below
 $filename1 = $_[$#_];
 $filename1 =~ s/_1.fq/_$sp/g;
 $file2 = $files{$file1};
 @_ = split(/\//, $file2);
 $filename2 = $_[$#_];
 
 if ($filename2 eq "") {
  $cmd = "bwa mem -t $threads $ref_fa $path_lib"."$file1 > $filename1.sam";
  print LOG "single\n\n";
 }
 else {
  $cmd = "bwa mem -t $threads $ref_fa $path_lib"."$file1 $path_lib"."$file2 > $filename1.sam";
  print LOG "paired\n\n";
 }
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

 push (@strains, $filename1);

 print LOG "\n### rest of pipe from sam\n";
 $cmd = "samtools view -bSh $filename1.sam | samtools sort -T $sp -O bam -o $filename1"."_s.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "rm $filename1.sam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "samtools index $filename1"."_s.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# filtering step, remove unmapped reads why 4?
 $cmd = "samtools view -bF 4 $filename1"."_s.bam > $filename1"."_s_f.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "rm $filename1"."_s.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# add read groups
 $cmd = "gatk AddOrReplaceReadGroups --INPUT $path_run"."$filename1"."_s_f.bam --OUTPUT $path_run"."$filename1"."_s_f_RGs.bam --RGID ".$filename1."_run --RGLB lib1 --RGPL Illumina --RGPU 1 --RGSM ".$filename1." --VALIDATION_STRINGENCY SILENT";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "samtools index $filename1"."_s_f_RGs.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "rm $filename1"."_s_f.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# mark and remove PCR duplicates
 $cmd = "gatk MarkDuplicates --INPUT $path_run"."$filename1"."_s_f_RGs.bam --OUTPUT $path_run"."$filename1"."_RD.bam --METRICS_FILE $path_run"."$filename1"."_dupl_metrics.txt --REMOVE_SEQUENCING_DUPLICATES true";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "samtools index $filename1"."_RD.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
 $cmd = "rm $filename1"."_s_f_RGs.bam";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# calculate intervals and recalibrate bam would go here but seems gatk4 does not need it

# generate interval list (no actual intervals, full seqs in fasta); avoid telos and rDNA here?
 $cmd = "awk \'BEGIN{RS=\">\"}{if (NR > 1) print \$1}\' $ref_fa > $ref_name.intervals";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# calculate coverage, ran here to avoid duplicates in final coverage
 $cmd = "gatk DepthOfCoverage --reference $ref_fa --input $filename1"."_RD.bam --output $filename1"."_coverage --intervals $ref_name.intervals --output-format TABLE";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# call SNPs, single sample VCF, emit & call ?
 $cmd = "gatk HaplotypeCaller --reference $ref_fa --input $filename1"."_RD.bam --standard-min-confidence-threshold-for-calling 30.0 --min-base-quality-score 20 --output $filename1.vcf";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";

# call SNPs, generate GVCFs 
 $cmd = "gatk HaplotypeCaller --reference $ref_fa --input $filename1"."_RD.bam --standard-min-confidence-threshold-for-calling 30.0 --min-base-quality-score 20 --emit-ref-confidence GVCF --output $filename1.g.vcf";
 if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
}
=pod
print LOG "\n\n\n-->joint genotyping\n";
$cmd = "java -jar /software/shared/apps/x86_64/GATK/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $ref_fa";
foreach (@strains){ $cmd .= " --variant $_".".g.vcf"; }
$cmd .= " -o joint$sp.vcf";
if($clean_run == 1) {system ($cmd)}; print LOG "$cmd\n";
=cut
close (LOG);
exit;
