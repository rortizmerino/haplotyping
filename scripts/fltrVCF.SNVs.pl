#!/usr/bin/perl
use strict;

my $usage = "\nUsage: fltrVCF.SNVs.pl <tag> <path>\n";
$usage .= " tag:  to distinguish project, used in input files e.g. Scer_v0\n";
$usage .= " path: to vcf file\n\n";

die $usage unless @ARGV > 0;

my $tag  = $ARGV[0];
my $path = $ARGV[1];

my $j = 0;
my (@strains, %strains_mincov, %var_pos);
my (%homoAltCounts, %heteroCounts, %twoAltCounts, %LocPos_Status);
my ($chr, $pos, $AF, $AD, $DP, $GQ, $TAF, $covcutoff);
my ($A, $B, $R, $AA, $AB, $BB);

my (%Loc_n_len);
my $genomeSize = 0;
my $window_size = 1000;
my $window_start = 1;

my %T_strts =(
"AP014606.1" => '1250',
"AP014605.1" => '1250',
"AP014604.1" => '2000',
"AP014603.1" => '2000',
"AP014602.1" => '2250',
"AP014601.1" => '2250',
"AP014600.1" => '5000',
"AP014599.1" => '5000');

my %T_stops =(
"AP014606.1" => '944000',
"AP014605.1" => '940100',
"AP014604.1" => '1211250',
"AP014603.1" => '1350000',
"AP014602.1" => '1420500',
"AP014601.1" => '1563500',
"AP014600.1" => '1688000',
"AP014599.1" => '1721250');

my $rDNA_loc = "AP014603.1";
my $rDNAstrt = 1160000;
my $rDNAstop = 1179500;

# change here when more strains are included, do dynamically
my $strainNO = 1;
my $fieldCounter = $strainNO + 8;
#print "$fieldCounter\n"; exit;

open (COV, '<', "$path"."$tag.coverage_per10kb.tab");
while (<COV>) {
 if ($_ !~ /^Strain/) { 
  @_ = split(/\t|,/, $_); my $val = $_[1]*0.1; $strains_mincov{$_[0]} = "$val";
 }
}
close(COV); 

#=pod
open (TXT, '>', "$path"."/"."$path.SNVs.txt");
open (VCF, '>', "$path"."/"."$path.SNVs.vcf");

# do this from ARGV
open (IN,  '<', "$path"."/"."$tag.vcf");
 while (<IN>){
  # header
  if ($_ =~ /^##/) {
   print VCF "$_" ;
   # get contig sizes
   if ($_ =~ /^##contig/) { 
    @_ = split(/=|,/, $_); $_[4] =~ s/\W//g; $Loc_n_len{"$_[2]"} = $_[4]; 
	#print "$_[2] $_[4]\n";
   }
  }
  if ($_ =~ /#CHROM/) {
   # get genome size
   while (my ($key, $val) = each %Loc_n_len) {$genomeSize = $genomeSize + $val;}
   chomp; @_ = split(/\t/, $_); $_[0] =~ s/#//g;
   print TXT "$_[0]\t$_[1]\t$_[3]\t$_[4]\t";
   print VCF "#$_[0]\t$_[1]\t$_[2]\t$_[3]\t$_[4]\t$_[5]\t$_[6]\t$_[7]\t$_[8]\t";
   for (my $i=9; $i<=$fieldCounter; $i++){
    # get strain names
    (my $ID = $_[$i]) =~ s/Kmar_//g; $strains_mincov{$ID} = ""; push @strains, $ID;
    if ($i < $fieldCounter) { print TXT "$ID\t"; print VCF "$ID\t"; }
    else { print TXT "$ID\n"; print VCF "$ID\n";}
   }
  }
  # body
  if ($_ !~ /#/) {
   chomp; @_ = split(/\t/, $_); $j=0;
   # filter out * deleted alleles
   if ($_[4] !~ /\*/) {
    my $altcount = () = $_[4] =~ /,/g; $altcount++;
#print "$_\n";
    # ignore rDNA
#    unless ($_[0] =~ $rDNA_loc && ($_[1] >= $rDNAstrt && $_[1] <= $rDNAstop)){
#print "$_\n";
    # avoid telomeres
#     if ($_[1] > $T_strts{$_[0]} && $_[1] < $T_stops{$_[0]}){
#print "$_\n";
      # filters out indels and multiallelic sites
      #unless ($_[3] =~ /\w{2,}/ || $_[4] =~ /\w{2,}/ || $_[4] =~ /,/ || $_[0] !~ /^gi/) {
	  #unless ($_[3] =~ /\w{2,}/ || $_[4] =~ /\w{2,}/ || $_[4] =~ /,/) {

      # filters out indels (2 or more nts)
	   unless ($_[3] =~ /\w{2,}/ || $_[4] =~ /\w{2,}/ ) {
       my $TAFrecord = ""; my $VCFrecord = "";

	   # strain iteration
#       for (my $i=9; $i<=$fieldCounter; $i++){
		my $i = 9;
	    #print "$strains[$j]\n";
        my $preTAFrecord = ""; 
        my @INFO = split(/:/, $_[$i]);
        my @DEPTHS = split(/,/, $INFO[1]);
        $DP = $INFO[2];
        $covcutoff = $strains_mincov{$strains[$j]};
        if ($INFO[3] !~ /""/) {$GQ = $INFO[3];}
        else {$GQ = 0;}
        # Absolute coverage >5 reads; coverage > 10% average ; GQ > 20
        if ($DP <= $covcutoff || $DP <= 5 || $GQ < 20 ) { $TAF = "NA"; $preTAFrecord .= "$TAF"; }
        else{
         #modified to count ref for AA calling
         for (my $k = 0; $k <= $altcount; $k++){
          $AD = $DEPTHS[$k];
          $TAF = $AD/$DP;
          if ($TAF <= 0.15){ $TAF = "NA"; }
          if ($TAF >= 0.85){ $TAF = 1; }
          $preTAFrecord .= "$TAF";
          if($k < $altcount){ $preTAFrecord .= ","; }
         }
        }
        # avoid multi NA records, set variants with NA TAFs to non-variant vcf record
        if ($preTAFrecord !~ /\d/) { $preTAFrecord = "NA"; }
        if ($preTAFrecord eq "NA") { 
		 $_[$i] = "./.:0,0:0"; $VCFrecord .= "$_[$i]";
        }
        if ($preTAFrecord ne "NA") { $VCFrecord .= "$_[$i]"; }
#=pod
		# preTAF check to print only alternative AFs, including biallelic sites
		my $chckdTAFrecord = "NA";
		my @TAFs = split(/,/, $preTAFrecord);
		my $TAFsNo = @TAFs;
		# not sure this 1st check is necessary
		if ($TAFsNo == 1 || $preTAFrecord eq "NA") {$chckdTAFrecord = "NA";}
		if ($TAFsNo == 2) {$chckdTAFrecord = $TAFs[1];}
		if ($TAFsNo  > 2) {
		 for ( my $cnt = 1; $cnt <= $TAFsNo; $cnt++ ) {
          unless ($TAFs[$cnt] eq "NA") {
		   if ($chckdTAFrecord eq "NA") {$chckdTAFrecord = "$TAFs[$cnt]";}
		   else {$chckdTAFrecord .= ",$TAFs[$cnt]";}
		  }
         }
		}
		$preTAFrecord = $chckdTAFrecord;
#=cut
        $TAFrecord .= $preTAFrecord;
        # count variant classes
        if ($_[$i] !~ /0\/0/ && $_[$i] !~ /\.\/\./) {
         if ($_[$i] =~ m{(\d+)/(\d+)}){
          if ($1 == $2) {$homoAltCounts{$strains[$j]}++;}
          else {$heteroCounts{$strains[$j]}++;}
         }
        }
        if ($i < $fieldCounter){ $TAFrecord .= "\t"; $VCFrecord .= "\t";}
        # check for variation in at least one strain
        if ($i == $fieldCounter && $TAFrecord =~ /\d/){ 
		 
		 # fixing weird formatting introduced during preTAF check
         if ($TAFrecord =~ /^\t/)  { $TAFrecord =~ s/^\t/NA\t/g; }
		 if ($TAFrecord =~ /\t\t/) { $TAFrecord =~ s/\t\t/\tNA\t/g; }
		 if ($TAFrecord =~ /\t\t/) { $TAFrecord =~ s/\t\t/\tNA\t/g; }
		 if ($TAFrecord =~ /\tNA\n/) { $TAFrecord =~ s/\tNA\n//g; }
		 $TAFrecord =~ s/,\t/\t/g;
		 
         print TXT "$_[0]\t$_[1]\t$_[3]\t$_[4]\t$TAFrecord\n"; 
         print VCF "$_[0]\t$_[1]\t$_[2]\t$_[3]\t";
		 print VCF "$_[4]\t$_[5]\t$_[6]\t$_[7]\t$_[8]\t$VCFrecord\n";
         #my @DipRecords = split(/\t/, $VCFrecord); 
         # do the labeling here
         # homozygous reference
         #if ($DipRecords[$strain] =~ /0\/0/) { $AA = 1; }
         # heterozygous alternative
         #if ($DipRecords[$strain] =~ /1\/1/) { $BB = 1; }
         # heterozygous
         #if ($DipRecords[$strain] =~ /0\/1/) { $AB = 1; }
        }
        $j++ ; 
#       } # end of # strain iteration

       #$LocPos_Status{"$_[0]\t$_[1]"} = "$AA\t$BB\t$AB";
       #$AA = $AB = $BB = 0;
      }
#     } #END# avoid telomeres
#    } #END# ignore rDNA

   }
  }
 }
 close(IN); close(TXT); close(VCF);
 
open (OUT3, '>', "$path"."/"."$path.SNVs.tab");
print OUT3 "strain\thetero\thomoAlt\n";
foreach (@strains){
 print OUT3 "$_\t$heteroCounts{$_}\t$homoAltCounts{$_}\n";
}
