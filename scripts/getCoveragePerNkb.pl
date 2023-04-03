#!/usr/bin/perl
use strict;

my $usage = "\nUsage: getCoveragePerNkb.pl <reference> <path> <window> <version>\n";
$usage .= " reference: fasta file to create indexes and use for mapping\n";
$usage .= " path:      path to coverage files\n";
$usage .= " window:    window size (e.g. 10000 for 10 kb)\n";
$usage .= " version:   e.g. Kmar_v0 \n\n";

die $usage unless @ARGV > 0;

my $refpath = $ARGV[0];
my $dirpath = $ARGV[1];
my $init_win_size = $ARGV[2];
my $ver = $ARGV[3];

my $kbsize = $init_win_size / 1000;

# changed expected for average, makes more sense
# calculates total average and by chromosome

my $filecount = 0;
my (%files, $path, $covPath, $IN, $OUT, $TAB);
my ($posIndex, $cov_in_chr, $cov_in_window, $Loc, $Pos, $totCov, $GenomeAverage);
my (%window_pos, %pos_cov, %Loc_n_len, @total_cov, @chr_cov, @window_cov);
my (%Loc_pos_n_cov);

my ($ID, $len);
my $genomeSize = 0;

open $TAB, '>', "$dirpath/"."$ver.coverage_per$kbsize"."kb.tab" or die "$!";
print $TAB "Strain\tGenome";

open my $ref, '<', $refpath or die "$!";
while (<$ref>){
 chomp;
 if ($_ =~ />/){
  @_ = split(/\s/, $_);
  ($ID = $_[0]) =~ s/>//g;
  print $TAB "\t$ID";
 }
 else { 
  $_ =~ s/\W//g;
  if   ($Loc_n_len{$ID}) { my $tmp0 = $Loc_n_len{$ID}; $Loc_n_len{$ID} = $tmp0 + length($_); }
  else {$Loc_n_len{$ID} = length($_); }
  my $tmp1 = $genomeSize; $genomeSize = $tmp1 + length($_); 
 }
}
print $TAB "\n";

my %tmp;
while ( my ($k,$v) = each %Loc_n_len ) { $tmp{"$v\t$k\t"} = $v; }
%Loc_n_len = %tmp;

my $rDNA_loc = "AP014603.1";
my $rDNAstrt = 1161094;
my $rDNAstop = 1180658;
my $rDNA_len = $rDNAstop - $rDNAstrt;
my $corrGenomeSize = $genomeSize - $rDNA_len;

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

#my $pattern = $ver."_coverage";
my $pattern = "_coverage";
opendir(my $dh, $dirpath) || die "Can't open $path: $!";
while (readdir $dh) {
 if ($_ =~ /$pattern$/){
  $path = "$dirpath/$_";
  print "\n parsing $path\n";
  (my $strain = $path) =~ s/$dirpath\///g; 
  @_ = split (/_/, $strain);
  $strain = $_[0];
  print $TAB "$strain\t";
 
  # info from sample_summary
  open $IN, '<', "$path.sample_summary" or die "$!";
  #while (<$IN>){ if ($_ =~ /^Total/){@_ = split(/\s/, $_); $totCov = $_[1];}} close($IN);
  while (<$IN>){ if ($_ =~ /^Total/){@_ = split(/,/, $_); $totCov = $_[1];}} close($IN);
 
  # coverage without rDNA 
  open $IN,  '<', $path or die "$!";
  $cov_in_window = 0;
  my $rDNA_cov = 0;
  while (<$IN>){
   if ($_ !~ /^Locus/){
    @_ = split(/\s/,$_);
    #@_ = split(/,/,$_);
    my @Loci_n_pos = split (/:/,$_[0]);
    $Loc = $Loci_n_pos[0]; $Pos = $Loci_n_pos[1];
    unless ($Loc =~ $rDNA_loc && ($Pos >= $rDNAstrt && $Pos <= $rDNAstop)){
     $cov_in_window = $cov_in_window + $_[1];
     push (@total_cov, $_[1]);
    }
   }
  }  
  close($IN);

  # this average is what some others call expected
  my $GenomeAverage = &average(@total_cov);
  my $GenomeStD = &stdev(@total_cov);
  print $TAB "$GenomeAverage,$GenomeStD";
  @total_cov=();
  
  # full thing
  open $IN,  '<', $path or die "$!";
  open my $OUT, '>', "$path"."_per".$kbsize."kb" or die "$!";
  print $OUT "Loc\tPos\tObserved\taverage\tDivision\n";
  $cov_in_window = 0; $cov_in_chr = 0; @chr_cov = ();
  my $window_size = $init_win_size;
  my $average_in_win;
  my ($mean_cov_in_chr, $mean_cov_in_window);
  while (<$IN>){
   if ($_ !~ /^Locus/){
    @_ = split(/\s/,$_);
    #@_ = split(/,/,$_);
    my @Loci_n_pos = split (/:/,$_[0]);
    $Loc = $Loci_n_pos[0]; $Pos = $Loci_n_pos[1];
#print "$Loc\n";
    push (@window_cov, $_[1]);
    my $Pos_n_Loc = "$Pos\t$Loc\t";
    $Loc_pos_n_cov{$Pos_n_Loc}=$_[1];
   
=pod
    # normalize rDNA coverage
    if ($Loc =~ $rDNA_loc && ($Pos >= $rDNAstrt && $Pos <= $rDNAstop)){
     $cov_in_window = $cov_in_window + $GenomeAverage;
     $cov_in_chr = $cov_in_chr + $GenomeAverage;
     push @chr_cov, $GenomeAverage;
    }
    # normalize telomeric coverage
    if ($Pos <= $T_strts{$Loc} || $Pos >= $T_stops{$Loc}){
     $cov_in_window = $cov_in_window + $GenomeAverage;
     $cov_in_chr = $cov_in_chr + $GenomeAverage;
     push @chr_cov, $GenomeAverage;
    }
    else {
     $cov_in_window = $cov_in_window + $_[1]; 
     $cov_in_chr = $cov_in_chr + $_[1];
     push @chr_cov, $_[1];
    }
=cut
    # no rDNA, nor telomeric "normalization"
    $cov_in_window = $cov_in_window + $_[1]; 
    $cov_in_chr = $cov_in_chr + $_[1];
    push @chr_cov, $_[1];
   
    if ($Pos_n_Loc == $Loc_n_len{$Pos_n_Loc} || $Pos == $window_size) {
     # why is it multiplied by window size?
     # shouldn't it be $cov_in_window/$win_average only?
     my $win_average = $GenomeAverage * scalar(@window_cov);
     my $division = $cov_in_window/$win_average;
     if ($Pos_n_Loc == $Loc_n_len{$Pos_n_Loc}) {
      print $OUT "$Loc\t$window_size\t$cov_in_window\t$win_average\t$division\n";
      $window_size = $init_win_size;
      #my $ChrAverage = $cov_in_chr/$Loc_n_len{$Pos_n_Loc};
      #print "\t"; print scalar @chr_cov; print "\n";
      my $ChrAverage = &average(@chr_cov);
      my $ChrStD = &stdev(@chr_cov);
      print $TAB "\t$ChrAverage,$ChrStD";
      @chr_cov = ();
      $cov_in_chr = 0;
     }
     if ($Pos == $window_size) { 
      print $OUT "$Loc\t$window_size\t$cov_in_window\t$win_average\t$division\n";
      $window_size = $window_size + $init_win_size;
     } 
     %Loc_pos_n_cov = (); @window_cov = (); $cov_in_window = 0;
    }
	# gaps ain't calculated, breaks count
    elsif ($Pos > $window_size) { 
     my $win_average = $GenomeAverage * scalar(@window_cov);
	 my $division;
	 if ($win_average == 0) { $division = 0; }
     else { $division = $cov_in_window/$win_average; }
     #print "$Loc\t$window_size\t$cov_in_window\t$win_average\t$division\n";
	 print $OUT "$Loc\t$window_size\t0\t0\t0\n";
     $window_size = $window_size + $init_win_size;
	 %Loc_pos_n_cov = (); @window_cov = (); $cov_in_window = 0;
    } 
	
   }
  }
  close($IN); close($OUT);
  print $TAB "\n";
 }
}

sub average{
 my @data = @_;
 if (not @data) { die("Empty array\n"); }
 my $total = 0;
 foreach (@data) { $total += $_; }
 my $average = $total / scalar @data;
 my $rounded = sprintf("%.3f", $average);
 return $rounded;
}

sub stdev{
 my @data = @_;
 if(@data == 1){ return 0; }
 my $average = &average(@data);
 my $sqtotal = 0;
 foreach(@data) { $sqtotal += ($average-$_) ** 2; }
 my $std = ($sqtotal / (scalar @data - 1)) ** 0.5;
 my $rounded = sprintf("%.3f", $std);
 return $rounded;
}

exit;
