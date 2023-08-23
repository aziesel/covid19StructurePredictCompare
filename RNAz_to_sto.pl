#!/usr/bin/perl -w

#read in the FASTA files and parse the FASTA header off of them
open(Bov, "<Bov");
@Bov = <Bov>;
close(Bov);
$Bov_seq = join('\n', @Bov);
$Bov_seq =~ s/\>(.*?)\n//;      #remove FASTA header in string
$Bov_seq =~ s/\n//gi;          #remove interstitial line breaks in string
$Bov_seq =~ s/\\n//gi;

open(COV19, "<COV19");
@COV19 = <COV19>;
close(COV19);
$COV19_seq = join('\n', @COV19);
$COV19_seq =~ s/\>(.*?)\n//;
$COV19_seq =~ s/\n//gi;
$COV19_seq =~ s/\\n//gi;

open(HKU4, "<HKU4");
@HKU4 = <HKU4>;
close(HKU4);
$HKU4_seq = join('\n', @HKU4);
$HKU4_seq =~ s/\>(.*?)\n//;
$HKU4_seq =~ s/\n//gi;
$HKU4_seq =~ s/\\n//gi;

open(HKU5, "<HKU5");
@HKU5 = <HKU5>;
close(HKU5);
$HKU5_seq = join('\n', @HKU5);
$HKU5_seq =~ s/\>(.*?)\n//;
$HKU5_seq =~ s/\n//gi;
$HKU5_seq =~ s/\\n//gi;

open(HKU9, "<HKU9");
@HKU9 = <HKU9>;
close(HKU9);
$HKU9_seq = join('\n', @HKU9);
$HKU9_seq =~ s/\>(.*?)\n//;
$HKU9_seq =~ s/\n//gi;
$HKU9_seq =~ s/\\n//gi;

open(Hu229E, "<Hu229E");
@Hu229E = <Hu229E>;
close(Hu229E);
$Hu229E_seq = join('\n', @Hu229E);
$Hu229E_seq =~ s/\>(.*?)\n//;
$Hu229E_seq =~ s/\n//gi;
$Hu229E_seq =~ s/\\n//gi;

open(HuHKU1, "<HuHKU1");
@HuHKU1 = <HuHKU1>;
close(HuHKU1);
$HuHKU1_seq = join('\n', @HuHKU1);
$HuHKU1_seq =~ s/\>(.*?)\n//;
$HuHKU1_seq =~ s/\n//gi;
$HuHKU1_seq =~ s/\\n//gi;

open(HuNL63, "<HuNL63");
@HuNL63 = <HuNL63>;
close(HuNL63);
$HuNL63_seq = join('\n', @HuNL63);
$HuNL63_seq =~ s/\>(.*?)\n//;
$HuNL63_seq =~ s/\n//gi;
$HuNL63_seq =~ s/\\n//gi;

open(HuOC43, "<HuOC43");
@HuOC43 = <HuOC43>;
close(HuOC43);
$HuOC43_seq = join('\n', @HuOC43);
$HuOC43_seq =~ s/\>(.*?)\n//;
$HuOC43_seq =~ s/\n//gi;
$HuOC43_seq =~ s/\\n//gi;

open(MERS, "<MERS");
@MERS = <MERS>;
close(MERS);
$MERS_seq = join('\n', @MERS);
$MERS_seq =~ s/\>(.*?)\n//;
$MERS_seq =~ s/\n//gi;
$MERS_seq =~ s/\\n//gi;

open(MHV, "<MHV");
@MHV = <MHV>;
close(MHV);
$MHV_seq = join('\n', @MHV);
$MHV_seq =~ s/\>(.*?)\n//;
$MHV_seq =~ s/\n//gi;
$MHV_seq =~ s/\\n//gi;

open(RaTG13, "<RaTG13");
@RaTG13 = <RaTG13>;
close(RaTG13);
$RaTG13_seq = join('\n', @RaTG13);
$RaTG13_seq =~ s/\>(.*?)\n//;
$RaTG13_seq =~ s/\n//gi;
$RaTG13_seq =~ s/\\n//gi;

open(SARS, "<SARS");
@SARS = <SARS>;
close(SARS);
$SARS_seq = join('\n', @SARS);
$SARS_seq =~ s/\>(.*?)\n//;
$SARS_seq =~ s/\n//gi;
$SARS_seq =~ s/\\n//gi;

#is there a way i could get this hash to accept command line arguments instead of
#hardcoding the genomes to be used?
%genome_vars = (Bov => 'Bov_seq',
  COV19 => 'COV19_seq',
  HKU4 => 'HKU4_seq',
  HKU5 => 'HKU5_seq',
  HKU9 => 'HKU9_seq',
  Hu229E => 'Hu229E_seq',
  HuHKU1 => 'HuHKU1_seq',
  HuNL63 => 'HuNL63_seq',
  HuOC43 => 'HuOC43_seq',
  MERS => 'MERS_seq',
  MHV => 'MHV_seq',
  RaTG13 => 'RaTG13_seq',
  SARS => 'SARS_seq'
);

#read in the RNAz output and parse out the 1) name, 2) length and 3) start position
#this really ought to be a command line argument
open(INPUT, "<corona_RNAz_out");
@RNAz = <INPUT>;
close(INPUT);


#collect the appropriate sequence and package it into individual FASTA format lines
mkdir(RNAz_flanked_sto);
chdir "RNAz_flanked_sto";
@output = ();
$genome = "";
$start = 0;
$len = 0;
$dir = "";

$outnum = 0;
for (my $i = 0; $i < @RNAz; $i++) {
  if ($RNAz[$i] =~ /RNAz/) {
    $filename = "RNAz_" . $outnum . ".sto";
    open(FILEOUT, ">$filename");

  }
  if ($RNAz[$i] =~ /^>(.*?)\s(\d*)\s(\d*)\s(.{1})\s(\d*)/) {
    if ($RNAz[$i] !~ /^>consensus/) {
      $genome = $1;
      $start = $2;
      $len = $3;
      $dir = $4;
      $genome_len = $5;
      $header = $genome . "_" . $start . " ";
      push(@output, $header);
      $genome = $genome_vars{$genome};
      $flankleft = "";
      $flankright = "";

      if ($dir eq "+") {
        $flanking = 20;
        $beg = $start - 20;
        if ($beg - 20 < 0) {
          $flanking = 20 - $start;
          $beg = 0;
        }
        if ($start == 0) {
          $flanking = 0;
          $beg = 0;
        }
        $flankleft = substr($$genome, $beg, $flanking);
        if ($flanking < 20) {
          $pad = '-' x (20 - $flanking);
          $flankleft = $pad . $flankleft;
        }
        if ($flanking == 0) {           #edgecase for 0
          $flankleft = '-' x 20;
        }
        $flankleft =~ tr/T/U/;
        $flankleft =~ s/\n//gi;

        $flanking = 20;
        $beg = $start + $len;
        if ($beg + $flanking > $genome_len) {
          $flanking = ($beg + $flanking - 1) - $genome_len;
        }
        $flankright = substr($$genome, $beg, $flanking);
        if ($genome_len - ($beg + $flanking) < 20) {
          $pad = '-' x (20 - ($genome_len - ($beg + $flanking)));
          $flankright = $flankright . $pad;
        }
        if ($genome_len - ($beg + $flanking) == 0) {
          $pad = '-' x 20;
          $flankright = $pad;
        }


        $flanking = 20;
        $beg = $start + $len;
        if ($beg + $flanking > $genome_len) {
          $flanking = ($beg + $flanking - 1) - $genome_len;
        }
        $flankright =~ tr/T/U/;
        $flankright =~ s/\n//gi;

      }
      if ($dir eq "-") {
        print("only works on forward right now\n");
      }
#     for reverse strand, need to: swap left and right, and take their reverse
#     complement (use reverse() and tr//). i had no reverse alignments to test on
      $middle = $RNAz[$i + 1];
      chomp($middle);
      $flanked = $flankleft . $middle . $flankright . "\n";
      push(@output, $flanked);

    }
  }
  if ($RNAz[$i] =~ /^>consensus/) {
    push(@output, "\n");
    print FILEOUT "# STOCKHOLM 1.0\n";
    print FILEOUT @output;
    print FILEOUT "//";
    close(FILEOUT);
    $outnum += 1;
    @output = ();
  }
}
