#!/bin/perl

my $snpfile = $ARGV[0];
my $genomepath = $ARGV[1];
my $readlength = $ARGV[2];
my $readquality = $ARGV[3];
my $sequencefile;
foreach my $chr (1..22,X,Y) {
  read_files();
  print_overlaps($chr);
}

sub read_files {
  open(SNP, '<', $snpfile) or die "cannot open snpfile";
  open(SEQUENCE, '<', "$genomepath/chr.$chr.fa") or die "cannot open sequencefile";
#DE: chr-formatted ref (UCSC):  open(SEQUENCE, "zcat $genomepath/chr$chr.fa.gz |") or die "cannot open sequencefile";
}

sub print_overlaps {

  $chr = shift @_;
  @sequence = <SEQUENCE>;
  chomp(@sequence);
  $sequence = join("",@sequence);
  $sequence =~ s/>chr$chr//g;
  #print $sequence;
  close(SEQUENCE);

  #DE: removes header: <SNP>;
  while ($line = <SNP>) {
    chomp($line);
    my @line = split(/ /, $line);
    my ($id, $thischr, $pos, $strand, $ref, $alt) = @line;
    #print("$pos\n");
    my $plusbase;
    #print "$thischr\tchr$chr\n";
#DE: chr-formatted ref:    if ($thischr eq "chr$chr") {
    if ($thischr eq $chr) {
      if ($strand eq "+") {
	$plusbase = $ref
      } elsif ($strand eq "-") {
	$plusbase = $ref;
	$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
      }

      for (my $i=0; $i<$readlength; $i++) {
	my $short = substr($sequence, ($pos-($readlength-1-$i)-1), $readlength-1-$i).$plusbase.substr($sequence, ($pos), $i);
	my $RC = $short;
	$RC =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
	$RC = reverse($RC);
	print('@', "$snpid_$thischr_$pos_REF_+_$plusbase", "_$i", "\n$short\n", '+', "$snpid_$thischr_$pos_REF_+_$plusbase", "_$i", "\n$readquality\n");
	print('@', "$snpid_$thischr_$pos_REF_-_$plusbase", "_$i", "\n$RC\n", '+', "$snpid_$thischr_$pos_REF_-_$plusbase", "_$i", "\n$readquality\n");
      }

      if ($strand eq "+") {
	$plusbase = $alt
      } elsif ($strand eq "-") {
	$plusbase = $alt; 
	$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
      }

      for (my $i=0; $i<$readlength; $i++) {
	my $short = substr($sequence, ($pos-($readlength-1-$i)-1), $readlength-1-$i).$plusbase.substr($sequence, ($pos), $i);
	my $RC = $short;
	$RC =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
	$RC = reverse($RC);
	print('@', "$snpid_$thischr_$pos_NONREF_+_$plusbase", "_$i", "\n$short\n", '+', "$snpid_$thischr_$pos_NONREF_+_$plusbase", "_$i", "\n$readquality\n");
	print('@', "$snpid_$thischr_$pos_NONREF_-_$plusbase", "_$i", "\n$RC\n", '+', "$snpid_$thischr_$pos_NONREF_-_$plusbase", "_$i", "\n$readquality\n");
      }
    }
  }
  close(SNP);
}


