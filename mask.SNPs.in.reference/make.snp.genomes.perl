#!/bin/perl
use strict;
use warnings;

my $snpfile = $ARGV[0];
my $sequencefile = $ARGV[1];
my $chromosomename=$ARGV[2];
#my $chr;

open(SNP, "$snpfile");
open(SEQUENCE, "zcat $sequencefile |");
open(NRSEQ, ">$chromosomename");

print_chromosome();

sub print_chromosome {
my @line;
my $line;
my $plusbase;
my $otherbase;
<SEQUENCE>;
my @NRsequence= <SEQUENCE>;
chomp(@NRsequence);
my $NRsequence=join("",@NRsequence);
my $test_string;
<SNP>;
while($line=<SNP>){
@line = split(/ /, $line);
#print @line, "\n";
if($line[3] eq "+"){
$plusbase=$line[4];
$otherbase=$line[5];
}
elsif($line[3] eq "-"){
$plusbase = $line[4];
$otherbase=$line[5];
$plusbase =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
$otherbase  =~ tr/ATGCYRKMBDHVatgcyrkmbdhv/TACGRYMKVHDBtacgrymkvhdb/;
}

print(substr($NRsequence, $line[2]-5, 10), "\n");
$test_string = substr($NRsequence, $line[2]-1, 1);
if( uc($plusbase) eq uc($test_string) ) {
print "OK\n";
if(uc($plusbase) ne "C" && uc($otherbase) ne "C"){
substr $NRsequence, $line[2]-1, 1, "C";
} elsif(uc($plusbase) ne "G" && uc($otherbase) ne "G"){
substr $NRsequence, $line[2]-1, 1, "G";
} elsif(uc($plusbase) ne "A" && uc($otherbase) ne "A"){
substr $NRsequence, $line[2]-1, 1, "A";
} elsif(uc($plusbase) ne "T" && uc($otherbase) ne "T"){
substr $NRsequence, $line[2]-1, 1, "T";
}
} elsif($test_string eq "N"){
print "Base was N in reference\n";
}
else {
die "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
print "THE BASE GIVEN WAS NOT THE REFERENCE BASE, ABORTING MISSION!";
}
print(substr($NRsequence, $line[2]-5, 10),"\n");
}
print NRSEQ (">$chromosomename\n$NRsequence\n");
}
