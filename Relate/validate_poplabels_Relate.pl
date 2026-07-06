#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);

my $pop; my $out; my $hap = 0;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "-pop"){
		$pop = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-o"){
		$out = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-hap"){
		$hap = 1;
	}
}

unless ($pop && -e $pop){
	die "Usage: perl validate_poplabels_Relate.pl -pop POPLABELS -o OUTPUT [-hap]\n";
}
unless ($out){
	die "-o is required.\n";
}

open(IN, "<$pop") || die "Cannot open $pop: $!\n";
my @lines = <IN>;
close(IN);
chomp(@lines);

my $header = shift(@lines);
my @header_eles = split(/\s+|\t/, $header);
unless (scalar(@header_eles) >= 4 && lc($header_eles[0]) eq "sample" && lc($header_eles[1]) eq "population" && lc($header_eles[2]) eq "group" && lc($header_eles[3]) eq "sex"){
	die BOLD "-pop: header should be: sample population group sex\n", RESET;
}

my $sex_value = $hap == 1 ? "1" : "NA";
open(OUT, ">$out") || die "Cannot write $out: $!\n";
print OUT "sample population group sex\n";
foreach my $line (@lines){
	next unless $line =~ /\S/;
	my @eles = split(/\s+|\t/, $line);
	die "-pop: each line should contain at least four columns: sample population group sex. Problem line: $line\n" unless scalar(@eles) >= 4;
	die "-pop: sample column is empty. Problem line: $line\n" unless defined $eles[0] && $eles[0] ne "";
	die "-pop: population column is empty for sample $eles[0].\n" unless defined $eles[1] && $eles[1] ne "";
	die "-pop: group column is empty for sample $eles[0].\n" unless defined $eles[2] && $eles[2] ne "";
	$eles[3] = $sex_value;
	print OUT join(" ", @eles[0..3]), "\n";
}
close(OUT);
