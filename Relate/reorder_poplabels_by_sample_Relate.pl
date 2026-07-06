#!/usr/bin/perl
use strict;
use warnings;
use Term::ANSIColor qw(:constants);

my $pop; my $sample; my $out;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "-pop"){
		$pop = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-sample"){
		$sample = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-o"){
		$out = $ARGV[$i+1];
	}
}

unless ($pop && -e $pop && $sample && -e $sample && $out){
	die "Usage: perl reorder_poplabels_by_sample_Relate.pl -pop POPLABELS -sample SAMPLE -o OUTPUT\n";
}

open(POP, "<$pop") || die "Cannot open $pop: $!\n";
my @pop_lines = <POP>;
close(POP);
chomp(@pop_lines);
my $pop_header = shift(@pop_lines);
my @pop_header = split(/\s+|\t/, $pop_header);
unless (scalar(@pop_header) >= 4 && lc($pop_header[0]) eq "sample" && lc($pop_header[1]) eq "population" && lc($pop_header[2]) eq "group" && lc($pop_header[3]) eq "sex"){
	die BOLD "-pop: header should be: sample population group sex\n", RESET;
}

my %poplabels;
foreach my $line (@pop_lines){
	next unless $line =~ /\S/;
	my @eles = split(/\s+|\t/, $line);
	die "-pop: each line should contain at least four columns. Problem line: $line\n" unless scalar(@eles) >= 4;
	$poplabels{$eles[0]} = join(" ", @eles[0..3]);
}

open(SAMPLE, "<$sample") || die "Cannot open $sample: $!\n";
my @sample_lines = <SAMPLE>;
close(SAMPLE);
chomp(@sample_lines);
shift(@sample_lines); # header
shift(@sample_lines); # 0 0 0

my @sample_order;
foreach my $line (@sample_lines){
	next unless $line =~ /\S/;
	my @eles = split(/\s+|\t/, $line);
	push(@sample_order, $eles[0]);
}

my @missing;
foreach my $sample_id (@sample_order){
	push(@missing, $sample_id) unless exists $poplabels{$sample_id};
}
if (@missing){
	die BOLD "The following samples are missing in the poplabels file $pop:\n", RESET, join("\n", @missing), "\n";
}

open(OUT, ">$out") || die "Cannot write $out: $!\n";
print OUT "sample population group sex\n";
foreach my $sample_id (@sample_order){
	print OUT "$poplabels{$sample_id}\n";
}
close(OUT);
