#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
#use PerlIO::gzip;
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
chomp($input);
if (-e $input){
	if ($input =~ /vcf$/ || $input =~ /vcf\.gz$/){}
	else {
		&usage;
		die BOLD "\[$time\]\: $input is not a vcf.", RESET, "\n";
	}
}
else {
	&usage;
	exit;
}
#my $output = $ARGV[1];
my $dp = 3;
my $out = $input;
if ($out =~ /\//){
	my @temp = split(/\//, $out);
	pop(@temp);
	$out = join("\/", @temp);
	$out .= "\/statistics_vcf\.txt";
}
else {
	$out = "statistics_vcf\.txt";
}

print "\[$time\]\: Counting start...\n";
if ($input =~ /\.gz$/){
	open(INPUT, "-|", "gzip -dc $input") || die BOLD "\[$time\]\: Cannot open $input: $!", RESET, "\n";
#	$out =~ s/\.gz//;
}
else {
	open (INPUT, "<$input") || die BOLD "\[$time\]\: Cannot open $input: $!", RESET, "\n";
}
my @sample_two_alleles;
my $total = 0; my @list; my @missing; my @hetero;
while (my $line = <INPUT>){
	my @ALT_alleles_noAster; my $which_allele_aster; my $where_in_FORMAT_is_DP;
	chomp $line;
	if ($line =~ /^#/){
		if ($line =~ /^#CHROM/){
			my @samples = split(/\t/, $line);
			for (my $x=9; $x<=$#samples; $x++){
				push(@list, $samples[$x]);
				push(@missing, 0);
				push(@hetero, 0);
			}
		}
	}
	else {
		$total += 1;
		my @line_vec = split('\t', $line);
		for (my $l=9; $l<=$#line_vec; $l++){
			my @sample_vec;
			@sample_vec = split(':', $line_vec[$l]);
			@sample_two_alleles = split(/\/|\|/, $sample_vec[0]);			
			if ($sample_two_alleles[0] eq "." && $sample_two_alleles[1] eq "."){
				$missing[$l-9] += 1;
			}
			elsif ($sample_two_alleles[0] ne $sample_two_alleles[1]){
				$hetero[$l-9] += 1;
			}
			elsif ($sample_two_alleles[0] != $sample_two_alleles[1]){
				$hetero[$l-9] += 1;
			}
		}
	}
}
#print "@missing\n";
#print "@hetero\n";
close(INPUT);
open (OUTPUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
my $m_rate; my $h_rate;
print OUTPUT "Total variants\t$total\n";
print OUTPUT "Sample\tMissing site\tMissing rate\tHetero site\tHetero rate\n";
for (my $i=0; $i<=$#list; $i++){
	if ($missing[$i] != 0){
		$m_rate = $missing[$i]/$total*100;
		$m_rate = sprintf "%.2f", $m_rate;
	}
	else {
		$m_rate = 0;
	}
	if ($hetero[$i] != 0){
		$h_rate = $hetero[$i]/$total*100;
		$h_rate = sprintf "%.2f", $h_rate;
	}
	else {
		$h_rate = 0;
	}
	print OUTPUT "$list[$i]\t$missing[$i]\t$m_rate\t$hetero[$i]\t$h_rate\n";
}

close(OUTPUT);
#system("gzip $out");
$time = scalar localtime();
print "\[$time\]\: Counting done.\n";

sub usage {
	print BOLD "Usage: perl count_vcf_missing_rate.pl INPUT_VCF_FILE\n", RESET;
	return 1;
} #print usage
