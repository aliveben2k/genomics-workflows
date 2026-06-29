#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
#use PerlIO::gzip;
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
#my $output = $ARGV[1];
my $min_depth = 3;
my $out = $input;

my @tmp_outs = split(/\//, $out);

if ($tmp_outs[-1] =~ /raw/){
	$tmp_outs[-1] =~ s/raw/star_filtered/;
}
else {
	$tmp_outs[-1] =~ s/\.vcf/\.star_filtered\.vcf/;
}
$out = join("\/", @tmp_outs);

print "\[$time\]\: Filtering start...\n";
if ($input =~ /\.gz$/){
	open(INPUT, "-|", "gzip -dc $input") || die BOLD "Cannot open $input: $!", RESET, "\n";
#	$out =~ s/\.gz//;
}
else {
	open (INPUT, "<$input") || die BOLD "Cannot open $input: $!", RESET, "\n";
	$out = $out."\.gz";
}
open (OUTPUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
#open (OUTPUT, ">$out") || die BOLD ("Cannot write $out: $!\n"), RESET;
my $num = 1; 
my @sample_two_alleles;
my $fixed;
my $fix_num = 0; my $fix_count = 0; my $fix_stat = 0; my $total = 0;
my $no_DP = 0; my $less_min = 0; my $has_aster = 0;
while (my $line = <INPUT>){
	$fixed = 0;
	my @ALT_alleles_noAster; my $which_allele_aster; my $where_in_FORMAT_is_DP;
	chomp $line;
	if ($line =~ /^#/){
		print OUTPUT "$line\n";
	}
	else {
		$total += 1;
		my @line_vec = split('\t', $line);
		my @ALT_alleles = split(',', $line_vec[4]); #original ALT column
		for (my $i=0; $i<=$#ALT_alleles; $i++){
			if ($ALT_alleles[$i] !~ /\*/){ #if the ALT allele is not *, put it in @ALT_alleles_noAster
				push(@ALT_alleles_noAster, $ALT_alleles[$i]);
			}
			else {
				$which_allele_aster = $i+1; #allele number of '*'
			}
		}
		if (@ALT_alleles_noAster eq () || @ALT_alleles_noAster == ()){ #if no ALT allele after filtering, skip the position
#			push(@ALT_alleles_noAster, ".");
			next;
		}
		my $ALT_join = join(',', @ALT_alleles_noAster); #new ALT column
		if ($line !~ /\*/){
			$which_allele_aster = 500;
		}
		for (my $j=0; $j<=8; $j++){ #print information columns of the position except INFO
			if ($j == 4){
				$line_vec[$j] = $ALT_join;
			}
			if ($j == 0){
				print OUTPUT "$line_vec[$j]";
			}
			else {
				print OUTPUT "\t$line_vec[$j]";
			}
		}
		my @format = split (':', $line_vec[8]);
		for (my $k=0; $k<=$#format; $k++){
			if ($format[$k] =~ /DP/){ #get DP index in the INFO column
				$where_in_FORMAT_is_DP = $k;
			}
		}
		for (my $l=9; $l<=$#line_vec; $l++){
			my @sample_vec;
			@sample_vec = split(':', $line_vec[$l]);
			@sample_two_alleles = split(/\/|\|/, $sample_vec[0]);
			if ($sample_two_alleles[0] eq "." && $sample_two_alleles[1] eq "."){}
			elsif ($sample_vec[$where_in_FORMAT_is_DP] eq () || $sample_vec[$where_in_FORMAT_is_DP] == ()){
				$sample_two_alleles[0] = ".";
				$sample_two_alleles[1] = ".";
				$fixed = 1;
				$no_DP += 1;
			}
			elsif ($sample_vec[$where_in_FORMAT_is_DP] eq "." || int($sample_vec[$where_in_FORMAT_is_DP]) < $min_depth){
				$sample_two_alleles[0] = ".";
				$sample_two_alleles[1] = ".";
				$fixed = 1;
				$less_min += 1;
			}
			elsif (int($sample_two_alleles[0]) eq $which_allele_aster || int($sample_two_alleles[1]) eq $which_allele_aster){
				$sample_two_alleles[0] = ".";
				$sample_two_alleles[1] = ".";
				$fixed = 1;
				$has_aster += 1;
			}
			foreach my $m (0..$#sample_two_alleles){
				if ($sample_two_alleles[$m] ne "."){
					if (int($sample_two_alleles[$m]) > $which_allele_aster){
						$sample_two_alleles[$m] = int($sample_two_alleles[$m]) - 1;
					}
				}
			}
			if ($sample_vec[0] =~ /\|/){
				if ($sample_two_alleles[0] ne "." && $sample_two_alleles[1] ne "."){
					$sample_vec[0] = join("\|", @sample_two_alleles);
				}
				else {
					$sample_vec[0] = join("\/", @sample_two_alleles);
				}
			}
			else {
				$sample_vec[0] = join("\/", @sample_two_alleles);
			}
			if ($sample_vec[0] eq '.'){
				$sample_vec[0] = './.';
			}
			my $join = join(':', @sample_vec);
			print OUTPUT "\t$join";	
		}
		print OUTPUT "\n";
	}
	if ($fixed == 1){
		$fix_num += 1;
		$fix_count = int($fix_num/5000);
		if ($fix_count > $fix_stat){
			$time = scalar localtime();
			print "\[$time\]\: $fix_num variants have been fixed.\n";
			$fix_stat = $fix_count;
		}
	}
}
close(INPUT);
close(OUTPUT);
#system("gzip $out");
$time = scalar localtime();
print "\[$time\]\: $fix_num out of $total variants have been fixed.\n";
print "\[$time\]\: Including $no_DP sample position\(s\) with no DP; $less_min sample position\(s\) with depth less than 3, and $has_aster sample position\(s\) with asterisk.\n";
print "\[$time\]\: Filter_vcf done.\n";
