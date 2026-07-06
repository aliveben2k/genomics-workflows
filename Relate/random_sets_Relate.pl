#!/usr/bin/perl
use Term::ANSIColor qw(:constants);

my $list = $ARGV[0];
my $repeats = $ARGV[1];
my $sample_num = $ARGV[2];
my $hap = $ARGV[3];
my $seed = $ARGV[4];
if ($ARGV[3] eq "hap"){
	$hap = 1;
}
elsif (!defined $hap || $hap eq ""){
	$hap = "NA";
}
if (defined $seed && $seed ne ""){
	srand($seed);
}

open(IN, "<$list") || die "Cannot open $list: $!\n";
my @lists = <IN>;
chomp(@lists);
close(IN);
shift(@lists);

my %pop;
my @sample_order;
foreach my $i (0..$#lists){
	my @eles = split(/\t|\s+/, $lists[$i]);
	push(@sample_order, $eles[0]);
	if (@{$pop{$eles[1]}}){
		push(@{$pop{$eles[1]}}, $eles[0]);
	}
	else {
		@{$pop{$eles[1]}} = $eles[0];
	}
}

foreach my $l (1..$repeats){
	my $out = $list;
	$out =~ s/poplabels$/$l.poplabels/;
	open(OUT, ">$out") || die "Cannot write $out: $!\n";
    print OUT "sample population group sex\n";
	my %picked;
	foreach my $key (sort keys %pop){
		my @samples = @{$pop{$key}};
		if (scalar(@samples) < $sample_num){
			die "Sample number is not enough for $key.\n";
		}
		my @picked_samples;
		foreach my $j (0..$sample_num-1){
			my $idx = int(rand(scalar(@samples)));
			push(@picked_samples, $samples[$idx]);
			splice(@samples, $idx, 1);
		}
		foreach my $sample (@picked_samples){
			$picked{$sample} = $key;
		}
	}
	foreach my $sample (@sample_order){
		next unless exists $picked{$sample};
		print OUT "$sample $picked{$sample} $picked{$sample} $hap\n";
	}
	close(OUT);
}
