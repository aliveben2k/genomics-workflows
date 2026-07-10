#!/usr/bin/perl
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
my $r_env = '$HOME/miniconda3/envs/R-4.1/bin';

chomp(@ARGV);
my $input; my $output; my $exc; my $mem = 16; my $local; my $pop;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-i" || $ARGV[$i] eq "\-I"){
		if (-e $ARGV[$i+1]){
			$input = $ARGV[$i+1];
		}
		else {
			print "Input file is not exist.\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-o" || $ARGV[$i] eq "\-O"){
		if ($ARGV[$i+1] =~ /^\-/){
			print "Output prefix is not acceptible.\n";
			exit;
		}
		$output = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-exc"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-mem"){
		$mem = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-local"){
		$local = "-cj_local ";
	}
	if ($ARGV[$i] eq "\-pop"){
		$pop = $ARGV[$i+1];
		unless (-e $pop){
			die "Cannot find the population file: $pop.\n";
		}
	}
}


if ($input !~ /\w/i || $output !~ /\w/i){
	print "no input or output directions.\n";
	&usage;
	exit;
}

unless (-e "vcf2table_missingNA_large.pl"){
	print "File \"vcf2table_missingNA_large.pl\" is required.\n";
}
unless (-e "Calculate_pairwise_dist_simple_large1.R"){
	print "File \"Calculate_pairwise_dist_simple_large1.R\" is required.\n";
}
unless (-e "Calculate_pairwise_dist_simple_large2.R"){
    print "File \"Calculate_pairwise_dist_simple_large2.R\" is required.\n";
}
unless (-e "Tree_PCoA.R"){
	print "File \"Tree_PCoA.R\" is required.\n";
}

if ($mem == 0){
	$mem = "";	
}
else {
	$mem = "-cj_mem $mem ";
}

RE:
unless ($ran){
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "qsub_files\/$ran\_phylo\.q" && $sn != 1){
	$ran = undef;
	goto RE;
}
my $out;
unless (-e "$output.1.tmp.txt.gz"){
	$out = "perl vcf2table_missingNA_large.pl $input $output\.txt $pop\\n";
    &pbs_setting("$exc$mem$local\-cj_quiet -cj_sn $ran -cj_qname vcf2table -cj_qout . $out");
    &status($ran);
}
unless (-e "$output\.rda"){
    my $path = $output;
    if ($path =~ /\//){
        my @tmp = split(/\//, $output);
        pop(@tmp);
        $path = join("\/", @tmp);
    }
    else {
        $path = '.';
    }
    my @files = <$path\/*.tmp.txt.gz>;
    unless (@files){
        die "Cannot find the files for processing.\n";
    }
    foreach my $i (0..$#files){
        my $outfile = $files[$i];
        $outfile =~ s/\.txt\.gz$//;
        $out = "Rscript Calculate_pairwise_dist_simple_large1.R $files[$i] $outfile\\n";
        &pbs_setting("$exc$mem$local\-cj_quiet -cj_env $r_env -cj_sn $ran -cj_qname partition_matrix_$i -cj_qout . $out");
    }
    &status($ran);
    $out = "Rscript Calculate_pairwise_dist_simple_large2.R $path $output\\n";
    &pbs_setting("$exc$mem$local\-cj_quiet -cj_env $r_env -cj_sn $ran -cj_qname merge_matrix -cj_qout . $out");
	#$out .= "Rscript Calculate_pairwise_dist_simple_mo.R $output\.txt.gz $output\\n";
    &status($ran);
    if (-e "$output.rda"){
        foreach (@files){
            system("rm $_");
            $rda_file = $_;
            $rda_file =~ s/txt\.gz$/rda/;
            system("rm $rda_file");
        }
    }
}
$out = "Rscript Tree_PCoA.R $output.rda $output\\n";
&pbs_setting("$exc$mem$local\-cj_quiet -cj_env $r_env -cj_qname phylo -cj_sn $ran -cj_qout . $out");
print "The qsub file is: qsub_files\/$ran\_phylo\.q\n";

sub usage {
	print BOLD "perl vcf2phylo.pl -i INPUT_FILE.vcf\(.gz\) -o OUTPUT_PREFIX_NAME [-pop POPULATION_INFO] [-local] [-mem INT] [-exc]\n", RESET;
	print BOLD "R library \"ape\" is required. Please make sure you have this package before execution.\n", RESET;
	return 1;
}
