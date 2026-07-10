#!/usr/bin/perl
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

my $plink = 'plink64';

chomp(@ARGV);
if ($ARGV[0] !~ /\w/){
	&usage;
	exit;
}
my $exc = 0; my $vcf; my $bed; my $ran; my $sn; my $kv; my $pru = 1; my $rep = 1; my $mem = 0;
my $cj_exc; my $local; my $ow;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-exc"){
		$exc = 1;
		$cj_exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-ow"){
		$ow = 1;
	}
	if ($ARGV[$i] eq "\-vcf"){
		if (-e $ARGV[$i+1] && $ARGV[$i+1] =~ /vcf$|vcf\.gz$/){
			$vcf = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-bed"){
		if (-e $ARGV[$i+1] && $ARGV[$i+1] =~ /bed$/){
			$bed = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-k"){
		if ($ARGV[$i+1] =~ /[0-9]/ && $ARGV[$i+1] > 0 && $ARGV[$i+1] <= 30){
			$kv = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-no_prune"){
		$pru = 0;
	}
	if ($ARGV[$i] eq "\-mem"){
		$mem = $ARGV[$i+1];
		unless ($mem !~ /[^0-9]/){
            die "-mem should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-sn"){
		$ran = $ARGV[$i+1];
		$sn = 1;
	}
	if ($ARGV[$i] eq "\-rep"){
		if ($ARGV[$i+1] !~ /[^0-9]/ && $ARGV[$i+1] <= 10){
			$rep = $ARGV[$i+1];
		}
		else {
			die "-rep should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-h" || $ARGV[$i] eq "\-\-help"){
		&usage;
		exit;
	}
	if ($ARGV[$i] eq "\-local"){
		$local = "-cj_local ";
	}
}
if ($local){
	$exc = 0;
}

RE:
if ($ran){}
else {
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "qsub_files\/$ran\_admixture_1.q" && $sn != 1){
	$ran = undef;
	goto RE;
}
if ($mem == 0){
	$mem = 6;
}

print "The q files are: qsub_files\/$ran\_\*\n";
my $out; my $vo;
my $p_mem = $mem * 1024;
if ($bed && $ow == 0){}
elsif ($vcf) {
	$vo = $vcf;
	if ($vo =~ /vcf$/){
		$vo =~ s/\.vcf$//g;
	}
	if ($vo =~ /vcf\.gz$/){
		$vo =~ s/\.vcf\.gz$//g;
	}
	if ($pru != 1){
		$out = "$plink \-\-vcf $vcf \-\-double\-id \-\-allow\-extra\-chr \-\-set\-missing\-var\-ids \@\:\# \-\-make\-bed \-\-threads 6 \-\-memory $p_mem \-\-pca \-\-out $vo\\n";
	}
	else {
		$out = "$plink \-\-vcf $vcf \-\-double\-id \-\-allow\-extra\-chr \-\-set\-missing\-var\-ids \@\:\# \-\-threads 6 \-\-memory $p_mem \-\-indep\-pairwise 50 10 0.1 \-\-out $vo\\n";
		$out .= "$plink \-\-vcf $vcf \-\-double\-id \-\-allow\-extra\-chr \-\-set\-missing\-var\-ids \@\:\# \-\-threads 6 \-\-memory $p_mem \-\-extract $vo\.prune\.in \-\-make\-bed \-\-pca \-\-out $vo\\n";
	}
	$bed = "$vo"."\.bed";
}
else {
	print "Cannot find proper vcf or bed file.\n";
	exit;
}

my @paths; my $path; my $name;
if ($bed =~ /\//){
	@paths = split(/\//, $bed);
	chomp(@paths);
	$name = pop(@paths);
	$path = join("\/", @paths);
}
else {
	$name = $bed;
	$path = ".";
}
$name =~ s/\.bed$//;
if ($out){
	unless (-e $bed && $ow == 0){
		if ($mem/6 > 6){
			$resource = "-cj_mem $mem ";
		}
		else {
			$resource = "-cj_ppn 6 -cj_mem 36 ";
		}
        &pbs_setting("$cj_exc$local$resource\-cj_quiet -cj_qname plink -cj_sn $ran -cj_qout . $out");
	}
}
if ($exc == 1){
	&status($ran);
}

#my $has_job = 0;
if ($mem != 0){
	$mem = "-cj_mem $mem ";
}
else {
    $mem = "";
}
for (my $i=1; $i<=$kv; $i++){
	$out = "";
	unless ($local){
		open(BASH, ">my_bash_admixture_$ran.sh") || print "Cannot write my_bash_admixture_$ran.sh: $!\n";
	}
	my $repeat_cnt = 1;
	for (my $j=1; $j<=$rep; $j++){
		unless (-e "$path\/$j\_$name\.$i\.Q" && $ow == 0){
			$out = "admixture \-\-cv $bed $i \-j1 \| tee $path\/$j\_log$i\.out\\n";
			$out .= "mv $name\.$i\.P $path\/$j\_$name\.$i\.P\\n";
			$out .= "mv $name\.$i\.Q $path\/$j\_$name\.$i\.Q\\n";
			$repeat_cnt += 1;
		}
		my $return = &pbs_setting("$cj_exc$mem$local\-cj_quiet -cj_qname admixture_$i\_rep$j -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
}
unless ($local){
	close(BASH);
}
if ($exc == 1){
	&status($ran);
}

$out = "grep \-h CV $path\/\*log\*\.out > $path\/K_values.txt\n";
&pbs_setting("$cj_exc$local\-cj_quiet -cj_qname admixture_cal -cj_sn $ran -cj_qout . $out");
if ($exc == 1){
	&status($ran);
}
print "Finished.\n";

sub usage {
	print BOLD "Usage: perl admixture_qsub.pl -vcf|-bed VCF_FILE|BED_FILE -k INT [-rep INT][-sn SN][-no_prune][-mem INT][-exc][-ow][-h\|--help]\n", RESET;
	return 1;
} #print usage
