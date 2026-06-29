#!/usr/local/bin/perl

#Usage: perl GATK_pipeline.pl (gatk3) -p PATH -r REFERENCE_FILE -g Group_name [-sp START_STEP][-esp END_STEP][-adp ADAPTER_FILE][-db DB_PATH][--list LIST_FILE][-f FILTERED_VCF_FOLDER_PATH][-sn SERIAL_NUMBER][-prx SELECTED_VCF_FOLDER][-pox LIST.args][--ps LIST.args][-mmc INT][-mm FLOAT][-mq INT][-pfdp PREFIX][-kp LIST.txt][-x][-s][-a][-q][-c][-gv][-gvc][-ow][-as][-pf][-exc][-rad][-bao][-wi][-lc][-res FILE][-vqsr][-an]
#[Required]:
#-p: path to raw sequence data directory
#-r: path to reference file
#-g for group name of bam file (eg. Vigna or angularis (any group name that you want is acceptable))
#[Optional]: 
#gatk3: use GATK3 pipeline, if defined, this argument must be put at first place
#-sp: step of pipeline you want to proceed (0s, 0p, 1s, 1p, ....., 6s)
#-esp: end step of pipeline (1 to 6)
#-x for alternative indexing program: tabix
#-s HaplotypeCaller for Spark multithreads program (beta, results might be different)
#-a: alternative pipeline using FastqToSam-->BwaSpark-->MergeBamAlignment instead of bwa mem-->samtools view-->samtools sort
#-q: quality check of sequence file
#-adp: adapter file (if provided, the adapter trimming pipe will be turn on)
#-db: path to GenomicsDB
#--list: list file for importing samples to GenomicsDB
#-ow: over-write GenomicsDB
#-as: use -all-site option to generate a vcf file including all non-variant sites.
#-pf: prefix of chromosome name in vcf file. (Eg. Chr01-->the prefix would be "chr")
#-exc: Execute qsub query.
#-f: output folder path of the filtered vcf file.
#-c: using CombineGVCFs instead of GenomicsDBImport
#-gv: gathering vcfs of each chromosome file into one vcf file.
#-gvc: gathering vcfs of each chromosome file into one vcf file. this is multi-thread mode.
#-sn: Serial number of each batch run
#-prx: exclude sample from SelectVariants
#-pox: exclude sample from SelectVariants
#-ps: pre-select vcf to run GenotypeGVCFs
#-rad: skip MarkDuplicates step
#-bao: vcftools: turn off bi-allele selection
#-ri: vcftools: remove-indels
#-mm: vcftools: --max-missing
#-mmc: vcftools: --max-missing-count
#-mq: vcftools: --minQ
#-pfdp: prefix of chromosome name in vcf file for step 7. (Eg. Chr01-->the prefix would be "chr")
#-sf: skip star_filtering in step 5
#-dstr: enable DRAGEN-GATK pipeline in HaplotypeCaller
#-ng: do not generate gvcf files in step 2, please use -sp 2s or -esp 2 at the same time.
#-kp: keep user defined samples only when processing step 6.
#-lc: run the pipeline locally
#-ignc: ignore checking if input file numbers and gvcf file numbers are matched.
#-gpu: use GPU version of gatk4


use Cwd qw(getcwd);
use FindBin;
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
else {
	die "Cannot find required subroutine file: qsub_subroutine.pl\n";
}

my $p_ver = "3.4"; #this perl script version

print "GATK_pipeline v.$p_ver\nThe script is written by Ben Chien. 2024.11.\n";
print BOLD "Usage: perl GATK_pipeline_v$p_ver\_gpu.pl \(gatk3\) -p PATH -r REFERENCE_FILE [-g GROUP_NAME][-rad][-wes][-cnv][-sp START_STEP][-esp END_STEP][-adp ADAPTER_FILE][-db DB_PATH][--list LIST_FILE][-f FILTERED_VCF_FOLDER_PATH][-sn SERIAL_NUMBER][-prx SELECTED_VCF_FOLDER][-pox LIST.args][-ps LIST.args][-mmc INT][-mm FLOAT][-maf FLOAT][-mq INT][-pfdp PREFIX][-kp LIST.txt][-l FILE][-ip INT][-proj PROJECT_ID][-res FILE][-an NAME][-vqsr][-x][-s][-a][-q][-c][-gv][-gvc][-ow][-as][-pf][-sf][-exc][-bao][-wi][-lc][-gpu]\n", RESET;
print "For detail functions, please read \"Manual of GATK_pipeline.\"\n";
print RED "IMPORTANT: \*.g.vcf.gz files are incompatible between GATK3 and GATK4. However, \*.bam files are compatible.\n\n", RESET;
#print "GATK4 is highly recommend, many bugs are fixed with this pipeline.\n\n";

chomp(@ARGV);
my ($v3, $adp, $db, $ow, $as, $esp);
my $p_dir = $FindBin::Bin;
unless ($p_dir){
	$p_dir = ".";
}
my $gk3 = "0";
my $exc; my $sp = "0p"; my $spark = "0"; my $alter = "0"; my $qua="0"; my $ow = "0"; my $db = "GenomicsDB"; my $sf = "0"; my $ignc = "0"; my $cnv = 0;
my $pre = "0"; my $ref = "0"; my $path = "0"; my $fo = "03\-filter_vcf"; my $c = "1"; my $gv = "0"; my $sn; my $dbset; my $gname = "0"; my $dstr = "0"; my $ng = "0";
my $dblist = "0"; my @return; my $xlsn = "0"; my $pxlsn = "0"; my $prsn = "0"; my $posn = "0"; my $prese = "0"; my $rad = "0"; my $pre_DP = "chr"; my $ns = "0"; my $nlc = "0";
my $mmc = "0"; my $mm = "\-\-max\-missing 0\.9 "; my $ba = "\-\-min\-alleles 2 \-\-max\-alleles 2 "; my $ri = "\-\-remove\-indels "; my $minQ = "\-\-minQ 30 "; my $maf; my $wes; my $ip;
my $keep = 0; my $local; my $interval; $ip = 150; my $proj; my $trimmo = 0; my @reses; my @resis; my $vqsr = 0; my $as_name = "GRCh38"; my $gpu = 0; #my $clean = 0;
my $time = scalar localtime();
my $vep_env = '$HOME/miniconda3/envs/vep/bin';
my $r_env = '$HOME/miniconda3/envs/R-4.1/bin';

print "Input command line:\n";
print "perl GATK_pipeline_v$p_ver\_gpu.pl @ARGV\n\n";

if ($gk3 eq "1"){
	$v3 = `java \-jar \$GATKFILE \-version`;
	print "The GATK3 version: $v3\n";
}
else {
	my @v4s = `\$gatk2 \-version 2>&1`;
	foreach (@v4s){
		if ($_ =~ /The Genome Analysis Toolkit/){
			my @v4s_t = split(/\s/, $_);
			print "The GATK4 version: $v4s_t[-1]\n";
		}
	}
}

#check input arguments
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "gatk3"){
		$gk3 = "1";
		print "\[$time\]\: GATK3 pipeline will be used.\n";
	}
	if ($ARGV[$i] eq "\-sp" || $ARGV[$i] eq "\-\-step"){
		$sp = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-esp" || $ARGV[$i] eq "\-\-end\-step"){
		if ($esp !~ /[^0-9]/){
			$esp = int($ARGV[$i+1]);
		}
		else {
			print BOLD "\[$time\]\: ERROR\: -esp is set, but the input value is illegal. Please check!", RESET, "\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-exc" || $ARGV[$i] eq "\-\-execute"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-sn" || $ARGV[$i] eq "\-\-serial\-number"){
		$sn = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-ns" || $ARGV[$i] eq "\-\-no\-separation"){
		$ns = "1";
	}
	if ($ARGV[$i] eq "\-r" || $ARGV[$i] eq "\-\-reference"){ #step_0
		if (-e $ARGV[$i+1]){
			$ref = $ARGV[$i+1];
			chomp($ref);
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find the reference file!", RESET, "\n";
			exit;
		}		
	}
	if ($ARGV[$i] eq "\-p" || $ARGV[$i] eq "\-\-path"){ #step_1
		$time = scalar localtime();
		if ($ARGV[$i+1] =~ /\/$/){
			$ARGV[$i+1] =~ s/\/$//;
		}
		if (-d $ARGV[$i+1]){
			$path = $ARGV[$i+1];
			print "\[$time\]\: Sequence directory loaded.\n";
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find the path!", RESET, "\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-q" || $ARGV[$i] eq "\-\-quality"){ #step_1
		$qua = "1";
	}
	if ($ARGV[$i] eq "\-trimmo" || $ARGV[$i] eq "\-\-trimmomatic"){ #step_1
		$trimmo = 1;
	}
	if ($ARGV[$i] eq "\-adp" || $ARGV[$i] eq "\-\-adapter\-file"){ #step_1
		if ($gk3 eq "1"){
			$adp = "1";
		}
		else {
			$time = scalar localtime();
			if (-e $ARGV[$i+1]){
				$adp = $ARGV[$i+1];
				print "\[$time\]\: Adapter file loaded! $adp\n";
			}
			elsif ($i == $#ARGV || $ARGV[$i+1] =~ /^-/){
				$adp = "";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find the adapter file!", RESET, "\n";
				exit;
			}
		}	
	}
	if ($ARGV[$i] eq "\-rad" || $ARGV[$i] eq "\-\-RAD"){ #step_2
		$rad = "1";
	}
	if ($ARGV[$i] eq "\-g" || $ARGV[$i] eq "\-\-group\-name"){ #step_2
		if ($ARGV[$i+1] =~ /\W+/i){
			print BOLD "\[$time\]\: ERROR\: Only letters, numbers and \"_\" are acceptable for the group name!", RESET, "\n";
			exit;
		}
		if ($ARGV[$i+1] !~ /^\-/){
			$gname = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-ng" || $ARGV[$i] eq "\-\-no-gvcf"){ #step_2
		$ng = "1";
	}
	if ($ARGV[$i] eq "\-dstr" || $ARGV[$i] eq "\-\-dragSTR"){ #step_2
		$dstr = "1";
	}
	if ($ARGV[$i] eq "\-s" || $ARGV[$i] eq "\-\-spark"){ #step_2
		$spark = "1";
	}
	if ($ARGV[$i] eq "\-a" || $ARGV[$i] eq "\-\-alternative\-pipeline"){ #step_0 & 2
		$alter = "1";
	}
	if ($ARGV[$i] eq "\-db" || $ARGV[$i] eq "\-\-database\-path"){ #step_3
		if ($ARGV[$i+1] =~ /\/$/){
			$ARGV[$i+1] =~ s/\/$//;
		}
		if ($ARGV[$i+1] !~ /^\-/){
			$db = $ARGV[$i+1];
		}
		$dbset = "1";
	}
	if ($ARGV[$i] eq "\-\-list"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			$dblist = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-ow" || $ARGV[$i] eq "\-\-over\-write"){ #step_3
		$ow = "1";
	}
	if ($ARGV[$i] eq "\-nlc" || $ARGV[$i] eq "\-\-no-list-checking"){ #step_3
		$nlc = "1";
	}
	if ($ARGV[$i] eq "\-as" || $ARGV[$i] eq "\-\-all\-sites"){ #step_3
		$as = "\-all\-sites";
	}
	if ($ARGV[$i] eq "\-pf" || $ARGV[$i] eq "\-\-prefix"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			$pre = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-f" || $ARGV[$i] eq "\-\-folder"){ #step_3
		if ($ARGV[$i+1] =~ /\/$/){
			$ARGV[$i+1] =~ s/\/$//;
		}
		if ($ARGV[$i+1] !~ /^\-/){
			$fo = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-d" || $ARGV[$i] eq "\-\-use\-database"){ #step_3
		$c = "0";
	}
	if ($ARGV[$i] eq "\-\-pre\-xlsn" || $ARGV[$i] eq "\-prx"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			if (-e $ARGV[$i+1] && $ARGV[$i+1] !~ /\.args$/){
				print BOLD "\[$time\]\: ERROR\: The list file for \-pre\-xlsn must uses the extension \.args!", RESET, "\n";
				exit;
			}
			elsif (-e $ARGV[$i+1]) {
				$pxlsn = $ARGV[$i+1];
			}
			else {
				print BOLD "\[$time\]\: ERROR\: \-\-pre\-xlsn is set, but no file is found.", RESET, "\n";
				exit;
			}
		}
	}
	if ($ARGV[$i] eq "\-\-pre\-sn" || $ARGV[$i] eq "\-prn"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			if (-e $ARGV[$i+1] && $ARGV[$i+1] !~ /\.args$/){
				print BOLD "\[$time\]\: ERROR\: The list file for \-pre\-sn must uses the extension \.args!", RESET, "\n";
				exit;
			}
			elsif (-e $ARGV[$i+1]) {
				$prsn = $ARGV[$i+1];
			}
			else {
				print BOLD "\[$time\]\: ERROR\: \-\-pre\-sn is set, but no file is found.", RESET, "\n";
				exit;
			}
		}
	}	
	if ($ARGV[$i] eq "\-\-post\-xlsn" || $ARGV[$i] eq "\-pox"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			if (-e $ARGV[$i+1] && $ARGV[$i+1] !~ /\.args$/){
				print BOLD "\[$time\]\: ERROR\: The list file for \-post\-xlsn must uses the extension \.args!", RESET, "\n";
				exit;
			}
			elsif (-e $ARGV[$i+1]) {
				$xlsn = $ARGV[$i+1];
			}
			else {
				print BOLD "\[$time\]\: ERROR\: \-\-post\-xlsn is set, but no file is found.", RESET, "\n";
				exit;
			}
		}
	}
	if ($ARGV[$i] eq "\-\-post\-sn" || $ARGV[$i] eq "\-pon"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			if (-e $ARGV[$i+1] && $ARGV[$i+1] !~ /\.args$/){
				print BOLD "\[$time\]\: ERROR\: The list file for \-post\-sn must uses the extension \.args!", RESET, "\n";
				exit;
			}
			elsif (-e $ARGV[$i+1]) {
				$posn = $ARGV[$i+1];
			}
			else {
				print BOLD "\[$time\]\: ERROR\: \-\-post\-sn is set, but no file is found.", RESET, "\n";
				exit;
			}
		}
	}
	if ($ARGV[$i] eq "\-\-pre\-selected" || $ARGV[$i] eq "\-ps"){ #step_3
		if ($ARGV[$i+1] !~ /^\-/){
			if (-d $ARGV[$i+1]){
				print BOLD "\[$time\]\: ERROR\: Cannot find the pre-selected gvcf file folder!", RESET, "\n";
				exit;
			}
			else {
				$prese = $ARGV[$i+1];
				$prese =~ s/\/$//;
			}
		}
	}	
	if ($ARGV[$i] eq "\-gv" || $ARGV[$i] eq "\-\-gather\-vcfs"){ #step_gv
		$ARGV[$i+1] = int($ARGV[$i+1]);
		if ($ARGV[$i+1] !~ /[^0-9]/ && $ARGV[$i+1] > 3 && $ARGV[$i+1] < 7){
			$gv = "$ARGV[$i+1]";
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Only 4, 5, 6 are allowed for -gv\/\-\-gather\-vcfs.", RESET, "\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-gvc" || $ARGV[$i] eq "\-\-gather\-vcfs\-cloud"){ #step_gv
		if ($gv ne "4" || $gv ne "5" || $gv ne "6"){
			$ARGV[$i+1] = int($ARGV[$i+1]);
			if ($ARGV[$i+1] !~ /[^0-9]/ && $ARGV[$i+1] > 3 && $ARGV[$i+1] < 7){
				$gv = "2"."$ARGV[$i+1]";
			}
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Only 4, 5, 6 are allowed for -gv\/\-\-gather\-vcfs.", RESET, "\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-sf" || $ARGV[$i] eq "\-\-skip-filtering"){ #step_5
		$sf = "1";
	}
	if ($ARGV[$i] eq "\-bao" || $ARGV[$i] eq "\-\-bi\-allele\-off"){ #step_6
		$ba = "0";
	}
	if ($ARGV[$i] eq "\-wi" || $ARGV[$i] eq "\-\-with\-indels"){ #step_6
		$ri = "0";
	}
	if ($ARGV[$i] eq "\-mmc" || $ARGV[$i] eq "\-\-max\-missing\-count"){ #step_6
		$mmc = $ARGV[$i+1];
		if ($mmc =~ /[^0-9]/){
			print BOLD "\[$time\]\: ERROR\: \-mmc value should be an integer.\n", RESET;
			exit;
		}
		$mmc = int($ARGV[$i+1]);
		if ($mmc == 0){
			$mmc = "0";
		}
	}
	if ($ARGV[$i] eq "\-mm" || $ARGV[$i] eq "\-\-max\-missing"){ #step_6
		$mm = $ARGV[$i+1];
		if ($mm > 0 && $mm <= 1){
			$mm = "\-\-max\-missing $mm ";
		}
		elsif ($mm == 0 || $mm eq "0"){
			$mm = "0";		
		}
		else {
			print BOLD "\[$time\]\: ERROR\: \-mm value should be between 0 to 1\n", RESET;
			exit;			
		}
	}
	if ($ARGV[$i] eq "\-maf" || $ARGV[$i] eq "\-\-maf"){ #step_6
		$maf = $ARGV[$i+1];
		if ($maf > 0 && $maf <= 1){
			$maf = "\-\-maf $maf ";
		}
		elsif ($maf == 0 || $maf eq "0"){
			$maf = "";		
		}
		else {
			print BOLD "\[$time\]\: ERROR\: \-maf value should be between 0 to 1\n", RESET;
			exit;			
		}
	}
	if ($ARGV[$i] eq "\-mq" || $ARGV[$i] eq "\-\-minQ"){ #step_6
		$minQ = $ARGV[$i+1];
		if ($minQ =~ /[^0-9]/){
			print BOLD "\[$time\]\: ERROR\: \-minQ value should be an integer.\n", RESET;	
			exit;		
		}
		$minQ = int($minQ);
		if ($minQ == 0){
			$minQ = "0";
		}
		else {
			$minQ = "\-\-minQ $minQ ";
		}
	}
	if ($ARGV[$i] eq "\-kp" || $ARGV[$i] eq "\-\-keep"){ #step_6
		$keep = $ARGV[$i+1];
		if (-e $keep){
			$keep = "\-\-keep $keep ";
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find the file for -kp\/--keep.\n", RESET;		
			exit;	
		}
	}
	if ($ARGV[$i] eq "\-pfdp" || $ARGV[$i] eq "\-\-prefix\-DP"){ #step_7
		if ($ARGV[$i+1] !~ /^\-/){
			$pre_DP = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-lc" || $ARGV[$i] eq "\-\-local"){ #step_3
		$local = "-cj_local ";
	}
	if ($ARGV[$i] eq "\-ignc" || $ARGV[$i] eq "\-\-ignore-gvcf-number-checking"){ #step_3
		$ignc = "1";
	}
	if ($ARGV[$i] eq "\-wes" || $ARGV[$i] eq "\-\-WES"){ #for wes
		$wes = 1;
	}
	if ($ARGV[$i] eq "\-ip" || $ARGV[$i] eq "\-\-interval-padding"){ #for wes
		$ip = $ARGV[$i+1];
		if ($ip =~ /[^0-9]/){
			print BOLD "\[$time\]\: ERROR\: Only digitals are allowed for -ip\/--interval-padding.\n", RESET;		
			exit;			
		}
	}
	if ($ARGV[$i] eq "\-l" || $ARGV[$i] eq "\-\-interval"){ #for wes
		$interval = $ARGV[$i+1];
		unless (-e $interval){
			print BOLD "\[$time\]\: ERROR\: Cannot find the file for -l\/--interval.\n", RESET;		
			exit;		
		}
	}
	if ($ARGV[$i] eq "\-proj" || $ARGV[$i] eq "\-\-project"){ #for wes
		$proj = "-cj_proj $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-res" || $ARGV[$i] eq "\-\-resource"){ #for vcf annotation
		my $res = $ARGV[$i+1];
		if (-e $res){
			open(RES, "<$res") || die BOLD "\[$time\]\: ERROR\: Cannot find the file for -res\/--resource.", RESET, "\n";
			@reses = <RES>;
			chomp(@reses);
			close(RES);
			foreach my $j (0..$#reses){
				$reses[$j] = "--resource:$reses[$j]";
			}
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find the file for -res\/--resource.\n", RESET;
			exit;	
		}
	}
	if ($ARGV[$i] eq "\-resi" || $ARGV[$i] eq "\-\-resource-indel"){ #for vcf annotation
		my $resi = $ARGV[$i+1];
		if (-e $resi){
			open(RES, "<$resi") || die BOLD "\[$time\]\: ERROR\: Cannot find the file for -resi\/--resource-indel.", RESET, "\n";
			@resis = <RES>;
			chomp(@resis);
			close(RES);
			foreach my $j (0..$#resis){
				$resis[$j] = "--resource:$resis[$j]";
			}
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find the file for -resi\/--resource-indel.\n", RESET;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-vqsr" || $ARGV[$i] eq "\-\-VQSR"){ #for vcf annotation
		$vqsr = 1;
	}
	if ($ARGV[$i] eq "\-an" || $ARGV[$i] eq "\-\-assembly-name"){ #for vcf annotation
		$as_name = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-cnv" || $ARGV[$i] eq "--copy-number-variant"){ #step_6
		$cnv = 1;
	}
	if ($ARGV[$i] eq "\-gpu" || $ARGV[$i] eq "--GPU"){
	    $gpu = 1;
	}
=start
	if ($ARGV[$i] eq "\--clean"){ #clean the unnecessary files
		$clean = 1;
	}
=cut
}

#check dependencies
@return = &check_dep($exc, $gk3, $qua, $p_dir, $trimmo, $gpu);
my $bwa;
if ($return[0] eq "2"){
	print "\[$time\]\: Pipeline shotting down.\n";
	exit;
}
else {
	$bwa = $return[0];
}

#detect server
my @server = `ip addr`;
chomp(@server);
my $check_serv = 0;
foreach (@server){
	if ($_ =~ /140.110.148/){
		$check_serv = 1;
	}
}
if ($check_serv == 0){
	$proj = "";
}

#check if -p and -r are set
if ($path eq "0" || $ref eq "0"){
	$time = scalar localtime();
	if (($sp eq "0s" && $path eq "0") || ($sp eq "1s" && $ref eq "0")){}
	else {
		print BOLD "\[$time\]\: ERROR\: \-p or \-r are not defined.", RESET, "\n";
		exit;
	}
}

#check steps to run
my @stps = &step($sp, $esp);
$sp = $stps[0];
$esp = $stps[1];
$time = scalar localtime();
if ($sp eq "exit"){
	print "\[$time\]\: \"p\" or \"s\" should be supplied after step number.\n";
	exit;
}
print "\[$time\]\: \-sp set to $sp.\n";
if ($sp =~ /p/){
	print "\[$time\]\: \-esp set to $esp.\n";
}

#check -wes, -gv, -d, and -ns option
if ($wes == 1){
	$ns = "1";
	$alter = "0";
	$c = "1"; #use combineGVCFs
	unless ($interval){
		print BOLD "\[$time\]\: ERROR\: -l should be set when using -wes.", RESET, "\n";
		exit;
	}
	if ($gk3 == "1"){
		print BOLD "\[$time\]\: ERROR\: GATK3 for WES analysis is not supported in this pipeline.", RESET, "\n";
		exit;	
	}
	if ($ip == 0){
		$ip = "";
	}
} 
if ($ns eq "1"){
	$gv = "0";
}
if ($vqsr == 1){
	$gv = "4";
	unless (@reses){
		print BOLD "\[$time\]\: ERROR\: -res must be set when using -vqsr option.", RESET, "\n";
		exit;
	}
	unless (@resis){
		@resis = @reses;
	}
}

#check -d argument if -db is set
=start
if ($c eq "1"){
	if ($dbseq eq "1" || $wes != 1){
		print RED "\[$time\]\: WARNING\: \-db is set, but \-d is not set.\n", RESET;
		print RED "\[$time\]\: Switch to \-d mode.\n", RESET;
		$c = "0";
	}
}
=cut

#check coexistance of -prx/-prn or -pox/-pon
my $svsn1 = "\-xl\-sn";
my $svsn2 = "\-xl\-sn";
if ($pxlsn ne "0"){}
elsif ($prsn ne "0"){
	$pxlsn = $prsn;
	$svsn1 = "\-sn";
}
if ($xlsn ne "0"){}
elsif ($posn ne "0"){
	$xlsn = $posn;
	$svsn2 = "\-sn";
}

#check if -gpu is set
=start
if ($gpu == 1){
    $c = "1";
    $pxlsn = "0";
    $xlsn = "0";
    $dstr = "0";
}
=cut

#check step 6 arguments
if ($mmc ne "0"){
	$mm = "\-\-max\-missing\-count $mmc ";
}

#check existing SN and decide SN to use
my $dir = getcwd;
opendir DIR, $dir  || die BOLD "\[$time\]\: ERROR\: Cannot open $dir: $!", RESET, "\n";
my @sub_d = readdir DIR;
chomp(@sub_d);
my $sub = join("\t", @sub_d);
RAND:
my $ran = &rnd_str(4, "A".."Z", 0..9);
if ($sn){
	if ($sub !~ /$sn/){
		print RED "\[$time\]\: WARNING: Cannot find existing serial number.\n", RESET;
		$ran = $sn;
	}
	elsif ($sub =~ /$sn/){
		$ran = $sn;
	}
}
else {
	if ($sub =~ /$ran/){
		goto RAND;
	}
}
if ($fo eq "03\-filter_vcf"){
	$fo = "03\-filter_vcf_$ran";
}
print "\[$time\]\: Serial Number \(SN\) of this run: $ran\n";

#make a folder to store qsub files
if (-d "qsub_files"){
	$time = scalar localtime();
	print "\[$time\]\: The qsub will be stored at qsub_files\n";
}
else{
	system("mkdir qsub_files");
	system("mkdir qsub_files\/out");
	$time = scalar localtime();
	print "\[$time\]\: Make a folder\: qsub_files.\n[$time\]\: The gvcf file(s) will be stored there.\n";
}

#check paired-end fastq file names
unless ($sp =~ /0s|0p/){
	print "\[$time\]\: Check name of fastq files.\n";
	&pair($path);
}

#execute step 0: indexing reference file
if ($sp =~ /0s|0p/){
	if (-e "my_bash_00_$ran\.sh"){
		system("rm my_bash_00_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 0.\n";
	$ref = &ref_ind($proj, $bwa, $ran, $gk3, $ref, $alter, $exc, $gpu);
	if ($ref eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	if ($sp =~ /0p/){
		$sp = "1p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 0 is done.\n";
}
if ($sp eq "0s"){$time = scalar localtime(); print BOLD "\[$time\]\: Pipline running is done!", RESET,"\n"; exit;}

#check if the reference is indexed. If not, index it.
$time = scalar localtime();
my $ver;
my $dict = $ref;
$dict =~ s/fasta$|fas$|fa$/dict/g;
if ($sp =~ /2s|2p|3s|3p|4s|4p/){
	unless (-e "$ref\.fai" && -e $dict){
		$time = scalar localtime();
		$ref = &ref_ind($proj, $bwa, $ran, $gk3, $ref, $alter, $exc, $gpu);
		if ($exc){
			unless ($local){
				&status($ran);
			}
			print "\[$time\]\: A copy of the reference file and the indexed file are stored at $ref.\n";			
		}
	}
}

#execute step 1: quality check and trimming of fastq files
if ($sp =~ /1s|1p/){
	if (-e "my_bash_01_$ran\.sh"){
		system("rm my_bash_01_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 1.\n";
	if ($qua eq "1"){
		@return = &quality($proj, $ran, $exc, $path, $local);
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		if ($exc){
			$time = scalar localtime();
			unless ($local){
				&status($ran);
			}
		}
	}
	if ($gk3 eq "1"){
		@return = &trim_3($proj, $ran, $exc, $path, $adp, $local);
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	}
	else{
		if ($trimmo == 1){
			@return = &trim_4($proj, $ran, $exc, $path, $adp, $local);
		}
		else {
			@return = &trim_fastp($proj, $ran, $exc, $path, $local);
		}
		if ($return[0] eq "2"){
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	}
	$time = scalar localtime();
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	if ($sp =~ /1p/){
		$sp = "2p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 1 is done.\n";
}
if ($sp eq "1s" || $esp <= 1){print BOLD ("\[$time\]\: Pipline running is done!"), RESET, "\n"; exit;}

#execute step 2: generate raw gvcfs
if ($sp =~ /2s|2p/){
	if (-e "my_bash_02_$ran\.sh"){
		system("rm my_bash_02_$ran\.sh");
	}
	if ($gname eq "0" || $gname =~ /\-/){
		$time = scalar localtime();
		print BOLD "\[$time\]\: ERROR\: Group name must be defined.", RESET, "\n";
		print "\[$time\]\: Pipeline shotting down.\n";
		exit;
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 2.\n";
	if ($gpu == 1){
	    my @returns = &gpu_gvcf($cnv, $proj, $wes, $interval, $ip, $ng, $rad, $ran, $exc, $path, $ref, $local, $gname);
	    $path = $returns[0];
	    if ($path eq "2"){
			    $time = scalar localtime();
			    print "\[$time\]\: Pipeline shotting down.\n";
			    exit;
		    }
	    if ($exc){
		    $time = scalar localtime();
		    unless ($local){
			    &status($ran);
		    }
	    }
	    open (BASH2, ">my_bash_02_stat_$ran\.sh") || die BOLD "Cannot write my_bash_02_stat_$ran\.sh: $!", RESET, "\n";
	    @stats = @{$returns[-1]};
	    foreach my $cnt (0..$#stats){
	        my $shown_cnt = $cnt+1;
	        my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 2 -cj_time 168\:0\:0 -cj_qname gatk_02_stat_$shown_cnt -cj_sn $ran -cj_qout . $stats[$cnt]");
	        print BASH2 "$return\n";
	    }
	    close(BASH2);
	    if ($exc){
		    $time = scalar localtime();
		    unless ($local){
			    &status($ran);
		    }
	    }
	}
	else {
	    my @returns = &raw_gvcf($cnv, $proj, $wes, $interval, $ip, $ng, $dstr, $bwa, $rad, $ran, $exc, $path, $gk3, $ref, $gname, $spark, $alter, $local);
        $path = $returns[0];
        if ($path eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	    if ($exc){
		    $time = scalar localtime();
		    unless ($local){
			    &status($ran);
		    }
	    }
	}
	my $file_cnt = $returns[1];
	if ($exc || $local){
		&bam_stat($ran);
		if ($ignc eq "1"){
			my $check_file_num = `ls 02-get_gvcf_$ran \| grep \"vcf.gz\$\" \| wc`;
			unless ($check_file_num =~ /\b$file_cnt\b/){
				my @tmp = split(/\t+|\s+/, $check_file_num);
				$time = scalar localtime();
				print BOLD "\[$time\]\: ERROR\: Some gvcf file are not generated. Input: $file_cnt files. Output: $tmp[1] files.", RESET, "\n";
				exit;
			}
		}
	}
	if ($sp =~ /2p/){
		$sp = "3p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 2 is done.\n";
}
if ($sp eq "2s" || $esp <= 2){$time = scalar localtime(); print BOLD ("\[$time\]\: Pipline running is done!"), RESET, "\n"; exit}

#execute step 3: import gvcf into database or combine individual sample gvcf into combined interval (chromosome) gvcf
if ($sp =~ /3s|3p/){
	if (-e "my_bash_03_$ran\.sh"){
		system("rm my_bash_03_$ran\.sh");
	}
	if ($gk3 eq "1"){
		$sp =~ s/3/4/;
		print "\[$time\]\: Step 3 is skipped using GATK3 pipeline.\n";
		goto STEP4;
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 3.\n";
	if ($c eq "0"){
		@return = &db_import_4($proj, $wes, $interval, $ip, $ns, $fo, $ran, $exc, $path, $ref, $pre, $db, $dblist, $ow, $local);
		$ns = $return[5];
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		if ($exc){
			$time = scalar localtime();
			if ($return[4] > 1){
				for (my $i=1; $i<=$return[4]; $i++){
		            unless ($local){
			            &status($ran);
		            }
	            }
				$time = scalar localtime();
				print "\[$time\]\: Writing database sample list in $return[0]\$return[2]\_samples\.list\n";
				if ($return[1] =~ /0/){
					open(LIST, ">$return[0]\$return[2]\_samples\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $return[0]\$return[2]\_samples\.list: $!", RESET, "\n";
				}
				elsif ($return[1] =~ /1/){
					open(LIST, ">>$return[0]\$return[2]\_samples\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $return[0]\$return[2]\_samples\.list: $!", RESET, "\n";
				}
				print LIST "$return[3]\t";
				close(LIST);
			}
			elsif ($local !~ /[a-z]/i){
				&status($ran);
				$time = scalar localtime();
				print "\[$time\]\: Writing database sample list in $return[0]\$return[2]\_samples\.list\n";
				if ($return[1] =~ /0/){
					open(LIST, ">$return[0]\$return[2]\_samples\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $return[0]\$return[2]\_samples\.list: $!", RESET, "\n";
				}
				elsif ($return[1] =~ /1/){
					open(LIST, ">>$return[0]\$return[2]\_samples\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $return[0]\$return[2]\_samples\.list: $!", RESET, "\n";
				}
				print LIST "$return[3]\t";
				close(LIST);
			}
		}
		@return = ();
	}
	elsif ($c eq "1"){
		@return = &CombineGVCFs_4($proj, $nlc, $ns, $ran, $fo, $exc, $path, $ref, $ow, $pre, $dblist, $local);
		$ns = $return[1];
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		if ($exc){
			$time = scalar localtime();
			unless ($local){
				&status($ran);
				$time = scalar localtime();
				print "\[$time\]\: Writing sample list in $fo\/c_vcf\.list\n";
				if ($ow eq "1"){
					open(LIST, ">$fo\/c_vcf\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $fo\/c_vcf\.list: $!", RESET, "\n";
				}
				elsif ($ow eq "0"){
					open(LIST, ">>$fo\/c_vcf\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $fo\/c_vcf\.list: $!", RESET, "\n";
				}
				print LIST "$return[0]\t";
				close(LIST);
			}
		}
		@return = ();
	}
	if ($sp =~ /3p/){
		$sp = "4p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 3 is done.\n";
}
if ($sp eq "3s" || $esp <= 3){$time = scalar localtime(); print BOLD "\[$time\]\: Pipline running is done!", RESET, "\n"; exit;}

#execute step 4: compile gvcfs
STEP4:
if ($sp =~ /4s|4p/){
	if (-e "my_bash_04_$ran\.sh"){
		system("rm my_bash_04_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 4.\n";
	if ($gk3 eq "1"){
		@return = &GenotypeGVCFs_3($proj, $nlc, $ns, $ow, $ran, $fo, $exc, $path, $ref, $pre, $dblist, $as, $local);
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	}
	else {
		@return = &GenotypeGVCFs_4($proj, $wes, $interval, $ip, $ns, $prese, $pxlsn, $xlsn, $svsn1, $svsn2, $ran, $c, $fo, $exc, $path, $ref, $pre, $db, $as, $local, $gpu);
		if ($return[0] eq "2"){
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	}
	if ($gv =~ /4/ && $ns eq "0"){
		my $err_4gv = &GatherVcfs($proj, $xlsn, $gv, $ran, $fo, $exc, $path, $ref, $pre, $local);
		if ($err_4gv eq "2"){
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		$ns = "1";
	}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
			my $ree; my $try_n;
			do {
				$ree = &check_err($ran, $fo, $sp);
				if ($ree eq "2"){
					$ow = "0";
					if ($gk3 eq "1"){
						&GenotypeGVCFs_3($proj, $nlc, $ns, $ow, $ran, $fo, $exc, $path, $ref, $pre, $dblist, $as);
					}
					elsif ($c eq "1"){
						&GenotypeGVCFs_4($proj, $wes, $interval, $ip, $ns, $prese, $pxlsn, $xlsn, $svsn1, $svsn2, $ran, $c, $fo, $exc, $path, $ref, $pre, $db, $as, $local, $gpu);
					}
					$try_n += 1;
				}
			} until ($ree != 2 || $try_n == 3);
			$ree = &check_err($ran, $fo);
			if ($ree eq "2"){
				$time = scalar localtime();
				print BOLD "\[$time\]\: ERROR\: Some raw vcf file cannot finish.\n", RESET;
				print BOLD "\[$time\]\: Please use \"check_log.pl\" to see the problem.", RESET, "\n";
				print "\[$time\]\: Pipeline shotting down.\n";
				exit;
			}
			if ($gk3 eq "1"){
				$time = scalar localtime();
				print "\[$time\]\: Writing sample list in $fo\/c_vcf\.list\n";
				open(LIST, ">$fo\/c_vcf\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $fo\/c_vcf\.list: $!", RESET, "\n";
				if ($return[0] ne "1"){
					print LIST "$return[0]\t";
				}
				close(LIST);
			}
		}
	}
	if ($sp =~ /4p/){
		$sp = "5p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 4 is done.\n";
}
if ($sp eq "4s" || $esp <= 4){$time = scalar localtime(); print BOLD "\[$time\]\: Pipline running is done!", RESET, "\n"; exit;}

#execute step 5: filtering vcfs by perl script / VQSR filtering
if ($sp =~ /5s|5p/ && $sf eq "0"){
	if (-e "my_bash_05_$ran\.sh"){
		system("rm my_bash_05_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 5.\n";
	#print "debug: -ns: $ns\n";
	if ($vqsr == 1){
		$sp = "5s";
		@return = &VQSR($as_name, $r_env, $vep_env, $proj, $ran, $fo, $exc, $path, $ref, $pre, $local, \@reses, \@resis);
	}
	else {
		@return = &filter_vcf($proj, $gv, $ns, $ran, $ref, $path, $fo, $exc, $pre, $xlsn, $local, $p_dir);
		if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		if ($gv =~ /5/ && $ns eq "0"){
			my $err_5gv = &GatherVcfs($proj, $xlsn, $gv, $ran, $fo, $exc, $path, $ref, $pre, $local);
			if ($err_5gv eq "2"){
				print "\[$time\]\: Pipeline shotting down.\n";
				exit;
			}
			$ns = "1";
		}
	}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	if ($sp =~ /5p/){
		$sp = "6p";
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 5 is done.\n";
}
if ($sf eq "1" && $sp =~ /5s|5p/){
	$time = scalar localtime();
	print "\[$time\]\: \-sf is set. Skip step 5.\n";
	if ($sp eq "5p"){
		$sp = "6p";
	}
}
if ($sp eq "5s" || $esp <= 5){$time = scalar localtime(); print BOLD "\[$time\]\: Pipline running is done!", RESET, "\n"; exit;}

#execute step 6: get bi-allelic SNPs and further filtering
if ($sp =~ /6s|6p/){
	if (-e "my_bash_06_$ran\.sh"){
		system("rm my_bash_06_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 6.\n";
	@return = &bi_allele($proj, $keep, $sf, $xlsn, $minQ, $mm, $ri, $ba, $ns, $ran, $gv, $fo, $exc, $path, $ref, $pre, $local, $maf);
	if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	if ($gv =~ /6/ && $ns eq "0"){
		my $err_5gv = &GatherVcfs($proj, $xlsn, $gv, $ran, $fo, $exc, $path, $ref, $pre, $local);
		if ($err_5gv eq "2"){
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
		$ns = "1";
	}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 6 is done.\n";
	if ($sp =~ /6p/ && $esp eq "7"){
		$sp = "7s";
	}
}
if ($sp eq "6s" || $esp <= 6){$time = scalar localtime(); print BOLD "\[$time\]\: Pipline running is done!", RESET, "\n"; exit;}

#execute step 7: grep DP column from vcfs
if ($sp =~ /7|7s/){
	if (-e "my_bash_07_$ran\.sh"){
		system("rm my_bash_07_$ran\.sh");
	}
	$time = scalar localtime();
	print "\[$time\]\: Execute step 7.\n";
	@return = &Grep_DP($proj, $xlsn, $ns, $ran, $gv, $fo, $exc, $path, $ref, $pre, $pre_DP, $local);
	if ($return[0] eq "2"){
			$time = scalar localtime();
			print "\[$time\]\: Pipeline shotting down.\n";
			exit;
		}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	$time = scalar localtime();
	print "\[$time\]\: Step 7 is done.\n";
}

$time = scalar localtime();
print BOLD "\[$time\]\: Pipline running is done!\n", RESET;


#TOOLS
sub check_dep {
	my $exc = shift; my $gk3 = shift; my $qua = shift; my $p_dir = shift; my $trimmo_check = shift; my $gpu =shift;
	$time = scalar localtime();
	my $java = `java \-version 2>&1`; # check java version
	if ($java !~ /1.8/ && $gpu == 0){
		print BOLD "\[$time\]\: ERROR\: Java 8 JDK, OpenJDK 8 or OracleJDK 8 is required for GATK pipeline.\n", RESET;
		print "System java version:\n$java\n";
		return 2;
	}
	if ($gk3 eq "1" || $qua eq "1"){
		my $so = `SolexaQA\+\+ 2>&1`; # check trimmomatic software
		if ($so !~ /Usage/){
			print BOLD "\[$time\]\: ERROR\: Cannot find the dependency: SolexaQA\+\+.\n", RESET;
			return 2;
		}
	}
	my $fastp = `fastp 2>&1`;
	if ($fastp !~ /fastp/){
			print BOLD "\[$time\]\: ERROR\: Cannot find the dependency: fastp.\n", RESET;
			return 2;		
	}
	if ($trimmo_check == 1){
		my $trimmo = `java -jar $TRIMMO 2>&1`; # check trimmomatic software
		if ($trimmo !~ /Usage/){
			print BOLD "\[$time\]\: ERROR\: Cannot find the dependency: trimmomatic.\n", RESET;
			return 2;
		}
	}
	my $bwa = `bwa-mem2 2>&1`; # check bwa software
	if ($bwa !~ /bwa\-mem2 \<command\> \<arguments\>/ || $gpu == 1){
		$bwa = `bwa 2>&1`;
		if ($bwa !~ /bwa \<command\>/){
			print BOLD "\[$time\]\: ERROR\: Cannot find the dependency: bwa or bwa\-mem2.\n", RESET;
			return 2;
		}
		else {
			$bwa = "bwa";
		}
	}
	else {
		$bwa = "bwa-mem2";
	}
	my $sam = `samtools 2>&1`; # check samtools software
	if ($sam !~ /samtools \<command\>/){
		print BOLD "\[$time\]\: ERROR\: Cannot find the dependency: samtools.\n", RESET;
		return 2;
	}
	unless (-e "$p_dir\/filter_vcf_2.4\.pl"){ # check filter_vcf_2.4.pl file
		if ($exc){
			print "\[$time\]\: WARNING\: filter_vcf_2.4\.pl file is required!\n";
		}
		else {
			print "\[$time\]\: WARNING\: filter_vcf_2.4\.pl file is required! Make sure you have put the file at the same root with pipeline.pl before sending qsub job(s).\n";
		}
	}
	return ($bwa);
} #check dependencies
sub step {
	my $time = scalar localtime();
	my $steps = shift; my $e_stp = shift;
	if ($steps =~ /0s|0p|1s|1p|2s|2p|3s|3p|4s|4p|5s|5p|6s|6p|7s|7p|7/){
		if ($steps =~ /7p|7/){
			$steps = "7s";
			print "\[$time\]\: There is no \"7p\" step.\n\\[$time\]\: \-sp Set to \"7s\".\n";
		}
	}
	elsif ($steps =~ /0|1|2|3|4|5|6/){
		return "exit";
	}
	else {
		print "\[$time\]\: Cannot recognize the input value of -sp.\nSet to the default value \"0p\".\n";
		$steps = "0p";
	}
	if ($e_stp =~ /[0-7]/ && $e_stp !~ /[a-z]|\-/gi){
		$e_stp = int($e_stp);
	}
	else {
		$e_stp = 6;
	}
	if ($steps =~ /p/){
		my $s_stp = $steps;
		$s_stp =~ s/p//i;
		$s_stp = int($s_stp);
		if ($s_stp > $e_stp){
			$e_stp = 7;
			$steps =~ s/p/s/;
		}
	}
	return ($steps, $e_stp);
} #check steps to run
sub pair {
	my $path = shift;
	my $time; my @names; my $mo;
	opendir DIR, $path || die BOLD "ERROR\: Cannot open $path: $!", RESET, "\n";
	my @names = readdir DIR;
	chomp(@names);
	foreach my $name (@names){
		$time = scalar localtime();
		my $name_m = $name;
		$mo = "0";
		if ($name_m =~ /\.fq|\.fastq/gi){
			if ($name_m =~ /_1\.|_2\./){
				$name_m =~ s/_1\./_R1\./;
				$name_m =~ s/_2\./_R2\./;
				$mo = "1";
			}
			if ($name_m =~ /\.fq/){
				$name_m =~ s/\.fq/\.fastq/;
				$mo = "1";
			}
			
		}
		if ($mo eq "1"){			
			my $re = `mv $path\/$name $path\/$name_m 2>&1`;
			if ($re =~ /[a-z]/i){
				print BOLD  "\[$time\]\: ERROR\: $re", RESET, "\n";
				print BOLD "\[$time\]\: ERROR\: file name of $path\/$name cannot be modified for the pipeline. It may not process.", RESET, "\n";
			}
			else {
				print "\[$time\]\: $path\/$name has been modified to $path\/$name_m.\n";
			}
		}
	}
	return 1;	
} #check file name of paired-end files and modify them if necessary
sub check_chrs {
	my $time = scalar localtime();
	my $path_o = shift; my $ref = shift; my $pre = shift;
	my @Ls;
	if (-e $ref){
		$time = scalar localtime();
		if ($ref =~ /\.gz$/){
			open(FILE, "-|", "gzip -dc $ref") || die BOLD "\[$time\]\: ERROR\: Failed to check chromosome number: $!", RESET, "\n";
		}
		else {
			open(FILE, "<$ref") || die BOLD "\[$time\]\: ERROR\: Failed to check chromosome number: $!", RESET, "\n";
		}
	}
	else{
		my @ref_t = split("\/", $ref);
		my $ref2 = "$path_o\/$ref_t[-1]";
		$time = scalar localtime();
		print "\[$time\]\: WARNING: Make sure you have executed step 0 before executing this step.\n";
		if (-e $ref2){
			if ($ref2 =~ /\.gz$/){
				open(FILE, "-|", "gzip -dc $ref2") || die BOLD "\[$time\]\: ERROR\: Failed to check chromosome number: $!", RESET, "\n";
			}
			else {
				open(FILE, "<$ref2") || die BOLD "\[$time\]\: ERROR\: Failed to check chromosome number: $!", RESET, "\n";
			}
		}
		else{
			die BOLD "\[$time\]\: ERROR\: Failed to check chromosome number.", RESET, "\n";
		}
	}
	while (my $line = <FILE>){
		chomp($line);
		if ($line =~ /\>/){
			$line =~ s/\>//;
			if ($pre eq "0"){
				if ($line =~ /\s+|\t+/){
					my @temp = split(/\s+|\t+/, $line);
					push(@Ls, $temp[0]);
				}
				else {
					push(@Ls, $line);
				}
			}
			else {
				push(@Ls, $line) if ($line =~ /^$pre/i);
			}
		}
	}
	close(FILE);
	$time = scalar localtime();
	my @unique;
	if (@Ls eq ()){
		die BOLD "\[$time\]\: ERROR\: Cannot find prefix \"$pre\" in the file.", RESET, "\n";
	}
	else {
		my %seen;
		@unique = do { %seen; grep { !$seen{$_}++ } @Ls };
	}
	return (@unique);
} #check interval (chromosome) names and numbers
sub chr_length {
	$time = scalar localtime();
	my $ran = shift; my $folder = shift; my $path_o = shift; my $pre = shift;
	my @gvcf = <$folder\/*.vcf>;
#	print "$folder\n";
	if ($gvcf[0] !~ /\.vcf/){
		@gvcf = <$folder\/*.vcf.gz>;
		if ($gvcf[0] !~ /\.vcf\.gz/){
			@gvcf = <$path_o\/*.vcf>;
			if ($gvcf[0] !~ /\.vcf/){
				@gvcf = <$path_o\/*.vcf.gz>;
				if ($gvcf[0] !~ /\.vcf\.gz/){
					@gvcf = <02\-get_gvcf_$ran\/*.vcf.gz>;
				}
			}
		}
	}
	chomp(@gvcf);
	my @content; my @line; my @id; my @len; my $un; my $no_gz;
	if ($gvcf[0] =~ /\.vcf\.gz/){
		@content = `gzip \-cd $gvcf[0] \| head \-n 10000`;
	}
	elsif ($gvcf[0] =~ /\.vcf$/){
		@content = `head -n 10000 $gvcf[0]`;
	}
	else {
		print "\[$time\]\: WARNING\: Need a vcf file for \"\-si\" argument. Cannot find any vcf file.\n";
		print "\[$time\]\: WARNING\: \"\-si\" argument is skipped.\n";
		print "\[$time\]\: WARNING\: If interval length is too long, step 4 might be disrupted by the GATK3 bug.\n";
		print "\[$time\]\: WARNING\: Or it will excess walltime during GenomicsDBImport step.\n";
		return ("no_sp");
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
			if ($pre ne "0"){
				if ($line[3] =~ /$pre/i){
					push(@id, $line[3]);
					push(@len, $line[5]);
				}
			}
			elsif ($pre eq "0"){
				push(@id, $line[3]);
				push(@len, $line[5]);			
			}
		}
	}
	my $ids = join("\t", @id);
	unshift(@len, $ids);
	return (@len);
} #get interval length from vcf header
sub check_err {
	my $time = scalar localtime();
	my $ran = shift; my $folder = shift; my $sp = shift;
	my $out_d = "qsub_files\/out";
	if ($sp){
		$sp =~ s/s|p//;
	}
	opendir DIR, $out_d || die BOLD "\[$time\]\: ERROR\: Cannot open \/$out_d: $!", RESET, "\n";
	my @out_files = readdir DIR; chomp(@out_files); 
	my $out_contant; my $q_file; my $q_contant; my @q_lines; my $path;
	foreach (@out_files){
		if ($sp){
			if ($_ !~ /$ran_gatk_0$sp/){
				next();
			}
		}
		if ($_ =~ /$ran/){
			$out_contant = `cat qsub_files\/out\/$_`;
			$q_file = $_;
			$q_file =~ s/out/q/;
			if (-e "qsub_files\/$q_file") {
				$q_contant = `cat qsub_files\/$q_file`;
			}
			else {
				print "\[$time\]\: WARNING\: There is no qsub file: qsub_files\/$q_file, which is required for checking. Skip.\n";
				next(); 
			}	
		}
		else {
			next();
		}
		if ($_ =~ /gatk_04/){
			@q_lines = split("\n", $q_contant);
			foreach my $q_line (@q_lines){
				if ($q_line =~ /GenotypeGVCFs/){
					$path = &check_file($q_line);
				}
			}
			if ($q_contant =~ /GenotypeGVCFs/){
				if ($out_contant !~ /GenotypeGVCFs done|Done\./){
					print BOLD "\[$time\]\: WARNING\: Step 4 \(GenotypeGVCFs\) is incomplete: $_\n", RESET;
					system("rm $path");
					system("rm $path\.tbi");
					print "\[$time\]\: $path has been deleted.\n";
					return 2;
				}
			}
		}	
	}
	return 1;
} #check qsub job results if there is any error.
sub check_file {
	my $q_line = shift;
	my $path; my $path2;
	my @q_elements = split(/\s/, $q_line);
	chomp(@q_elements);
	for (my $i=0; $i<=$#q_elements; $i++) {
		if ($q_line =~ /SolexaQA\+\+ analysis/ && $q_elements[$i] eq "analysis"){
			$path = $q_elements[$i+1];
		}
		if ($q_line =~ /trimmomatic\.jar|\$TRIMMO/ && $q_elements[$i] eq "\-threads"){
			if ($q_elements[$i+1] =~ /_R1\./){
				$path = $q_elements[$i+4];
			}
			else {
				$path = $q_elements[$i+3];
			}
		}
		elsif ($q_line =~ /SolexaQA\+\+ dynamictrim/ && $q_elements[$i] eq "dynamictrim"){
			$path = $q_elements[$i+1];
			$path =~ s/\.fastq/\.fastq\.trimmed/;
		}
		if ($q_line =~ /bwa|bwa-mem2/){
			$path = $q_elements[-1];
		}
		if ($q_line =~ /IndexFeatureFile|tabix/ && $q_elements[$i] eq "\-I"){
			$path = $q_elements[$i+1];
		}
		elsif ($q_line =~ /IndexFeatureFile|tabix/ && $q_elements[$i] eq "\-p"){
			$path = $q_elements[$i+2];
		}
		if ($q_line =~ /HaplotypeCaller/ && $q_elements[$i] =~ "\-O"){
			$path = $q_elements[$i+1];
		}
		if ($q_line =~ /MarkDuplicatesSpark/){
			if ($q_elements[$i] =~ /\-O/){
				$path = $q_elements[$i];
				$path =~ s/\-O //;
			}
			if ($q_elements[$i] =~ /\-M/){
				$path2 = $q_elements[$i];
				$path =~ s/\-M //;
			}
		}
		if ($q_line =~ /CombineGVCFs/ && $q_elements[$i] eq "\-O"){
			$path = $q_elements[$i+1];
		}
		if ($q_line =~ /GenotypeGVCFs/ && ($q_elements[$i] eq "\-O" || $q_elements[$i] eq "\-o")){
			$path = $q_elements[$i+1];
		}
	}
	if ($path2){
		$path = $path." $path2";
	}
	return $path;
} #sub script for "check_err"
sub bam_stat {
	my $time = scalar localtime();
	my $ran = shift;
	my @bam = <01\-fq_trim_$ran\/*.samtools.stats>;
	chomp(@bam);
	if ($bam[0] =~ /stats$/){}
	else {
		print RED "\[$time\]\: WARNING\: Skip extracting *.samtools.stats info to bam_alignment_report.txt\n", RESET;
		return 1;
	}
	open (OUTPUT, ">01\-fq_trim_$ran\/bam_alignment_report\.txt") || die "Cannot write 01\-fq_trim_$ran\/bam_alignment_report\.txt: $!\n";
	print OUTPUT "Sample\tTotal_reads\tTotal_mapped\tMapping_rate\n";
	my @temp; my $total; my $mapped; my $rate; my @samples; my $count = 0;
	foreach (@bam){
		$count += 1;
		my $sample = $_;
		$sample =~ s/01\-fq_trim_$ran\///;
		$sample =~ s/_aln_sort.samtools.stats//;
		push(@samples, $sample);
		$time = scalar localtime();
		print "\[$time\]\: Counting $_ file for reporting...\n";
		my @contents = `head -n 30 $_ 2>&1`;
		foreach my $line (@contents){
			if ($line =~ /raw total sequences\:/){
				@temp = split(/\s+|\t+/, $line);
				$total = $temp[4];
			}
			if ($line =~ /reads mapped\:/){
				@temp = split(/\s+|\t+/, $line);
				$mapped = $temp[3];
			}
			if ($total != 0){
				$rate = $mapped/$total*100;
				$rate = sprintf "%.2f", $rate;
			}
			else {
				$rate = 0;
			}
		}
		print OUTPUT "$sample\t$total\t$mapped\t$rate\n";
		$total = 0;
		$mapped = 0;
		$rate = 0;
	}
	print OUTPUT "Counted samples: $count\nSample list:\n@samples\n";
	@samples = ();
	$count = "";
	$time = scalar localtime();
	print "\[$time\]\: Report generated.\n";
	close(OUTPUT);
} #generate bam alignment report
sub bgzip {
#	$ran\_gatk_04_gz_$_\.q
	my $time = scalar localtime();
	my $proj = shift; my $file = shift; my $chr = shift, my $exc = shift; my $ran = shift; my $local = shift;
	if (-e "$file.idx"){
		if (-e $file){
			my $out = "bgzip \-c \-\@ 8 $file\\n";
			$out .= "\$gatk2 IndexFeatureFile \-I $file\.gz\\n";
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 8 -cj_qname gatk_04_gz_$chr -cj_sn $ran -cj_qout . $out");
			return $return;
		}
		else {
			print BOLD "\[$time\]\: ERROR: Cannot find $file. Skip compressing.", RESET, "\n";
			return 2;
		}
	}
	else {
		print BOLD "\[$time\]\: ERROR: $file might be truncated. Skip compressing.", RESET, "\n";
		return 2;
	}	
} #a tool for compress g.vcf.gz before GenotypeGVCF in GATK4

#GATK3 & GATK4
sub ref_ind {
	my $proj = shift; my $bwa = shift; my $ran = shift; my $gk3 = shift; my $ref = shift; my $alter = shift; my $local = shift; my $exc = shift; my $gpu = shift;
	my $time = scalar localtime();
	my $out; my $ver;
	if ($gk3 eq "1"){
		$ver = "gatk3";
	}
	elsif ($gk3 eq "0"){
		$ver = "gatk4";
	}
	if (-e $ref){
		my $dict = $ref;
		$dict =~ s/fasta$|fas$|fa$/dict/g;
        if (-e "$ref\.fai" && -e "$ref\.amb" && -e "$ref\.ann" && (-e "$ref.bwt\.2bit\.64" || -e "$ref.bwt") && -e "$ref\.pac" && -e "$ref\.sa" && -e "$dict"){
			if ($gk3 eq "0" && $alter eq "1"){
				if (-e "$ref\.img"){
					return $ref;
				}
			}
			else {
				return $ref;
			}
		}
	}
	else {
		print BOLD "\[$time\]\: ERROR\: Could not find $ref\n", RESET;
		return 2;
	}
	if (-d "00\-ref_$ver"){
		print "\[$time\]\: The reference file is stored at 00\-ref_$ver\n";		
	}
	else {
		print "\[$time\]\: Make a folder\: 00\-ref_$ver\n\[$time\]\: The reference file will be stored there.\n";
		$out .= "mkdir 00\-ref_$ver\\n";
	}
	if ($ref =~ /\//){
		my @temp = split(/\//, $ref);
		chomp(@temp);
		my $ref2 = $temp[-1];
		$time = scalar localtime();
		unless (-e "00\-ref_$ver\/$ref2"){
			$out .= "cp $ref 00\-ref_$ver\/$ref2\\n";
			print "\[$time\]\: $ref will be copied to \.\/00\-ref_$ver\.\n";
		}
		$ref = "00\-ref_$ver\/$ref2";
	}
	else {
		$out .= "cp $ref 00\-ref_$ver\/$ref\\n";
		$ref = "00\-ref_$ver\/$ref";
	}
	#check again with new path
	my $dict = $ref;
	$dict =~ s/fasta$|fas$|fa$/dict/g;
    if (-e "$ref\.fai" && -e "$ref\.amb" && -e "$ref\.ann" && (-e "$ref.bwt\.2bit\.64" || -e "$ref.bwt") && -e "$ref\.pac" && -e "$ref\.sa" && -e "$dict"){
		if ($gk3 eq "0" && $alter eq "1"){
			if (-e "$ref\.img"){
				return $ref;
			}
		}
		else {
			return $ref;
		}
	}
	if ($ref =~ /fasta\.gz$|fa\.gz$|fas\.gz$/i){
		$ref =~ s/\.gz//;
		$time = scalar localtime();
		unless (-e $ref){
			$out .= "gunzip $ref\.gz\\n";
			print "\[$time\]\: $ref\.gz will be unzipped.\n";
		}
	}
	if ($ref !~ /fasta$/i && $ref !~ /fa$/i && $ref !~ /fas$/i){
		$time = scalar localtime();
		print BOLD "\[$time\]\: ERROR\: Could not find the fasta file!\n", RESET;
		return 2;
	}
	my $skip = 0;
	if ((-e "$ref\.img" && -e "$ref") || (-e "$ref\.bwt\.2bit\.64" && -e "$ref") || (-e "$ref\.bwt" && -e "$ref")){
		$time = scalar localtime();
		if ($alter eq "0" && -e "$ref\.bwt\.2bit\.64"){
			print "\[$time\]\: Indexing file of $ref exists!\n";
			$skip = 1;
		}
		elsif ($alter eq "0" && -e "$ref\.bwt"){
			print "\[$time\]\: Indexing file of $ref exists!\n";
			$skip = 1;
		}
		elsif ($alter eq "1" && -e "$ref\.img"){
            print "\[$time\]\: Image file of $ref exists!\n";
            $skip = 1;
		}
	}
    if ($gk3 eq "0" && $alter eq "1"){
        $time = scalar localtime();
        unless ($skip == 1){
        	unless (-e "$ref\.img"){
            	$out .= "\$gatk2 BwaMemIndexImageCreator \-I $ref \-O $ref\.img\\n";
            }
        }
    }
    else {
        if ($skip == 0){
        	unless (-e "$ref.bwt" || -e "$ref\.bwt\.2bit\.64"){
            	$out .= "$bwa index $ref\\n";
            }
        }
        if ($gpu == 1){
            unless (-e "$ref.bwt"){
                $out .= "bwa index $ref\\n";
            }
        }
    }
    my $ref_o = $ref;
    $ref_o =~ s/\.fasta$|\.fa$|\.fas$//i;
    unless (-e "$ref_o\.dict"){
    	if ($gk3 eq "1"){
        	$out .= "java \-jar \$PICARDFILE CreateSequenceDictionary R\=$ref O\=$ref_o\.dict\\n";
    	}
    	else {
        	$out .= "\$gatk2 CreateSequenceDictionary -R $ref -O $ref_o\.dict\\n";
    	}
    }
    unless (-e "$ref\.fai"){
        $out .= "samtools faidx $ref\\n";
    }
    &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 32 -cj_qname gatk_00 -cj_sn $ran -cj_qout . $out");
	return $ref;
} #indexing reference file
sub quality {
	my $proj = shift; my $ran = shift; my $exc = shift; my $path = shift; my $local = shift;
	my $time = scalar localtime();
	my @files = <$path\/*.fastq.gz>;
	my $gz = "\.gz"; my $file_2;
	if (@files eq ()){
		@files = <$path\/*.fastq>;
		$gz = ();
	}
	if (-d "01\-QA_report_$ran"){
		print "\[$time\]\: The trimmed file(s) are stored at 01\-QA_report_$ran\n";
	}
	else {
		system("mkdir 01\-QA_report_$ran");
		print "\[$time\]\: Make a folder\: 01\-QA_report_$ran\n\[$time\]\: The report file(s) will be stored there.\n";
	}
	my @temp;
	my $cnt = 0;
	my ($str1, $str2, $c_a, $type);
	unless ($local){
		open (BASH, ">my_bash_QA_$ran\.sh") || die BOLD ("Cannot write my_bash_QA_$ran\.sh: $!\n"), RESET;
	}
	foreach my $file (@files){
		next if ($file =~ /_R2/);
		my $out;
		if ($file =~ /_R1/){
			@temp = split(/_R1/, $file);
			$temp[0] = s/$path\///;
		}
		else {
			@temp = split(/\./, $file);
			$temp[0] =~ s/$path\///;
		}
		$time = scalar localtime();
		if (-e "01\-QA_report_$ran\/$temp[0]\.fastq\.quality\.pdf" || -e "01\-QA_report_$ran\/$temp[0]\_R1\.fastq\.quality\.pdf"){
			print "\[$time\]\: Quality files of $temp[0] exist.\n";
			next();
		}
		if ($file =~ /_R1/){
			$file_2 = $file;
			$file_2 =~ s/_R1/_R2/;
			unless (-e "$file_2" || $exc ne "-cj_exc "){
				print RED "\[$time\]\: WARNING\: It seems that the paired file of $file_2 is missing!\n", RESET;
				print RED "\[$time\]\: Skip quality check of $file.\n", RESET;
				next();
			}			
			$file_2 = " $file_2";
		}
		else {$file_2 = ();}
		print "\[$time\]\: Quality check of $file$file_2: Yes.\n";
		$out .= "SolexaQA\+\+ analysis $file$file_2 \-d 01\-QA_report_$ran\\n";
		if ($out =~ /[a-z]/i){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 4 -cj_qname gatk_QA_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
	}
	unless ($local){
		close(BASH);
	}
	return 1;
} #check quality of fastq files
sub raw_gvcf {
	my $cnv = shift; my $proj = shift; my $wes = shift; my $interval = shift; my $ip = shift; my $ng = shift; my $dstr = shift; my $bwa = shift; my $rad = shift; my $ran = shift; my $exc = shift; my $path_o = shift; my $gk3 = shift; my $ref = shift; my $gname = shift; my $spark = shift; my $alter = shift; my $local = shift;
	my $time = scalar localtime();
	my @files; my $stp; my $path = "01\-fq_trim_$ran"; my @che; my $j_che; my $drag; my $file_cnt = 0;
	if ($gk3 eq "1"){
		$alter = "0"; $spark = "0";
	}
	my $j_che; my @che;
	if (-d $path){
		opendir DIR, $path || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /trimmed|.bam/){
			print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path. Using path: $path_o\n";
			$path = $path_o;
			opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
			@che = readdir DIR; chomp(@che);
			close(DIR);
			$j_che = join(" ", @che);
			if ($j_che !~ /trimmed|.bam/){
				unless ($exc){
					print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.\n";
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.", RESET, "\n";
					exit;
				}
			}
		}
	}
	else {
		$path = $path_o;
		opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /trimmed|.bam/){
			unless ($exc){
				print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.", RESET, "\n";
				exit;
			}
		}
	}
	if ($ref =~ /gz$/){
		print BOLD "\[$time\]\: ERROR\: GATK does not support gzipped reference file, please unzip first.", RESET, "\n";
		return 2;
	}
	if (-e "$ref"){
		if (-e "$ref\.fai"){
			print "\[$time\]\: Reference $ref is verified.\n";
		}
	}
	elsif ($exc ne "-cj_exc "){
		print "\[$time\]\: WARNING: Cannot find reference file. Make sure reference is available before sending jobs.\n";
	}
	else {
		print BOLD "\[$time\]\: ERROR\: Cannot find reference file. Please check again.", RESET, "\n";
		return 2;
	}
	if (-d "02\-get_gvcf_$ran"){
		print "\[$time\]\: The gvcf file(s) will be stored at 02\-get_gvcf_$ran\n";
	}
	else {
		system("mkdir 02\-get_gvcf_$ran");
		print "\[$time\]\: Make a folder\: 02\-get_gvcf_$ran. The gvcf file(s) will be stored there.\n";
	}
	my @temp; my $q = "'"; my $cnt = 0; my ($f_l, $f_l_a); my $cnt = 0;
	my $in; my $list; my $file_2; my $gz = "\.gz";
	unless ($local){
		open (BASH, ">my_bash_02_$ran\.sh") || die BOLD "Cannot write my_bash_02_$ran\.sh: $!", RESET, "\n";
	}
	foreach my $file (@che){
		next if ($file =~ /_R2|unpaired|txt$|list$|bai$|stats$/);
		next if ($file !~ /[a-z]/i);
		my $out;
		$in = "0";
		$time = scalar localtime();
		@temp = split(/\./, $file);
		$temp[0] =~ s/_R1|_aln_sort_MD|_aln_sort|_metrics|_unaln|_aln//gi;
		my @tmp_list = split(/\t/, $list);
		my $exist = 0;
		foreach my $x (0..$#tmp_list){
			if ($tmp_list[$x] eq $temp[0]){
				$exist = 1;
				last;
			}
		}
		next if ($exist == 1);
		#next if ($list =~ /\b$temp[0]\b/);
		$list .= "$temp[0]\t";
		if (-e "$path\/$temp[0]_R1\.fastq\.trimmed" || -e "$path\/$temp[0]\.fastq\.trimmed"){
			if (-e "$path\/$temp[0]_R1\.fastq\.trimmed"){
				if (-e "$path\/$temp[0]_R2\.fastq\.trimmed"){
					$out .= "bgzip $path\/$temp[0]_R1\.fastq\.trimmed\\n";
					$out .= "bgzip $path\/$temp[0]_R2\.fastq\.trimmed\\n";
				}
			}
			else{
				$out .= "bgzip $path\/$temp[0]\.fastq\.trimmed\\n";
			}
		}		
		if (-e "$path\/$temp[0]_R1\.fastq\.trimmed\.gz" || -e "$path\/$temp[0]_R1\.fastq\.trimmed"){
			$stp = "0";
			$file = "$path\/$temp[0]_R1\.fastq\.trimmed\.gz";
			$in = "$bwa mem\" or \"FastqToSam";
			
		}
		if (-e "$path\/$temp[0]\.fastq\.trimmed\.gz" || -e "$path\/$temp[0]\.fastq\.trimmed"){
			$stp = "0";
			$file = "$path\/$temp[0]\.fastq\.trimmed\.gz";
			$in = "$bwa mem or FastqToSam";	
		}
		if (-e "$path\/$temp[0]\_aln_sort\.bam" && $rad eq "0" && $alter eq "0"){
			$stp = "1";
			$file = "$path\/$temp[0]\_aln_sort\.bam";
			$in = "MarkDuplicatesSpark";
		}
		if (-e "$path\/$temp[0]\_aln_pre_sort\.bam" && $rad eq "0" && $alter eq "1"){
			$stp = "1";
			$file = "$path\/$temp[0]\_aln_pre_sort\.bam";
			$in = "MarkDuplicatesSpark";
		}
		if (-e "$path\/$temp[0]\_aln_sort\.bam" && $rad eq "1" && $alter eq "0"){
			$stp = "2";
			$file = "$path\/$temp[0]\_aln_sort\.bam";
			$in = "samtools index";
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf"){
				system("rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf");
			}
		}
		if (-e "$path\/$temp[0]\_aln_pre_sort\.bam" && $rad eq "1" && $alter eq "1"){
			$stp = "2";
			$file = "$path\/$temp[0]\_aln_pre_sort\.bam";
			$in = "samtools index";
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf"){
				system("rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf");
			}
		}
		if (-e "$path\/$temp[0]\_aln_sort_MD\.bam" && $rad eq "0"){
			$stp = "2";
			$file = "$path\/$temp[0]\_aln_sort_MD\.bam";
			$in = "samtools index";
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf"){
				system("rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf");
			}
		}
		if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz"){
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi"){
				print "\[$time\]\: The raw vcf file of $temp[0] is ready. Nothing to do with this sample.\n";
				$file_cnt++;
				next();
			}
		}
		if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf" && $gk3 eq "1"){
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi"){
				print "\[$time\]\: The raw vcf file of $temp[0] is ready. Nothing to do with this sample.\n";
				$file_cnt++;
				next();
			}
			else {
				$stp = "3";
				$file = "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf";
				$in = "tabix";
			}
		}
		if ($exc ne "-cj_exc " && $in eq "0" && -e "$path_o\/$file" && $file =~ /fastq/i){
			$stp = "0";
			if ($file =~ /_R1/){
				$file = "$path\/$temp[0]\_R1.fastq\.trimmed\.gz";
			}
			else {
				$file = "$path\/$temp[0]\.fastq\.trimmed\.gz";
			}
		}
		elsif ($in eq "0") {next;}
		$time = scalar localtime();
		if ($in ne "0"){
			print "\[$time\]\: Validation\: $file exists.\n\[$time\]\: This file will start from \"$in\" check point.\n";
			$file_cnt++;
		}
		else {
			print "\[$time\]\: WARNING\: $file is simulated. Make sure you have executed step 1 before sending job\(s\).\n";
		}
		#WES								
		my $L; my $IP;
		if ($wes){
			$L = "-L $interval ";
			if ($ip){
				$IP = "-ip $ip ";
			}
		}
		#WES
		if ($stp eq "0"){
			if ($file =~ /_R1/){
				$file_2 = $file;
				$file_2 =~ s/_R1/_R2/gi;
				unless (-e $file_2 || $exc ne "-cj_exc "){
					print RED "\[$time\]\: WARNING\:\ It seems that the paired file of $file_2 is missing!\n", RESET;
					print RED "\[$time\]\: Skip generating the qsub job for $file.\n", RESET;
					next();
				}
				$f_l = "$file $file_2";
				$f_l_a = "-F1 $file -F2 $file_2";
			}
			else {
				$f_l = "$file";
				$f_l_a = "-F1 $file";
			}
			if ($alter eq "0"){ #non-gatk pipeline
				if (-e "$path\/$temp[0]\_aln_sort\.bam"){
					print "\[$time\]\: $path\/$temp[0]\_aln_sort\.bam exists.\n";
				}
				else {
					if (-e "$path\/$temp[0]\_aln\.sam"){
						$out .= "rm $path\/$temp[0]\_aln\.sam\\n";
					}
					$out .= "$bwa mem \-t 2 \-M \-R $q\@RG\\tID\:$gname\\tSM\:$temp[0]$q $ref $f_l > $path\/$temp[0]\_aln\.sam\\n";
					$out .= "samtools view \-\@ 2 \-bu $path\/$temp[0]\_aln.sam \| samtools sort \-\@ 2 \| samtools view \-\@ 2 $L\-b > $path\/$temp[0]\_aln_sort.bam\\n";
					$out .= "rm $path\/$temp[0]\_aln\.sam\\n";
					$out .= "samtools index \-\@ 2 $path\/$temp[0]\_aln_sort.bam\\n";
				}
			}
			elsif ($alter eq "1"){ #gatk pipeline
				if (-e "$path\/$temp[0]\_aln_pre_sort\.bam"){
					print "\[$time\]\: $path\/$temp[0]\_aln_pre_sort\.bam exists.\n";
				}
				else {
					if (-e "$path\/$temp[0]\_aln\.bam"){
						print "\[$time\]\: $path\/$temp[0]\_aln\.bam exists.\n";
					}
					else {
						if (-e "$path\/$temp[0]\_unaln\.bam"){
							print "\[$time\]\: $path\/$temp[0]\_unaln\.bam exists.\n";
						}
						else {
							$out .= "\$gatk2 FastqToSam $f_l_a -O $path\/$temp[0]\_unaln\.bam -SM $temp[0] -RG $gname\\n";
							$out .= "samtools index \-\@ 2 $path\/$temp[0]\_unaln\.bam\\n";
						}		
						if (-d "$path\/$temp[0]\_aln\.bam\.parts"){
							$out .= "rm \-r $path\/$temp[0]\_aln\.bam\.parts\\n";
						}
						if ($file =~ /_R1/){
							$out .= "\$gatk2 BwaSpark -I $path\/$temp[0]\_unaln\.bam \-O $path\/$temp[0]\_aln\.bam \-R $ref\\n";
							$out .= "samtools index \-\@ 2 $path\/$temp[0]\_aln\.bam\\n";
						}
						else {
							$out .= "\$gatk2 BwaSpark -I $path\/$temp[0]\_unaln\.bam \-O $path\/$temp[0]\_aln\.bam \-R $ref \-se\\n";
							$out .= "samtools index \-\@ 2 $path\/$temp[0]\_aln\.bam\\n";
						}							
					}
					$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g\" MergeBamAlignment -ALIGNED $path\/$temp[0]\_aln\.bam -UNMAPPED $path\/$temp[0]\_unaln\.bam -O $path\/$temp[0]\_aln_pre_sort\.bam -R $ref --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant -MAX_GAPS -1\\n";
					$out .= "samtools index \-\@ 2 $path\/$temp[0]\_aln_pre_sort\.bam\\n";
					$out .= "rm $path\/$temp[0]\_unaln.bam\*\\n";
					$out .= "rm $path\/$temp[0]\_aln\.bam\*\\n";
				}
			}
			if ($rad eq "0"){
				$stp = "1";
			}
			elsif ($rad eq "1"){
				$stp = "2";
			}
		}
		$time = scalar localtime();
		if ($stp eq "1" && $rad eq "0"){
			if (-e "$path\/$temp[0]_aln_sort_MD.bam"){
				print "\[$time\]\: $path\/$temp[0]_aln_sort_MD\.bam exists.\n";
			}
			else {
				if ($gk3 eq "1"){
					$out .= "java \-jar \$PICARDFILE MarkDuplicates -I $path\/$temp[0]_aln_sort.bam -O $path\/$temp[0]_aln_sort_MD.bam -M $path\/$temp[0]_aln_sort_MD_metrics\.txt\\n";
					$out .= "samtools index \-\@ 2 $path\/$temp[0]_aln_sort_MD.bam\\n";
					$out .= "rm $path\/$temp[0]_aln_sort.bam\*\\n";
				}
				elsif ($alter eq "0"){ #non-gatk pipeline
					$out .= "\$gatk2 MarkDuplicates -I $path\/$temp[0]_aln_sort.bam -O $path\/$temp[0]_aln_sort_MD.bam -M $path\/$temp[0]_aln_sort_MD_metrics\.txt\\n";
					#$out .= "samtools index \-\@ 2 $path\/$temp[0]_aln_sort_MD.bam\\n";
					$out .= "rm $path\/$temp[0]_aln_sort.bam\*\\n";
				}
				else { #gatk pipeline
					if (-e "$path\/$temp[0]_aln_sort_MD.bam"){
						print "\[$time\]\: $path\/$temp[0]_aln_sort_MD.bam exists.\n";
					}
					else {
						if (-e "$path\/$temp[0]_aln_MD.bam"){
							print "\[$time\]\: $path\/$temp[0]_aln_MD.bam exists.\n";
						}
						else {
							$out .= "\$gatk2 MarkDuplicatesSpark $L$IP\-I $path\/$temp[0]_aln_pre_sort.bam -O $path\/$temp[0]_aln_MD.bam -M $path\/$temp[0]_aln_sort_MD_metrics\.txt\\n";
							$out .= "samtools index \-\@ 2 $path\/$temp[0]_aln_MD.bam\\n";
							$out .= "rm $path\/$temp[0]_aln_pre_sort.bam\*\\n";
						}
						$out .= "\$gatk2 SortSamSpark $L$IP\-I $path\/$temp[0]_aln_MD.bam -O $path\/$temp[0]_aln_sort_MD.bam -SO coordinate\\n";
						$out .= "samtools index \-\@ 2 $path\/$temp[0]_aln_sort_MD.bam\\n";
						$out .= "rm $path\/$temp[0]_aln_MD.bam\*\\n";
					}
				}
			}
			$stp = "2";
		}
		$time = scalar localtime();
		if ($stp eq "2"){
			my $sp2_rad; my $sp2_in; my $stat;
			if ($rad eq "0"){
				$sp2_in = "$path\/$temp[0]_aln_sort_MD.bam";
			}
			elsif ($rad eq "1"){
				if ($alter eq "0"){}
				elsif ($alter eq "1"){ #gatk pipeline
					if (-e "$path\/$temp[0]_aln_sort.bam"){}
					else{
						$out .= "\$gatk2 SortSam $L$IP\-I $path\/$temp[0]_aln_pre_sort.bam -O $path\/$temp[0]_aln_sort.bam -SO coordinate\\n";
						$out .= "samtools index \-\@ 2 $path\/$temp[0]_aln_sort.bam\\n";
						$out .= "rm $path\/$temp[0]_aln_pre_sort.bam\*\\n";
					}
				}
				$sp2_in = "$path\/$temp[0]_aln_sort.bam";
			}
			if (-e "$sp2_in\.bai"){
				print "\[$time\]\: $sp2_in exists.\n";
			}
			else {
				$out .= "samtools index \-\@ 2 $sp2_in\\n";
			}
			$stat = $sp2_in;
			$stat =~ s/bam$/samtools\.stats/;
			if (-e $stat){
				print "\[$time\]\: $stat exists.\n";			
			}
			else {
				$out .= "samtools stats \-\@ 2 $sp2_in > $stat\\n";
			}
			if ($ng eq "1"){
				goto STEP2END;
			}
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz"){
				print "\[$time\]\: 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz exists.\n";
			}
			else {
				if ($gk3 eq "1"){
					$out .= "java \-jar \$GATKFILE \-T HaplotypeCaller \-\-emitRefConfidence GVCF \-nct 2 \-R $ref \-I $sp2_in \-o 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\\n";
				}
				else {
					if ($dstr eq "1"){ #GRAGEN-GATK pipeline
						$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g\" ComposeSTRTableFile $L$IP\-R $ref \-I $sp2_in \-O 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.zip\\n";
						$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g\" CalibrateDragstrModel $L$IP --threads 3 \-R $ref \-I $sp2_in \-str 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.zip \-O 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.table\\n";
						$drag = "\-\-dragen\-mode \-\-dragstr\-params\-path 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.table ";
					}
					my $anno;
					$anno = "-A AS_QualByDepth -A AS_FisherStrand -A AS_StrandOddsRatio -A AS_RMSMappingQuality -A AS_MappingQualityRankSumTest -A AS_ReadPosRankSumTest ";
					if ($spark eq "0"){
						$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g\" HaplotypeCaller $L$IP$anno\-ERC GVCF $drag\-R $ref \-I $sp2_in \-O 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\\n";
					}
					elsif ($spark eq "1"){
						$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g\" HaplotypeCallerSpark $L$IP$anno\-ERC GVCF $drag\-R $ref \-I $sp2_in \-O 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\\n";
					}
				}
			}
			if ($gk3 eq "1"){
				$out .= "bgzip 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\\n";
				$out .= "rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.idx\\n";
				$stp = "3";
			}
			if ($dstr eq "1" && $gk3 eq "0"){
				$out .= "rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.zip\\n";
				$out .= "rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.table\\n";
			}
		}
		$time = scalar localtime();
		if ($stp eq "3"){
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi"){
				print "\[$time\]\: 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi exists.\n";
				next();
			}
			else {
				$out .= "tabix \-p vcf 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\\n";
			}
		}
		STEP2END:
		if ($out =~ /[a-z]/ig){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 2 -cj_mem 60 -cj_qname gatk_02_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
		#print "debug: $cnt\n";
		$stp = "0"; @temp = ();
	}
	unless ($local){
		close(BASH);
	}
	#print BOLD "debug: $file_cnt", RESET, "\n";
	return ("02\-get_gvcf_$ran", $file_cnt);
} #generate raw gvcf files
sub gpu_gvcf {
    my $cnv = shift; my $proj= shift; my $wes = shift; my $interval = shift; my $ip = shift; my $ng = shift;
    my $rad = shift; my $ran = shift; my $exc = shift; my $path = shift; my $ref = shift; my $local = shift; my $gname = shift;
    my $time = scalar localtime();
	my @files; my $path = "01\-fq_trim_$ran"; my @che; my $j_che; my $drag; my $file_cnt = 0;
	if (-d $path){
		opendir DIR, $path || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /trimmed|.bam/){
			print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path. Using path: $path_o\n";
			$path = $path_o;
			opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
			@che = readdir DIR; chomp(@che);
			close(DIR);
			$j_che = join(" ", @che);
			if ($j_che !~ /trimmed|.bam/){
				unless ($exc){
					print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.\n";
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.", RESET, "\n";
					exit;
				}
			}
		}
	}
	else {
		$path = $path_o;
		opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /trimmed|.bam/){
			unless ($exc){
				print "\[$time\]\: WARNING\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find any fastq.trimmed\(.gz\) or .bam file in $path.", RESET, "\n";
				exit;
			}
		}
	}
    if ($ref =~ /gz$/){
		print BOLD "\[$time\]\: ERROR\: GATK does not support gzipped reference file, please unzip first.", RESET, "\n";
		return 2;
	}
	if (-e "$ref"){
		if (-e "$ref\.fai"){
			print "\[$time\]\: Reference $ref is verified.\n";
		}
	}
	elsif ($exc ne "-cj_exc "){
		print "\[$time\]\: WARNING: Cannot find reference file. Make sure reference is available before sending jobs.\n";
	}
	else {
		print BOLD "\[$time\]\: ERROR\: Cannot find reference file. Please check again.", RESET, "\n";
		return 2;
	}
	if (-d "02\-get_gvcf_$ran"){
		print "\[$time\]\: The gvcf file(s) will be stored at 02\-get_gvcf_$ran\n";
	}
	else {
		system("mkdir 02\-get_gvcf_$ran");
		print "\[$time\]\: Make a folder\: 02\-get_gvcf_$ran. The gvcf file(s) will be stored there.\n";
	}
	my @temp; my $q = "'"; my $cnt = 0; my $f_l; my @stats;
	my $list; my $file_2; my $gz = "\.gz";
	unless ($local){
		open (BASH, ">my_bash_02_$ran\.sh") || die BOLD "Cannot write my_bash_02_$ran\.sh: $!", RESET, "\n";
	}
	foreach my $file (@che){
		next unless ($file =~ /trimmed|bam$/);
		next if ($file =~ /_R2/);
		next if ($file !~ /[a-z]/i);
		my $out;
		$time = scalar localtime();
		@temp = split(/\./, $file);
		$temp[0] =~ s/_R1|_aln_sort_MD|_aln_sort|_metrics|_unaln|_aln//gi;
		my @tmp_list = split(/\t/, $list);
		my $exist = 0;
		foreach my $x (0..$#tmp_list){
		    #print "debug: $tmp_list $temp[0]\n";
			if ($tmp_list[$x] eq $temp[0]){
				$exist = 1;
				last;
			}
		}
		next if ($exist == 1);
		#next if ($list =~ /\b$temp[0]\b/);
		$list .= "$temp[0]\t";
		#print "debug: $list\n";
		if (-e "$path\/$temp[0]_R1\.fastq\.trimmed" || -e "$path\/$temp[0]\.fastq\.trimmed"){
			if (-e "$path\/$temp[0]_R1\.fastq\.trimmed"){
				if (-e "$path\/$temp[0]_R2\.fastq\.trimmed"){
					$out .= "bgzip $path\/$temp[0]_R1\.fastq\.trimmed\\n";
					$out .= "bgzip $path\/$temp[0]_R2\.fastq\.trimmed\\n";
				}
			}
			else{
				$out .= "bgzip $path\/$temp[0]\.fastq\.trimmed\\n";
			}
		}
		if (-e "$path\/$temp[0]_R1\.fastq\.trimmed.gz"){
		    $file = "$path\/$temp[0]_R1\.fastq\.trimmed.gz"
		}
		else {
		    $file = "$path\/$temp[0]\.fastq\.trimmed.gz";
		}
		$time = scalar localtime();	
		#WES								
		my $L; my $IP;
		if ($wes){
			$L = "--interval-file $interval ";
			if ($ip){
				$IP = "-ip $ip ";
			}
		}
		#WES
		$out = "";
		if ($file =~ /_R1/){
			$file_2 = $file;
			$file_2 =~ s/_R1/_R2/gi;
			unless (-e $file_2 || $exc ne "-cj_exc ") {
				print RED "\[$time\]\: WARNING\:\ It seems that the paired file of $file_2 is missing!\n", RESET;
				print RED "\[$time\]\: Skip generating the qsub job for $file.\n", RESET;
				next();
			}
			$f_l = "--in-fq $file $file_2 $q\@RG\\tID\:$gname\\tLB:lib1\\tPL:bar\\tSM\:$temp[0]\\tPU\:unit1$q ";
		}
		else {
			$f_l = "--in-se-fq $file $q\@RG\\tID\:$gname\\tLB:lib1\\tPL:bar\\tSM\:$temp[0]\\tPU\:unit1$q ";
		}
		if ($rad eq "0"){
		    unless (-e "$path\/$temp[0]_aln_sort_MD.bam"){
			    $out = "pbrun fq2bam --ref $ref $L$f_l$IP\--out-bam $path\/$temp[0]_aln_sort_MD.bam --num-gpus 4 --low-memory\\n";
			}
			else {
			    if ($ng eq "0"){
			        print "\[$time\]\: Validation\: $path\/$temp[0]_aln_sort_MD.bam exists.\n\[$time\]\: This file will start from \"haplotypecaller\" check point.\n";
			    }
			    else {
			        print "\[$time\]\: The bam file, $path\/$temp[0]_aln_sort_MD.bam, is ready. Nothing to do with this sample.\n";
			    }
			    $file_cnt++;
			}
		}
		elsif ($rad eq "1"){
		    unless (-e "$path\/$temp[0]_aln_sort.bam"){
			    $out = "pbrun fq2bam --ref $ref --no-markdups $L$f_l$IP\--out-bam $path\/$temp[0]_aln_sort.bam --num-gpus 4 --low-memory\\n";
			}
			else {
			    if ($ng eq "0"){
			        print "\[$time\]\: Validation\: $path\/$temp[0]_aln_sort.bam exists.\n\[$time\]\: This file will start from \"haplotypecaller\" check point.\n";
			    }
			    else {
			        print "\[$time\]\: The bam file, $path\/$temp[0]_aln_sort.bam, is ready. Nothing to do with this sample.\n";
			    }
			    $file_cnt++;			
			}
		}
		unless ($ng eq "1"){
		    $time = scalar localtime();
		    my $sp1_in;
			if ($rad eq "0"){
				$sp1_in = "$path\/$temp[0]_aln_sort_MD.bam";
			}
			else {
			    $sp1_in = "$path\/$temp[0]_aln_sort.bam";
			}
			unless (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz"){
			    $out .= "pbrun haplotypecaller --ref $ref --in-bam $sp1_in --gvcf --out-variants 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz --num-gpus 4\\n";
			    #$out .= "bgzip 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\\n";
		    }
		    else {
		        print "\[$time\]\: 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz exists.\n";
		    }
		    $time = scalar localtime();
			if (-e "02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi"){
				print "\[$time\]\: 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\.tbi exists.\n";
				next();
			}
			else {
				$out .= "pbrun indexgvcf --input 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.gz\\n";
				#$out .= "rm 02\-get_gvcf_$ran\/$temp[0]\.g\.vcf\.idx\\n";
			}
		}
		if ($out =~ /[a-z]/ig){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_docker nvcr-clara-parabricks-4401 -cj_gpu 4 -cj_time 168\:0\:0 -cj_qname gatk_02_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
		unless ($ng eq "1"){
		    $time = scalar localtime();
			my $sp1_in; my $stat;
			$out = "";
			if ($rad eq "0"){
				$sp1_in = "$path\/$temp[0]_aln_sort_MD.bam";
			}
			else {
			    $sp1_in = "$path\/$temp[0]_aln_sort.bam";
			}
			if (-e "$sp1_in\.bai"){
				print "\[$time\]\: $sp1_in exists.\n";
			}
			else {
				$out .= "samtools index \-\@ 2 $sp1_in\\n";
			}
			$stat = $sp1_in;
			$stat =~ s/bam$/samtools\.stats/;
			if (-e $stat){
			    my $check_s = `head -n 1 $stat`;
			    if ($check_s !~ /\w/){
			        $out .= "samtools stats \-\@ 2 $sp1_in > $stat\\n";
			    }
			}
			else {
			    $out .= "samtools stats \-\@ 2 $sp1_in > $stat\\n";
			}
			my $s_exc = $exc;
			unless ($local){
			    $s_exc = "";
			}
		    if ($out =~ /[a-z]/ig){
			    push(@stats, $out);
		    }
		}
		@temp = ();
	}
	unless ($local){
		close(BASH);
	}
	return ("02\-get_gvcf_$ran", $file_cnt, \@stats);
} #generate raw gvcf file using GPU-based GATK4
sub filter_vcf {
	my $proj = shift; my $gv = shift; my $ns = shift; my $ran = shift; my $ref = shift; my $path_o = shift; my $folder = shift; my $exc = shift; my $pre = shift; my $xlsn = shift; my $local = shift; my $p_dir = shift;
	my $time = scalar localtime();
	my @chrs;
	if ($ns eq "0"){ 
		@chrs = &check_chrs($path_o, $ref, $pre);
	}
	elsif ($ns eq "1" || $gv =~ /4/){
		@chrs = ("all");
	}
	unless ($local){
		open (BASH, ">my_bash_05_$ran\.sh") || die BOLD "Cannot write my_bash_05_$ran\.sh: $!", RESET, "\n";
	}
	foreach my $chr (@chrs){
		my $gz = "\.gz";
		my $out;
		if (-e "$folder\/select_vcf_star_filtered_$chr\.vcf\.gz"){
			print "\[$time\]\: select_vcf_star_filtered_$chr\.vcf\.gz file exists.\n";
			next();
		}
		elsif (-e "$folder\/vcf_star_filtered_$chr\.vcf\.gz"){
			print "\[$time\]\: vcf_star_filtered_$chr\.vcf\.gz file exists.\n";
			next();
		}
		elsif (-e "$folder\/select_vcf_raw_$chr\.vcf\.gz"){
			if (-e "$folder\/select_vcf_raw_$chr\.vcf\.gz\.tbi"){}
			else {
				print "\[$time\]\: WARNING\: Index file of $folder\/select_vcf_raw_$chr\.vcf\.gz doesn\'t exist. Try select_$folder\/vcf_raw_$chr\.vcf\n";
				if (-e "$folder\/select_vcf_raw_$chr\.vcf"){
					if (-e "$folder\/select_vcf_raw_$chr\.vcf\.idx"){
						$gz = "";
					}
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Index file of $folder\/select_vcf_raw_$chr\.vcf\.gz doesn\'t exist.", RESET, "\n";
					return 2;
				}
			}
		}
		elsif (-e "$folder\/vcf_raw_$chr\.vcf\.gz"){
			if (-e "$folder\/vcf_raw_$chr\.vcf\.gz\.tbi"){
			}
			else {
				print "\[$time\]\: WARNING\: Index file of $folder\/vcf_raw_$chr\.vcf\.gz doesn\'t exist. Try $folder\/vcf_raw_$chr\.vcf\n";
				if (-e "$folder\/vcf_raw_$chr\.vcf"){
					if (-e "$folder\/vcf_raw_$chr\.vcf\.idx"){
						$gz = "";
					}
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Index file of $folder\/vcf_raw_$chr\.vcf\.gz doesn\'t exist.", RESET, "\n";
					return 2;
				}					
			}
		}
		elsif ($exc){
			print RED "\[$time\]\: WARNING\: Cannot find $folder\/vcf_raw_$chr\.vcf\(\.gz\) or $folder\/select_vcf_raw_$chr\.vcf\(\.gz\)\.", RESET, "\n";
			print "\[$time\]\: Skip processing $folder\/vcf_raw_$chr\.vcf\(\.gz\) or $folder\/select_vcf_raw_$chr\.vcf\(\.gz\)\.\n";
			next();
		}
		else {
			print "\[$time\]\: $folder\/vcf_raw_$chr\.vcf\.gz has been simulated.\n";
			print "\[$time\]\: Please make sure you have execute step 4 before sending job\(s\).\n";
		}
		if (-e "$folder\/select_vcf_raw_$chr\.vcf$gz" || ($exc ne "-cj_exc " && $xlsn ne "0")){
			$out .= "perl $p_dir\/filter_vcf_2.4\.pl $folder\/select_vcf_raw_$chr\.vcf$gz\\n";
		}
		else {
			$out .= "perl $p_dir\/filter_vcf_2.4\.pl $folder\/vcf_raw_$chr\.vcf$gz\\n";
		}
		my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 24 -cj_qname gatk_05_$chr -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return 1;
} #filtering vcf files
sub GatherVcfs {
	my $dir = getcwd;
	my $proj = shift; my $xlsn = shift; my $gv = shift; my $ran = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $local = shift;
	my $time = scalar localtime();
	my @l_chrs; my @chrs; my $ppn;
	my $sel = "select_"; my $stpn;
	@chrs = &check_chrs($path_o, $ref, $pre);
	if ($gv =~ /4/){	
		$stpn = "vcf_raw_";
	}
	elsif ($gv =~ /5/){	
		$stpn = "vcf_star_filtered_";
	}
	elsif ($gv =~ /6/){	
		$stpn = "vcf_passed_";
	}
	my @l_chrs_tmp = <$folder\/$sel$stpn\*.vcf.gz>;
	if ($l_chrs_tmp[0] !~ /[a-z]/i && $exc){
		@l_chrs_tmp = <$folder\/$stpn\*.vcf.gz>;
		$sel = "";
	}
	elsif ($xlsn eq "0"){
		$sel = "";
	}
	elsif ($xlsn ne "0"){
		$sel = "select_";
	}
	foreach (@chrs){
		foreach $l_chr_tmp (@l_chrs_tmp){
			if ($l_chr_tmp =~ /$_\.vcf/){
				push(@l_chrs, $l_chr_tmp);
			}
		}
	}
	if ($l_chrs[0] !~ /\.vcf\.gz/ && $l_chrs[0] !~ /$stpn/){
		@l_chrs = ();
		@chrs = &check_chrs($path_o, $ref, $pre);
		foreach (@chrs){
			push(@l_chrs, "$folder\/$sel$stpn$_\.vcf\.gz");
		}
	}
	$time = scalar localtime();
	if (-e "$folder\/$sel$stpn\all\.vcf\.gz"){
		print "\[$time\]\: $folder\/$sel$stpn".'all'."\.vcf already exists\n";
		return 1;
	}
	chomp(@l_chrs);
	unless ($exc){
		unless (-d "$folder"){
			system("mkdir $folder");
			print "\[$time\]\: $sel\vcf_chr.list will be stored at $folder\.\n";
		}
	}
	unless (@l_chrs){
		print BOLD "\[$time\]\: ERROR\: Cannot find proper file\(s\) for gathering!", RESET, "\n";
		return 2;
	}
	open(INPUT, ">$folder\/$sel\vcf_chr.list") || die BOLD "\[$time\]\: ERROR\: Cannot generate $sel\vcf_chr\.list: $!", RESET, "\n";
	foreach my $l_chr (@l_chrs){
		print INPUT "$l_chr\n";
	}
	close(INPUT);
	$time = scalar localtime();
	print "\[$time\]\: GatherVcfs merges multiple vcfs into one vcf.\n";
	if ($gv =~ /2/){
		$gv = "Cloud";
		$ppn = 12;
	}
	else {
		$gv = "";
		$ppn = 3;
	}
	my $out = "\$gatk2 GatherVcfs$gv \-I $folder\/$sel\vcf_chr.list \-O $folder\/$sel$stpn".'all'."\.vcf\.gz\\n";
	$out .= "tabix $folder\/$sel$stpn".'all'."\.vcf\.gz\\n";
	$time = scalar localtime();
	unless ($local){
		&status($ran);			
	}
	&pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn $ppn -cj_mem 32 -cj_qname gatk_gv -cj_sn $ran -cj_qout . $out");
	return 1;
} #combine chrmosome vcfs into one vcf file
sub Grep_DP {
	my $time = scalar localtime();
	my $proj = shift; my $xlsn = shift; my $ns = shift; my $ran = shift; my $gv = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $pre_DP = shift; my $local = shift;
	my $sel = "select_";
	if ($files[0] !~ /[a-z]/i && $exc){
		@files = <$folder\/vcf_passed_*>;
		$sel = "";
	}
	elsif ($exc ne "-cj_exc " && $xlsn eq "0"){
		$sel = "";
	}
	elsif ($exc ne "-cj_exc " && $xlsn ne "0"){
		$sel = "select_";
	}	

	chomp(@files);
	my @chrs;
	if ($gv =~ /4|5|6/ || $ns eq "1"){
		if (-e "$folder\/$sel\vcf_passed_all\.vcf\.gz" || $exc ne "-cj_exc"){
			@chrs = ("all");
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find $folder\/$sel\vcf_passed_all\.vcf\.gz file.", RESET, "\n";
			return 2;
		}
	}
	else {
		if ($files[0] !~ /vcf_passed_/){
			if ($exc){
				print BOLD "\[$time\]\: ERROR\: Cannot find vcf_passed vcf file\(s\).", RESET, "\n";
				return 2;
			}
			else {
				@chrs = &check_chrs($path_o, $ref, $pre);
			}
		}
		elsif ($files[0] =~ /vcf_passed_/) {
			foreach my $file (@files){
				if ($file =~ /$sel\vcf_passed_all\.vcf/){
					next();
				}
				$file =~ s/$folder//g;
				$file =~ s/\/+//g;
				$file =~ s/$sel//g;
				$file =~ s/vcf_passed_//g;
				$file =~ s/\.vcf\.gz$//;
				push(@chrs, $file);
			}
		}
	}
	@files = ();
	if ($chrs[0] =~ /$pre_DP/i){
		$pre_DP = substr $chrs[0], 0, 3;
		$pre_DP =~ s/[0-9]//;
	}
	unless ($local){
		open (BASH, ">my_bash_07_$ran\.sh") || die BOLD "\[$time\]\: ERROR\: Cannot write my_bash_07_$ran\.sh: $!", RESET, "\n";
	}
	foreach my $chr (@chrs){
		my $out;
		if (-e "$folder\/$sel\vcf_passed_DP_$chr\.txt"){
			if ($ow eq "0"){
				print "\[$time\]\: WARNING\: $folder\/$sel\vcf_passed_DP_$chr\.txt exists. Skip grepping \"DP\" in the vcf file.\n";
				next();
			}
			elsif ($ow eq "1"){
				$out = "rm $folder\/$sel\vcf_passed_DP_$chr\.txt\\n";
			}
		}
		$out .= "zgrep \"\^$pre_DP\" $folder\/$sel\vcf_passed_$chr\.vcf\.gz \| cut \-f 8 \-d\$\'\\t\' \| zgrep \-P \'\(\?\<\=DP\=\)\\d\+\' \-o \> $folder\/$sel\vcf_passed_DP_$chr\.txt\\n";
		my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 8 -cj_qname gatk_07_$chr -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return 1;
} #grep DP field of vcf files
sub bi_allele {
	my $time = scalar localtime();
	my $proj = shift; my $keep = shift; my $sf = shift; my $xlsn = shift; my $minQ = shift; my $mm = shift; my $ri = shift; my $ba = shift; my $ns = shift; my $ran = shift; my $gv = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $local = shift; my $maf = shift;
	my $filtered;
	if ($sf eq "1"){
		$filtered = "raw";
	}
	else {
		$filtered = "star_filtered";
	}	
	my @files = <$folder\/select_vcf_$filtered\_*.gz>;
	if ($files[0] !~ /[a-z]/i && $exc){
		@files = <$path_o\/select_vcf_$filtered\_*.gz>;
	}
	my $sel = "select_";
	#check vcftools
	my @vcft = `vcftools`; my $check = 0;
	chomp(@vcft);
	foreach (@vcft){
		if ($_ =~ /^VCFtools/){
			$check = 1;
		}
	}
	if ($check == 0){
		print BOLD "\[$time\]\: ERROR\: Cannot find vcftools. Please install vcftools first or set vcftools path to the environment.", RESET, "\n";
		return 2;
	}
	
	if ($files[0] !~ /[a-z]/i && $exc){
		@files = <$folder\/vcf_$filtered\_*.gz>;
		if ($files[0] !~ /[a-z]/i && $exc){
			@files = <$path_o\/vcf_$filtered\_*.gz>;
		}
		$sel = "";
	}
	elsif ($xlsn ne "0" && $exc ne "-cj_exc "){
		$sel = "select_";
	}
	elsif ($xlsn eq "0" && $exc ne "-cj_exc "){
		$sel = "";
	}
	chomp(@files);
	my @chrs;
	if ($gv =~ /4|5/ || $ns eq "1"){
		if (-e "$folder\/$sel\vcf_$filtered\_all\.vcf\.gz" || $exc ne "-cj_exc "){
			@chrs = ("all");
		}
		else {
			print BOLD "\[$time\]\: ERROR\: Cannot find $folder\/$sel\vcf_$filtered\_all\.vcf\.gz file.", RESET, "\n";
			return 2;
		}
	}
	else {
		if ($files[0] !~ /vcf_$filtered\_/){
			if ($exc){
				print BOLD "\[$time\]\: ERROR\: Cannot find $filtered vcf\.gz file\(s\).", RESET, "\n";
				return 2;
			}
			else {
				@chrs = &check_chrs($path_o, $ref, $pre);
			}
		}
		else {
			foreach my $file (@files){
				if ($file =~ /$sel\$filtered\_all\.vcf\.gz/){
					next();
				}
				if ($file =~ /$folder/){
					$file =~ s/$folder//g;
				}
				if ($file =~ /$path_o/){
					$file =~ s/$path_o//g;
				}				
				$file =~ s/\/+//g;
				$file =~ s/$sel//g;
				$file =~ s/vcf_$filtered\_//g;
				$file =~ s/\.vcf\.gz$//;
				push(@chrs, $file);
			}
		}
	}
	@files = ();
	unless ($local){
		open (BASH, ">my_bash_06_$ran\.sh") || die BOLD "\[$time\]\: ERROR\: Cannot write my_bash_06_$ran\.sh: $!", RESET, "\n";
	}
	foreach my $chr (@chrs){
		my $out;
		if (-e "$folder\/$sel\vcf_passed_$chr\.vcf\.gz"){
			if ($ow eq "0"){
				print "\[$time\]\: WARNING\: $folder\/$sel\vcf_passed_$chr\.vcf\.gz exists. Skip further filtering by vcftools.\n";
				next();
			}
			elsif ($ow eq "1"){
				$out = "rm $folder\/$sel\vcf_passed_$chr\.vcf\.gz\\n";
			}
		}
		if ($ri eq "0"){
			$ri = "";
		}
		if ($ba eq "0"){
			$ba = "";
		}
		if ($mm eq "0"){
			$mm = "";
		}
		if ($minQ eq "0"){
			$minQ = "";
		}
		if ($keep eq "0"){
			$keep = "";
		}
		$out .= "vcftools \-\-gzvcf $folder\/$sel\vcf_$filtered\_$chr\.vcf\.gz $ri$ba$mm$minQ$maf$keep\-\-recode \-\-recode\-INFO\-all \-\-stdout \| bgzip -c \> $folder\/$sel\vcf_passed_$chr\.vcf\.gz\\n";
		my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 8 -cj_qname gatk_06_$chr -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return 1;
} #get bi-allelic SNPs and further filtering

#GATK3
sub trim_3 {
	my $proj = shift; my $ran = shift; my $exc = shift; my $path = shift; my $adp = shift; my $local = shift;
	my $time = scalar localtime();
	my @files = <$path\/*.fastq.gz>;
	my $gz = "\.gz";
	if ($files[0] !~ /fastq\.gz/){
		@files = <$path\/*.fastq>;
		$gz = ();
	}
	my $cnt = 0;
	my @temp;
	chomp(@files);
	if ($files[0] !~ /fastq/){
		print RED "\[$time\]\: WARNING\: Cannot find appropriate file to proceed step 1.\n", RESET;
		return 1;
	}
	if (-d "01\-fq_trim_$ran"){
		print "\[$time\]\: The trimmed file(s) are stored at 01\-fq_trim_$ran\n";
	}
	else{
		system("mkdir 01\-fq_trim_$ran");
		print "\[$time\]\: Make a folder\: 01\-fq_trim_$ran\n\[$time\]\: The trimmed file(s) will be stored there.\n";
	}
	my @temp; my $file_2;
	my $cnt = 0;
	my ($str, $type);
	unless ($local){
		open (BASH, ">my_bash_01_$ran\.sh") || die BOLD "Cannot write my_bash_01_$ran\.sh: $!", RESET, "\n";
	}
	my @re_list;
	foreach my $file (@files){
		my $out;
		$time = scalar localtime();
		next if ($file =~ /_R2/);
		if ($file =~ /_R1/){
			@temp = split(/_R1/, $file);
			$temp[0] =~ s/$path\///;
		}
		else {
			@temp = split(/\./, $file);
			$temp[0] =~ s/$path\///;
		}
		REDO:
		if (-e "01\-fq_trim_$ran\/$temp[0]\_aln_sort\.bam" || -e "01\-fq_trim_$ran\/$temp[0]\_aln_sort_MD\.bam"){
			print "\[$time\]\: step 2 file of $temp[0] exists. Skip this sample from step 1.\n";
		}		
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz"){
			print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz exists.\n";
		}
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz"){
			if (-e "01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed$gz"){
				print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz and 01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed$gz exist.\n";
			}
			else {
				system("rm 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz");
				goto REDO;
			}
		}
		else{
			if ($file =~ /_R1/){
				print "\[$time\]\: $file is a paired-end file.\n";
				$file_2 = $file;
				$file_2 =~ s/_R1/_R2/;
#				print "$file_2\n";
				unless (-e "$file_2" || $exc ne "-cj_exc "){
					print RED "\[$time\]\: WARNING\: It seems that the paired file of $file_2 is missing!\n", RESET;
					print RED "\[$time\]\: Skip trimming of $file.\n", RESET;
					next();
				}			
				$file_2 = " $file_2";
				$str = "01\-fq_trim_$ran\/$temp[0]\_R1\.cutadapt\.fastq$gz \-p 01\-fq_trim_$ran\/$temp[0]\_R2\.cutadapt\.fastq$gz 01\-fq_trim_$ran\/$temp[0]\_R1\.fastq\.trimmed$gz 01\-fq_trim_$ran\/$temp[0]\_R2\.fastq\.trimmed$gz";
			}
			else{
				print "\[$time\]\: $file is a single-end file.\n";
				$file_2 = ();
				$str = "01\-fq_trim_$ran\/$temp[0]\.cutadapt\.fastq$gz \-p 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz";
			}
			$out .= "SolexaQA\+\+ dynamictrim $file$file_2 -h 10 -d 01\-fq_trim_$ran\\n";
			if ($adp eq "1"){
				$out .= "cutadapt \-m 30 \-f fastq \-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \-o $str\\n";
				if ($file =~ /_R1/){
					$out .= "rm 01\-fq_trim_$ran\/$temp[0]\*\.fastq\.trimmed$gz\\n";
					$out .= "mv 01\-fq_trim_$ran\/$temp[0]\_R1\.cutadapt\.fastq$gz 01\-fq_trim_$ran\/$temp[0]\_R1\.fastq\.trimmed$gz\\n";
					$out .= "mv 01\-fq_trim_$ran\/$temp[0]\_R2\.cutadapt\.fastq$gz 01\-fq_trim_$ran\/$temp[0]\_R2\.fastq\.trimmed$gz\\n";
				}
				else{
					$out .= "rm 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz\\n";
					$out .= "mv 01\-fq_trim_$ran\/$temp[0]\.cutadapt\.fastq$gz 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz\\n";
				}
			}
		}
		if ($out =~ /[a-z]/i){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 4 -cj_qname gatk_01_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
		push(@re_list, $temp[0]);
	}
	unless ($local){
		close(BASH);
	}
	return @re_list;
} #trim fastq files
sub GenotypeGVCFs_3 {
	my $proj = shift; my $nlc = shift; my $ns = shift; my $ow = shift; my $ran = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $dblist = shift; my $as = shift; my $local = shift;
	my $dir = getcwd;
	my $path = "02\-get_gvcf_$ran";
	my $time = scalar localtime();
	my $j_che; my @che;
	if (-d $path){
		opendir DIR, $path || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path. Using path: $path_o\n";
			$path = $path_o;
			opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
			@che = readdir DIR; chomp(@che);
			close(DIR);
			$j_che = join(" ", @che);
			if ($j_che !~ /g.vcf.gz/){
				unless ($exc){
					print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
					exit;
				}
			}
		}
	}
	else {
		$path = $path_o;
		opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			unless ($exc){
				print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
				exit;
			}
		}
	}
	if (-d $path){
		unless ($exc){
			print "\[$time\]\: The list file will be stored at $path.\n";
		}
		else {
			print BOLD "\[$time\]\: WARNING: Make sure you have executed step 2.", RESET, "\n";
		}
	}
	if ($as eq "\-all\-sites"){
		$as = "\-allSites";
	}
	my @lists;
	if (-e $dblist){
		print "\[$time\]\: gvcf list is loaded.\n";
		open(DBLIST, "<$dblist") || die BOLD "\[$time\]\: ERROR\: Cannot open $dblist: $!", RESET, "\n";
		my @db_line = <DBLIST>;
		chomp(@db_line);
		my @valis;
		foreach (@db_line){
			@valis = split("\t", $_);
			chomp(@valis);
			if (-e $valis[1]){
				push(@lists, $valis[1]);
			}
			else{
				print RED "\[$time\]\: WARNING\: Cannot find the path of $valis[0]. Skip importing this sample.\n", RESET;
			}
		}
		close(DBLIST);
	}
	else{
		print "\[$time\]\: gvcf list is not defined.\n\[$time\]\: Using defined folder files as a gvcf list.\n";
		@lists = <$path\/*.vcf.gz>;
		chomp(@lists);
		if ($lists[0] !~ /vcf\.gz/){
			@lists = <$path_o\/*.vcf.gz>;
			if ($lists[0] !~ /vcf\.gz/){
				if ($exc){
					print BOLD "\[$time\]\: Cannot find \*\.vcf\.gz file!", RESET, "\n";
					return 2;
				}
				else {
					print "\[$time\]\: WARNING\: Cannot find \*\.vcf\.gz file\(s\)\! Make sure you have executed previous step before sending job\(s\).\n";
					@lists = <$path_o\/*.fastq.gz>;
					chomp(@lists);
					for (my $z=0; $z<=$#lists; $z++){
						if ($lists[$z] =~ /_R2/){
							$lists[$z] = undef;
							next();
						}
						if ($lists[$z] =~ /_R1/){
							$lists[$z] =~ s/_R1//gi;
						}
						$lists[$z] =~ s/fastq\.gz/g\.vcf\.gz/;
						$lists[$z] =~ s/$path_o\///;
						$lists[$z] = "$path\/$lists[$z]";
					}
				}
			}
			$path = $path_o;
		}
	}
	my @temp; my @temp_l;  my @fi_lists; my $fi_tb; my $out_list;
	open(LIST, ">$path\/gvcf_list\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $path\/gvcf_list\.list: $!", RESET, "\n";
	open(PATH, ">$path\/sample_path\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $path\/sample_path\.list: $!", RESET, "\n";
	foreach my $list (@lists){
		$time = scalar localtime();
		@temp = split (/\//, $list);
		@temp_l = $temp[-1];
		$temp[0] =~ s/\.gz$//;
		$temp[0] =~ s/\.g\.vcf$//;
		if ($list !~ /[a-z]/i){
			next();
		}
		if (-e "$folder\/c_vcf\.list" && $nlc eq "0"){
			open(FILE, "<$folder\/c_vcf\.list") || die BOLD "\[$time\]\: ERROR\: Cannot open $folder\/c_vcf\.list: $!", RESET, "\n";
			my $cc = <FILE>;
			chomp($cc);
			$time = scalar localtime();
			my @ccs = split(/\t/, $cc);
			my $exist = 0;
			foreach my $x (0..$#ccs){
				if ($ccs[$x] =~ /\b$temp_l[0]\b/g && $ow eq "0"){
					print "\[$time\]\: WARNING: Sample $temp_l[0] already in the gvcf file\n";
					print "\[$time\]\: WARNING: Skip combining $temp_l[0] into gvcf file.\n";
					$exist = 1;
					last;
				}
			}
			if ($exist == 1){
				next;
			}
		}
		$fi_tb = 0;
		foreach my $x (0..$#fi_lists){
			my @spl = split(/\//, $fi_lists[$x]);
			@spl = $spl[-1];
			$spl[0] =~ s/\.gz$//;
			$spl[0] =~ s/\.g\.vcf$//;
			if ($spl[0] eq $temp_l[0]){
				print "\[$time\]\: WARNING: Sample $temp_l[0] has been detected before.\n";
				print "\[$time\]\: WARNING: Skip combining the duplicated $temp_l[0] into gvcf file.\n";
				$fi_tb = 1;
				last;
			}
		}
		if ($if_tb == 1){
			next;
		}
		if ($list !~ /^\./ && $dblist eq "0"){
			$list = "$dir\/$list";
		}
		print LIST "$temp_l[0]\t$list\n";
		print PATH "$list\n";
		$out_list = "$out_list"."$temp_l[0]\t";	
	}
	close(LIST);
	close(PATH);
	open (PATH, "<$path\/sample_path\.list") || die BOLD "Cannot open $path\/sample_path\.list: $!", RESET, "\n";

	my $path_con = <PATH>;
	close(PATH);
	if ($path_con !~ /[a-z]/i){
		print RED "\[$time\]\: WARNING\: No sample to process!", RESET, "\n";
		system("rm $path\/sample_path\.list");
		return 1;
	}
	@lists = (); @temp = (); @temp_l = ();
	$time = scalar localtime();
	if (-d $folder){
		print "\[$time\]\: The gvcf file(s) will be stored at $folder\n";
	}
	else{
		print "\[$time\]\: Make a folder\: $folder.\n\[$time\]\: The star_filtered vcf file will be stored there.\n";
		system("mkdir $folder");
	}
	my @chrs;
	if ($ns eq "0"){
		@chrs = &check_chrs($path_o, $ref, $pre);
	}
	else {
		@chrs = ("all");
	}
	unless ($local){
		open (BASH, ">my_bash_04_$ran\.sh") || die BOLD "Cannot write my_bash_04_$ran\.sh: $!", RESET, "\n";
	}
	my $L;
	foreach my $chr (@chrs){
		if ($ns eq "0"){
			$L = "\-L $chr ";
		}
		elsif ($ns eq "1"){
			$L = "";
		}
		if (-e "$folder\/vcf_raw_$chr\.vcf.gz"){
			print "\[$time\]\: File $folder\/vcf_raw_$chr\.vcf\.gz exists. Skip the step.\n";
			next();
		}
		if (-e "$ref\.fai"){}
		else {
			unless ($exc){
				print "\[$time\]\: WARNING: Cannot find indexed reference file at 00\-ref folder. Make sure you have execute step 0 before sending the job.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR: Cannot find the indexed reference file.", RESET, "\n";
				return 2;
			}
		}
		my $out = "java \-jar \$GATKFILE \-T GenotypeGVCFs \-nt 4 \-R $ref $L\-\-variant $path\/sample_path\.list \-o $folder\/vcf_raw_$chr\.vcf\.gz $as\\n";
		my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 4 -cj_mem 24 -cj_qname gatk_04_$chr -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return ($out_list);
} #compile gvcfs

#GATK4
sub trim_4 {
	my $proj = shift; my $ran = shift; my $exc = shift; my $path = shift; my $adp = shift; my $local = shift;
	my $time = scalar localtime();
	my @files = <$path\/*.fastq.gz>;
	my $gz = "\.gz";
	my $cnt = 0;
	my @temp;
	my ($file_2, $str1, $str2, $c_a, $type);
	
	if ($files[0] !~ /fastq\.gz/){
		@files = <$path\/*.fastq>;
		$gz = ();
	}
	chomp(@files);
	if ($files[0] !~ /fastq/){
		print RED "\[$time\]\: WARNING\: Cannot find appropriate file to proceed step 1.\n";
		return 1;
	}
	if (-d "01\-fq_trim_$ran"){
		print "\[$time\]\: The trimmed file(s) are stored at 01\-fq_trim_$ran\n";
	}
	else{
		system("mkdir 01\-fq_trim_$ran");
		print "\[$time\]\: Make a folder\: 01\-fq_trim_$ran\n\[$time\]\: The trimmed file(s) will be stored there.\n";
	}
	unless ($local){
		open (BASH, ">my_bash_01_$ran\.sh") || die BOLD "Cannot write my_bash_01_$ran\.sh: $!", RESET, "\n";
	}
	my @re_list;
	foreach my $file (@files){
		$time = scalar localtime();
		next if ($file =~ /_R2/);
		my $out;
		if ($file =~ /_R1/){
			@temp = split(/_R1/, $file);
			$temp[0] =~ s/$path\///;
		}
		else {
			@temp = split(/\./, $file);
			$temp[0] =~ s/$path\///;
		}
		REDO:
		if (-e "01\-fq_trim_$ran\/$temp[0]\_aln_sort\.bam" || -e "01\-fq_trim_$ran\/$temp[0]\_aln_sort_MD\.bam"){
			print "\[$time\]\: step 2 file of $temp[0] exists. Skip this sample from step 1.\n";
		}
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz"){
			print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz exists.\n";
		}
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz"){
			if (-e "01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed$gz"){
				print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz and 01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed$gz exist.\n";
			}
			else {
				system("rm 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed$gz");
				goto REDO;
			}
		}
		else{
			if ($file =~ /_R1/){
				$file_2 = $file;
				$file_2 =~ s/_R1/_R2/;
				unless (-e $file_2 || $exc ne "-cj_exc "){
					print RED "\[$time\]\: WARNING\: It seems that the paired file of $file_2 is missing!\n", RESET;
					print RED "\[$time\]\: Skip trimming of $file.\n", RESET;
					next();
				}			
				$file_2 = " $file_2";
				$type = "PE";
				$str1 = "01\-fq_trim_$ran\/$temp[0]\_R1\.fastq\.trimmed$gz 01\-fq_trim_$ran\/$temp[0]\_R1\.fastq\.unpaired$gz";
				$str2 = " 01\-fq_trim_$ran\/$temp[0]\_R2\.fastq\.trimmed$gz 01\-fq_trim_$ran\/$temp[0]\_R2\.fastq\.unpaired$gz";
			}
			else{
				$type = "SE";
				$file_2 = ();
				$str1 = "01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed$gz";
				$str2 = "";
			}
			if ($adp ne ""){
				$c_a = " ILLUMINACLIP\:$adp\:2\:30\:10";
			}
			elsif ($adp eq ""){
				$c_a = "";
			}
			$out .= "java -jar \$TRIMMO $type \-threads 2 $file$file_2 $str1$str2$c_a LEADING\:3 TRAILING\:3 SLIDINGWINDOW\:4\:10 MINLEN\:30\\n";
		}
		if ($out =~ /[a-z]/i){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 2 -cj_mem 48 -cj_qname gatk_01_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
		push(@re_list, "$temp[0]");
	}
	unless ($local){
		close(BASH);
	}
	return @re_list;
} #trim fastq files
sub trim_fastp {
	my $proj = shift; my $ran = shift; my $exc = shift; my $path = shift; my $local = shift;
	my $time = scalar localtime();
	my $gz = "\.gz";
	my @files = <$path\/*.fastq.gz>;
	unless (@files){
		@files = <$path\/*.fastq.bz2>;
		$gz = "\.bz2";
		unless (@files){
			@files = <$path\/*.fastq>;
			$gz = "";
		}
	}
	my $cnt = 0;
	my @temp;
	my ($file_2, $str1, $str2, $c_a, $type);
	chomp(@files);
	if ($files[0] !~ /fastq/){
		print RED "\[$time\]\: WARNING\: Cannot find appropriate file to proceed step 1.\n";
		return 1;
	}
	if (-d "01\-fq_trim_$ran"){
		print "\[$time\]\: The trimmed file(s) are stored at 01\-fq_trim_$ran\n";
	}
	else {
		system("mkdir 01\-fq_trim_$ran");
		print "\[$time\]\: Make a folder\: 01\-fq_trim_$ran\n\[$time\]\: The trimmed file(s) will be stored there.\n";
	}
	unless ($local){
		open (BASH, ">my_bash_01_$ran\.sh") || die BOLD "Cannot write my_bash_01_$ran\.sh: $!", RESET, "\n";
	}
	my @re_list;
	foreach my $file (@files){
		$time = scalar localtime();
		next if ($file =~ /_R2/);
		my $out;
		if ($file =~ /_R1/){
			@temp = split(/_R1/, $file);
			$temp[0] =~ s/$path\///;
		}
		else {
			@temp = split(/\./, $file);
			$temp[0] =~ s/$path\///;
		}
		REDO:
		if (-e "01\-fq_trim_$ran\/$temp[0]\_aln_sort\.bam" || -e "01\-fq_trim_$ran\/$temp[0]\_aln_sort_MD\.bam"){
			print "\[$time\]\: step 2 file of $temp[0] exists. Skip this sample from step 1.\n";
		}
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed.gz"){
			print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed.gz exists.\n";
		}
		elsif (-e "01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed.gz"){
			if (-e "01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed.gz"){
				print "\[$time\]\: 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed.gz and 01\-fq_trim_$ran\/$temp[0]\_R2.fastq\.trimmed.gz exist.\n";
			}
			else {
				system("rm 01\-fq_trim_$ran\/$temp[0]\_R1.fastq\.trimmed.gz");
				goto REDO;
			}
		}
		else{
			if ($file =~ /_R1/){
				$file_2 = $file;
				$file_2 =~ s/_R1/_R2/;
				unless (-e $file_2 || $exc ne "-cj_exc "){
					print RED "\[$time\]\: WARNING\: It seems that the paired file of $file_2 is missing!\n", RESET;
					print RED "\[$time\]\: Skip trimming of $file.\n", RESET;
					next();
				}
				$file = "--in1 $file ";
				$file_2 = "--in2 $file_2 ";
				$str1 = "--out1 01\-fq_trim_$ran\/$temp[0]\_R1\.fastq\.trimmed.gz ";
				$str2 = "--out2 01\-fq_trim_$ran\/$temp[0]\_R2\.fastq\.trimmed.gz ";
			}
			else{
				$file = "--in1 $file ";
				$file_2 = ();
				$str1 = "--out1 01\-fq_trim_$ran\/$temp[0]\.fastq\.trimmed.gz ";
				$str2 = ();
			}
			$out .= "fastp $file$file_2$str1$str2\-l 50 -g -h 01\-fq_trim_$ran\/$temp[0]\.html \&\> 01\-fq_trim_$ran\/$temp[0]\.log\\n";
		}
		if ($out =~ /[a-z]/i){
			$cnt += 1;
			my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 2 -cj_mem 48 -cj_qname gatk_01_$cnt -cj_sn $ran -cj_qout . $out");
			print BASH "$return\n";
		}
		push(@re_list, "$temp[0]");
	}
	unless ($local){
		close(BASH);
	}
	return @re_list;		
}
sub db_import_4 {
#	print "debug: db_import_4\n";
	my $dir = getcwd;
	my $proj = shift; my $wes = shift; my $interval = shift; my $ip = shift; my $ns = shift; my $folder = shift; my $ran = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $DB_path = shift; my $dblist = shift; my $ow = shift; my $local = shift;
	my $time = scalar localtime();
	my $DB_c = "0"; my @dir; my $dir_j; my $ow_db;
	my $cnt_int = 0;
	my $path = "02\-get_gvcf_$ran";
	my $j_che; my @che;
	if (-d $path){
		opendir DIR, $path || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path. Using path: $path_o\n";
			$path = $path_o;
			opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
			@che = readdir DIR; chomp(@che);
			close(DIR);
			$j_che = join(" ", @che);
			if ($j_che !~ /g.vcf.gz/){
				unless ($exc){
					print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
					exit;
				}
			}
		}
	}
	else {
		$path = $path_o;
		opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			unless ($exc){
				print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
				exit;
			}
		}
	}
	my @chrs = &check_chrs($path, $ref, $pre); 
	if ($ns eq "1"){
		@chrs = ("all");
	}
	unless (-d $path){
		unless ($exc){
			print "\[$time\]\: The list file will be stored at $path.\n";
		}
		else {
			print "\[$time\]\: WARNING: Make sure you have executed step 2.\n";
		}
	}
	my @DB_ps = split(/\//, $DB_path);
	chomp(@DB_ps);
	my $DB_name = pop(@DB_ps);
	my $DB_p = join("\/", @DB_ps);
	if ($DB_path !~ /\//){
		$DB_p = "$dir\/";
	}
	else{
		$DB_p = $DB_p."\/";
	}
	foreach my $chr (@chrs){
		if (-d "$DB_path\_$chr"){
			print "\[$time\]\: $chr database path set to \"$DB_path\_$chr\"\n";
			opendir DIR, "$DB_path\_$chr" || die BOLD "\[$time\]\: ERROR\: Cannot open $DB_path\_$chr: $!", RESET, "\n";
			@dir = readdir DIR;
			chomp(@dir);
			$dir_j = join("", @dir);
			$time = scalar localtime();
			if ($dir_j !~ /vidmap\.json/gi){
				print BOLD "\[$time\]\: ERROR\: $DB_path\_$chr has some files that are not database files.\n\[$time\]\: Database path should contain no file or database file\(s\)\.", RESET, "\n";
				return 2;
			}
			else{
				if ($ow == "1"){
					print RED "\[$time\]\: WARNING: Over-write $DB_path\_$chr. Exist files will be deleted.\n", RESET;
					$ow_db = "\-\-overwrite\-existing\-genomicsdb\-workspace ";
					$DB_c = "0";
				}
				elsif ($ow == "0"){
					print RED "\[$time\]\: WARNING: The folder has some files!\n", RESET;
					$DB_c = "1";
				}
			}
		}
		else{
			$time = scalar localtime();
			print "\[$time\]\: $chr database path set to \"$DB_path\_$chr\"\n";
			$DB_c = "0";
			if ("$DB_path\_$chr" eq "GenomicsDB\_$chr"){
				print "\[$time\]\: WARNING: DB folder is not specified or the same with default.\n\[$time\]\: DB files will be stored at GenomicsDB_$chr.\n";
			}
			my $err = `mkdir $DB_path\_$chr 2>&1`;
			if ($err =~ /[a-z]/i){
				print BOLD "\[$time\]\: ERROR\: Cannot make a folder \"$DB_path\_$chr\" for database storage.", RESET, "\n";
				return 2;
			}
			if (-d "$DB_path\_$chr"){
				system("rm -r $DB_path\_$chr");
			}
		}		
	}

	my @lists;
	$time = scalar localtime();
	if (-e $dblist){
		print "\[$time\]\: gvcf list is loaded.\n";
		open(DBLIST, "<$dblist") || die BOLD "\[$time\]\: ERROR\: Cannot open $dblist: $!", RESET, "\n";
		my @db_line = <DBLIST>;
		chomp(@db_line);
		my @valis;
		foreach (@db_line){
			@valis = split("\t", $_);
			chomp(@valis);
			if (-e $valis[1]){
				push(@lists, $valis[1]);
			}
			else{
				print RED "\[$time\]\: WARNING\: Cannot find the path of $valis[0]. Skip importing this sample.\n", RESET;
			}
		}
		close(DBLIST);
	}
	else{
		print "\[$time\]\: gvcf list is not defined.\n\[$time\]\: Using files in the defined path as a gvcf list.\n";
		@lists = <$path\/*.vcf.gz>;
		chomp(@lists);	
		if ($lists[0] !~ /vcf\.gz/){
			if ($exc){	
				print BOLD "\[$time\]\: ERROR\: Cannot find \*\.vcf\.gz file!", RESET, "\n";
				return 2;
			}
			else {
				print "\[$time\]\: WARNING\: Cannot find \*\.vcf\.gz file\(s\)\! Make sure you have \*\.vcf\.gz file before execute qsub jobs.\n";
				@lists = <01\-fq_trim_$ran\/*.fastq.gz>;
				unless (@lists){
					print "\[$time\]\: ERROR: Cannot find the substituted path to generate the command line. Please run the pipeline with previous step first.\n";
					exit;
				}
				for (my $z=0; $z<=$#lists; $z++){
					if ($lists[$z] =~ /_R2/){
						$lists[$z] = undef;
						next();
					}
					if ($lists[$z] =~ /_R1/){
						$lists[$z] =~ s/_R1//gi;
					}
					$lists[$z] =~ s/fastq\.gz/g\.vcf\.gz/;
					$lists[$z] =~ s/$path_o\///;
					$lists[$z] = "$path\/$lists[$z]";
				}
			}
			$path = $path_o;
			print "\[$time\]\: Using vcfs in $path_o.\n";
		}
	}
	my @temp; my @temp_l; my @im_lists; my $e_list; my $im_tb;  my @o_lists;
	$time = scalar localtime();
	if (-e "$DB_p\$DB_name\_samples\.list"){
		open(IN, "<$DB_p\$DB_name\_samples\.list") || die BOLD "\[$time\]\: ERROR\: Failed to open $DB_p\$DB_name\_samples\.list: $!", RESET, "\n";
		$e_list = <IN>;
		chomp($e_list);
		close(IN);
	}
	foreach my $list (@lists){
		$time = scalar localtime();
		@temp = split (/\//, $list);
		@temp_l = split (/\./, $temp[-1]);
		if ($e_list){
			if ($e_list =~ /\b$temp_l[0]\b/ && $ow == "0"){
				print "\[$time\]\: WARNING: $temp_l[0] already in the database. Skip importing this sample.\n";
				next();
			}
		}
		$im_tb = join("\t", @im_lists);
		if ($im_tb =~ /\b$temp_l[0]\b/g){
			print "\[$time\]\: WARNING: Sample $temp_l[0] has been detected before.\n";
			print "\[$time\]\: WARNING: Skip importing the duplicated $temp_l[0] into database.\n";
			next();
		}
		if ($list !~ /[a-z]/i){
			next();
		}
		if ($list !~ /^\./ && $dblist eq "0"){
			$list = "$dir\/$list";
		}
		push(@o_lists, "$temp_l[0]\t$list");
		push(@im_lists, $temp_l[0]);
	}
	@temp = ();
	@temp_l = ();
	my $im_list = join("\t", @im_lists);

	my $k; my $ct = 0; my $start = 0;
	my $list_num = int(($#o_lists+1)/50);
	$list_num++ if ((($#o_lists+1)/50) > int(($#o_lists+1)/50));
#	print "Number: $list_num\n";
	my $L; my $IP;
	for (my $i=1; $i<=$list_num; $i++){
		if ($list_num == 1){
			$k = "";
		}
		else {
			$k = "_$i";
		}
		open(LIST, ">$path\/gvcf_list$k\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $path\/gvcf_list$k\.list: $!", RESET, "\n";
		if ($i < $list_num){
			for (my $n=$start; $n<($start+50); $n++){
				print LIST "$o_lists[$n]\n";
			}
		}
		else {
			for (my $n=$start; $n<$#o_lists; $n++){
				print LIST "$o_lists[$n]\n";
			}		
		}
		close(LIST);
		my $db_eq; my $dict; my $mcinp;
		if ($ns eq "1"){
			if ($cnt_int >= 100){
				$mcinp = "\-\-merge\-contigs\-into\-num\-partitions 15";
			}
			$L = "\-\-merge\-input\-intervals \-\-max\-num\-intervals\-to\-import\-in\-parallel 4 $mcinp";
			@chrs = ("all");
		}
		$time = scalar localtime();
		if (-e "$path\/gvcf_list$k\.list"){
			unless ($local){
				open (BASH, ">my_bash_03$k\_$ran\.sh") || die BOLD "\[$time\]\: ERROR\: Cannot write my_bash_03$k\_$ran\.sh: $!", RESET, "\n";
			}
			foreach my $chr (@chrs){
				my $out;
				if ($wes){
					$L = "-L $interval";
					if ($ip){
						$IP = "-ip $ip ";
					}
				}
				else {
					$L = "\-L $chr";
				}
				if ($exc){
					print "\[$time\]\: Create\\Update $DB_name\_$chr$k.\n";
				}
				if ($DB_c == "0"){
					if ($i == 1){
						$db_eq = "\-\-genomicsdb\-workspace\-path";
					}
					else {
						$db_eq = "\-\-genomicsdb\-update\-workspace\-path";
						$ow_db = "";
					}
				}
				elsif ($DB_c == "1"){
					$db_eq = "\-\-genomicsdb\-update\-workspace\-path";
				}
				$dict = $ref;
				$dict =~ s/fasta$|fas$|fa$/dict/;
				$out = "TILEDB_DISABLE_FILE_LOCKING\=1 \$gatk2 \-\-java\-options \"\-Xmx60g\" GenomicsDBImport $db_eq $DB_path\_$chr \-\-sample\-name\-map $path\/gvcf_list$k\.list \-\-sequence\-dictionary $dict $ow_db\-\-reader\-threads 12 $IP$L\\n";
				my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 6 -cj_mem 60 -cj_qname gatk_03_$chr$k -cj_sn $ran -cj_qout . $out");
				print BASH "$return\n";
				#close(INPUT);
			}
			unless ($local){
				close(BASH);
			}
		}
		else{
			die BOLD "\[$time\]\: ERROR\: gvcf_list$k\.list is missing!", RESET, "\n";
		}
		$start += 50;
	}
	print "[$time\]\: Please see the imported sample list in $DB_p\$DB_name\_samples\.list after sending job\(s\)\.\n";
	return ($DB_p, $DB_c, $DB_name, $im_list, $list_num, $ns);
} #import gvcf files into database
sub CombineGVCFs_4 {
	#print "debug: CombineGVCFs_4\n";
	my $dir = getcwd;
	my $proj = shift; my $nlc = shift; my $ns = shift; my $ran = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; $ow = shift; $pre = shift; my $dblist = shift; my $local = shift;
	my $time = scalar localtime();
	my $cnt_int = 0;
	my $path = "02\-get_gvcf_$ran";
	my $j_che; my @che;
	if (-d $path){
		opendir DIR, $path || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path. Using path: $path_o\n";
			$path = $path_o;
			opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
			@che = readdir DIR; chomp(@che);
			close(DIR);
			$j_che = join(" ", @che);
			if ($j_che !~ /g.vcf.gz/){
				unless ($exc){
					print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
				}
				else {
					print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
					exit;
				}
			}
		}
	}
	else {
		$path = $path_o;
		opendir DIR, $path_o || die BOLD "\[$time\]\: ERROR\: Cannot open \/$path_o: $!", RESET, "\n";
		@che = readdir DIR; chomp(@che);
		close(DIR);
		$j_che = join(" ", @che);
		if ($j_che !~ /g.vcf.gz/){
			unless ($exc){
				print "\[$time\]\: WARNING\: Cannot find any g.vcf.gz file in $path.\n";
			}
			else {
				print BOLD "\[$time\]\: ERROR\: Cannot find any g.vcf.gz file in $path.", RESET, "\n";
				exit;
			}
		}
	}
	if (-d $path){
		unless ($exc){
			print "\[$time\]\: The list file will be stored at $path.\n";
		}
		else {
			print "\[$time\]\: WARNING: Make sure you have executed step 2.\n";
		}
	}
	if (-e $dblist){
		print "\[$time\]\: gvcf list is loaded.\n";
		open(DBLIST, "<$dblist") || die BOLD "\[$time\]\: ERROR\: Cannot open $dblist: $!", RESET, "\n";
		my @db_line = <DBLIST>;
		chomp(@db_line);
		my @valis;
		foreach (@db_line){
			@valis = split("\t", $_);
			chomp(@valis);
			if (-e $valis[1]){
				push(@lists, $valis[1]);
			}
			else{
				print RED "\[$time\]\: WARNING\: Cannot find the path of $valis[0]. Skip importing this sample.\n", RESET;
			}
		}
		close(DBLIST);
	}	
	else{
		print "\[$time\]\: gvcf list is not defined.\n\[$time\]\: Using files in $path as a gvcf list.\n";
		@lists = <$path\/*.vcf.gz>;
		chomp(@lists);	
		if ($lists[0] !~ /vcf\.gz/){
			@lists = <$path_o\/*.vcf.gz>;
			if ($lists[0] !~ /vcf\.gz/ && $exc){
				if ($lists[0] !~ /vcf\.gz/){
					print BOLD "\[$time\]\: ERROR\: Cannot find \*\.vcf\.gz file!", RESET, "\n";
					return (2, $ns);
				}
			}
			if ($lists[0] !~ /vcf\.gz/ && $exc ne "-cj_exc "){
				print "\[$time\]\: WARNING\: Cannot find \*\.vcf\.gz file\(s\)\! Make sure you have \*\.vcf\.gz file before execute qsub jobs.\n";
				@lists = <$path_o\/*.fastq.gz>;
				for (my $z=0; $z<=$#lists; $z++){
					if ($lists[$z] =~ /_R2/){
						$lists[$z] = undef;
						next();
					}
					if ($lists[$z] =~ /_R1/){
						$lists[$z] =~ s/_R1//gi;
					}
					$lists[$z] =~ s/fastq\.gz/g\.vcf\.gz/;
					$lists[$z] =~ s/$path_o\///;
					$lists[$z] = "$path\/$lists[$z]";
				}
			}
		}
	}
	my @temp; my @temp_l; my @fi_lists; my $fi_tb; my $out_list;
	open(LIST, ">$path\/gvcf_list\.list") || die BOLD "\[$time\]\: ERROR\: Cannot write $path\/gvcf_list\.list: $!", RESET, "\n";
	$time = scalar localtime();	
	#print "debug: \n", join("\n", @lists), "\n";
	foreach my $list (@lists){
		$seen = 0;
		$time = scalar localtime();
		@temp = split (/\//, $list);
		@temp_l = $temp[-1];
		$temp[0] =~ s/\.gz$//;
		$temp[0] =~ s/\.g\.vcf$//;
		if ($list !~ /[a-z]/i){
			next();
		}
		if (-e "$folder\/c_vcf\.list" && ($nlc eq "0" || $ow eq "0")){
			open(FILE, "<$folder\/c_vcf\.list") || die BOLD "\[$time\]\: ERROR\: Cannot open $folder\/c_vcf\.list: $!", RESET, "\n";
			my $cc = <FILE>;
			chomp($cc);
			$time = scalar localtime();
			my @ccs = split(/\t/, $cc);
			my $exist = 0;
			foreach $x (0..$#ccs){
				if ($ccs[$x] =~ /\b$temp_l[0]\b/g){
					print "\[$time\]\: WARNING: Sample $temp_l[0] already in the gvcf file.\n";
					print "\[$time\]\: WARNING: Skip combining $temp_l[0] into gvcf file.\n";
					$exist = 1;
					last;
				}
			}
			if ($exist == 1){
				next;
			}
		}
		$fi_tb = 0;
		foreach my $x (0..$#fi_lists){
			my @spl = split(/\//, $fi_lists[$x]);
			@spl = $spl[-1];
			$spl[0] =~ s/\.gz$//;
			$spl[0] =~ s/\.g\.vcf$//;
			if ($spl[0] eq $temp_l[0]){
				print "\[$time\]\: WARNING: Sample $temp_l[0] has been detected before.\n";
				print "\[$time\]\: WARNING: Skip combining the duplicated $temp_l[0] into gvcf file.\n";
				$fi_tb = 1;
				last;
			}
		}
		if ($if_tb == 1){
			next;
		}
		if ($list !~ /^\./ && $dblist eq "0"){
			$list = "$dir\/$list";
		}
		print LIST "$temp_l[0]\t$list\n";
		$out_list = "$out_list"."$temp_l[0]\t";
		push(@fi_lists, $list);
	}
	close(LIST);
	if ($fi_lists[0] !~ /\w/){
		print RED "\[$time\]\: WARNING\: No sample to combine!", RESET, "\n";
		return (1, $ns);
	}
	@lists = (); @temp = (); @temp_l = ();	
	$time = scalar localtime();
	if (-e "$folder"){
		print "\[$time\]\: The combined g\.vcf file is stored at $folder\n";
	}
	else {
		system ("mkdir $folder");
		print "\[$time\]\: Make a folder\: $folder\n\[$time\]\: The combined g\.vcf file will be stored there.\n";
	}
	my @chrs = &check_chrs($path_o, $ref, $pre);
	if ($ns eq "1"){
		@chrs = ("all");
	}
	unless ($local){
		open (BASH, ">my_bash_03_$ran\.sh") || die BOLD "\[$time\]\: ERROR\: Cannot write my_bash_03_$ran\.sh: $!", RESET, "\n";
	}
	my $L; my $IP;
	foreach my $chr (@chrs){
		my $out;
		if (-e "$folder\/combined_vcf_$chr\.g\.vcf\.gz"){
			if ($nlc eq "1"){
				next();
			}
			elsif ($ow eq "0"){
				print "\[$time\]\: $folder\/combined_vcf_$chr\.g\.vcf\.gz exist. Updating $folder\/combined_vcf_$chr\.g\.vcf\.gz file.\n";
				$out = "mv $folder\/combined_vcf_$chr\.g\.vcf\.gz $folder\/combined_vcf_$chr\.t\.g\.vcf\.gz\nmv $folder\/combined_vcf_$chr\.g\.vcf\.gz\.tbi $folder\/combined_vcf_$chr\.t\.g\.vcf\.gz\.tbi\\n";
			}
			elsif ($ow eq "1"){
				print "\[$time\]\: $folder\/combined_vcf_$chr\.g\.vcf\.gz exist. Over-write the file.\n";
				$out = "rm $folder\/combined_vcf_$chr\.g\.vcf\.gz\\n";
			}
		}
		if ($ns eq "0"){
			$L = "\-L $chr ";
		}
		else {
			$L = "";
		}
		$out .= "\$gatk2 CombineGVCFs \-R $ref $IP$L";
		foreach (@fi_lists){
			chomp($_);
			$out .= "\-V $_ ";
		}
		if ($out =~ /\.t\.g\.vcf\.gz/){
			$out .= "\-V $folder\/combined_vcf_$chr\.t\.g\.vcf\.gz ";
		}
		$out .= "\-O $folder\/combined_vcf_$chr\.g\.vcf\.gz\\n";
		my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_mem 48 -cj_qname gatk_03_$chr -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return ($out_list, $ns);
} #combine individual sample gvcf into merged gvcf
sub GenotypeGVCFs_4 {
	my $proj = shift; my $wes = shift; my $interval = shift; my $ip = shift; my $ns = shift; my $prese = shift; my $pxlsn = shift; my $xlsn = shift; my $svsn1 = shift; my $svsn2 = shift; my $ran = shift; my $c = shift; 
	my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $DB_path = shift; my $as = shift; my $local = shift; my $gpu = shift;
	my $time = scalar localtime();
	my @chrs; my $v; my $folder_out; my $gz_back;
	#it's not working for gpu right now, because the genotypgvcf in parabricks is not supported yet.
	$gpu = 0;
	if ($gpu == 1){
	    if ($pxlsn ne "0" || $xlsn ne "0"){
            print RED "\[$time\]\: WARNING\: -prx or -pox cannot be used when -gpu is used.", RESET, "\n";
            $pxlsn = "0";
            $xlsn = "0";
	    }
	}
	#print "debug: -ns $ns\n";
	if ($as =~ /[a-z]/i){
		$as = "$as ";
		if ($gpu == 1){
		    $as = 1;
		}
	}
	my $dict;
	$dict = $ref;
	$dict =~ s/fa$|fas$|fasta$/dict/i;
	@chrs = &check_chrs($path_o, $ref, $pre);
	if ($wes || $ns eq "1"){
		@chrs = "all";
	}
	unless ($local){
		open (BASH, ">my_bash_04_gz_$ran\.sh") || die BOLD "Cannot write my_bash_04_gz_$ran\.sh: $!", RESET, "\n";
	}
	foreach (@chrs){
		if ($c eq "0"){
			if (-d "$DB_path\_$_"){}
			else {
				if ($exc){
					$time = scalar localtime();
					print BOLD "\[$time\]\: ERROR\: Cannot find $DB_path\_$_.", RESET, "\n";
					return 2;
				}
			}
		}
		elsif ($c eq "1"){
			$time = scalar localtime();
			if (-e "$folder\/combined_vcf_$_\.g\.vcf\.gz"){}
			elsif (-e "$path_o\/combined_vcf_$_\.g\.vcf\.gz"){
				$folder_out = $folder;
				$folder = $path_o;
			}
			elsif (-e "$folder\/combined_vcf_$_\.g\.vcf"){
				$gz_back = &bgzip($proj, "$folder\/combined_vcf_$_\.g\.vcf", $_, $exc, $ran, $local);
				if ($gz_back eq "2"){
					return 2;
				}
				else {
					print BASH "$gz_back\n";
				}
			}
			elsif (-e "$path_o\/combined_vcf_$_\.g\.vcf"){
				$folder_out = $folder;
				$folder = $path_o;	
				$gz_back = &bgzip($proj, "$folder\/combined_vcf_$_\.g\.vcf", $_, $exc, $ran, $local);
				if ($gz_back eq "2"){
					return 2;
				}
				else {
					print BASH "$gz_back\n";
				}
			}
			elsif ($exc){
				print BOLD "\[$time\]\: ERROR: Cannot find combined_vcf_$_\.g\.vcf\.gz file.", RESET, "\n";
				return 2;
			}	
		}
	}
	unless ($local){
		close(BASH);
		if ($exc){
			&status($ran);
		}
	}
	foreach (@chrs){
		if (-e "$folder\/combined_vcf_$_\.g\.vcf"){
			if (-e "$folder\/combined_vcf_$_\.g\.vcf\.gz" && -e "$folder\/combined_vcf_$_\.g\.vcf\.gz\.tbi"){
				system("rm $folder\/combined_vcf_$_\.g\.vcf");
				if (-e "$folder\/combined_vcf_$_\.g\.vcf\.idx"){
					system("rm $folder\/combined_vcf_$_\.g\.vcf\.idx");
				}
			}
		}
	}	
	$time = scalar localtime();
	if ($folder_out){
		if (-e $folder_out){
			print "\[$time\]\: The gvcf file(s) will be stored at $folder_out\n";
		}
		else {
			print "\[$time\]\: Make a folder\: $folder_out.\n\[$time\]\: The star_filtered vcf file will be stored there.\n";
			my $err = `mkdir $folder_out 2>&1`;
			if ($err =~ /[a-z]/i){
				print BOLD "\[$time\]\: ERROR: Cannot make folder \"$folder_out\".", RESET, "\n";
				return 2;
			}		
		}
	}
	elsif (-d $folder){
		print "\[$time\]\: The gvcf file(s) will be stored at $folder\n";
	}
	else{
		print "\[$time\]\: Make a folder\: $folder.\n\[$time\]\: The star_filtered vcf file will be stored there.\n";
		my $err = `mkdir $folder 2>&1`;
		if ($err =~ /[a-z]/i){
			print BOLD "\[$time\]\: ERROR: Cannot make folder \"$folder\".", RESET, "\n";
			return 2;
		}
	}
	unless ($local){
		open (BASH, ">my_bash_04_$ran\.sh") || die BOLD "Cannot write my_bash_04_$ran\.sh: $!", RESET, "\n";
	}
	my $pre_t; my $L; my $pre_select; my $IP;
	foreach my $chr (@chrs){
		my $out;
		$time = scalar localtime();
		if ($c eq "0"){
			$pre_t = "TILEDB_DISABLE_FILE_LOCKING\=1 ";
			$v = "gendb\:\/\/$DB_path\_$chr";
			if ($ns eq "0"){
				$L = "\-L $chr";
			}
			elsif ($wes){
				$L = "\-L $interval";
				$IP = "\-ip $ip ";
			}
			else {
				$L = "";
			}
			unless ($exc){
				unless (-d "$DB_path\_$chr"){
					print RED "\[$time\]\: WARNING\: Cannot find database $DB_path\_$chr.", RESET, "\n";
					print "\[$time\]\: Skip processing $DB_path\_$chr.\n";
				}
			}
		}
		elsif ($c eq "1"){
			$v = "$folder\/combined_vcf_$chr\.g\.vcf\.gz";
			$L = "";
			unless (-e "$folder\/combined_vcf_$chr\.g\.vcf\.gz\.tbi"){
			    if ($exc){
				    print RED "\[$time\]\: WARNING\: Cannot find index file of combined_vcf_$chr\.g\.vcf\.gz\. This file might be truncated.", RESET, "\n";
				    print "\[$time\]\: Skip processing $folder\/combined_vcf_$chr\.g\.vcf\.gz\.\n";
				    next();
				}
			}
		}
		if ($pxlsn ne "0"){
			if ($folder_out){
				if (-e "$folder_out\/pre_select_vcf_$chr\.g\.vcf\.gz"){}
				else {
					print "\[$time\]\: Pre-selected list is loaded.\n";
					$pre_select = "$pre_t\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=2\" SelectVariants \-R $ref \-V $v \-O $folder_out\/pre_select_vcf_$chr\.g\.vcf\.gz \-\-sequence\-dictionary $dict $svsn1 $pxlsn $IP$L\\n";
					$v = "$folder_out\/pre_select_vcf_$chr\.g\.vcf\.gz";
				}
			}
			elsif (-e "$folder\/pre_select_vcf_$chr\.g\.vcf\.gz"){}
			else {
				print "\[$time\]\: Pre-selected list is loaded.\n";
				$pre_select = "$pre_t\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=2\" SelectVariants \-R $ref \-V $v \-O $folder\/pre_select_vcf_$chr\.g\.vcf\.gz \-\-sequence\-dictionary $dict $svsn1 $pxlsn $IP$L\\n";
				$v = "$folder\/pre_select_vcf_$chr\.g\.vcf\.gz";
			}
		}
		if ($prese ne "0"){
			if (-e "$prese\/pre_select_vcf_$chr\.g\.vcf\.gz"){
				print "\[$time\]\: Pre-selected files are loaded.\n";
				$v = "$prese\/pre_select_vcf_$chr\.g\.vcf\.gz";
			}
			else {
				print RED "\[$time\]\: WARNING\: $prese\/pre_select_vcf_$chr\.g\.vcf\.gz is not found.", RESET, "\n";
			}
		}
		if ($folder_out){
			if (-e "$folder_out\/vcf_raw_$chr\.vcf\.gz"){
				if ($xlsn eq "0"){
					print "\[$time\]\: File $folder_out\/vcf_raw_$chr\.vcf\.gz exists. Skip the step.\n";
					next();
				}			
			}
		}
		elsif (-e "$folder\/vcf_raw_$chr\.vcf\.gz"){
			if ($xlsn eq "0"){
				print "\[$time\]\: File $folder\/vcf_raw_$chr\.vcf\.gz exists. Skip the step.\n";
				next();
			}
		}
		if ($pxlsn ne "0" && $xlsn eq "0"){
			if ($folder_out){
				if (-e "$folder_out\/pre_select_vcf_$chr\.g\.vcf\.gz"){}
				else {
					$out .= "$pre_select";
				}					
			}
			elsif (-e "$folder\/pre_select_vcf_$chr\.g\.vcf\.gz"){}
			else {
				$out .= "$pre_select";
			}
		}
		if ($folder_out){
			if (-e "$folder_out\/vcf_raw_$chr\.vcf\.gz"){}
			else {
			    if ($gpu == 1){
			        $out .= "pbrun genotypegvcf \--ref $ref \--in-gvcf $v \--out-vcf $folder_out\/vcf_raw_$chr\.vcf --num-threads 32\\n";
			        $out .= "bgzip $folder_out\/vcf_raw_$chr\.vcf\\n";
			        if ($as == 1){
			            $out .= "awk -v OFS\=\"\\t\" \'\{print \$1\, 0\, \$2\}\' $ref.fai \> $ref.bed\\n";
			            $out .= "bcftools mpileup --threads 32 -f $ref --regions-file $ref.bed -Ou $folder_out\/vcf_raw_$chr\.vcf \| bcftools call --threads 32 --ploidy 2 -m -o $folder_out\/vcf_raw_$chr\.vcf\.gz\\n";
			            $out .= "bcftools merge -R $ref.bed --missing-to-ref -f PASS -Oz -o $folder_out\/vcf_raw_$chr\.allsites.vcf\.gz $folder_out\/vcf_raw_$chr\.vcf\.gz\\n";
			            $out .= "mv $folder_out\/vcf_raw_$chr\.allsites.vcf\.gz $folder_out\/vcf_raw_$chr\.vcf\.gz\\n";
			        }
			        $out .= "tabix $folder_out\/vcf_raw_$chr\.vcf\.gz\\n";
			    }
			    else {
				    $out .= "$pre_t\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=4\" GenotypeGVCFs \-R $ref \-V $v \-O $folder_out\/vcf_raw_$chr\.vcf\.gz \-\-sequence\-dictionary $dict $as$IP$L\\n";
			    }
			}				
		}
		elsif (-e "$folder\/vcf_raw_$chr\.vcf\.gz"){}
		else {
		    if ($gpu == 1){
		        $out .= "pbrun genotypegvcf \--ref $ref \--in-gvcf $v \--out-vcf $folder\/vcf_raw_$chr\.vcf --num-threads 32\\n";
		        $out .= "bgzip $folder\/vcf_raw_$chr\.vcf\\n";
			    if ($as == 1){
			        $out .= "awk -v OFS\=\"\\t\" \'\{print \$1\, 0\, \$2\}\' $ref.fai \> $ref.bed\\n";
			        $out .= "bcftools mpileup --threads 32 -f $ref --regions-file $ref.bed -Ou $folder\/vcf_raw_$chr\.vcf \| bcftools call --threads 32 --ploidy 2 -m -o $folder\/vcf_raw_$chr\.vcf\.gz\\n";
			        $out .= "bcftools merge -R $ref.bed --missing-to-ref -f PASS -Oz -o $folder\/vcf_raw_$chr\.allsites.vcf\.gz $folder\/vcf_raw_$chr\.vcf\.gz\\n";
			        $out .= "mv $folder\/vcf_raw_$chr\.allsites.vcf\.gz $folder\/vcf_raw_$chr\.vcf\.gz\\n";
			    }
		        $out .= "tabix $folder\/vcf_raw_$chr\.vcf\.gz\\n";
		    }
		    else {
			    $out .= "$pre_t\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=4\" GenotypeGVCFs \-R $ref \-V $v \-O $folder\/vcf_raw_$chr\.vcf\.gz \-\-sequence\-dictionary $dict $as$IP$L\\n";
		    }
		}
		my $nas;
		if ($xlsn ne "0" && $pxlsn eq "0"){
			if ($as !~ /[a-z]/i){
				$nas = "\-\-exclude\-non\-variants ";
			}
			else {
				$nas = ();
			}
			if ($folder_out){
				print "\[$time\]\: Post-selected list is loaded.\n";
		        $out .= "\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=4\" SelectVariants \-R $ref \-V $folder_out\/vcf_raw_$chr\.vcf\.gz \-O $folder_out\/select_vcf_raw_$chr\.vcf\.gz \-\-sequence\-dictionary $dict $nas$svsn2 $xlsn\\n";
			}
			else {
				print "\[$time\]\: Post-selected list is loaded.\n";
				$out .= "\$gatk2 \-\-java\-options \"\-Xmx60g \-XX\:ParallelGCThreads=4\" SelectVariants \-R $ref \-V $folder\/vcf_raw_$chr\.vcf\.gz \-O $folder\/select_vcf_raw_$chr\.vcf\.gz \-\-sequence\-dictionary $dict $nas$svsn2 $xlsn\\n";
			}
		}
		my $return;
		if ($gpu == 1){
		    $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_docker nvcr-clara-parabricks-4401 -cj_time 168\:0\:0 -cj_gpu 4 -cj_qname gatk_04_$chr -cj_sn $ran -cj_qout . $out");
		}
		else {
		    $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_ppn 4 -cj_mem 60 -cj_qname gatk_04_$chr -cj_sn $ran -cj_qout . $out");
		}
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	return 1;
} #compile gvcfs
sub VQSR {
	my $as_name = shift; my $r_env = shift; my $vep_env = shift; my $proj = shift; my $ran = shift; my $folder = shift; my $exc = shift; my $path_o = shift; my $ref = shift; my $pre = shift; my $local = shift; my @reses = @{$_[-2]}; my @resis = @{$_[-1]};
	my $time = scalar localtime(); my $out;
	unless ($folder){
		$folder = $path_o;
	}
	my $file = "$folder\/select_vcf_raw_all.vcf.gz";
	unless (-e $file){
		$file = "$path_o\/select_vcf_raw_all.vcf.gz";
	}
	my $sel = "select_";
	unless (-e $file){
		$file = "$folder\/vcf_raw_all.vcf.gz";
		unless (-e $file){
			$file = "$path_o\/vcf_raw_all.vcf.gz";
		}
		$sel = "";
	}
	unless ($file){
		$file = "$folder\/vcf_raw_all.vcf.gz";
	}
	unless (-e "$folder\/$sel\Qrecalibrate_SNP\E.recal"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" VariantRecalibrator -R $ref -V $file @reses -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -O $folder\/$sel\Qrecalibrate_SNP\E.recal --tranches-file $folder\/$sel\Qrecalibrate_SNP\E.tranches --rscript-file $folder\/$sel\Qrecalibrate_SNP_plots\E.R\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated_SNP\E.vcf"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" ApplyVQSR -R $ref -V $file --mode SNP --truth-sensitivity-filter-level 99.0 --recal-file $folder\/$sel\Qrecalibrate_SNP\E.recal --tranches-file $folder\/$sel\Qrecalibrate_SNP\E.tranches -O $folder\/$sel\Qrecalibrated_SNP\E.vcf\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrate_INDEL\E.recal"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" VariantRecalibrator -R $ref -V $folder\/$sel\Qrecalibrated_SNP\E.vcf @resis -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --mode INDEL --max-gaussians 6 -O $folder\/$sel\Qrecalibrate_INDEL\E.recal --tranches-file $folder\/$sel\Qrecalibrate_INDEL\E.tranches --rscript-file $folder\/$sel\Qrecalibrate_INDEL_plots\E.R\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated_INDEL\E.vcf"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" ApplyVQSR -R $ref -V $folder\/$sel\Qrecalibrated_SNP\E.vcf --mode INDEL --truth-sensitivity-filter-level 99.0 --recal-file $folder\/$ser\Qrecalibrate_INDEL\E.recal --tranches-file $folder\/$sel\Qrecalibrate_INDEL\E.tranches -O $folder\/$sel\Qrecalibrated_INDEL\E.vcf\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated\E.filtered.vcf"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" SelectVariants -R $ref -V $folder\/$sel\Qrecalibrated_INDEL\E.vcf -O $folder\/$sel\Qrecalibrated\E.filtered.vcf --exclude-non-variants --remove-unused-alternates\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated_SNP\E.filtered.vcf"){
		$out .= "\$gatk2 --java-options \"\-Xmx60g \-XX\:ParallelGCThreads=1\" SelectVariants -R $ref -V $folder\/$sel\Qrecalibrated_SNP\E.vcf -O $folder\/$sel\Qrecalibrated_SNP\E.filtered.vcf --exclude-non-variants --remove-unused-alternates\\n";
	}
	#unless (-e "$folder\/$sel\Qrecalibrated\E.filtered.vcf"){
	#	$out .= "bcftools concat -Ov -o $folder\/$sel\Qrecalibrated\E.filtered.vcf $folder\/$sel\Qrecalibrated_SNP\E.filtered.vcf $folder\/$sel\Qrecalibrated_INDEL\E.filtered.vcf\\n";
	#}
#vep annotation
	#check if VEP exist	
	my $ck_vep = `export PATH\=$vep_env\:\$PATH; vep`;
	if ($ck_vep !~ /ENSEMBL/){
			print BOLD "\[$time\]\: ERROR\: Cannot find the VEP program.", RESET, "\n";
			return 2;
	}
	unless (-d "vep_cache"){
		system("mkdir vep_cache");
		unless (-d "vep_cache"){
			print BOLD "\[$time\]\: ERROR\: Cannot create \"vep_cache\" for VEP annotation.", RESET, "\n";
			return 2;
		}
		system("mkdir vep_cache\/Plugins");
		unless (-d "vep_cache\/Plugins"){
			print BOLD "\[$time\]\: ERROR\: Cannot create \"vep_cache\/Plugins\" for VEP annotation.", RESET, "\n";
			return 2;
		}
	}
	unless ($local){
		open (BASH, ">my_bash_05-1_$ran\.sh") || die BOLD "Cannot write my_bash_05-1_$ran\.sh: $!", RESET, "\n";
	}
	my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_env $r_env -cj_mem 60 -cj_qname gatk_05-1_VQSR -cj_sn $ran -cj_qout . $out");
	unless ($local){
		print BASH "$return\n";
		close(BASH);
	}
	if ($exc){
		$time = scalar localtime();
		unless ($local){
			&status($ran);
		}
	}
	#$out .= "export PATH\=$vep_env\:\$PATH\\n";
	my $snp_only = "_SNP";
	$out = "";
	if (-e "$folder\/$sel\Qrecalibrated\E.filtered.vcf"){
		$snp_only = "";
		#$out .= "rm $folder\/\*INDEL\*\\n";
	}
	#rename ref. contig name as numeric
	my $num_ref = $ref;
	$num_ref =~ s/fasta$|fa$|fas$/num.fasta/;
	unless (-e $num_ref){
		$out .= "cp $ref $num_ref\\n";
		$out .= "sed -i \'s\/chrM\/MT\/g\' $num_ref\\n";
		$out .= "sed -i \'s\/chr\/\/g\' $num_ref\\n";
	}
	#rename contig name in vcf file
	my $num_vcf = "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vcf";
	$num_vcf =~ s/vcf$/num.vcf/;
	unless (-e "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.num.vcf"){
	my @chr_names = &chr_name("$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vcf");
		open(REN, ">$folder\/$sel\Qrecalibrated\E$snp_only.rename.list") || die "Cannot write $folder\/$sel\Qrecalibrated\E$snp_only.rename.list: $!\n";
		foreach my $x (0..$#chr_names){
			my $new = $chr_names[$x];
			$new =~ s/chrM/MT/g;
			$new =~ s/chr//ig;
			print REN "$chr_names[$x] $new\n";
		}
		close(REN);
		$out .= "bcftools annotate -Ov -o $num_vcf --rename-chrs $folder\/$sel\Qrecalibrated\E$snp_only.rename.list $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vcf\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.vcf.gz"){
		$out .= "vep --cache vep_cache --dir_cache vep_cache --fasta $num_ref --fork 8 --assembly\=$as_name --offline --vcf --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i $num_vcf -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.vcf --force_overwrite\\n";
		#$out .= "vep --fasta $ref --fork 8 --assembly\=$as_name --offline --vcf --everything --terms SO --pick --coding_only --transcript_version -i $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vcf.gz -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.vcf --force_overwrite\\n";
		$out .= "bgzip $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.vcf\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.tsv"){
		$out .= "vep --cache vep_cache --dir_cache vep_cache --fasta $num_ref --fork 8 --assembly\=$as_name --offline --tab --plugin Downstream --everything --terms SO --pick --coding_only --transcript_version -i $num_vcf -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.tsv --force_overwrite\\n";
		#$out .= "vep --fasta $ref --fork 8 --assembly\=$as_name --offline --tab --everything --terms SO --pick --coding_only --transcript_version -i $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vcf -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.tsv --force_overwrite\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.focusing.vcf.gz"){
		$out .= "filter_vep --format vcf -i $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.vcf.gz -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.focusing.vcf --filter \"\(MAX_AF \< 0.001 or not MAX_AF\) and \(\(IMPACT is HIGH\) or \(IMPACT is MODERATE and \(SIFT match deleterious or PolyPhen match damaging\)\)\)\" --force_overwrite\\n";
		$out .= "bgzip $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.focusing.vcf\\n";
	}
	unless (-e "$folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.focusing.tsv"){
		$out .= "filter_vep --format tab -i $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.tsv -o $folder\/$sel\Qrecalibrated\E$snp_only.filtered.vep.focusing.tsv --filter \"\(MAX_AF \< 0.001 or not MAX_AF\) and \(\(IMPACT is HIGH\) or \(IMPACT is MODERATE and \(SIFT match deleterious or PolyPhen match damaging\)\)\)\" --force_overwrite\\n";
	}
	unless ($local){
		open(BASH, ">my_bash_05-2_$ran\.sh") || die BOLD "Cannot write my_bash_05-2_$ran\.sh: $!", RESET, "\n";
	}
	my $return = &pbs_setting("$proj$exc$local\-cj_quiet -cj_env $vep_env -cj_mem 60 -cj_qname gatk_05-2_VQSR -cj_sn $ran -cj_qout . $out");
	unless ($local){
		print BASH "$return\n";
		close(BASH);
	}
	return 1;
} #WES filtering






