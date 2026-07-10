#!/usr/local/bin/perl
use Cwd qw(getcwd);
use FindBin;
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
else {
	die "Cannot find required subroutine file: qsub_subroutine.pl\n";
}
chomp(@ARGV);
print BOLD "Usage: perl merge_cnv_vcf.pl -p PATH -r REFERENCE_FILE  [-ped][-l INTERVALS][-o OUTPUT][-sn SN][-ow][-exc]\n", RESET;
print "Input command line:\n";
print "perl merge_cnv_vcf.pl @ARGV\n\n";

my $exc; my $sn; my $path; my $ow = 0; my $l_file; my $ref; my $output; my $ped;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "-exc" || $ARGV[$i] eq "--execute"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "-sn" || $ARGV[$i] eq "--serial-number"){
		$sn = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-l" || $ARGV[$i] eq "--intervals"){
		$l_file = $ARGV[$i+1];
		unless (-e $l_file){
            die "Cannot find the interval file: $l_file\n";
		}
	}
	if ($ARGV[$i] eq "-r" || $ARGV[$i] eq "--reference"){
		$ref = $ARGV[$i+1];
		unless (-e $ref){
            die "Cannot find the reference file: $ref\n";
		}
	}
	if ($ARGV[$i] eq "-ped" || $ARGV[$i] eq "--pedigree"){
		$ped = $ARGV[$i+1];
		unless (-e $ped){
            die "Cannot find the pedigree file: $ref\n";
		}
	}
	if ($ARGV[$i] eq "-p"){
		$path = $ARGV[$i+1];
		$path =~ s/\/$//;
		unless (-d $path){
            die "Cannot find the path: $path\n";
		}
	}
	if ($ARGV[$i] eq "-ow" || $ARGV[$i] eq "--overwrite"){
		$ow = 1;
	}
	if ($ARGV[$i] eq "-o" || $ARGV[$i] eq "--output"){
		$output = $ARGV[$i+1];
	}
}
unless ($path){
    die "-p is required.\n";
}
unless ($ref){
    die "-r is required.\n";
}

#merge all segments vcfs
##step 1: generate a pedigree file
my $new_sample = 0; my $out;
my @vcfs = <$path\/*segments*vcf.gz>;
@vcfs = grep(!/^all/, @vcfs);
unless (@vcfs){
    die "Cannot find any CNV vcfs.\n";
}
unless ($ped){
	if (-e "$path\/sample_pedigree.ped" && $ow == 0){ #check if all samples are in the list
    	open (PED, "<$path\/sample_pedigree.ped") || die BOLD "Cannot open $path\/sample_pedigree.ped: $!", RESET, "\n";
    	my @p_content = <PED>;
    	chomp(@p_content);
    	close(PED);
    	open (PED, ">>$path\/sample_pedigree.ped") || die BOLD "Cannot write $path\/sample_pedigree.ped: $!", RESET, "\n";
    	foreach my $i (0..$#vcfs){
        	my $sample_name = &cnv_sample_name($vcfs[$i]);
			if ($sample_name =~ /^all/i){
				next;
			}
        	if ($p_content !~ /\b$sample_name\b/){
            	$new_sample = 1;
            	print PED "NO_FAMILY_ID\t$sample_name\tNO_PARENTS\tNO_PARENTS\tNO_SEX\tNO_PHENOTYPE\n";
        	}
    	}
    	close(PED);
	}
	else {
    	open (PED, ">$path\/sample_pedigree.ped") || die BOLD "Cannot write $path\/sample_pedigree.ped: $!", RESET, "\n";
    	foreach my $i (0..$#vcfs){
        	my $sample_name = &cnv_sample_name($vcfs[$i]);
        	print PED "NO_FAMILY_ID\t$sample_name\tNO_PARENTS\tNO_PARENTS\tNO_SEX\tNO_PHENOTYPE\n";
    	}
    	close(PED);
    	$new_sample = 0;
	}
	$ped = "$path\/sample_pedigree.ped";
}
##step 2: do merging segments vcfs
my $out;
if (-e $l_file){
    $l_file = "--model-call-intervals $l_file ";
}
else {
    $l_file = "";
}
unless ($output){
    $output = "$path\/all.PostCNVCalls.segments.vcf.gz";
}
foreach my $i (0..$#vcfs){
	unless (-e "$vcfs[$i].tbi"){
		$out .= "tabix $vcfs[$i]\\n";
	}
}
unless (-e $output && $new_sample == 0 && $ow == 0){
    $out .= "\$gatk2 JointGermlineCNVSegmentation -R $ref $l_file\-ped $ped ";
    foreach my $i (0..$#vcfs){
        $out .= "-V $vcfs[$i] ";
    }
    $out .= "-O $path\/all.PostCNVCalls.segments.vcf.gz\\n";
}
my @return;
print "debug: $exc -cj_mem 60 -cj_qname merge_cnv_vcf -cj_sn $ran -cj_qout . $out\n";
@return = &pbs_setting("$exc\-cj_mem 60 -cj_qname merge_cnv_vcf -cj_qout . $out");
print "@return\n";
sub cnv_sample_name {
	my $sample = shift;
	my @tmps = split(/\//, $sample);
	$tmps[-1] =~ s/\.read\.counts\.hdf5$|\.PostCNVCalls\.copy_ratios\.tsv$|\.PostCNVCalls\.segments\.vcf\.gz$//;
	return $tmps[-1];
} #get the sample name from hdf5 files
