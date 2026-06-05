#!/usr/bin/perl

use Term::ANSIColor qw(:constants);
#use Time::HiRes qw(time);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}
my $conda = "R-4.4";

print "This script is writtern by Ben Chien. Sep. 2025\n";
print "Input command line:\n";
print "perl IntroPainter_v2\.pl @ARGV\n\n";

sub usage {
    print "Usage: perl IntroPainter_v2.pl -vcf PATH -list LIST_FILE [-win INT] [-step INT] [-pre PREFIX] [-p1c COLOR] [-p2c COLOR] [-ci] [-rcal] [-rplot] [-conda NAME] [-sn SN] [-m FLOAT] [-d FLOAT] [-n] [-mem INT] [-local] [-ow] [-exc]\n";
    print "-list: trios list, population_name follows by sample_names. 3 lines: 2 parents at first two lines, and the testing population in the 3rd line.\n";
    print "-win: window size (SNP number). Default: 1000\n";
    print "-step: step size (SNP number). Default: 500\n";
    print "-p1c: represent color of ancestor 1 in the figure. If you use color code, please add \'\'. Default: blue\n";
    print "-p2c: represent color of ancestor 2 in the figure. If you use color code, please add \'\'. Default: red\n";
    print "-ci: shown 95% confidence interval or not. Default: False\n";
    print "-m: missing rate of the sites. Default: 0.8\n";
    print "-d: allele differential rate of the sites between two parents. Default: 0.8\n";
    print "-rcal: re-calcutate introgression ratios only.\n";
    print "-rplot: re-plot introgression figures only.\n";
    print "-conda: conda environment name. Default: R-4.4\n";
    print "-n: threads to use. Default: 1\n";
    print "-mem: memory to use. Default: 12Gb\n";
    print "-ow: over-write the outputs.\n";
    print "-sn: serial number for the run. Required if -rcal, -rplot, or -ow is used.\n";
    print "-local: run the pipeline in the local machine. (Warning: It might be very slow)\n";
    print "-exc: send jobs to execute.\n\n";
    print "Dependencies: bcftools, conda environment for R, R packages: ggplot2, dplyr\n";
}

my $path_v; my $path; my @vcfs; my $list; my $out; my $ran; my $exc = 0; my $missing = 0.8; my $diff = 0.8;
my $mem = 12; my $thread_in = 1; my $local; my $pre; my $cj_exc; my $ci; my $rc = 0; my $rp = 0;
my $win = 1000; my $step = 500; my $ow = 0; my $sn = 0; my $p1c = "blue"; my $p2c = "red";
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-exc"){
		$exc = 1;
		$cj_exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-vcf"){
		if (-d $ARGV[$i+1]){
            $path_v = $ARGV[$i+1];
            if ($path_v =~ /\/$/){
                $path_v =~ s/\/$//;
            }
            @vcfs = <$path_v\/*.vcf.gz>;
            unless (@vcfs){
                @vcfs = <$path_v\/*.vcf>;
                unless(@vcfs){
                    die "Cannot find vcfs.\n";
                }
            }
		}
        elsif (-e $ARGV[$i+1]){
            @vcfs = $ARGV[$i+1];
        }
		else {
            die "-vcf: Cannot find the file\(s\).\n";
		}
	}
	if ($ARGV[$i] eq "\-list"){ #trio_list
        if (-e $ARGV[$i+1]){
            $list = $ARGV[$i+1];
        }
        else {
			die "-list: Cannot find the file.\n";
        }
        if ($list =~ /\//){
			my @tmps = split(/\//, $list);
			pop(@tmps);
			$path = join("\/", @tmps);
        }
		else {
			$path = ".";
		}
		$path = &check_path($path);
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]";
        $sn = 1;
	}
	if ($ARGV[$i] eq "\-n"){ #threads
        $thread_in = "$ARGV[$i+1]";
        if ($thread_in =~ /[^0-9]/){
        	die "-n: parameter should be an integer number.\n";
        }
	}
	if ($ARGV[$i] eq "\-mem"){
		$mem = $ARGV[$i+1];
		unless ($mem !~ /[^0-9]/){
            die "-mem: should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-win"){
		$win = $ARGV[$i+1];
		unless ($win !~ /[^0-9]/){
            die "-win: should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-step"){
		$step = $ARGV[$i+1];
		unless ($step !~ /[^0-9]/){
            die "-step: should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-m"){
		$missing = $ARGV[$i+1];
		unless ($step !~ /[^0-9\.]/ && $missing >= 0 && $missing <= 1){
            die "-m: should be a number between 0\~1.\n";
		}
	}
	if ($ARGV[$i] eq "\-d"){
		$diff = $ARGV[$i+1];
		unless ($diff !~ /[^0-9\.]/ && $diff >= 0 && $diff <= 1){
            die "-d: should be a number between 0\~1.\n";
		}
	}
	if ($ARGV[$i] eq "\-local"){
		$local = "-cj_local ";
	}
	if ($ARGV[$i] eq "\-pre"){
		$pre = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-ow"){
		$ow = 1;
	}
	if ($ARGV[$i] eq "\-ci"){
		$ci = "-ci ";
	}
	if ($ARGV[$i] eq "\-rcal"){
		$rc = 1;
	}
	if ($ARGV[$i] eq "\-rplot"){
		$rp = 1;
	}
	if ($ARGV[$i] eq "\-p1c"){
		$p1c = $ARGV[$i+1];
		if ($p1c =~ /\#/){
			$p1c = "\\$p1c";
		}
	}
	if ($ARGV[$i] eq "\-p2c"){
		$p2c = $ARGV[$i+1];
		if ($p2c =~ /\#/){
			$p2c = "\\$p2c";
		}
	}
	if ($ARGV[$i] eq "\-conda"){
		$conda = $ARGV[$i+1];
	}
}

unless (@vcfs && $list){
	die "Not enough parameter is provided.\n";
}
if ($rp == 1 || $rc == 1 || $ow == 1){
	unless($ran){
		die "-rcal, -rplot, -ow: -sn is required.\n";
	}
}

RE:
if ($ran){}
else {
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "my_bash_introgression_1_$ran.sh" && $sn != 1){
	$ran = undef;
	goto RE;
}

print "The serial number is: $ran\n";

#get genome information
my @chrs = &chr_name($vcfs[0], $pre);

#indexing the vcf if necessary
open (BASH, ">my_bash_introgression_1_$ran.sh") || print "Cannot write my_bash_introgression_1_$ran.sh: $!\n";
foreach my $i (0..$#vcfs){
	unless (-e "$vcfs[$i].tbi"){
		$out = "tabix $vcfs[$i]\\n";
		my $return = &pbs_setting("$cj_exc$local\-cj_quiet -cj_qname indexing_$i -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
}
unless ($local){
	close(BASH);
}
if ($exc == 1){
	&status($ran);
}
my @ck_vcfs;
foreach my $cnt (0..$#vcfs){
	my @ck_chrs;
	if (scalar(@vcfs) == 1){
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			push(@ck_chrs, $_);
    		}
    	}
	} 
	else {
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			push(@ck_chrs, $_);
    		}
    	}
   	}
   	if (scalar(@ck_chrs) > 1){
   		open (BASH, ">my_bash_seperate_chr_$ran.sh") || print "Cannot write my_bash_seperate_chr_$ran.sh: $!\n";
   		foreach my $i (0..$#ck_chrs){
   			my $name = $vcfs[$cnt];
   			$name =~ s/.gz$//;
   			$name =~ s/.vcf$/.$ck_chrs[$i].trio_filtered.vcf.gz/;
   			unless (-e $name && $ow == 0){
   				$out = "bcftools view -Oz -o $name -m2 -M2 -v snps -r $ck_chrs[$i] --threads $thread_in $vcfs[$cnt]\\n";
   				$out .= "tabix $name\\n";
   				my $return = &pbs_setting("$cj_exc$local\-cj_quiet -cj_ppn $thread_in -cj_qname seperate_chr_$ck_chrs[$i] -cj_sn $ran -cj_qout . $out");
   				print BASH "$return\n";
   			}
   			push(@ck_vcfs, $name);
   		}
   		close(BASH);
		unless ($local){
			close(BASH);
		}
		if ($exc == 1){
			unless ($local){
				&status($ran);
			}
		}
   	}
   	else {
		my $name = $vcfs[$cnt];
   		$name =~ s/.gz$//;
   		$name =~ s/.vcf$/.trio_filtered.vcf.gz/;
   		unless (-e $name && $ow == 0){
   			$out = "bcftools view -Oz -o $name -m2 -M2 -v snps --threads $thread_in $vcfs[$cnt]\\n";
   			$out .= "tabix $name\\n";
   			&pbs_setting("$cj_exc$local\-cj_quiet -cj_ppn $thread_in -cj_qname filtering_$cnt -cj_sn $ran -cj_qout . $out");
   		}
   		push(@ck_vcfs, $name);
   	}
}
if ($exc == 1){
	unless ($local){
		&status($ran);
	}
}
#update the vcf list
@vcfs = @ck_vcfs;

#convert vcf to trios format
my @trios_files;
open (BASH, ">my_bash_introgression_2_$ran.sh") || print "Cannot write my_bash_introgression_2_$ran.sh: $!\n";
foreach my $i (0..$#vcfs){
	my $trios_file_name = $vcfs[$i];
	$trios_file_name =~ s/\.vcf\./.$ran.trios./;
	unless (-e $trios_file_name && $ow == 0){
		$out = "perl vcf2trios_thread_v2.pl -vcf $vcfs[$i] -list $list -sn $ran -n $thread_in\\n";
		my $return = &pbs_setting("$cj_exc$local\-cj_quiet -cj_mem $mem -cj_ppn $thread_in -cj_qname convertion_$i -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	push(@trios_files, $trios_file_name);
}
unless ($local){
	close(BASH);
}
if ($exc == 1){
	unless ($local){
		&status($ran);
	}
}

#do introgression calculation
unless (-d "$path\/$ran\_output"){
	system("mkdir $path\/$ran\_output");
}
my $ck_files = `ls -l $path\/$ran\_output \| grep \"introgression_${diff}.rda\"`;
if ($thread_in > 4){
    $thread_in = 4;
}
if ($mem <= 24){
    if ($thread_in > 2){
        $thread_in = 2;
    }
}
unless ($ck_files && $ow == 0 && $rc == 0){
	open (BASH, ">my_bash_introgression_3_$ran.sh") || print "Cannot write my_bash_introgression_3_$ran.sh: $!\n";
	foreach my $i (0..$#trios_files){
		$out = "Rscript new_intro_count_2_majorAF.R -g $trios_files[$i] -t $list -gi $path\/$ran\_genome_info.txt -p $path\/$ran\_output -w $win -s $step -m $missing -d $diff -n $thread_in\\n";
		my $return = &pbs_setting("$cj_exc$local\-cj_quiet -cj_ppn $thread_in -cj_mem $mem -cj_conda $conda -cj_qname calculation_$i -cj_sn $ran -cj_qout . $out");
		print BASH "$return\n";
	}
	unless ($local){
		close(BASH);
	}
	if ($exc == 1){
		unless ($local){
			&status($ran);
		}
	}
}

#plotting
$ck_files = "";
$ck_files = `ls -l $path\/$ran\_output \| grep \"introgression_${diff}.pdf\"`;
unless ($ck_files && $ow == 0 && $rp == 0){
	$out = "Rscript intro_plot_2.R -p $path\/$ran\_output -d $diff -t $list -gi $path\/$ran\_genome_info.txt $ci\-p1c $p1c -p2c $p2c\\n";
	&pbs_setting("$cj_exc$local\-cj_quiet -cj_conda $conda -cj_qname plotting -cj_sn $ran -cj_qout . $out");
	if ($exc == 1){
		unless ($local){
			&status($ran);
		}
	}
}
