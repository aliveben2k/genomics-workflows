#!/usr/bin/perl
use Term::ANSIColor qw(:constants);
use Cwd qw(getcwd);
use FindBin;

my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}
my $p_dir = $FindBin::Bin;
unless ($p_dir){
	$p_dir = ".";
}

print "This script is written by Ben Chien. Oct.2023\n";
print "Input command line:\n";
print "perl geovar_pipeline.pl @ARGV\n\n";

my $vcf; my $pop_list_f; my @pop_list; my $rand;
my $proj; my $ran; my $ow; my $sn; my $exc; my $o_path;
my $top50 = 0; my $ratio = 0; my $private = 0; my @rns;
for (my $i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] eq "\-vcf"){
		if (-e $ARGV[$i+1]){
			$vcf = $ARGV[$i+1];
		}
		else {
			&usage;
			print BOLD "Cannot find the vcf file.\n", RESET;
			exit;
		}
	}
    if ($ARGV[$i] eq "\-list"){
		if (-e $ARGV[$i+1]){
			$pop_list_f = $ARGV[$i+1];
			open(LIST, "<$pop_list_f") || die BOLD "Cannot open $pop_list_f: $!.", RESET, "\n";
			@pop_list = <LIST>;
			chomp(@pop_list);
			close(LIST);
			if ($ARGV[$i+1] =~ /\//){
				my @tmp_paths = split(/\//, $ARGV[$i+1]);
				pop(@tmp_paths);
				$o_path = join("\/", @tmp_paths);
			}
			else {
				$o_path = ".";
			}
		}
		else {
			&usage;
			print BOLD "Cannot find the list file.\n", RESET;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-rand"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$rand = $ARGV[$i+1];
		}
		else {
			die "-rand value is wrong.\n";
		}
	}
=start
	if ($ARGV[$i] eq "\-top50"){
		$top50 = 1;
	}
	if ($ARGV[$i] eq "\-r"){
		$ratio = 1;
	}
=cut
	if ($ARGV[$i] eq "\-p"){
		$private = 1;
	}
	if ($ARGV[$i] eq "\-ow"){
		$ow = 1;
	}
	if ($ARGV[$i] eq "\-sn"){
		$ran = $ARGV[$i+1];
		$sn = 1;
	}
	if ($ARGV[$i] eq "\-exc"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "\-proj"){
		$proj = "-cj_proj $ARGV[$i+1] ";
	}
    if ($ARGV[$i] eq "\-rn"){ #specify which loop number to run
        if ($ARGV[$i+1] !~ /[^0-9\-\,]/){
            my $rn = $ARGV[$i+1];
            if ($rn =~ /\,/){
            	my @tmp = split(/\,/, $rn);
            	foreach (@tmp){
            		if ($_ =~ /\-/){
            			my @start_end = split(/\-/, $rn);
            			foreach my $n ($start_end[0]..$start_end[1]){
            				push(@rns, $n);
            			}
            		}
            		else {
            			push(@rns, $_);
            		}
            	}
            }
            else {
            	if ($rn =~ /\-/){
            		my @start_end = split(/\-/, $rn);
            		foreach my $n ($start_end[0]..$start_end[1]){
            			push(@rns, $n);
            		}
            	}
            	else {
            		push(@rns, $rn);
            	}
            }
        }
        else {
 			print BOLD "Invalid loop number value.\n", RESET;
			exit;       	
        }
    }
}

unless ($vcf){
	print BOLD "Cannot find the input vcf file.\n", RESET;
	exit;
}

unless (-d "qsub_files"){
	system ("mkdir qsub_files");
}
unless (-d "qsub_files\/out"){
	system ("mkdir qsub_files\/out");
}

RE:
unless ($ran){
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "my_bash_geovar_1_$ran\.sh" && $sn != 1){
	$ran = undef;
	goto RE;
}
print "The qsub SN is: $ran\n";

unless ($rand) {
	$rand = 1;
}

my $out;
if ($rand >= 1){
	my $tmp_file = $pop_list_f;
	$tmp_file =~ s/txt$|list$/$rand.txt/;
	unless (-e $tmp_file && $ow == 0){
        $out = "perl random_sets_geovar.pl $pop_list_f $rand\\n";
		&pbs_setting("$exc$proj\-cj_quiet -cj_qname geovar_random -cj_sn $ran -cj_qout . $out");
	}
	if ($exc){
		&status($ran);
	}
}

#step0: subset samples from the vcf file
my @vcfs;
foreach my $z (1..$rand){
	if (@rns){ #if loop number is specified, only run the specified loop numbers.
		my $run = 0;
		foreach (@rns){
			if ($z == $_){
				$run = 1;
			}
		}
		if ($run == 0){
			next;
		}
	}
	my $list_file = $pop_list_f;
	$list_file =~ s/txt$|list$/$z.txt/;
    open (LIST, "<$list_file") || die BOLD "Cannot read $list_file: $!", RESET, "\n";
    my @content = <LIST>;
    chomp(@content);
    shift(@content);
    my @samples;
    foreach my $i (0..$#content){
    	my @tmp = split(/\t|\s+/, $content[$i]);
    	push(@samples, $tmp[0]);
    }
    my $sample = join("\,", @samples);
    open (BASH, ">my_bash_geovar_0_$ran\.sh") || die BOLD "Cannot write my_bash_geovar_0_$ran\.sh: $!", RESET, "\n";
    my $outname = $list_file;
    $outname =~ s/txt$/vcf.gz/;
    my $cntname = $list_file;
    $cntname =~ s/txt$|list$/freq.csv.cnt/;
    unless (-e $outname && $ow == 0){
    	if (-e $cntname && $ow == 0){
    		next;
    	}
        print BASH "qsub \.\/qsub_files\/$ran\_geovar_0_$z.q\n";
        $out = "bcftools view -Oz -o $outname --threads 4 -s $sample $vcf\\n";
        &pbs_setting("$exc$proj\-cj_quiet -cj_ppn 4 -cj_qname geovar_0_$z -cj_sn $ran -cj_qout . $out");
    }
    close(BASH);
}
if ($exc){
	&status($ran);
}

#step1: generate csv table files from the vcf file
foreach my $z (1..$rand){
	if (@rns){ #if loop number is specified, only run the specified loop numbers.
		my $run = 0;
		foreach (@rns){
			if ($z == $_){
				$run = 1;
			}
		}
		if ($run == 0){
			next;
		}
	}
	my $list_file = $pop_list_f;
	$list_file =~ s/txt$|list$/$z.txt/;
	my $curr_vcf = $list_file;
	$curr_vcf =~ s/txt$/vcf.gz/;
    my $csv_file = $list_file;
    $csv_file =~ s/txt$/freq.csv/;
    my $cntname = $list_file;
    $cntname =~ s/txt$|list$/freq.csv.cnt/;
    open (BASH, ">my_bash_geovar_1_$ran\.sh") || die BOLD "Cannot write my_bash_geovar_1_$ran\.sh: $!", RESET, "\n";
    unless (-e $csv_file && $ow == 0){
    	if (-e $cntname && $ow == 0){
    		next;
    	}
        print BASH "qsub \.\/qsub_files\/$ran\_geovar_1_$z.q\n";
        $curr_vcf = &check_path($curr_vcf);
        $list_file = &check_path($list_file);
        $csv_file = &check_path($csv_file);
        $out = "python geovar_freq.py $curr_vcf $list_file $csv_file\\n"; #generate the frequency file
        $out .= "rm $curr_vcf\\n";
        &pbs_setting("$exc$proj\-cj_quiet -cj_mem 18 -cj_qname geovar_1_$z -cj_sn $ran -cj_qout . $out");
    }
    close(BASH);
}
if ($exc){
	&status($ran);
}

#step2: generate count files from the csv files
foreach my $z (1..$rand){
	if (@rns){ #if loop number is specified, only run the specified loop numbers.
		my $run = 0;
		foreach (@rns){
			if ($z == $_){
				$run = 1;
			}
		}
		if ($run == 0){
			next;
		}
	}
	my $csv_file = $pop_list_f;
	$csv_file =~ s/txt$|list$/$z.freq.csv/;
    my $cnt_file = "$csv_file.cnt";
    open (BASH, ">my_bash_geovar_2_$ran\.sh") || die BOLD "Cannot write my_bash_geovar_2_$ran\.sh: $!", RESET, "\n";
    unless (-e $cnt_file && $ow == 0){
        print BASH "qsub \.\/qsub_files\/$ran\_geovar_2_$z.q\n";
        $out = "perl geovar_count_UC.pl $csv_file\\n";
        $out .= "rm $csv_file\\n";
        &pbs_setting("$exc$proj\-cj_quiet -cj_mem 12 -cj_qname geovar_2_$z -cj_sn $ran -cj_qout . $out");
    }
    close(BASH);
}
if ($exc){
	&status($ran);
}

#step3: generate private SNP counting files from the count files
foreach my $z (1..$rand){
	if (@rns){ #if loop number is specified, only run the specified loop numbers.
		my $run = 0;
		foreach (@rns){
			if ($z == $_){
				$run = 1;
			}
		}
		if ($run == 0){
			next;
		}
	}
	my $cnt_file = $pop_list_f;
	$cnt_file =~ s/txt$|list$/$z.freq/;
    if ($private == 1){
        open (BASH, ">my_bash_geovar_3_$ran\.sh") || die BOLD "Cannot write my_bash_geovar_3_$ran\.sh: $!", RESET, "\n";
        my $pdf = $cnt_file;
        #$pdf =~ s/\.txt$|\.list$//;
        unless (-e "$pdf.csv.Rplot_private.C_only.pdf" && $ow == 0){
            print BASH "qsub \.\/qsub_files\/$ran\_geovar_3_$z.q\n";
            $out = "Rscript geovar_private_plot_UC.R $pdf.csv.Rplot_private.txt\\n";
            &pbs_setting("$exc$proj\-cj_quiet -cj_qname geovar_3_$z -cj_sn $ran -cj_qout . $out");
        }
        close(BASH);
    }
}

sub usage {
	print BOLD "Usage: perl perl geovar_pipeline.pl -vcf VCF_FILE -list POP_LIST [-rand INT] [-p] [-proj PROJECT_NAME] [-rn NUMBER] [-ow] [-sn] [-exc]\n", RESET;
	return 1;
} #print usage
