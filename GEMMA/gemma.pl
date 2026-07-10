#!/usr/bin/perl

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
my $r_env = '/home/hpc/crlee/miniconda3/envs/R-4.1/bin';
my $gemma = 'gemma-0.98.6';
my $check_gemma = `$gemma`;
if ($check_gemma !~ /GEMMA/){
	$gemma = 'gemma';
	$check_gemma = `$gemma`;
	if ($check_gemma !~ /GEMMA/){
		die "Cannot find the GEMMA program.\n";
	}
}

chomp(@ARGV);
print "The script is written by Ben Chien. Feb. 2022. Ver.2.0\n";
print "Input command line:\n";
print "perl gemma\.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;	
}

my $p_dir = $FindBin::Bin;
unless ($p_dir){
	$p_dir = ".";
}
my $exc; my $ran; my $sn; 
my $cov; my $anno; my $eigen = 0; my $bslmm; my $inverse = 0; my $ioff;
my $phenofile; my $predict = 0; my $ow = 0; my $bimbam; my $maf; my $lm = 0;
my $o_path; my @vcf; my $pre; my $ec = 0; my $mem; my $highlight; my $thr = 1;
my $local;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-exc"){
		$exc = "-cj_exc ";
	}
	if ($ARGV[$i] eq "-sn"){
		$ran = $ARGV[$i+1];
		$sn = 1;
	}
	if ($ARGV[$i] eq "\-h"){
		&usage;
		exit;
	}
	if ($ARGV[$i] eq "\-o"){
		$o_path = $ARGV[$i+1];
		if ($o_path =~ /\/$/){
			$o_path =~ s/\/$//;
		}
		unless (-d $o_path){
			my $return = `mkdir $o_path`;
			if ($return){
				die "Cannot make the directory $o_path: $return\n";
			}
		}
	}
    if ($ARGV[$i] eq "\-vcf"){
        if (-e $ARGV[$i+1] && ($ARGV[$i+1] =~ /vcf.gz$/ || $ARGV[$i+1] =~ /vcf$/)){
            @vcf = $ARGV[$i+1];
        }
        elsif (-e $ARGV[$i+1] && ($ARGV[$i+1] =~ /list$/ || $ARGV[$i+1] =~ /txt$/)){
        	open(VCFS, "<$ARGV[$i+1]") || die "Cannot find $ARGV[$i+1]: $!\n";
        	@vcf = <VCFS>;
        	chomp(@vcf);
        }
        else {
            die "Cannot find the vcf.\n";
        }
    }
    if ($ARGV[$i] eq "\-g"){
        if (-e $ARGV[$i+1] && ($ARGV[$i+1] =~ /bimbam.gz$/ || $ARGV[$i+1] =~ /bimbam$/)){
            $bimbam = $ARGV[$i+1];
        }
        else {
            die "Cannot find bimbam file.\n";
        }
    }
	if ($ARGV[$i] eq "\-p"){
		if (-e $ARGV[$i+1]){
			$phenofile = $ARGV[$i+1];
		}
		else {
			print "Cannot find the phenotype file.\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-a"){
		if (-e $ARGV[$i+1]){
			$anno = "-a $ARGV[$i+1] ";
		}
		else {
			print "Cannot find annotation file.\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-c"){
		if (-e $ARGV[$i+1]){
			$cov = $ARGV[$i+1];
		}
		else {
			print "Cannot find the covariants file.\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-ec"){
		if (-e $ARGV[$i+1]){
			$cov = $ARGV[$i+1];
			$ec = 1;
		}
		else {
			print "Cannot find the covariants file.\n";
			exit;
		}
	}
	if ($ARGV[$i] eq "\-predict"){
        $predict = 1;
	}
	if ($ARGV[$i] eq "\-nor"){
        $inverse = 1;
	}
	if ($ARGV[$i] eq "\-bslmm"){
        if ($ARGV[$i+1] eq 1 || $ARGV[$i+1] eq 2 || $ARGV[$i+1] eq 3){
        	$bslmm = $ARGV[$i+1];
        }
        else {
        	print "-bslmm is not set correctly, using default value: 1.\n";
        	$bslmm = 1;
        }
	}
	if ($ARGV[$i] eq "\-eigen"){
        $eigen = 1;
	}
	if ($ARGV[$i] eq "\-prefix"){
        $pre = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-maf"){
        $maf = $ARGV[$i+1];
        unless ($maf !~ /[^0-9\.]/){
        	die "-maf value is wrong, only numbers can be accepted.\n";
        }
	}
	if ($ARGV[$i] eq "\-ow"){
        $ow = 1;
	}
	if ($ARGV[$i] eq "\-local"){
        $local = "-cj_local ";
	}
	if ($ARGV[$i] eq "\-mem"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$mem = "-cj_mem $ARGV[$i+1] ";
		}
	}
	if ($ARGV[$i] eq "-ioff"){
		$ioff = "-ioff ";
	}
	if ($ARGV[$i] eq "-renv"){
		$r_env = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-hl"){
		$highlight = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-lm"){
		$lm = 1;
	}
	if ($ARGV[$i] eq "\-n"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$thr = $ARGV[$i+1];
		}
		else {
			print "-n value is wrong, using default value: 1.\n";
		}
	}
}
unless (@vcf || $bimbam){
    die "-vcf or -g is required.\n";
}
if ($bimbam){
    unless (@vcf){
        die "-vcf is required.\n";
    }
}
if ($predict == 1){
	unless ($bslmm){
		print "For prediction function, -bslmm is not set, using default value: 1.\n";
		$bslmm = 1;
	}
}

my @samples = &sample_name($vcf[0]);
if (scalar(@samples) < 58){
	my $s_cnt = scalar(@samples);
	print "Sample size is too small for Gemma analysis.\n";
	print "Required size: \>57\tActural size: $s_cnt\n";
	undef($s_cnt);
	exit;
}

unless ($phenofile){
    print "-p should be used to indicate phenotype file path.\n";
    exit;
}

RE:
if ($ran){}
else {
	$ran = &rnd_str(4, "A".."Z", 0..9);
}
if (-e "\.\/qsub_files\/$ran\_gemma_0_1.q" && $sn != 1){
	$ran = undef;
	goto RE;
}
print "The qsub SN is: $ran\n";
my $out;
unless ($bimbam){
    my $bout = $phenofile;
    $bout =~ s/txt$|list$/$ran.bimbam\.gz/;
    unless (-e $bout && $ow == 0){
        if (scalar(@vcf) == 1){
            unless (-e "$vcf[0].tbi" || -e "$vcf[0].idx"){
                die "$vcf[0] needs to be indexed first.\n";
            }
    	}
    	my @chrs = &chr_name($vcf[0], $pre);
    	my $c_cnt = 1;
    	my @tmp_outs;
    	if (scalar(@vcf) == 1){
    		foreach (@chrs){
    			$out = "";
        		my $tmp_vcf = $vcf[0];
        		$tmp_vcf =~ s/.gz$//;
        		$tmp_vcf =~ s/.vcf$/.$c_cnt.vcf.gz/;
        		$tmp_out = $tmp_vcf;
        		$tmp_out =~ s/.vcf.gz/.$ran.bimbam.gz/;;
        		unless (-e $tmp_vcf && $ow == 0){
                    $out = "tabix -h $vcf[0] $_ \| bgzip -c \> $tmp_vcf\\n";
        		}
        		unless (-e $tmp_out && $ow == 0){
                    $out .= "perl $p_dir\/vcf2bimbam_thread.pl -vcf $tmp_vcf -list $phenofile $ioff\-sn $ran -n $thr\\n";
        		}
        		&pbs_setting("$exc$local\-cj_quiet -cj_ppn $thr -cj_qname gemma_0_$c_cnt -cj_sn $ran -cj_qout . $out");
        		push(@tmp_outs, $tmp_out);
        		$c_cnt++;
        	}
        }
        else {
        	foreach (@vcf){
        		$out = "";
                my $tmp_out = $_;
        		$tmp_out =~ s/.gz$//;
        		$tmp_out =~ s/.vcf$/.$ran.bimbam.gz/;
        		unless (-e $tmp_out && $ow == 0){
        			if (-e "$p_dir\/vcf2bimbam_thread.pl"){
                    	$out = "perl $p_dir\/vcf2bimbam_thread.pl -vcf $_ -list $phenofile $ioff\-sn $ran -n $thr\\n";
                    }
                    elsif (-e "$p_dir\/vcf2bimbam.pl"){
                    	$out = "perl $p_dir\/vcf2bimbam.pl -vcf $_ -list $phenofile $ioff\-sn $ran\\n";
                    	if ($thr > 1){
                    		print "WARNING: Cannot find the vcf2bimbam_thread.pl, using the single thread script for bimbam format transformation.\n";
                    	}
                    }
                    else {
                    	die "ERROR: Cannot find the required vcf2bimbam.pl or vcf2bimbam_thread.pl file.\n";
                    }
        		}
        		&pbs_setting("$exc$local\-cj_quiet -cj_ppn $thr -cj_qname gemma_0_$c_cnt -cj_sn $ran -cj_qout . $out");
        		push(@tmp_outs, $tmp_out);  		
        		$c_cnt++;	
        	}
        }
        if ($exc eq "-cj_exc "){
            &status($ran);
            unless (-e $bout && $ow == 0){
                system("cat @tmp_outs \> $bout");
                foreach (@tmp_outs){
                    #system("rm $_");
                }
            }
        }
    }
    $bimbam = $bout;
}

#get the output path of bimbam path, if no output path defined 
unless ($o_path){
	$o_path = $bimbam;
	if ($o_path =~ /\//){
		my @o_tmp1 = split(/\//, $o_path);
		pop(@o_tmp1);
		$o_path = join("\/", @o_tmp1);
	}
	else {
		$o_path = ".";
	}
	$o_path .= "\/output_$ran";
	unless (-d $o_path){
		my $return = `mkdir $o_path`;
		if ($return){
			die "Cannot make the directory $o_path: $return\n";
		}
	}
}
#get the name of bimbam file
my @o_tmp2 = split(/\//, $bimbam);
my $o_prefix = pop(@o_tmp2);
$o_prefix =~ s/.bimbam.gz$//;

$bimbam = "-g $bimbam ";

#deal with covariant file
my @pfile; my @covs; my @covs_out; my $c_check; my $cov_out;
if ($cov){
	print "Modifying the covariant file.\n";
	open (COV, "<$cov") || die "Cannot open $cov: $!\n";
	@covs = <COV>;
	chomp(@covs);
	shift(@covs);
	close(COV);
	open (PHENO, "<$phenofile") || die "Cannot open $phenofile: $!\n";
	@pfile = <PHENO>;
	chomp(@pfile);
	close(PHENO);
	shift(@pfile);
	foreach my $j (0..$#pfile){
    	$pfile[$j] =~ s/[\x0A\x0D]//g;
    	my @tmp = split(/\t+|\s+/, $pfile[$j]);
    	$c_check = 0;
    	if ($ec == 0){
            foreach my $i (0..$#covs){
                $covs[$i] =~ s/[\x0A\x0D]//g;
                if ($covs[$i] =~ /^\b$tmp[0]\b/){
                    my @tmp2 = split(/\t+|\s+/, $covs[$i]);
                    shift(@tmp2);
                    unshift(@tmp2, "1");
                    push(@covs_out, join("\t", @tmp2));
                    $c_check = 1;
                    last;
                }
            }
    	}
    	else {
            $covs[$j] =~ s/[\x0A\x0D]//g;
            if ($covs[$j] =~ /^\b$tmp[0]\b/){
                my @tmp2 = split(/\t+|\s+/, $covs[$j]);
                shift(@tmp2);
                unshift(@tmp2, "1");
                push(@covs_out, join("\t", @tmp2));
                $c_check = 1;
            }
    	}
    	if ($c_check == 0){
            if ($ec == 0){
                die "Cannot find covariant values for $tmp[0]\n";
    		}
    		else {
                die "The sample order in the phenotype file does not match the covariant file.\n";
    		}
    	}
	}
	$cov_out = $cov;
	$cov_out =~ s/.txt$//;
	$cov_out .= ".cov";
	open (COVOUT, ">$cov_out") || die "Cannot write $cov_out: $!\n";
	print COVOUT join("\n", @covs_out), "\n";
	close(COVOUT);
	$cov = "-c $cov_out ";
}

#deal with phenotype file
my $line1; my @phenos; my $num;
print "Modifying the phenotype file.\n";
open (PHENO, "<$phenofile") || die "Cannot open $phenofile: $!\n";
@pfile = <PHENO>;
chomp(@pfile);
close(PHENO);
$line1 = shift(@pfile);
$line1 =~ s/[\x0A\x0D]//g;
@phenos = split(/\t+/, $line1); #1st element is "ID"
foreach my $i (0..$#phenos){
    if ($phenos[$i] =~ /\s+|\(|\)|\[|\]|\-/){
        $phenos[$i] =~ s/\s+|\(|\[/_/g;
        $phenos[$i] =~ s/\)|\]|\-//g;
        $phenos[$i] =~ s/_+/_/g;
    }
}
$num = $#phenos-1;
my $phem = $phenofile;
$phem =~ s/.txt$//;
$phem = $phem."\.pheno";
my $pheno_cnt = 0;
open (PHENOM, ">$phem") || die "Cannot write $phem: $!\n";
foreach (@pfile){ #phenotype file
    $_ =~ s/[\x0A\x0D]//g;
    my @tmp = split(/\t+|\s+/, $_);
    shift(@tmp); #remove "ID"
    $_ = join("\t", @tmp);
    print PHENOM "$_\n";
    $pheno_cnt++;
}
close(PHENOM);

#check if bimbam file exists
my $bimbam_check = $bimbam;
$bimbam_check =~ s/\-g|\s+//g;
unless (-e $bimbam_check){
	die "Please execute gemma.pl to generate *.bimbam file before continuing the pipeline.\n";
}
#check if sample sizes of the genotype and phenotype are the same
$bimbam_check = `zcat $bimbam_check \| head -n 1`;
@bimbam_eles = split("\, ", $bimbam_check);
if ($pheno_cnt != scalar(@bimbam_eles)-3){
	my $g_cnt = scalar(@bimbam_eles)-3;
	print "Phenotype sample size is not equal to genotype sample size.\n";
	print "Phenotype sample size: $pheno_cnt\n";
	print "Genotype sample size: $g_cnt\n";
	exit;
}
$phenofile = "-p $phem ";

#deal with covariants file used for BSLMM
my $bimbam_cov = 0; my $bimbam_ori; my $cov_tmp;
if ($predict == 1 || $bslmm){
	if ($cov_out){
		$bimbam_ori = $bimbam;
		$bimbam_ori =~ s/\-g|\s+//g;;
		$bimbam_cov = $bimbam_ori;
		$bimbam_cov =~ s/.gz$//;
		$bimbam_cov =~ s/txt$|bimbam$/cov.bimbam.gz/;
		$cov_tmp = $cov_out;
		$cov_tmp .= ".tmp.gz";
		if (-e $bimbam_cov && $ow == 0){
			goto PREDICT;
		}
		open(COV, "<$cov_out") || die "Cannot open $cov_out file: $!\n";
		@covs = <COV>;
		chomp(@covs);
		close(COV);
		my $num = scalar(split(/\t|\s/, $covs[0])) - 1;
		my %idx;
		foreach (@covs){
			my @line_eles = split(/\t|\s/, $_);
			foreach my $i (1..$num){
				my $j = $i;
				if ($j < 10){
					$j = "00$j";
				}
				if ($j >= 10 && $j < 100){
					$j = "0$j";
				}
				push(@{$idx{$j}}, $line_eles[$i]);
			}
		}
		open(COVOUT, "|-", "gzip \> $cov_tmp") || die "Cannot write $cov_tmp file: $!\n";
		foreach my $k (sort keys %idx){
			my $line = join("\, ", @{$idx{$k}});
			print COVOUT "Covariant_$k\, X\, Y\, $line\n";
		}
		close(COVOUT);
		unless (-e $bimbam_cov && $ow == 0){
    		$out = "cat $bimbam_ori $cov_tmp \> $bimbam_cov\\n";
    		$out .= "rm $cov_tmp\\n";
    		&pbs_setting("$exc$local\-cj_quiet -cj_qname gemma_bimbam_cov -cj_sn $ran -cj_qout . $out");
			if ($exc eq "-cj_exc "){
				&status($ran);
			}    		
    	}
	}
}

#deal with predict function
PREDICT:
my $prdt_out;
my $redo_kinship = 0;
$out = "";
if ($predict == 1){
	#deal with kinship matrix (temp, for prediction)
	my @matrixes = &kinship($local, $node, $ppn, $ow, $exc, $ran, $mem, $gemma, $o_path, $o_prefix, $eigen, $ow, $bimbam, $phenofile, $redo_kinship, \@phenos);
	&bslmm($local, $exc, $predict, $ran, $o_path, $o_prefix, $ow, $gemma, $bimbam, $bimbam_cov, $phenofile, $anno, $bslmm, \@phenos, \@matrixes);
    $prdt_out = $phem;
    $prdt_out =~ s/pheno$/prdt.pheno/;
    unless (-e $prdt_out && $ow == 0){
    	$out = "perl $p_dir\/combine_prdt_pheno.pl $phem $o_path\/$o_prefix\\n";
    }
    $phenofile = "-p $prdt_out ";
    &pbs_setting("$exc$local\-cj_quiet -cj_qname gemma_combine_prdt -cj_sn $ran -cj_qout . $out");
	if ($exc == 1){
		&status($ran);
	}    
}
#deal with BSLMM function
elsif ($bslmm){
	&bslmm($local, $exc, $predict, $ran, $o_path, $o_prefix, $ow, $gemma, $bimbam, $bimbam_cov, $phenofile, $anno, $bslmm, \@phenos, \@matrixes);
}

#deal with inverse normal transformation of the phenotype file
my $invs_out;
if ($inverse == 1){
	$out = "";
	if ($predict == 1){
		$invs_out = $prdt_out;
	}
	else {
		$invs_out = $phem;
	}
    $invs_out =~ s/.pheno$/.inverse.pheno/;
    unless (-e $invs_out && $ow == 0){
		if ($predict == 1){
			$out = "Rscript $p_dir\/inverse_normal_phenotype_gemma.R $prdt_out\\n";
    	}
    	else {
			$out = "Rscript $p_dir\/inverse_normal_phenotype_gemma.R $phem\\n";
    	}
    }
    $phenofile = "-p $invs_out ";
    &pbs_setting("$exc$local\-cj_quiet -cj_env $r_env -cj_qname gemma_inverse -cj_sn $ran -cj_qout . $out");
	if ($exc eq "-cj_exc "){
		&status($ran);
	}
	$redo_kinship = 1;
}

#deal with kinship matrix (final)
my @matrixes;
if ($lm == 0){
	@matrixes = &kinship($local, $node, $ppn, $ow, $exc, $ran, $mem, $gemma, $o_path, $o_prefix, $eigen, $bimbam, $phenofile, $redo_kinship, \@phenos);
}
#generating association files for each phenotype
my $cnt = 0; my $rep = 0; my $use_tmp_pheno = 0;
do {
    for (my $i=1; $i<=$#phenos; $i++){ #skip "ID" at $phenos[0]
    	$out = "";
        if ($maf !~ /[0-9]/){
            if (-e "$o_path\/trait_$phenos[$i]\.assoc" && -e "$o_path\/trait_$phenos[$i]\.man_plot.tiff" && $ow == 0){
                $cnt++;
                next;
            }
        }
        else {
            if (-e "$o_path\/trait_$phenos[$i].maf_$maf.assoc" && -e "$o_path\/trait_$phenos[$i].maf_$maf.man_plot.tiff" && $ow == 0){
                $cnt++;
                next;
            }
        }
        if ($rep == 5 && $use_tmp_pheno == 0){
        	$use_tmp_pheno = 1;
        	$rep = 0;
        	unless (-e "$o_path\/trait_$phenos[$i]\.assoc" || -e "$o_path\/trait_$phenos[$i].maf_$maf.assoc"){
        		if ($phenofile !~ /\.tmp\./){
        			if ($inverse == 1){
        				$phenofile =~ s/inverse.pheno$/tmp.inverse.pheno/;
        			}
        			else {
        				$phenofile =~ s/pheno$/tmp.pheno/;
        			}
        		}
        	}
        }
        if (-e "$o_path\/trait_$phenos[$i]\.assoc" && $ow == 0){
            if ($maf =~ /[0-9]/){
                unless (-e "$o_path\/trait_$phenos[$i].maf_$maf.assoc" && $ow == 0){
                    $out = "perl $p_dir\/modify_assoc.pl $o_path\/trait_$phenos[$i]\.assoc $maf\\n";
                }
            }
            else {
                $out = "perl $p_dir\/modify_assoc.pl $o_path\/trait_$phenos[$i]\.assoc\\n";
            }
        }
        else {
        	#remove singular column in the cov file
        	my $check_pheno = $phenofile;
        	$check_pheno =~ s/^\-p|\s+//g;
        	my $check_cov = $cov;
        	$check_cov =~ s/^\-c|\s+//g;
        	my $current_cov = $cov;
        	if (-e "remove_singular.R" && -e $check_cov){
				my $c_return = `Rscript $p_dir\/remove_singular.R $i $check_pheno $check_cov`;
				if (-e "$check_cov.tmp$i"){
					$current_cov = "-c $check_cov.tmp$i ";
				}
				elsif ($c_return =~ /F/){
					$current_cov = "";
				}
				else {
					$current_cov = "-c $current_cov ";
				}
        	}
        	elsif (-e $check_cov){
				print "WARNING: Cannot find the remove_singular.R file. The covariants file is not checked. It might cause an singular maxtrix error.\n";
        	}
        	if ($lm == 0){
            	$out = "$gemma $bimbam$phenofile$anno$matrixes[$i-1]$current_cov\-lmm 4 -n $i -outdir $o_path -o trait_$phenos[$i]\\n";
            }
            else {
            	open(TMP, "<$check_pheno") || die "Cannot open $check_pheno file: $!\n";
            	my @tmp_pheno = <TMP>;
            	chomp(@tmp_pheno);
            	close(TMP);
            	open(POUT, ">$check_pheno.tmp$i") || die "Cannot write $check_pheno.tmp$i file: $!\n";
            	foreach (@tmp_pheno){
            		my @ptmp_eles = split(/\t+|\s+/, $_);
            		print POUT "$ptmp_eles[$i-1]\n";
            	}
            	close(POUT);
            	my $phenofile2 = "-p $check_pheno.tmp$i ";
            	$out =  "$gemma $bimbam$phenofile2$anno$current_cov\-lm 4 -outdir $o_path -o trait_$phenos[$i]\\n";
            	$out .= "rm $check_pheno.tmp$i\\n";
            }
            $out .= "mv $o_path\/trait_$phenos[$i]\.assoc.txt $o_path\/trait_$phenos[$i]\.assoc\\n";
            $out .= "perl $p_dir\/modify_assoc.pl $o_path\/trait_$phenos[$i]\.assoc $maf\\n";
        }
        my $maf_out;
        if ($maf =~ /[0-9]/){
            $maf_out = ".maf_$maf";
        }
        unless (-e "$o_path\/trait_$phenos[$i]$maf_out.man_plot.tiff" && $ow == 0){
            $o_path = &check_path($o_path);
            if ($highlight){
            	$highlight = "-hl $highlight";
            }
            my $lm_on;
            if ($lm == 1){
            	$lm_on = "-lm ";
            }
            $out .= "Rscript $p_dir\/qqman_v2.R -f $o_path\/trait_$phenos[$i]$maf_out.assoc -o $o_path\/trait_$phenos[$i]$maf_out $lm_on$highlight\\n";
        }
        &pbs_setting("$exc$mem$local\-cj_env $r_env -cj_quiet -cj_qname gemma_main_$i -cj_sn $ran -cj_qout . $out");
    }
    if ($exc eq "-cj_exc "){
        &status($ran);
    }
	my $check_cov = $cov;
	$check_cov =~ s/^\-c|\s+//g;
	$check_cov =~ s/.txt$//;
    $rm_tmps = `echo $check_cov.tmp\*`;
    if ($rm_tmps !~ /\*/){
		system("rm $check_cov.tmp\*");
    }
    $rep++;
} until ($cnt == $num || $exc == 0 || $rep > 5);
if (-e "$o_path\/$o_prefix.prdt.txt"){
    system("rm $o_path\/$o_prefix.prdt.txt");
}
#system("rm $o_path\/*.tmp.*");

sub usage {
	print BOLD "Usage: perl gemma.pl -vcf VCF_FILE_PATH -p PHENOTYPE_FILE [-g BIMBAM_FILE] [-ioff] [-c\|-ec COVARIANTS_FILE] [-a ANNOTATION_FILE] [-o OUTPUT_PATH] [-eigen] [-lm] [-bslmm INT] [-predict] [-nor] [-prefix PREFIX_OF_THE_CONTIG_NAME] [-maf MAF_VALUE] [-mem MEMORY_USE_IN_GB] [-n INT] [-ow] [-sn SERIAL_NUMBER] [-local] [-exc] [-h]\n\n", RESET;
	return;
}
sub kinship {
    my $local = shift; my $node = shift; my $ppn = shift; my $ow = shift; my $exc = shift; my $ran = shift; my $mem = shift; my $gemma = shift; my $o_path = shift; my $o_prefix = shift; my $eigen = shift; my $bimbam = shift; my $phenofile = shift; my $redo_kinship = shift; my @phenos = @{$_[-1]};
    my @matrixes;
    for my $i (1..$#phenos){
    	my $out;
        my $matrix;
        unless (-e "$o_path\/$o_prefix.$i\.cXX.txt" && $ow == 0 && $redo_kinship == 0){
            $out = "$gemma $bimbam$phenofile\-gk 1 -n $i -outdir $o_path -o $o_prefix.$i\\n";
        }
        $matrix = "-k $o_path\/$o_prefix.$i\.cXX.txt ";
        if ($eigen == 1){
            unless (-e "$o_path\/$o_prefix.$i.eigenD.txt" && -e "$o_path\/$o_prefix.$i.eigenU.txt" && $ow == 0 && $redo_kinship == 0){
                $out .= "$gemma $bimbam$phenofile$matrix\-eigen -n $i -outdir $o_path -o $o_prefix.$i\\n";
            }
            $matrix ="-d $o_path\/$o_prefix.$i.eigenD.txt -u $o_path\/$o_prefix.$i.eigenU.txt ";
        }
        push(@matrixes, $matrix);
        &pbs_setting("$exc$mem$local\-cj_quiet -cj_qname gemma_kinship_$i -cj_sn $ran -cj_qout . $out");
    }
    if ($exc =~ /-cj_exc/){
        &status($ran);
    }
    return @matrixes;
}
sub bslmm {
	my $local = shift; my $exc = shift; my $predict = shift; my $ran = shift; my $o_path = shift; my $o_prefix = shift; my $ow = shift; my $gemma = shift;
	my $bimbam = shift; my $bimbam_cov = shift; my $phenofile = shift; my $anno = shift; my $bslmm = shift; my @phenos = @{$_[-2]}; my @matrixes = @{$_[-1]};
	for my $i (1..$#phenos){
		my $out;
    	unless (-e "$o_path\/$o_prefix.$i.param.txt" && -e "$o_path\/$o_prefix.$i.log.txt" && -e "$o_path\/$o_prefix.$i.bv.txt" && $ow == 0){
        	if ($bimbam_cov != 0){
        		my $bimbam_cov_out = "-g $bimbam_cov ";
        		$out = "$gemma $bimbam_cov_out$phenofile$anno\-bslmm $bslmm -n $i -outdir $o_path -notsnp -o $o_prefix.$i\\n";
        	}
        	else {
        		$out = "$gemma $bimbam$phenofile$anno\-bslmm $bslmm -n $i -outdir $o_path -o $o_prefix.$i\\n";
    		}
    	}
    	if ($predict == 1){
    		unless (-e "$o_path\/$o_prefix.$i.prdt.txt" && $ow == 0){
        		$out .= "$gemma $bimbam$phenofile$matrixes[$i-1]\-epm $o_path\/$o_prefix.$i.param.txt -emu $o_path\/$o_prefix.$i.log.txt -ebv $o_path\/$o_prefix.$i.bv.txt -predict 1 -n $i -outdir $o_path -o $o_prefix.$i\\n";
    		}
    	}
    	&pbs_setting("$exc$local\-cj_quiet -cj_qname gemma_bslmm_$i -cj_sn $ran -cj_qout . $out");
	}
	if ($exc =~ /-cj_exc/){
		&status($ran);
	}
	return 1;
}
