#!/usr/bin/perl

#if -spl is used, recomb_spline_v3.R is required
#vcftools, bcftools and tabix are required
#keep_only_GT.pl, add_contig_info.pl are required

use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
my @tmp = split(/\//, $home);
my $uid = $tmp[-1];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
else {
	die "Cannot find required subroutine file: qsub_subroutine.pl\n";
}
my $r_env = '/home/hpc/crlee/miniconda3/envs/ben/bin';
my $dir = getcwd;

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;	
}

print "Input command line:\n";
print "perl beagle5_v4\.pl @ARGV\n\n";

my $exc; my $ran; my $sn; my $pre; my @recomb; my $spl = 0; my $target; my $pre; my $ref;
my $r_pre; my $t_pre; my $m_path; my $t_path; my $mem = "-cj_mem 16 "; my $list; my $repeat;
my $ibd = 0; my $win = 40; my $lod = 3; my $trim = 0.15; my $len = 1.5; my $scale; my $ribd = 0;
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
	if ($ARGV[$i] eq "\-ref"){
		if (-e $ARGV[$i+1]){
			$ARGV[$i+1] = &check_path($ARGV[$i+1]);
			$ref = $ARGV[$i+1];
			$r_pre = $ref;
			$r_pre =~ s/.gz$//; $r_pre =~ s/.vcf$//;	
		}
	}
	if ($ARGV[$i] eq "\-g"){
		if (-e $ARGV[$i+1]){
			$ARGV[$i+1] = &check_path($ARGV[$i+1]);
			$target = $ARGV[$i+1];
			$t_pre = $target;
			$t_pre =~ s/.gz$//; $t_pre =~ s/.vcf$//;
			my @tmp = split(/\//, $t_pre);
			pop(@tmp);
			$t_path = join("\/", @tmp);
		}
	}
	if ($ARGV[$i] eq "\-pre"){
		$pre = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-ibd"){
		$ibd = 1;
	}
	if ($ARGV[$i] eq "\-ribd"){
		$ribd = 1;
	}
	if ($ARGV[$i] eq "\-w"){
		$win = $ARGV[$i+1];
		unless ($win !~ /[^0-9\.]/){
			die "-w value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-lod"){
		$lod = $ARGV[$i+1];
		unless ($lod !~ /[^0-9\.]/){
			die "-lod value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-trim"){
		$trim = $ARGV[$i+1];
		unless ($trim !~ /[^0-9\.]/){
			die "-trim value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-length"){
		$len = $ARGV[$i+1];
		unless ($len !~ /[^0-9\.]/){
			die "-length value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-scale"){
		$scale = $ARGV[$i+1];
		unless ($scale !~ /[^0-9\.]/){
			die "-scale value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-m"){
		if (-d $ARGV[$i+1]){
			if ($ARGV[$i+1] =~ /\/$/){
				$ARGV[$i+1] =~ s/\/$//;
			}
			$m_path = $ARGV[$i+1];
			$m_path = &check_path($m_path);
			@recomb = <$ARGV[$i+1]\/*.bmap>;
			unless (@recomb){
				die BOLD "Cannot find the \*.bmap files.", RESET, "\n";
			}
			foreach my $j (0..$#recomb){
				$recomb[$j] = &check_path($recomb[$j]);
			}
		}
		elsif (-e $ARGV[$i+1]){
			$ARGV[$i+1] = &check_path($ARGV[$i+1]);
			push(@recomb, $ARGV[$i+1]);
			if ($ARGV[$i+1] =~ /\//){
				my @tmp = split(/\//, $ARGV[$i+1]);
				pop(@tmp);
				$m_path = join("\/", @tmp);
			}
			else {
				$m_path = ".";
				$m_path = &check_path($m_path);
			}
		}
		else {
			&usage;
			print BOLD "ERROR: Cannot find the recombination rate file.\n", RESET;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-spl"){
		if ($ARGV[$i+1] =~ /^\-/){
			$spl = 25;
		}
		elsif ($ARGV[$i+1] =~ /[^0-9]/){
			&usage;
			print BOLD "ERROR: -spl value should be an integer.\n", RESET;
			exit;
		}
		else {
			$spl = int($ARGV[$i+1]);
		}
	}
	if ($ARGV[$i] eq "\-mem"){
		unless ($ARGV[$i+1] !~ /[^0-9]/){
			print BOLD "ERROR: -mem value is wrong. Use default value.\n", RESET;
			$mem = "-cj_mem 16 ";
		}
		else {
			$mem = "-cj_mem $ARGV[$i+1] ";
		}
	}
	if ($ARGV[$i] eq "\-list"){
		unless (-e $ARGV[$i+1]){
			print BOLD "ERROR: Cannot find the list file.\n", RESET;
			exit;
		}
		$list = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-rep"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$repeat = $ARGV[$i+1];	
		} 
		else {
			print BOLD "ERROR: -rep value is wrong. Skip\n", RESET;
		}
	}
}
unless ($target){
	&usage;
	print BOLD "ERROR: -g is required.\n", RESET;
	exit;
}
if ($repeat){
	unless ($list){
		print BOLD "ERROR: if -rep is used, -list must be provided.\n", RESET;
		exit;
	}
}
if ($ibd == 1 || $ribd == 1){
	if ($win <= $len*3){
		print BOLD "ERROR: window: $win\, length: $len. -w must be at least 3 times the -length parameter.\n", RESET;
		exit;
	}
}

unless (-e "keep_only_GT.pl"){
	die "keep_only_GT.pl file is required.\n";
}
unless (-e "add_contig_info.pl"){
	die "add_contig_info.pl file is required.\n";
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
if (-e "my_bash_beagle5_1_$ran\.sh" && $sn != 1){
	$ran = undef;
	goto RE;
}
print "The qsub SN is: $ran\n";

@chrs = &chr_name($target, $pre);
if ($chrs[0] eq "no"){
	print BOLD "ERROR: Cannot detect contig/chromosome name from vcf file.\n", RESET;
	exit;
}

#convert -m and -ref format
open (BASH, ">my_bash_beagle5_1_$ran\.sh") || die BOLD "ERROR: Cannot write my_bash_beagle5_1_$ran\.sh: $!", RESET, "\n";
my $num = 0; my @brefs; my $out;
foreach (@chrs){
	$out = "";
	print BASH "qsub \.\/qsub_files\/$ran\_beagle5_1_$_\.q\n";
	if (scalar(@recomb) == 1 && $recomb[0] !~ /bmap$/ && $num == 0){
		my $check = `head -n 1 $recomb[0]`;
		chomp($check);
		my @ele_tmp = split (/\s+|\t+/, $check);
		unless (-d "$m_path\/recomb_map_$ran"){
			system("mkdir $m_path\/recomb_map_$ran");
		}
		if ($ele_tmp[0] =~ /^chr|^contig/i && $ele_tmp[1] =~ /^pos|^position|^ps/i && $ele_tmp[2] =~ /rate|cM\/Mb|/i && $ele_tmp[3] =~ /mp|map|\bcM\b/i){}
		else {
			print BOLD "ERROR: It seems that $recomb[0] is not in correct format.\n", RESET;
			exit;				
		}
		if ($spl > 0){
			unless (-e "$m_path\/recomb_map_$ran\/$_.spline.bmap"){
				$out .= "Rscript recomb_spline_v3.R $recomb[0] $m_path\/recomb_map_$ran 0 $spl $pre\\n";
				print "$recomb[0] modified. Spline function is used. df \= $spl.\n";
			}
		}
		if ($spl == 0){
			unless (-e "$m_path\/recomb_map_$ran\/$_.bmap"){
				$out .= "Rscript recomb_spline_v3.R $recomb[0] $m_path\/recomb_map_$ran 1 0 $pre\\n";
				print "$recomb[0] modified. Spline function is NOT used.\n";
			}
		}
	}
	if ($ref){
		my $noN_ref = $ref;
		if (-e $ref){
			$noN_ref =~ s/\.vcf\./\.$_\.noN\.vcf\./i;
			my $noN_ref_GT = $noN_ref;
			if ($noN_ref_GT =~ /GTonly/){}
			elsif ($noN_ref_GT =~ /raw/){
				$noN_ref_GT =~ s/raw/GTonly/;
			}
			elsif ($noN_ref_GT =~ /filtered/){
    			$noN_ref_GT =~ s/filtered/GTonly/;
			}
			else {
				$noN_ref_GT =~ s/\.vcf/\.GTonly\.vcf/;
			}
			unless (-e $noN_ref_GT){
				unless (-e $noN_ref){
					$out .= "bcftools view --threads 2 -Oz -o $noN_ref -M2 -r $_ $ref\\n";
					if ($ref =~ /gz$/){
						#$out .= "vcftools --gzvcf $ref --recode --recode-INFO-all --max-missing 1 --max-alleles 2 --chr $_ --stdout \| bgzip -c \> $noN_ref\\n";	
					}
					else {
						#$out .= "vcftools --vcf $ref --recode --recode-INFO-all --max-missing 1 --max-alleles 2 --chr $_ --stdout \| bgzip -c \> $noN_ref\\n";	
					}
				}
			}
			unless ($noN_ref =~ /GTonly/){
				unless (-e $noN_ref_GT){
					$out .= "perl keep_only_GT.pl $noN_ref\\n";
					$out .= "rm $noN_ref\\n";
				}
				$noN_ref = $noN_ref_GT;
			}
			unless (-e "$noN_ref\.tbi"){
				$out .= "tabix $noN_ref\\n";
			}
			#print "WARNING: Reference_vcf must be phased and have no missing allele. Otherwise, the program will stop.\n";
		}
		else {
			die "Cannot find the reference file.\n";
		}
		my $bref3f = $noN_ref;
		$bref3f =~ s/vcf\.gz$/bref3/;
		unless (-e $bref3f){
			$out .= "java -jar \$bref3 $noN_ref \> $bref3f\\n";
			push(@brefs, $bref3f);
		}
	}
	&pbs_setting("$exc$mem\-cj_quiet -cj_ppn 2 -cj_env $r_env -cj_qname beagle5_1_$_ -cj_sn $ran -cj_qout . $out");
	$num++;
}
close(BASH);
if ($exc){
	&status($ran, $uid);
}
if (scalar(@recomb) == 1 && $recomb[0] !~ /bmap$/){
	@recomb = <$m_path\/recomb_map_$ran\/*.bmap>;
}

#deal with sample list
my @list_lines; my @rand_lists;
if ($list){
	open (LIST, "<$list") || die BOLD "Cannot open $list: $!", RESET, "\n";
	@list_lines = <LIST>;
	close(LIST);
	my @samples = &sample_name($target);
	foreach (@list_lines){
		my @list_tmp = split(/\t|\s+/, $_);
		foreach my $i (0..$#samples){
			if ($samples[$i] eq $list_tmp[0]){
				splice(@samples, $i, 1);
			}
		}
	}
	if (scalar(@samples) < 1 && $repeat eq ""){
		print BOLD "WARNING: All samples in the list are used.\n", RESET;
	}
	else {
		my $out_file = $list;
		$out_file =~ s/txt$|list$/1.txt/;
		open (LIST, ">$out_file") || die BOLD "Cannot write $out_file: $!", RESET, "\n";
		print LIST join("\n", @samples), "\n";
		close(LIST);
		push(@rand_lists, $out_file);
	}
}
if ($repeat){
	if (@rand_lists){
		@rand_lists = ();
	}
	my %pop_list;
	foreach (@list_lines){
		my @list_tmp = split(/\t|\s+/, $_);
		unless ($list_tmp[1]){
			print BOLD "ERROR: The list format is wrong.\n", RESET;
			exit;			
		}
		unless (@{$pop_list{$list_tmp[1]}}){
			@{$pop_list{$list_tmp[1]}} = $list_tmp[0];
		}
		else {
			push(@{$pop_list{$list_tmp[1]}}, $list_tmp[0]);
		}
	}
	foreach my $i (1..$repeat){
		my @samples = &sample_name($target);
		my $out_file = $list;
		$out_file =~ s/txt$|list$/$i.txt/;
		my $remained_file = $out_file;
		$remained_file =~ s/txt$/remained.txt/;
		open (LIST, ">$out_file") || die BOLD "Cannot write $out_file: $!", RESET, "\n";
		open (LIST2, ">$remained_file") || die BOLD "Cannot write $remained_file: $!", RESET, "\n";
		foreach my $key (sort keys %pop_list){
			my @key_array = @{$pop_list{$key}};
			my $random_idx = int(rand($#key_array));
			my $include = $key_array[$random_idx];
			foreach my $j (0..$#samples){
				if ($samples[$j] eq $include){
					splice(@samples, $j, 1);
					print LIST2 "$samples[$j]\t$key\n";
					last;
				}
			}			
		}
		print LIST join("\n", @samples), "\n";
		close(LIST);
		close(LIST2);
		push(@rand_lists, $out_file);
	}
}

#main beagle5 program
my $rep_num = 1;
if (@rand_lists){
	$rep_num = scalar(@rand_lists);
}

my @ori_bfiles;
foreach my $x (1..$rep_num){
	my $current_list;
	open (BASH, ">my_bash_beagle5_2_$ran\.sh") || die BOLD "Cannot write my_bash_beagle5_2_$ran\.sh: $!", RESET, "\n";
	$num = 0; my $bref; my $map;
	if (@rand_lists){
		$current_list = "excludesamples\=$rand_lists[$x-1] ";
	}
	foreach (@chrs){
		$out = "";
		print BASH "qsub \.\/qsub_files\/$ran\_beagle5_2_$_\_$x\.q\n";
		if (@brefs){
			$bref = "ref\=$brefs[$num] ";
		}
		if ($x > 1){
			if ($ibd == 1){
				$ribd = 1;
			}
		}
		if (@recomb){
			foreach my $chr_recomb (@recomb){
				my $check = `head -n 1 $chr_recomb`;
				if ($check =~ /^\b$_\b/){
					$map = "map\=$chr_recomb ";
				}
			}
		}
		unless (-e "$t_pre\_$_\.bimputed.vcf.gz"){
			my $chr_target = $target;
			$chr_target =~ s/\.vcf/\.$_.vcf/i;
			my $chr_target_GT = $chr_target;
			if ($chr_target_GT =~ /GTonly/){}
			elsif ($chr_target_GT =~ /raw/){
				$chr_target_GT =~ s/raw/GTonly/;
			}
			elsif ($chr_target_GT =~ /filtered/){
    			$chr_target_GT =~ s/filtered/GTonly/;
			}
			else {
				$chr_target_GT =~ s/\.vcf/\.GTonly\.vcf/;
			}
			unless (-e $chr_target_GT){
				unless (-e $chr_target){
					$out .= "bcftools view --threads 2 -M2 -Oz -o $chr_target -r $_ $target\\n";
					if ($chr_target =~ /gz$/){
						#$out .= "vcftools --gzvcf $target --recode --recode-INFO-all --max-alleles 2 --chr $_ --stdout \| bgzip -c \> $chr_target\\n";
					}
					else {
						#$out .= "vcftools --vcf $target --recode --recode-INFO-all --max-alleles 2 --chr $_ --stdout \| bgzip -c \> $chr_target\\n";
					}
				}
				$out .= "perl keep_only_GT.pl $chr_target\\n";
				$out .= "rm $chr_target\\n";
			}
			$chr_target = $chr_target_GT;			
			$out .= "java -Xmx16g -jar \$beagle gt\=$chr_target $bref";
			$out .= "out\=$t_pre\_$_\.bimputed $map";
			$out .= "chrom\=$_ nthreads\=2\\n";
		}
		if ($ibd == 1 || $ribd == 1){
			unless (-e "$t_pre\_$_\.bimputed.$x.ibd" && $ribd == 0){
				if (@recomb){
					foreach my $chr_recomb (@recomb){
						my $check = `head -n 1 $chr_recomb`;
						if ($check =~ /^\b$_\b/){
							$map = "map\=$chr_recomb ";
						}
					}
				}
				$out .= "java -Xmx$mem\g -jar \$ribd gt\=$t_pre\_$_\.bimputed.vcf.gz ";
				$out .= "out\=$t_pre\_$_\.bimputed.$x $map";
				$out .= "chrom\=$_ nthreads\=2 ";
				$out .= "window\=$win lod\=$lod trim\=$trim length\=$len $current_list";
				if ($scale){
					$out .= "scale\=$scale";
				}
				$out .= "\\n";
			}
		}
		unless (-e "$t_pre\_$_\.bimputed\.m\.vcf\.gz"){
			$out .= "perl add_contig_info.pl -d $t_pre\_$_\.bimputed.vcf.gz -v $target\\n";
		}
		push(@ori_bfiles, "$t_pre\_$_\.bimputed.vcf.gz");
		&pbs_setting("$exc$mem\-cj_quiet -cj_ppn 2 -cj_env $r_env -cj_qname beagle5_2_$_\_$x -cj_sn $ran -cj_qout . $out");
		$num++;
	}
	close(BASH);
	if ($exc){
		&status($ran, $uid);
		if ($ibd == 1 || $ribd == 1){
			unless (-d "$t_path\/beagle_IBD_$ran"){
				system("mkdir $t_path\/beagle_IBD_$ran");
			}
			system("mv $t_pre\_\*\.bimputed\*bd.gz $t_path\/beagle_imputed_$ran");
		}
	}
}

my @m_files;
foreach (@ori_bfiles){
	my $m_file = $_;
	$m_file =~ s/bimputed\.vcf\.gz$/bimputed\.m\.vcf\.gz/;
	push(@m_files, $m_file);
	if (-e $m_file){
		system("rm $_");
	}
}
if ($exc){
	my $merged_file = "$t_pre\.all\.bimputed.vcf.gz";
	unless (-e $merged_file){
		$out = "bcftools concat -Oz -o $merged_file @m_files\\n";
		&pbs_setting("$exc$mem\-cj_quiet -cj_qname beagle5_concat -cj_sn $ran -cj_qout . $out");
	}
	&status($ran, $uid);
	unless (-e "$merged_file.tbi"){
		$out = "tabix $merged_file\\n";
		&pbs_setting("$exc$mem\-cj_quiet -cj_qname beagle5_tabix -cj_sn $ran -cj_qout . $out");
	}
	&status($ran, $uid);
}

if (@rand_lists){
	foreach (@rand_lists){
		system("rm $_");
	}
}

sub usage {
	print BOLD "Usage: perl beagle5_v4.pl -g target_vcf [-ref reference_vcf] [-m recombination_map_file] [-pre prefix_of_contig_name] [-spl degree_of_freedom] [-ibd] [-ribd] [-w window_size_for_ibd] [-lod lod_value_for_ibd] [-trim trim_value_for_ibd] [-length length_value_for_ibd] [-scale scale_value_for_ibd] [-sn serial_number] [-mem memory_assigned] [-list sample_pop_list_file] [-rep repeat_number_for_ibd] [-exc] [-h]\n", RESET;
	print "-g\tTarget vcf that you want to impute. \(eg. vcf from RAD \+ WGS\)\n";
	print "-ref\t\[optional\] Reference vcf contained complete SNPs. \(eg. vcf from WGS\)\n";
	print "-m\t\[optional\] Unfomated recombination map file. If you point to a folder, it will assume you have right recombinant \n\tmap files for each contig.\n";
	print "-pre\t\[optional\] Limit the imputation to contig/chromosome that have specific prefix of name.\n";
	print "-spl\t\[optional\] Turn on the spline function and define the degree of freedom for creating\/adjusting recombination\n\trate at each position.\n";
	print "\tThe value must be an integer. If don't define the value, it will use default value\: 25.\n";
	print "-list\t\[optional\] sample_population list for IBD analysis.\n";
	print "-ibd\t\[optional\] Do additional IBD analysis.\n";
	print "-ribd\t\[optional\] Re-do additional IBD analysis.\n";
	print "-w\t\[optional\] Set window size for IBD analysis.\n";
	print "-lod\t\[optional\] Set LOD value for IBD analysis.\n";
	print "-trim\t\[optional\] Set trimming value for IBD analysis.\n";
	print "-length\t\[optional\] Set length value for IBD analysis.\n";
	print "-scale\t\[optional\] Set scale value for IBD analysis.\n";
	print "-rep\t\[optional\] Set repeat times for IBD analysis. -list must be used.\n";
	print "-mem\t\[optional\] Set memory \(Gb\) used for the job. Default value: 16\n";
	print "-sn\t\[optional\] Serial number generated in each run.\n";
	print "-exc\tSend the job\(s\) for execution.\n";
	print "-h\tHelp.\n";
	return;
}

