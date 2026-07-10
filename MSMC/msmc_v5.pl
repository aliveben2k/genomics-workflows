#!/usr/bin/perl
use Cwd qw(getcwd);
use FindBin;
use Term::ANSIColor qw(:constants);
use Time::HiRes qw(gettimeofday);

my $start_time = gettimeofday();
my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}

#please set the python3 and the bedtools environments path here
#otherwise the server cannot find the program
my $env = '-cj_conda R-4.4 ';
my $bedtool = '-cj_conda bedtools ';

my @server = `ip addr`;
chomp(@server);
my $serv = 0; my $thread = 6;
foreach (@server){
	if ($_ =~ /140.110.148.11/ || $_ =~ /140.110.148.12/){
		$env = '-cj_env $HOME/miniconda3/envs/ben/bin ';
		$thread = 20;
	}
}
my $mu = '1e-8';
my $gen = '1';

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}

print "This script is written by Ben Chien. Sep.2023\n";
print "Input command line:\n";
print "perl msmc_v5.pl @ARGV\n\n";

my @vcfs; my @pop_list; my $ran; my $exc; my $mask; my $ow = 0; my $sn; my $debug = 0;
my $o_path; my $syn; my $psa; my $rand; my $pop_list_f; my $synnum = 2; my @rns;
my $p_dir = $FindBin::Bin; my $out;
unless ($p_dir){
	$p_dir = ".";
}
my $pre; my $proj;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-vcf"){
		if (-d $ARGV[$i+1]){
			$ARGV[$i+1] =~ s/\/$//;
			@vcfs = <$ARGV[$i+1]\/*.vcf.gz>;
			unless (@vcfs){
				@vcfs = <$ARGV[$i+1]\/*.vcf>;
				unless (@vcfs){
					&usage;
					print BOLD "Cannot find the vcf file.\n", RESET;
					exit;					
				}
			}
		}
		elsif (-e $ARGV[$i+1]){
			push(@vcfs, $ARGV[$i+1]);
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
	if ($ARGV[$i] eq "\-random"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$rand = $ARGV[$i+1];
		}
		else {
			die "-random value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-t"){ #thread used
        $thread = $ARGV[$i+1];
		if ($thread =~ /[^0-9]/){
			print BOLD "-t only accept 0-9.\n", RESET;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-debug"){
		$debug = 1;
	}
	if ($ARGV[$i] eq "\-pre"){
		$pre = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-proj"){
		$proj = "-cj_proj $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-p"){ #population separation analysis
		if ($ARGV[$i+1] =~ /[^0-9\*\+]/){
			print BOLD "-p only accept 0-9, \* and \+.\n", RESET;
			exit;
		}
        $psa = "\-p $ARGV[$i+1] ";
	}
	if ($ARGV[$i] eq "\-mask" || $ARGV[$i] eq "\-mk"){
		$mask = $ARGV[$i+1];
		unless (-e $mask) {
			print BOLD "Cannot find the mask file.\n", RESET;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-syn"){
		$syn = "\.syn";
	}
	if ($ARGV[$i] eq "\-s_num"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$synnum = $ARGV[$i+1];
		}
		else {
			die "-syn_num value is wrong.\n";
		}
	}
	if ($ARGV[$i] eq "\-mu"){
		if ($ARGV[$i+1] !~ /[^0-9e\-\.]/){
			$mu = $ARGV[$i+1];
		}
		else {
			print BOLD "Invalid mutation rate value, set to default 1e\-8.\n";
		}
	}
	if ($ARGV[$i] eq "\-gen"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$gen = $ARGV[$i+1];
		}
		else {
			print BOLD "Invalid generation value, set to default 3.\n";
		}
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

unless (@vcfs){
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
if (-e "my_bash_msmc_0_$ran\.sh" && $sn != 1){
	$ran = undef;
	goto RE;
}
elsif (-e "my_bash_msmc_1_$ran\.sh" && $sn != 1){
	$ran = undef;
	goto RE;	
}
print "The qsub SN is: $ran\n";

my @chr_n = &chr_name($vcfs[0], $pre); #get chromosome name
#handling mask files
if (-e $mask){
	my $mask_out;
	if ($mask =~ /gff$|gff3$|gff\.gz$|gff3\.gz$/){
		my $out_gff = $mask;
		if ($out_gff !~ /gz$/){
			$out_gff .= ".gz";
		}
		$out_gff =~ s/\.gff\b|\.gff3\b/.sorted.gff/;
		$mask_out = $out_gff;
		$mask_out =~ s/\.gff/.bed/;
		unless (-e $mask_out && $ow == 0){
			$out = "bedtools sort -i $mask \| bgzip -\@ $thread -c \> $out_gff\\n";
			$out .= "bedtools merge -i $out_gff \| bgzip -\@ $thread -c \> $mask_out\\n";
			$out .= "tabix $mask_out\\n";
			&pbs_setting("$exc$proj$bedtool\-cj_quiet -cj_qname msmc_mask -cj_sn $ran -cj_qout . $out");
			$out = "";
		}
	}
	elsif ($mask =~ /bed$|bed\.gz$/){
		$mask_out = $mask;
		$mask_out =~ s/\bbed/sorted.bed/;
		unless (-e $mask_out && $ow == 0){
			my $tmp_bed = $mask_out;
			$tmp_bed =~ s/\bbed/bed.tmp/;
			$out = "bedtools sort -i $mask \| bgzip -\@ $thread -c \> $tmp_bed\\n";
			$out .= "bedtools merge -i $tmp_bed \| bgzip -\@ $thread -c \> $mask_out\\n";
			$out .= "rm $tmp_bed\\n";
			$out .= "tabix $mask_out\\n";
			&pbs_setting("$exc$proj$bedtool\-cj_quiet -cj_ppn $thread -cj_qname msmc_mask -cj_sn $ran -cj_qout . $out");
			$out = "";
		}
	}
	$mask = $mask_out;
	if ($exc){
		&status($ran);	
	}
}
foreach (@chr_n){
	if ($_ eq "0"){
		next;
	}
	my $chr_mask;
	if ($mask){
		$chr_mask = $mask;
		$chr_mask =~ s/\.bed\.gz/\.$_\.bed/;
		my $qout;
		unless (-e "$chr_mask\.gz" && -e "$chr_mask\.gz\.tbi" && $ow == 0){
			unless (-e "$chr_mask\.gz" && $ow == 0){
				$qout .= "zcat $mask \| grep \"$_\\b\" \> $chr_mask\\n";
				$qout .= "bgzip -\@ $thread $chr_mask\\n";
			}
			$qout .= "tabix $chr_mask\.gz\\n";
		}
		if ($qout){
			&pbs_setting("$exc$proj\-cj_quiet -cj_ppn $thread -cj_qname msmc_mask_comp_$_ -cj_sn $ran -cj_qout . $qout");
		}
	}
}
if ($exc){
	&status($ran);
}

my @new_vcfs;
if (scalar(@vcfs) == 1){
	open (BASH, ">my_bash_msmc_1_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_1_$ran\.sh: $!", RESET, "\n";
	print "Preprocessing the vcf file...\n";
	my @sep_vcfs;
	if ($vcfs[0] !~ /gz$/){
		$out = "bgzip -\@ $thread $vcf\\n";
		$out .= "tabix $vcf.gz\\n";
		print BASH "qsub \.\/qsub_files\/$ran\_msmc_1.1\_$cnt\.q\n";
		&pbs_setting("$exc$proj\-cj_quiet -cj_ppn $thread -cj_qname msmc_1.1 -cj_sn $ran -cj_qout . $out");
		$vcfs[0] = "$vcfs[0].gz";
		if ($exc){
			&status($ran);	
		}
	}
	unless (-e "$vcfs[0].tbi"){
		$out = "tabix $vcfs[0]\\n";
		print BASH "qsub \.\/qsub_files\/$ran\_msmc_1.2\_$cnt\.q\n";
		&pbs_setting("$exc$proj\-cj_quiet -cj_qname msmc_1.2 -cj_sn $ran -cj_qout . tabix $vcfs[0]");
		if ($exc){
			&status($ran);	
		}
	}
	my $exist = 0;
	foreach my $ck (0..$#chr_n){
		if ($vcfs[0] =~ /$chr_n[$ck]\b/){
			$exist = 1;
		}
	}
	if ($exist == 0){
		foreach my $n (1..scalar(@chr_n)){
			my $sep_out = $vcfs[0];
			$sep_out =~ s/vcf\.gz$/$chr_n[$n-1].vcf.gz/;
			$out = "";
			unless (-e $sep_out && $ow == 0){
				$out = "bcftools view --threads $thread -r $chr_n[$n-1] -Oz -o $sep_out $vcfs[0]\\n";
			}
			unless (-e "$sep_out.tbi" && $ow == 0){
				$out .= "tabix $sep_out\\n";
			}
			push(@sep_vcfs, $sep_out);
			if ($out){
				print BASH "qsub \.\/qsub_files\/$ran\_msmc_1.3\_$cnt\.q\n";
				&pbs_setting("$exc$proj\-cj_quiet -cj_ppn $thread -cj_qname msmc_1.3\_$n -cj_sn $ran -cj_qout . $out");
			}
		}
		if ($exc){
			&status($ran);	
		}
		@vcfs = @sep_vcfs;
	}
	close(BASH);
}

unless ($rand) {
	$rand = 1;
}

if ($rand >= 1){
	my $tmp_file = $pop_list_f;
	$tmp_file =~ s/txt$|list$/$rand.txt/;
	unless (-e $tmp_file && $ow == 0){
		my $is_syn = 0;
		if ($syn){
			$is_syn = 1;
		}
		system("perl random_sets.pl $pop_list_f $rand $is_syn $synnum");
	}
}

#loops from here
for my $z (1..$rand){
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
	if ($rand >= 1){
		if (-e "$home\/msmc\/population_plot_$mu\_$gen\_$z.pdf" || -e "$o_path\/population_plot_$mu\_$gen\_$z.pdf"){
			if ($ow == 0){
				print "population_plot_$mu\_$gen\_$z.pdf file eixsts. Skip the run.\n";
				next;
			}
		}
		print "Repeated run $z start...\n";
		my $tmp_file = $pop_list_f;
		$tmp_file =~ s/txt$|list$/$z.txt/;
		open(LIST, "<$tmp_file") || die "Cannot open $tmp_file: $!\n";
		@pop_list = <LIST>;
		chomp(@pop_list);
		close(LIST);	
	}

	#supposedly every contig are separated 
	open (BASH, ">my_bash_msmc_2_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_2_$ran\.sh: $!", RESET, "\n";
	print "Generating synthetic F1 vcf and bed file\(s\)...\n";
	$cnt = 1;
	my @in_vcfs; my @in_beds;
	foreach my $vcf (@vcfs){
		my $chr_name;
		foreach (@chr_n){
			$_ =~ s/\s+//g;
			if ($vcf =~ /$_\b/){
				$chr_name = $_;
			}
		}
		foreach my $list (@pop_list){
			my @syn_eles = split(/\t/, $list);
			my $s_name = $syn_eles[0]; #synthetic F1 ID
			my $s_pop;
			if ($#syn_eles == 3){
				$s_pop = $syn_eles[3];
			}
			elsif ($#syn_eles == 1){
				$s_pop = $syn_eles[1];
			}
			else {
				print BOLD "\-list should be a 2-column or 4-column table.\n", RESET;
				exit;			
			}
			if ($syn){
				if ($#syn_eles < 2){
					print BOLD "If you need synthetic F1, you need a 4-column table in the list file.\n", RESET;
					exit;
				}
				$s_name = $syn_eles[2];
			}
			#print "Generating command line for $s_name...\n";
			$out = "";
			if ($syn){
				unless (-e "$o_path\/$s_name\.$s_pop\.msmc\.syn\.$chr_name\.vcf\.gz" && $ow == 0){
					$out .= "perl $p_dir\/synthetic_f1_msmc_single_v4.pl -t $thread -chr $chr_name -i $vcf -syn -o $o_path -list @syn_eles\\n";
				}
			}
			else {
				unless (-e "$o_path\/$s_name\.$s_pop\.msmc\.$chr_name\.vcf\.gz" && $ow == 0){
					$out .= "perl $p_dir\/synthetic_f1_msmc_single_v4.pl -t $thread -chr $chr_name -i $vcf -o $o_path -list @syn_eles\\n";
				}		
			}
			unless (-e "$o_path\/$s_name\.$s_pop\.msmc\.syn\.$chr_name\.sorted\.bed\.gz" && $ow == 0){
				$out .= "bedtools sort -i $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.bed\.gz \| bgzip -\@ $thread -c \> $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.sorted\.tmp\.bed\.gz\\n";
				$out .= "bedtools merge -i $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.sorted\.tmp\.bed\.gz \| bgzip -\@ $thread -c \> $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.sorted\.bed\.gz\\n";
				$out .= "rm $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.sorted\.tmp\.bed\.gz\\n";
			}
			push(@in_vcfs, "$o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.vcf\.gz");
			push(@in_beds, "\-\-mask $o_path\/$s_name\.$s_pop\.msmc$syn\.$chr_name\.sorted\.bed\.gz");
			if ($out){
				print BASH "qsub \.\/qsub_files\/$ran\_msmc_2_$cnt\.q\n";
				&pbs_setting("$exc$proj$bedtool\-cj_quiet -cj_ppn $thread -cj_mem 128 -cj_qname msmc_2_$cnt -cj_sn $ran -cj_qout . $out");
				$out = "";
			}
			$cnt++;
		}
	}
	close(BASH);
	if ($exc){
		&status($ran);
	}

	my $sample; my $pop; my @pops;
	foreach my $list (@pop_list){
		my @tmp_eles = split(/\t/, $list);
		if ($syn) {
			$sample = $tmp_eles[2];
			$pop = $tmp_eles[-1];
			push(@pops, $pop);
		}
		else {
			$sample = $tmp_eles[0];
			$pop = $tmp_eles[-1];		
			push(@pops, $pop);
		}
	}

	my %seen;
	@pops = do { %seen; grep { !$seen{$_}++ } @pops };
	my @pops_sp4 = @pops;
	$cnt = 1;
	open (BASH3, ">my_bash_msmc_3_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_3_$ran\.sh: $!", RESET, "\n";
	print "Generating MSMC input files...\n";
	my @m3_files; my @m3_pop1_num; my @m3_pop2_num; my @m4_pop1; my @m4_pop2;
	do {
		for (my $x=1; $x<=$#pops; $x++){
			my $pop1_alleles; my $pop2_alleles;
			#print "debug: $pops[$x]\n";
			foreach (@chr_n){
				if ($_ eq "0"){
					next;
				}
				#print "debug: chr: $_\n";
				$out = "";
				my @comm_beds1; my @comm_vcfs1; my @comm_beds2; my @comm_vcfs2;
				foreach my $in_bed (@in_beds){
					if ($in_bed =~ /\.$_\./ && $in_bed =~ /\b$pops[0]\b/){
						push(@comm_beds1, $in_bed);
						#print "debug: bed: $in_bed\n";
					}
					if ($in_bed =~ /\.$_\./ && $in_bed =~ /\b$pops[$x]\b/){
						push(@comm_beds2, $in_bed);
						#print "debug: bed: $in_bed\n";
					}				
				}
				foreach my $in_vcf (@in_vcfs){
					if ($in_vcf =~ /\.$_\./ && $in_vcf =~ /\b$pops[0]\b/){
						push(@comm_vcfs1, $in_vcf);
					}
					if ($in_vcf =~ /\.$_\./ && $in_vcf =~ /\b$pops[$x]\b/){
						push(@comm_vcfs2, $in_vcf);
					}
				}
				#print "debug: $mask\n";
				my $chr_mask;
				if ($mask){
					$chr_mask = $mask;
					$chr_mask =~ s/bed\.gz$/$_.bed.gz/;
					$chr_mask = "\-\-negative_mask $chr_mask ";
				}
				#print "debug: $o_path\/$pops[0]\_$pops[$x]\.$_\.msmc_input\.txt\n";
				unless (-e "$o_path\/$pops[0]\_$pops[$x]\.$_\.msmc_input\.txt" && $ow == 0){
					my @comCCs = `find \$pwd \-iname \'generate_multihetsep.py\' \-type f`;
					chomp(@comCCs);
					my $comCC;
					foreach (@comCCs){
						if ($_ =~ /msmc-tools/){
							$comCC = $_;
							if ($comCC =~ /^\.\//){
								$comCC =~ s/^\.\///;
							}
						}
					}
					unless ($comCC){
						print "Cannot find generate_multihetsep.py file.\n";
						exit;
					}
					#print "debug: $chr_mask\n";
					$out .= "python3 $comCC \-\-chr $_ $chr_mask@comm_beds1 @comm_beds2 @comm_vcfs1 @comm_vcfs2 \> $o_path\/$pops[0]\_$pops[$x]\.$_\.msmc_input\.txt\\n";
				}
				$pop1_alleles = ($#comm_vcfs1 + 1) * 2 - 1;
				$pop2_alleles = ($#comm_vcfs2 + 1) * 2 - 1;
				if ($out){
					print BASH3 "qsub \.\/qsub_files\/$ran\_msmc_3\_$cnt\.q\n";
					&pbs_setting("$exc$proj\-cj_quiet -cj_mem 128 -cj_qname msmc_3\_$cnt -cj_sn $ran -cj_qout . $out");
				}
				$cnt++;
			}
			$cnt2++;
			push(@m3_files, "$o_path\/$pops[0]\_$pops[$x]\.\*\.msmc_input\.txt");
			push(@m4_pop1, $pops[0]);
			push(@m4_pop2, $pops[$x]);
			push(@m3_pop1_num, $pop1_alleles);
			push(@m3_pop2_num, $pop2_alleles);
		}
		shift(@pops);
	} until ($#pops == 0);
	close(BASH3);
	if ($exc){
		if (-e "my_bash_msmc_3_$ran\.sh"){
			&status($ran);
		}
	}

	$cnt = 1;
	my @folders;
	open (BASH4, ">my_bash_msmc_4_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_4_$ran\.sh: $!", RESET, "\n";
	print "Generating command line for population size estimation...\n";
	my @sp5_files; my @s_pop_files;
	for (my $i=0; $i<=$#m3_files; $i++){
		my $pop1 = $m4_pop1[$i];
		my $pop2 = $m4_pop2[$i];
		unless (-d "$o_path\/$pop1\_$pop2"){
			system("mkdir $o_path\/$pop1\_$pop2");
		}
		push(@folders, "$o_path\/$pop1\_$pop2");
		$m3_pop2_num[$i] = $m3_pop2_num[$i] + $m3_pop1_num[$i] + 1;
		my $pop2_start = $m3_pop1_num[$i] + 1;
		#$out = "chmod \+x \/home\/hpc\/crlee\/softwares\/msmc2\\n";
		$out = "";
		if (-e "$o_path\/$pop1\_$pop2\/$pop1\_msmc2\.$z.final\.txt" && $ow == 0){}
		else{
			$out .= "msmc2 -t $thread $psa\-o $o_path\/$pop1\_$pop2\/$pop1\_msmc2.$z \-I ";
			for (my $j=0; $j<=$m3_pop1_num[$i]; $j++){
				if ($j==0){
					$out .= "$j";
				}
				else {
					$out .= "\,$j";
				}
			}
			foreach (@chr_n){
				if ($_ eq "0"){
					next;
				}
				my $file = $m3_files[$i];
				$file =~ s/\*/$_/;
				$out .= " $file";
			}
			$out .= "\\n";
		}
		push(@sp5_files, "$o_path\/$pop1\_$pop2\/$pop1\_msmc2\.$z.final\.txt ");
		push(@s_pop_files, "$o_path\/$pop1\_$pop2\/$pop1\_msmc2\.$z.final\.txt");
		if (-e "$o_path\/$pop1\_$pop2\/$pop2\_msmc2\.$z.final\.txt" && $ow == 0){}
		else {
			$out .= "msmc2 -t $thread $psa\-o $o_path\/$pop1\_$pop2\/$pop2\_msmc2.$z \-I ";
			for (my $k=$pop2_start; $k<=$m3_pop2_num[$i]; $k++){
				if ($k==$pop2_start){
					$out .= "$k";
				}
				else {
					$out .= "\,$k";
				}
			}
			foreach (@chr_n){
				if ($_ eq "0"){
					next;
				}
				my $file = $m3_files[$i];
				$file =~ s/\*/$_/;
				$out .= " $file";
			}
			$out .= "\\n";
		}
		$sp5_files[-1] .= "$o_path\/$pop1\_$pop2\/$pop2\_msmc2\.$z.final\.txt ";
		push(@s_pop_files, "$o_path\/$pop1\_$pop2\/$pop2\_msmc2\.$z.final\.txt");
		if (-e "$o_path\/$pop1\_$pop2\/$pop1\_$pop2\_msmc2\.$z.final\.txt" && $ow == 0){}
		else {
			$out .= "msmc2 -t $thread -s $psa\-o $o_path\/$pop1\_$pop2\/$pop1\_$pop2\_msmc2.$z \-I ";
			for (my $l=0; $l<=$m3_pop1_num[$i]; $l++){
				for (my $m=$pop2_start; $m<=$m3_pop2_num[$i]; $m++){
					if ($l == 0 && $m == $pop2_start){
						$out .= "$l\-$m";
					}
					else {
						$out .= "\,$l\-$m";
					}
				}
			}
			foreach (@chr_n){
				if ($_ eq "0"){
					next;
				}
				my $file = $m3_files[$i];
				$file =~ s/\*/$_/;
				$out .= " $file";
			}
			$out .= "\\n";
		}
		$sp5_files[-1] .= "$o_path\/$pop1\_$pop2\/$pop1\_$pop2\_msmc2\.$z.final\.txt";
		if ($out){
			print BASH4 "qsub \.\/qsub_files\/$ran\_msmc_4\_$cnt\.q\n";
			&pbs_setting("$exc$proj\-cj_quiet -cj_qname msmc_4\_$cnt -cj_sn $ran -cj_qout . -cj_ppn $thread -cj_mem 128 $out");
		}
		$cnt++;
	}
	close(BASH4);
	if ($exc){
		if (-e "my_bash_msmc_4_$ran\.sh"){
			&status($ran);
		}
	}

	$cnt = 1;
	open (BASH5, ">my_bash_msmc_5_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_5_$ran\.sh: $!", RESET, "\n";
	print "Combining MSMC population results...\n";

	my $skip;
	for (my $i=0; $i<= $#m4_pop1; $i++){
		$skip = 0;
		my $pop1 = $m4_pop1[$i];
		my $pop2 = $m4_pop2[$i];
		my @three_files = split(/\s/, $sp5_files[$i]);
		foreach my $file (@three_files){
			if ($exc == 1){
				unless (-e $file){
					print BOLD "Cannot find $file, which is required for combining results. Skip.\n", RESET;
					$skip = 1;
				}
			}
		}
		if ($skip == 1){
			next;
		}
		$out = "";
		my @comCCs = `find \$pwd \-iname \'combineCrossCoal.py\' \-type f`;
		chomp(@comCCs);
		my $comCC;
		foreach (@comCCs){
			if ($_ =~ /msmc-tools/){
				$comCC = $_;
				if ($comCC =~ /^\.\//){
					$comCC =~ s/^\.\///;
				}
			}
		}
		unless ($comCC){
			print "Cannot find combineCrossCoal.py file.\n";
			exit;
		}
		unless (-e "$o_path\/$pop1\_$pop2\/$pop1\_$pop2\_msmc2\.$z.combined\.final\.txt" && $ow == 0){
			@three_files = ($three_files[2], $three_files[0], $three_files[1]);
			$out .= "python3 $comCC @three_files \> $o_path\/$pop1\_$pop2\/$pop1\_$pop2\_msmc2\.$z.combined\.final\.txt\\n";
		}
		if ($out){
			print BASH5 "qsub \.\/qsub_files\/$ran\_msmc_5\_$cnt\.q\n";
			&pbs_setting("$exc$env$proj\-cj_quiet -cj_qname msmc_5\_$cnt -cj_sn $ran -cj_qout . $out");
		}
		$cnt ++;
	}
	close(BASH5);
	if ($exc){
		if (-e "my_bash_msmc_5_$ran\.sh"){
			&status($ran);
		}
		if ($debug != 1){
			system("rm $o_path\/*.msmc.*");
			system("rm $o_path\/*.msmc_input.*");
		}
	}
	if ($exc){
		my @server = `ip addr`;
		chomp(@server);
		foreach (@server){
			if ($_ =~ /140.110.148.11/ || $_ =~ /140.110.148.12/){
				my $mv_home = $o_path;
				if ($mv_home =~ /home/){
					last;
				}
				else {
					$mv_home =~ s/work1/home/;
				}
				unless (-d $mv_home){
					system("mkdir $mv_home");
				}
				foreach my $folder (@folders){
					my @tmps;
					if ($folder =~ /\//){
						@tmps = split(/\//, $folder);
					}
					else {
						@tmps = $folder;
					}
					unless (-d "$mv_home\/$tmps[-1]"){
						system("mkdir $mv_home\/$tmps[-1]");
					}
					system("mv $folder\/\*.txt $mv_home\/$tmps[-1]");
					system("mv $o_path\/\*.pdf $mv_home");
					system("rm $o_path\/\*.log");
				}
			}	
		}
	}
	if ($rand > 1){
		print "Repeated run $z is done.\n";
	}
	if ($debug == 1){
		exit;
	}
}
#loops end here

if ($rand > 1){
	my @server = `ip addr`;
	chomp(@server);
	foreach (@server){
    	if ($_ =~ /140.110.148.11/ || $_ =~ /140.110.148.12/){
			unless ($o_path =~ /home/){
				$o_path =~ s/work1/home/;
			}
		}
	}
	my $check_comb = `echo $o_path\/\*all.final.out`;
	unless ($check_comb !~ /\*/ && $ow == 0){
		open (BASH6, ">my_bash_msmc_6_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_6_$ran\.sh: $!", RESET, "\n";
		print "Averaging lambda results...\n";
		$out = "";
		$out = "Rscript $p_dir\/averaging_lambda_new.R $o_path\\n";
		$out .= "Rscript $p_dir\/averaging_lambda_sep_new.R $mu $gen $o_path\\n";
		print BASH6 "qsub \.\/qsub_files\/$ran\_msmc_6\_avg\.q\n";
		&pbs_setting("$exc$env$proj\-cj_quiet -cj_qname msmc_6\_avg -cj_sn $ran -cj_qout . $out");
		if ($exc){
			if (-e "\.\/qsub_files\/$ran\_msmc_6\_avg\.q"){
				&status($ran);
			}
		}
	}
}

foreach (@server){
	if ($_ =~ /140.110.148.11/ || $_ =~ /140.110.148.12/){
		unless ($o_path =~ /home/){
			$o_path =~ s/work1/home/;
		}
	}
}
unless (-e "$o_path\/population_plot_$mu\_$gen\_$ran.pdf" && $ow == 0){
	open (BASH7, ">my_bash_msmc_7_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_7_$ran\.sh: $!", RESET, "\n";
	print "Plotting results...\n";
	$out = "Rscript $p_dir\/population_plot_new.R $mu $gen $ran $o_path\\n";
	print BASH7 "qsub \.\/qsub_files\/$ran\_msmc_7\_plot\.q\n";
	&pbs_setting("$exc$env$proj\-cj_quiet -cj_qname msmc_7\_plot -cj_sn $ran -cj_qout . $out");
}

my $check_sep_plot = `echo $o_path\/\*_rccr_plot_$mu\_$gen\_$ran.pdf`;
unless ($check_sep_plot !~ /\*/ && $ow == 0){
	open (BASH8, ">my_bash_msmc_8_$ran\.sh") || die BOLD "Cannot write my_bash_msmc_8_$ran\.sh: $!", RESET, "\n";
	print "Plotting rCCR results...\n";
	$out = "Rscript $p_dir\/population_plot_sep_new.R $mu $gen $ran $o_path\\n";
	print BASH8 "qsub \.\/qsub_files\/$ran\_msmc_8\_rCCR_plot\.q\n";
	&pbs_setting("$exc$env$proj\-cj_quiet -cj_qname msmc_8\_rCCR_plot -cj_sn $ran -cj_qout . $out");
}
my $end_time = gettimeofday();
my $elapsed_time = $end_time - $start_time;
print "Runtime: $elapsed_time s\n";

sub usage {
	print BOLD "Usage: perl msmc_v5.pl -vcf VCF_FILE -list SYNTHETIC_LIST -mask MASK_FILE [-syn] [-p] [-ow] [-pre PREFIX_OF_THE_CONTIG_NAME] [-mu NUM] [-gen INT] [-random INT] [-rn NUMBER] [-s_num INT] [-t INT] [-debug] [-proj PROJECT_NAME] [-sn] [-exc]\n", RESET;
	return 1;
} #print usage
