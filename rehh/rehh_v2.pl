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
my $r_env = '$HOME/miniconda3/envs/R-4.1/bin';

chomp(@ARGV);
print "The script is written by Ben Chien. May. 2023.\n";
print "Input command line:\n";
print "perl rehh_v2.pl @ARGV\n";

if ($#ARGV == -1){
	&usage;
	exit;
}

my $p_dir = $FindBin::Bin;
my $exc; my $ran; my $sn; my $ow = 0; my $local; my $pre;
my @target_chr; my @pos; my @vcfs; my $path; my $ebp = 200000; my $aid;
my $list; my $win = 10000; my $keep; my $maf = 0.05; my $bi; my $syn = 0; my $lim_ehh = 0.01;
my $o_path; my $method = "ehh"; my @popis; my @popis2; my $popi;
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
	if ($ARGV[$i] eq "\-ow"){
        $ow = 1;
	}
	if ($ARGV[$i] eq "\-local"){
        $local = "-cj_local ";
	}
	if ($ARGV[$i] eq "\-ka"){
        $keep = 1;
	}
	if ($ARGV[$i] eq "-pre"){
		$pre = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-vcf"){
        if (-d $ARGV[$i+1]){
            if ($ARGV[$i+1] =~ /\/$/){
                $ARGV[$i+1] =~ s/\/$//;
                $path = $ARGV[$i+1];
            }
            $path = &check_path($path);
            @vcfs = <$path\/*.vcf>;
            unless (@vcfs){
                @vcfs = <$path\/*.vcf.gz>;
                unless (@vcfs){
                    die "Cannot find vcf\(s\).\n";
                }
            }
            @vcfs = grep(!/anc\.vcf\.gz$|filtered\.vcf\.gz$|region\.vcf\.gz$/, @vcfs);
        }
        elsif (-e $ARGV[$i+1] && $ARGV[$i+1] =~ /txt$|list$/){
            open(LIST, "<$ARGV[$i+1]") || die "Cannot open $ARGV[$i+1]: $!\n";
            @vcfs = <LIST>;
            chomp(@vcfs);
            foreach my $i (0..$#vcfs){
                $vcfs[$i] =~ s/[\x0A\x0D]//g;
            }
            if ($vcfs[0] =~ /\//){
                my @tmp = split(/\//, $vcfs[0]);
                pop(@tmp);
                $path = join("\/", @tmp);
            }
            else {
                $path = ".";
            }
			$path = &check_path($path);
        }
        elsif (-e $ARGV[$i+1] && $ARGV[$i+1] =~ /vcf$|vcf\.gz$/){
        	@vcfs = $ARGV[$i+1];
            if ($vcfs[0] =~ /\//){
                my @tmp = split(/\//, $vcfs[0]);
                pop(@tmp);
                $path = join("\/", @tmp);
            }
            else {
                $path = ".";
            }
            $path = &check_path($path);        	
        }
        else {
            die "Cannot find the vcf\(s\). The input must be a folder contains seperated vcfs or a list file.\n";
        }
        $inprefix = "Relate_output";
    }
	if ($ARGV[$i] eq "\-bp"){
		if ($ARGV[$i+1] !~ /\:/){
			die "-bp: Chromosome\/contig name is required.\n";
		}
		my @bps;
		if ($ARGV[$i+1] =~ /\,/){
			@bps = split(/\,/, $ARGV[$i+1]);
		}
		else {
			@bps = $ARGV[$i+1];
		}
		foreach (@bps){
			my @eles = split(/\:|\-/, $_);
			if (scalar(@eles) != 2){
				die "-bp: Format is wrong. Please check.\n";
			}
			push(@target_chr, $eles[0]);
			if ($eles[1] !~ /[^0-9]/){
				push(@pos, $eles[1]);
			}
			else {
				die "-bp: Format is wrong. Please check.\n";
			}
		}
		undef(@bps);
	}
	if ($ARGV[$i] eq "\-ebp"){
		if ($ebp !~ /[^0-9]/){
			$ebp = $ARGV[$i+1];
		}
		else {
			die "-ebp: It should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-chr"){
		push(@target_chr, $ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "\-aid"){
		if (-e $ARGV[$i+1]){
			$aid = $ARGV[$i+1];
		}
		elsif ($ARGV[$i+1] =~ /\w/) {
			$aid = $ARGV[$i+1];
		}
		else {
			die "-aid: Cannot find the aid file.\n";
		}
	}
	if ($ARGV[$i] eq "\-list"){
		if (-e $ARGV[$i+1]){
			$list = $ARGV[$i+1];
		}
		else {
			die "-list: Cannot find the list file.\n";
		}
	}
	if ($ARGV[$i] eq "\-w"){ #window size for rehh
		if ($win !~ /[^0-9]/){
			$win = $ARGV[$i+1];
		}
		else {
			die "-w: It should be an integer.\n";
		}
	}
	if ($ARGV[$i] eq "\-bi"){
		$bi = 1;
	}
	if ($ARGV[$i] eq "\-syn"){
		$syn = 1;
	}
	if ($ARGV[$i] eq "\-maf"){
		if ($maf !~ /[^0-9\.]/){
			$maf = $ARGV[$i+1];
		}
		else {
			die "-maf: It should be a float number.\n";
		}
	}
	if ($ARGV[$i] eq "\-l"){
		if ($lim_ehh !~ /[^0-9\.]/){
			$lim_ehh = $ARGV[$i+1];
		}
		else {
			die "-l: It should be a float number.\n";
		}
	}
	if ($ARGV[$i] eq "\-o"){
        if (-d $ARGV[$i+1]){
			$o_path = $ARGV[$i+1];
            if ($o_path =~ /\/$/){
                $o_path =~ s/\/$//;
            }
            $o_path = &check_path($o_path);
		}
		else {
			my $return = `mkdir $o_path`;
			if ($return =~ /\w/){
				die "-o: Cannot create the output path.\n";
			}
		}
	}
	if ($ARGV[$i] eq "\-m"){
		$method = $ARGV[$i+1];
		if ($method !~ /ehh|ehhs|xpehh|ihs/i){
			die "-m: the possible values are: ehh, ehhs, xpehh, ihs. default: ehh\n";
		}
		else {
			lc $method;
		}
	}
	if ($ARGV[$i] eq "\-popi"){ #interest pop(s)
		$popi = $ARGV[$i+1];
	}
}

unless (-d "qsub_files"){
	system ("mkdir qsub_files");
}
unless (-d "qsub_files\/out"){
	system ("mkdir qsub_files\/out");
}

#RE:
unless ($ran){
	$ran = &rnd_str(4, "A".."Z", 0..9);
}

unless($o_path){
	$o_path = $path;
}
unless (-d $o_path){
	system("mkdir $o_path");
}

my @chrs = &chr_name($vcfs[0], $pre);
my @chr_len = &chr_lengths($vcfs[0], $pre);

if ($method eq "ihs" || $method eq "xpehh"){
	if (scalar(@pos) > 0){
		die "-bp: This option can only be used by ehh and ehhs analyses.\n";
	}
	unless (@target_chr){
		@target_chr = @chrs;
	}
}

#handling interested populations
if ($popi =~ /\,/){
	if ($method ne "xpehh"){
		if ($popi =~ /\(|\)/){
			die "-popi: \(\) can only be used for xpehh analysis.\n";
		}
		unless ($popi =~ /\[|\]/){ #all single
			my @tmp_popis = split(/\,/, $popi);
			push(@popis, @tmp_popis);
		}
		else { # [#,#,#],#,[#,#,#]; [#,#,#],[#,#,#],#; [#,#,#],#,#......
			my @tmp_popis = split(/\[|\]/, $popi);
			foreach (@tmp_popis){
				if ($_ eq ',' || $_ == ''){
					next;
				}
				$_ =~ s/^\,|\,$//g;
				push(@popis, $popi);
			}
		}
	}
	else { #XP-EHH
		my @tmp_pairs = split(/\(|\)/, $popi);
		print "debug1: ", join(" ", @tmp_pairs), "\n";
		foreach my $i (0..$#tmp_pairs){
			if ($tmp_pairs[$i] eq ','){
				next;
			}
			$tmp_pairs[$i] =~ s/^\,|\,$//g;
			unless ($tmp_pairs[$i] =~ /\[|\]/){ #single vs single
				my @tmp_popis = split(/\,/, $tmp_pairs[$i]);
				push(@popis, $tmp_popis[0]);
				push(@popis2, $tmp_popis[1]);
			}
			else { # #,[#,#,#]; [#,#,#],#; [#,#,#],[#,#,#]   
				my @multi_pops_tmp = split(/\[|\]/, $tmp_pairs[$i]);
				@multi_pops_tmp = grep {defined($_) && $_ ne ''} @multi_pops_tmp; #remove empty elements
				if (scalar(@multi_pops_tmp) == 3){ #remove comma
					splice @multi_pops_tmp, 1, 1;
				}
				foreach my $j (0..1){
					$multi_pops_tmp[$j] =~ s/^\,|\,$//g;
				}
				push(@popis, $multi_pops_tmp[0]);
				push(@popis2, $multi_pops_tmp[1]);
			}
		}
		print "debug: ", join(" ", @popis), "\n";
		print "debug: ", join(" ", @popis2), "\n";
	}
}
else {
	if ($method ne "xpehh"){
		$ARGV[$i+1] =~ s/\[|\]//g;
		@popis = "$ARGV[$i+1]";
	}
	else {
		die "-popi: cannot find pair(s) for XP-EHH analysis.\n";
	}
}

my $syn_list;
if ($syn == 1){ #make a synthetic F1 list
	my @aid_list; my @vcf_ids;
	if (-e $aid){
		open(AID, "<$aid") || die "Cannot open $aid: $!\n";
		@aid_list = <AID>;
		chomp(@aid_list);
		close(AID);
	}
	else {
		@aid_list = split(/\,/, $aid);
		chomp(@aid_list);
	}
	my $vcf_id;
	if ($vcfs[0] =~ /gz$/){
		$vcf_id = `zcat $vcfs[0] \| head -n 5000 \| grep \"CHROM\"`;
	}
	else {
		$vcf_id = `cat $vcfs[0] \| head -n 5000 \| grep \"CHROM\"`;
	}
	@vcf_ids = split(/\t+|\s+/, $vcf_id);
	@vcf_ids = splice @vcf_ids, 9;
	unless ($keep == 1){ #remove ancestral ID from the syn_F1 list if -ka is undefined
		my @new_vcf_ids;
		foreach my $i (0..$#vcf_ids){
			my $exist = 0;
			foreach my $j (0..$#aid_list){
				if ($vcf_ids[$i] eq $aid_list[$j]){
					$exist = 1;
				}
			}
			if ($exist == 0){
				push(@new_vcf_ids, $vcf_ids[$i]);
			}
		}
		@vcf_ids = @new_vcf_ids;
		undef(@new_vcf_ids);
	}
	my @lists; #final synthetic F1 list
	if ($list){ #keep only listed samples
		my @new_vcf_ids; my @new_lists;
		open(LIST, "<$list") || die "Cannot open $list: $!\n";
		@lists = <LIST>;
		chomp(@lists);
		close(LIST);
		foreach my $i (0..$#vcf_ids){
			foreach (@lists){
				if ($_ =~ /\b$vcf_ids[$i]\b/){
					push(@new_vcf_ids, $vcf_ids[$i]);
					push(@new_lists, $_);
				}
			}
		}
		@vcf_ids = @new_vcf_ids;
		undef(@new_vcf_ids);
		@lists = @new_lists;
		undef(@new_lists);
	}
	unless (@lists){
		@lists = @vcf_ids;
	}
	open (SYNLIST, ">$o_path\/$ran\_synthetic_f1.list") || die BOLD "Cannot write $ran\_synthetic_f1.list: $!", RESET, "\n";
	my @ck_list = split(/\t+|\s+/, $lists[0]);
	if (scalar(@ck_lists) == 1){ #no population info
		if (scalar(@lists)/2 > int(scalar(@lists)/2)){
			pop(@lists);
		}
		for (my $i=0; $i<=$#lists; $i+=2){
			my $syn_cnt = $i+1;
			print SYNLIST "$lists[$i]\t$lists[$i+1]\tsyn_$syn_cnt\n";
		}
	}
	else { #with population info
		my %pop_syn; my $syn_cnt = 1;
		foreach (@lists){
			my @line = split(/\t+|\s+/, $_);
			if (@{$pop_syn{$line[1]}}){
				push(@{$pop_syn{$line[1]}}, $line[0]);
			}
			else {
				@{$pop_syn{$line[1]}} = $line[0];
			}
		}
		foreach my $key (keys %pop_syn){
			my @pop_list = @{$pop_syn{$key}};
			if (scalar(@pop_list)/2 > int(scalar(@pop_list)/2)){
				pop(@pop_list);
			}
			for (my $i=0; $i<=$#pop_list; $i+=2){
				print SYNLIST "$pop_list[$i]\t$pop_list[$i+1]\tsyn_$syn_cnt\t$key\n";
				$syn_cnt++;
			}			
		}
	}
	close(SYNLIST);
	$syn_list = "$o_path\/$ran\_synthetic_f1.list";
}

foreach my $i (0..$#vcfs){
	$vcfs[$i] = &check_path($vcfs[$i]);
}

my $out;
foreach my $i (0..$#vcfs){
	unless ($vcfs[$i] =~ /gz$/){
		$out .= "bgzip -c $vcfs[$i]\\n";
		$vcfs[$i] = "vcfs[$i].gz";
	}
	unless (-e "$vcfs[$i].tbi"){
		$out .= "tabix $vcfs[$i]\\n";
	}
	&pbs_setting("$exc$local\-cj_quiet -cj_qname tabix_$i -cj_sn $ran -cj_qout . $out");
}
if ($exc =~ /-cj_exc/){
	&status($ran);
}

unless (-e $aid){
	open(AID, ">$o_path\/$ran\_aid_list.txt") || die "Cannot write aid list: $!\n";
	my @aid_tmp = split(/\,/, $aid);
	foreach (@aid_tmp){
		print AID "$_\n";
	}
	close(AID);
	$aid = "$o_path\/$ran\_aid_list.txt";
}

#for target region scan (ehh) or ehhs
if ($method eq "ehh" || $method eq "ehhs"){
	my $max_len;
	foreach my $i (0..$#pos){
		foreach my $j (0..$#chrs){
			if ($target_chr[$i] eq $chrs[$j]){
				$max_len = $chr_len[$j];
			}
		}
		my $min_exp = $pos[$i] - $ebp;
		my $max_exp = $pos[$i] + $ebp;
		if ($min_exp < 1){
			$min_exp = 1;
		}
		if ($max_exp > $max_len){
			$max_exp = $max_len;
		}
		$out = "";
		if (scalar(@vcfs) == 1){
			unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.region.vcf.gz" && $ow == 0){
				$out = "tabix $vcfs[0] $target_chr[$i]\:$min_exp\-$max_exp -h \| bgzip -c \> $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.region.vcf.gz\\n";
			}
		}
		else {
			foreach my $k (0..$#vcfs){
				my $return_ck = `tabix $vcfs[$k] $target_chr[$i] \| head -n 1`;
				if ($return_ck =~ /\w/){
					unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.region.vcf.gz" && $ow == 0){
						$out = "tabix $vcfs[$k] $target_chr[$i]\:$min_exp\-$max_exp -h \| bgzip -c \> $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.region.vcf.gz\\n";
						last;
					}
				}
			}
		}
		unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.filtered.vcf.gz" && $ow == 0){
			if ($bi == 1){
				$bi = "--min-alleles 2 --max-alleles 2 ";
			}
			$out .= "vcftools --gzvcf $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.region.vcf.gz --remove-indels --recode --recode-INFO-all $bi\--stdout \| bgzip -c \> $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.filtered.vcf.gz\\n";
		}
		if ($syn == 0){
			my $a_list; my $a_keep;
			if ($list){
				$a_list = "-list $list ";
			}
			if ($keep){
				$a_keep = "-keep ";
			}
			unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.thap" && $ow == 0){
				$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.filtered.vcf.gz $a_keep$a_list\-aid $aid -thap -o $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc\\n";
			}
		}
		else {
			unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn.thap" && $ow == 0){
				unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.vcf.gz" && $ow == 0){
					$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.filtered.vcf.gz -keep -aid $aid -o $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc\\n";
				}
				unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn.vcf.gz" && $ow == 0){
					$out .= "perl synthetic_f1_rehh.pl $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.vcf.gz $syn_list\\n";
				}
				unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn.thap" && $ow == 0){
					$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn.vcf.gz -thap -o $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn\\n";
					$out .= "rm $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.syn.vcf.gz\\n";
					$out .= "rm $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc.vcf.gz\\n";
				}
			}
		}
		my $plot_syn;
		if ($syn == 1){
			$plot_syn = "\.syn";
		}
		if ($method eq "ehh"){
			unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.$lim_ehh.$pos[$i].tiff" && $ow == 0){
				my $maf_out;
				if ($maf){
					$maf_out = "-maf $maf";
				}
				$out .= "Rscript rehh_calc.R -m ehh -thap $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.thap -map $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.map -chr $target_chr[$i] -min $min_exp -max $max_exp -pos $pos[$i] -l $lim_ehh $maf_out\\n";
			}
			&pbs_setting("$exc$local\-cj_quiet -cj_env $r_env -cj_qname region_ehh_$i -cj_sn $ran -cj_qout . $out");
		}
		else {
			my @popis_ehhs;
			unless (@popis){
				@popis_ehhs = 'all';
			}
			else {
				@popis_ehhs = @popis;
				
			}
			my $maf_out;
			if ($maf){
				$maf_out = "-maf $maf";
			}
			my $ehhs_list; my $ehhs_syn_on;
			if ($syn_list){
				$ehhs_list = "-pinfo $syn_list ";
				$ehhs_syn_on = "-syn ";
			}
			elsif ($list){
				$ehhs_list = "-pinfo $list ";
			}
			foreach my $j (0..$#popis_ehhs){
				my $popis_ehhs_name = $popis_ehhs[$j];
				$popis_ehhs_name =~ s/\,/_/g;
				unless (-e "$o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.$lim_ehh.$pos[$i].$popis_ehhs_name.tiff" && $ow == 0){				
					$out .= "Rscript rehh_calc.R -m ehhs -thap $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.thap -map $o_path\/$ran\_$target_chr[$i]\_$min_exp\_$max_exp.anc$plot_syn.map -chr $target_chr[$i] -min $min_exp -max $max_exp -pos $pos[$i] -l $lim_ehh -popi $popis_ehhs[$j] $ehhs_list$ehhs_syn_on$maf_out\\n";
				}
			}
			&pbs_setting("$exc$local\-cj_quiet -cj_env $r_env -cj_qname region_ehhs_$i -cj_sn $ran -cj_qout . $out");
		}
	}
	if ($exc =~ /-cj_exc/){
		&status($ran);
	}
}
#for target genome scan xpehh, ihs
if ($method eq "xpehh" || $method eq "ihs"){
	my $max_len;
	foreach my $i (0..$#target_chr){
		foreach my $j (0..$#chrs){
			if ($target_chr[$i] eq $chrs[$j]){
				$max_len = $chr_len[$j];
			}
		}
		$out = "";
		if (scalar(@vcfs) == 1){
			unless (-e "$o_path\/$ran\_$target_chr[$i].vcf.gz" && $ow == 0){
				$out = "tabix $vcfs[0] $target_chr[$i] -h \| bgzip -c \> $o_path\/$ran\_$target_chr[$i].vcf.gz\\n";
			}
		}
		else {
			foreach my $k (0..$#vcfs){
				my $return_ck = `tabix $vcfs[$k] $target_chr[$i] \| head -n 1`;
				if ($return_ck =~ /\w/){
					unless (-e "$o_path\/$ran\_$target_chr[$i].vcf.gz" && $ow == 0){
						$out = "tabix $vcfs[$k] $target_chr[$i] -h \| bgzip -c \> $o_path\/$ran\_$target_chr[$i].vcf.gz\\n";
						last;
					}
				}
			}
		}
		unless (-e "$o_path\/$ran\_$target_chr[$i].filtered.vcf.gz" && $ow == 0){
			if ($bi == 1){
				$bi = "--min-alleles 2 --max-alleles 2 ";
			}
			$out .= "vcftools --gzvcf $o_path\/$ran\_$target_chr[$i].vcf.gz --remove-indels --recode --recode-INFO-all $bi\--stdout \| bgzip -c \> $o_path\/$ran\_$target_chr[$i].filtered.vcf.gz\\n";
		}
		if ($syn == 0){
			my $a_list; my $a_keep;
			if ($list){
				$a_list = "-list $list ";
			}
			if ($keep){
				$a_keep = "-keep ";
			}
			unless (-e "$o_path\/$ran\_$target_chr[$i].anc.thap" && $ow == 0){
				$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i].filtered.vcf.gz $a_keep$a_list\-aid $aid -thap -o $o_path\/$ran\_$target_chr[$i].anc\\n";
			}
		}
		else {
			unless (-e "$o_path\/$ran\_$target_chr[$i].anc.syn.thap" && $ow == 0){
				unless (-e "$o_path\/$ran\_$target_chr[$i].anc.vcf.gz" && $ow == 0){
					$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i].filtered.vcf.gz -keep -aid $aid -o $o_path\/$ran\_$target_chr[$i].anc\\n";
				}
				unless (-e "$o_path\/$ran\_$target_chr[$i].anc.syn.vcf.gz" && $ow == 0){
					$out .= "perl synthetic_f1_rehh.pl $o_path\/$ran\_$target_chr[$i].anc.vcf.gz $syn_list\\n";
				}
				unless (-e "$o_path\/$ran\_$target_chr[$i].anc.syn.thap" && $ow == 0){
					$out .= "perl vcf2anc_vcf.pl -vcf $o_path\/$ran\_$target_chr[$i].anc.syn.vcf.gz -thap -o $o_path\/$ran\_$target_chr[$i].anc.syn\\n";
					$out .= "rm $o_path\/$ran\_$target_chr[$i].anc.syn.vcf.gz\\n";
					$out .= "rm $o_path\/$ran\_$target_chr[$i].anc.vcf.gz\\n";
				}
			}
		}
		my $plot_syn;
		if ($syn == 1){
			$plot_syn = "\.syn";
		}
		if ($method eq "ihs"){
			unless (-e "$o_path\/$ran\_$target_chr[$i].anc$plot_syn.$lim_ehh.ihs.tiff" && $ow == 0){
				my $maf_out;
				if ($maf){
					$maf_out = "-maf $maf";
				}
				$out .= "Rscript rehh_calc.R -m ihs -w $win -thap $o_path\/$ran\_$target_chr[$i].anc$plot_syn.thap -map $o_path\/$ran\_$target_chr[$i].anc$plot_syn.map -chr $target_chr[$i] -l $lim_ehh $maf_out\\n";
			}
			&pbs_setting("$exc$local\-cj_quiet -cj_env $r_env -cj_qname ihs_$i -cj_sn $ran -cj_qout . $out");
		}
		else { #XP-EHH
			my @popis_ehhs; my @popis_ehhs2;
			unless (@popis){
				@popis_ehhs = 'all';
				@popis_ehhs2 = 'all';
			}
			else {
				@popis_ehhs = @popis;
				@popis_ehhs2 = @popis2;
			}
			if (scalar(@popis_ehhs) != scalar(@popis_ehhs2)){
				die "XP-EHH: the population pair(s) seems impaired.\n";
			}
			my $maf_out;
			if ($maf){
				$maf_out = "-maf $maf";
			}
			my $ehhs_list; my $ehhs_syn_on;
			if ($syn_list){
				$ehhs_list = "-pinfo $syn_list ";
				$ehhs_syn_on = "-syn ";
			}
			elsif ($list){
				$ehhs_list = "-pinfo $list ";
			}
			foreach my $j (0..$#popis_ehhs){ #for each pair
				if ($popis_ehhs[$j] !~ /[0-9a-z]/i){
					next;
				}
				my $popis_ehhs_name = $popis_ehhs[$j];
				$popis_ehhs_name =~ s/\,/_/g;
				my $popis_ehhs_name2 = $popis_ehhs2[$j];
				$popis_ehhs_name2 =~ s/\,/_/g;
				unless (-e "$o_path\/$ran\_$target_chr[$i].anc$plot_syn.$lim_ehh.$popis_ehhs_name\_$popis_ehhs_name2\_xpehh.rda" && $ow == 0){				
					$out .= "Rscript rehh_calc.R -m xpehh -w $win -thap $o_path\/$ran\_$target_chr[$i].anc$plot_syn.thap -map $o_path\/$ran\_$target_chr[$i].anc$plot_syn.map -chr $target_chr[$i] -l $lim_ehh -popi $popis_ehhs[$j] -popi2 $popis_ehhs2[$j] $ehhs_list$ehhs_syn_on$maf_out\\n";
				}
			}
			&pbs_setting("$exc$local\-cj_quiet -cj_env $r_env -cj_qname xpehh_$i -cj_sn $ran -cj_qout . $out");
		}
	}
	if ($exc =~ /-cj_exc/){
		&status($ran);
	}
}
print "The qsub SN is: $ran\n";

sub usage {
	print BOLD "perl rehh_v2.pl -vcf VCF_FILE\/FOLDER -aid SAMPLE_LIST_FILE [-bp CHR\:POS][-chr CHR][-m METHOD][-syn][-list SAMPLE_LIST_FILE][-popi POP_NAMES][-ebp EXTEND_BP][-ka][-w WINDOW_SIZE][-maf FLOAT][-bi][-pre PREFIX][-o OUTPUT_PATH][-sn SN][-ow][-exc][-local][-h]\n", RESET;
	print "-bp: target positions. Multiple positions can be seperated by comma \(no space\). Only required for ehh and ehhs\n";
	print "-chr: assign target chromosome to analyze. The option only works for iHS and XP-EHH\n";
	print "-aid: the ancestral ID list file. One sample ID per line.\n";
	print "-m: method to use. possible value: ehh,ehhs,xpehh,ihs. default: ehh\n";
	print "-syn: use synthetic F1.\n";
	print "-list: the sample IDs of interests. One sample ID per line.\n";
	print "-ka: use to keep ancestral samples, otherwise, the ancestral samples will be removed from the analysis.\n";
	print "-ebp: extend size by bp. default: 200000\n";
	print "-maf: MAF filtering. default: 0.05\n";
	print "-pre: prefix of the chromosome name that you want do the analysis.\n";
	print "-bi: only keep bi-allele sites.\n";
	print "-l: a value of limehh in rehh package. default: 0.01\n";
	print "-o: output folder path\n";
	print "-popi: multiple populations as an group could be indicated by \[population_1,population_2\]. Multiple independent runs can be indicated by comma as population_1,population2.\nYou can combine these two functions as [population_1,population_2],population_3\n";
	print "\tThe first run will be population_1\+population_2, and the second run will be population_3.\n";
	print "-w: window size for calculation \(bp\). default: 10000 \(only for iHS, XP-EHH\)\n\n";
}
