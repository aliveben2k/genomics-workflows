#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
#use Time::HiRes qw(gettimeofday);
my $time = scalar localtime();

chomp(@ARGV);
if ($#ARGV == -1){
	exit;
}

my $thread; my $syn = 0; my $input; my $out_path; my $syn_list_s1; my $syn_list_s2; my $syn_list_name; my $syn_list_pop;
for (my $i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] eq "\-t"){ #thread
        $thread = $ARGV[$i+1];
		if ($thread =~ /[^0-9]/){
			print BOLD "-t only accept 0-9.\n", RESET;
			exit;
		}        
    }
    if ($ARGV[$i] eq "\-chr"){ #chromosome/contig name
        $chr_n = $ARGV[$i+1];      
    }
    if ($ARGV[$i] eq "\-i"){ #input vcf
        $input = $ARGV[$i+1];
        unless (-e $input){
            die "Cannot find $input.\n";
        }
    }
    if ($ARGV[$i] eq "\-syn"){ #use synthetic F1 or not
        $syn = 1;      
    }
    if ($ARGV[$i] eq "\-o"){ #output path
        if (-d $ARGV[$i+1]){
            $out_path = $ARGV[$i+1];
            if ($out_path =~ /\/$/){
                $out_path =~ s/\/$//;
            }
        }
        else {
            $out_path = $ARGV[$i+1];
            die "Cannot find $out_path.\n";
        }
    }
    if ($ARGV[$i] eq "\-list"){ #sample info list
        $syn_list_s1 = $ARGV[$i+1];
        $syn_list_s2 = $ARGV[$i+2];
        $syn_list_name = $ARGV[$i+3];
        $syn_list_pop = $ARGV[$i+4];
        if ($syn_list_name =~ /^\-/ || $syn_list_pop eq ""){
            $syn_list_pop = $ARGV[$i+2];
            $syn_list_s2 = "";
            $syn_list_name = "";
        }
    }
}

#print "debug:\ns1:$syn_list_s1\ns2:$syn_list_s2\nlist_name:$syn_list_name\npop:$syn_list_pop\n";

=skip
my $thread = $ARGV[0];
my $chr_n = $ARGV[1];
my $input = $ARGV[2];
my $syn = $ARGV[3];
my $out_path = $ARGV[4];
my $syn_list_s1 = $ARGV[5];
my $syn_list_s2 = $ARGV[6];
my $syn_list_name = $ARGV[7];
my $syn_list_pop = $ARGV[8];

if ($#ARGV == 6){
	$syn_list_pop = $ARGV[6];
}
=cut

if ($syn == 1){
	$syn = "\.syn";
}
else {
	$syn = "";
}

if ($syn =~ /syn/){}
else {
	$syn_list_s2 = $syn_list_s1;
	$syn_list_name = $syn_list_s1;
}
my $list = "$syn_list_s1\t$syn_list_s2\t$syn_list_name\t$syn_list_pop";

my $out; my $bed;
my @eles = split(/\t/, $list);
$out = "$out_path\/$eles[-2]\.$eles[-1]\.msmc$syn.$chr_n\.vcf.gz";
$bed = "$out_path\/$eles[-2]\.$eles[-1]\.msmc$syn.$chr_n\.bed.gz";
if (-e "$out_path\/$eles[-2]\.$eles[-1]\.msmc$syn.$chr_n\.vcf\.gz"){
	system("rm $out_path\/$eles[-2]\.$eles[-1]\.msmc$syn.$chr_n\.vcf\.gz");
}
if (-e "$out_path\/$eles[-2]\.$eles[-1]\.msmc$syn$chr_n\.bed\.gz"){
	system("rm $out_path\/$eles[-2]\.$eles[-1]\.msmc$syn.$chr_n\.bed\.gz");
}
my @chrs = &chr_name($input, $chr_n);
my @chr_lengss = &chr_lengths($input, $chr_n);
my $chr = $chrs[0];
my $chr_lengs = $chr_lengss[0];
#debug--
#$chr_lengs = 13857329;
#debug--
my $division; my $results_v; my $results_b;
if ($thread > 1){
	if ($input =~ /\.gz$/){
   		unless (-e "$input.tbi"){
   			system("tabix $input");
   		}
	}
	else {
   		system("bgzip -\@ $thread $input");
   		system("tabix $input.gz");
   		$vcf = "$input.gz";
	}
   	$division = $chr_lengs / $thread;
   	if ($division - int($division) > 0.5){
   		$thread = $thread - 1;
   	}
   	my @thr = (); my @joinable;
   	#my $start_time = gettimeofday();
   	my $cnt = 0;
   	do {
   		foreach my $d (1..$thread){
   			my $region_start; my $region_end;
   			if ($d < $thread){
   				if ($d == 1){
   					$region_start = 1;
   				}
   				else {
   					$region_start = int($chr_lengs/$thread*($d-1))+1;
   				}
   				$region_end = int($chr_lengs/$thread*($d));
   			}
   			else {
   				$region_start = int($chr_lengs/$thread*($d-1))+1;
   				$region_end = $chr_lengs;
   			}
   			$thr[$d] = threads->create({'context' => 'list'},'f1', ($input,$chr_lengs,$chr,$region_start,$region_end,$list));
   		}
   		my @running;
   		do {
   			sleep(30);
   			@running = threads->list(threads::running);
   			#my $end_time = gettimeofday();
   			#my $elapsed_time = $end_time - $start_time;
   			#if ($elapsed_time >= 3600){ # if the thread has run for an hour, detach all results.
   			#	foreach my $d (1..$thread){
   			#		$thr[$d]->detach();
   			#	}
   			#}
   		} until (scalar(@running) == 0);
   		@joinable = threads->list(threads::joinable);
   		if (scalar(@joinable) > 0){
   			foreach my $d (1..$thread){
   				my @results = $thr[$d]->join();
   				$results_v .= $results[0];
				$results_b .= $results[1];
   			}
   		}
   		undef(@thr);
   		$cnt++;
   	} until ($results_v =~ /[a-z]/i || $cnt == 3);
}
else {
	my $region_start = 0; my $region_end = 0;
	my @results = &f1($input,$chr_lengs,$chr,$region_start,$region_end,$list);
   	$results_v .= $results[0];
	$results_b .= $results[1];	
}


open(OUTVCF, '|-', "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
print OUTVCF $results_v;
close(OUTVCF);
open(OUTBED, '|-', "bgzip \> $bed") || die BOLD "Cannot write $bed: $!", RESET, "\n";
print OUTBED $results_b;
close(OUTBED);


print "Synthetic_F1 done.\n";

sub f1 {
my $input = shift; my $chr_lengs = shift; my $chr = shift; my $region_start = shift; my $region_end = shift; my $list = shift;
my @sample_line; my @samples; my $pop_sample = 0; my @pos;
my @last_mask;
if ($input =~ /\.gz$/){
   	unless (-e "$input.tbi"){
   		system("tabix $input");
   	}
   	if ($region_start != 0){
       	open(INPUT, "-|", "tabix -h $input $chr\:$region_start\-$region_end") || die BOLD "Cannot open $input: $!", RESET, "\n";
   	}
   	else {
   		open(INPUT, "-|", "bgzip -dc $input") || die BOLD "Cannot open $input: $!", RESET, "\n";
   	}
}
else {
    die BOLD "$input must be compressed first.", RESET, "\n";
}

my $first_line = 0; my @vcf_outs;
my $bed_out; my $vcf_out; my $cnt = 1; my $record = 2;
while (my $line = <INPUT>){
	$fixed = 0;
	chomp $line;
	my $fst_idx = 0;
	my $snd_idx = 1;
		if ($line =~ /^#/){
			if ($line =~ /^#CHROM/){
				@sample_line = split(/\t/, $line);
				for (my $i=0; $i<=8; $i++){
					if ($region_start == 1){
						push(@vcf_outs, "$sample_line[$i]\t");
						#print OUTVCF "$sample_line[$i]\t";
					}
				}
				for (my $x=0; $x<=$#sample_line; $x++){		
					if ($list =~ /$sample_line[$x]/){
						push(@pos, $x);
					}
				}
				my @list_tmps = split(/\t/, $list);
				if ($region_start == 1){
					push(@vcf_outs, "$list_tmps[2]\n");
					#print OUTVCF "$list_tmps[2]\n";
				}
			}
			else {
				if ($region_start == 1){
					push(@vcf_outs, "$line\n");
					#print OUTVCF "$line\n";
				}
			}
		}
		else {
			my @line_vec = split(/\t/, $line);
			my $first_8;
			for (my $j=0; $j<=8; $j++){
				if ($j <= 7){
					$first_8 .= "$line_vec[$j]\t";
				}
				elsif ($j == 8){
					if ($line_vec[$j] =~ /\:/){
						my @formats = split(/\:/, $line_vec[$j]);
						$line_vec[$j] = $formats[0];
					}
					$first_8 .= "$line_vec[$j]\t";
				}
			}
			my $fst_pos = int($pos[$fst_idx]);
			my $snd_pos = int($pos[$snd_idx]);
			my @fst_sample = split(/\/|\|/, $line_vec[$fst_pos]);
			my @snd_sample = split(/\/|\|/, $line_vec[$snd_pos]);
			my @returns; my $end; my $start;
			if ($syn !~ /syn/){
				$snd_sample[0] = $fst_sample[1];
			}
			if ($fst_sample[0] =~ /\./ || $snd_sample[0] =~ /\./){
				@returns = &beds("mask", $line_vec[0], $line_vec[1], $last_mask[$n], $chr, $record);
				$last_mask[$n] = $returns[0];
				$bed_out .= $returns[1];
				$record = $returns[2];
				if ($cnt == 1){
					$start = $region_start-1;
					$end = $line_vec[1]-2;
					if ($end > $start){
						$bed_out = "$chr\t$start\t$end\n".$bed_out;
					}
				}
				$cnt++;
			}
			elsif ($line_vec[6] =~ /LowQual/ || $line_vec[4] =~ /\,/){
				@returns = &beds("mask", $line_vec[0], $line_vec[1], $last_mask[$n], $chr, $record);
				$last_mask[$n] = $returns[0];
				$bed_out .= $returns[1];
				$record = $returns[2];
				if ($cnt == 1){
					$start = $region_start-1;
					$end = $line_vec[1]-2;
					if ($end > $start){
						$bed_out = "$chr\t$start\t$end\n".$bed_out;
					}
				}
				$cnt++;
			}
			elsif ($line_vec[4] =~ /\*/){
				@returns = &beds("mask", $line_vec[0], $line_vec[1], $last_mask[$n], $chr, $record);
				$last_mask[$n] = $returns[0];
				$bed_out .= $returns[1];
				$record = $returns[2];
				if ($cnt == 1){
					$start = $region_start-1;
					$end = $line_vec[1]-2;
					if ($end > $start){
						$bed_out = "$chr\t$start\t$end\n".$bed_out;
					}
				}
				$cnt++;
			}
			elsif ($fst_sample[0] eq '1' || $snd_sample[0] eq '1'){
				push(@vcf_outs, "$first_8$fst_sample[0]\|$snd_sample[0]\n");
				@returns = &beds("unmask", $line_vec[0], $line_vec[1], $last_mask[$n], $chr, $record);
				$last_mask[$n] = $returns[0];
				$bed_out .= $returns[1];
				$record = $returns[2];
				if ($cnt == 1){
					$start = $region_start-1;
					$end = $line_vec[1]-2;
					if ($end > $start){
						$bed_out = "$chr\t$start\t$end\n".$bed_out;
					}
				}
				$cnt++;
			}
			else {
				@returns = &beds("unmask", $line_vec[0], $line_vec[1], $last_mask[$n], $chr, $record);
				$last_mask[$n] = $returns[0];
				$bed_out .= $returns[1];
				$record = $returns[2];
				if ($cnt == 1){
					$start = $region_start-1;
					$end = $line_vec[1]-2;
					if ($end > $start){
						$bed_out = "$chr\t$start\t$end\n".$bed_out;
					}
				}
				$cnt++;
			}
			$fst_idx += 2;
			$snd_idx += 2;
			$first_8 = "";
		}
}
close(INPUT);
if ($record == 1){
	$bed_out .= "$region_end\n";
}
$vcf_out = join("", @vcf_outs);
return ($vcf_out, $bed_out);
}

sub beds { #bed file is 0-based, but vcf file is 1-based, so the position should minus 1
	my $mask = shift; my $chr = shift; my $start_pos = shift;
	my $last_mask = shift; my $chr = shift; my $record = shift;
	my @bed_out;
	if ($mask ne $last_mask || $last_mask eq "0"){
		if ($mask eq "unmask"){
			$start_pos -= 1;
			if ($record == 2){
				push(@bed_out, "$chr\t$start_pos\t");
				$record = 1;
			}
		}
		elsif ($mask eq "mask" && $last_mask ne "0"){
			$start_pos -= 2;
			#$start_pos -= 1;
			if ($record == 1){
				push(@bed_out, "$start_pos\n");
				$record = 2;
			}
		}
		$last_mask = $mask;
	}
	my $bedout = join("", @bed_out);
	return ($last_mask , $bedout, $record);
}

sub chr_name {
	$time = scalar localtime();
	my $file = shift; my $pre = shift;
	my @content; my @line; my @id;
	if (-e $file){}
	else {
		return "no";
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `bgzip \-cd $file \| head \-n 5000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 5000 $file`;
	}
	else {
		return "no";
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
			if ($pre){
                if ($line[3] =~ /\b$pre\b/i){
                    push(@id, $line[3]);
                }
			}
			else {
                push(@id, $line[3]);
			}
		}
	}
	return (@id);
} #get chromosome name

sub chr_lengths {
	$time = scalar localtime();
	my $vcf = shift; my $pre = shift;
	my @content; my @line; my @len;
	if ($vcf =~ /\.vcf\.gz$/){
		@content = `bgzip \-cd $vcf \| head \-n 5000`;
	}
	elsif ($vcf =~ /\.vcf$/){
		@content = `head -n 5000 $vcf`;
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
			if ($pre){
                if ($line[3] =~ /^\b$pre\b/i){
                    push(@len, $line[5]);
                }
			}
			else {
                push(@len, $line[5]);
			}
		}
	}
	return (@len);
} #get interval length from vcf header


