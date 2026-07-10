#!/usr/bin/perl

use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
#use Time::HiRes qw(time);
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

print "Input command line:\n";
print "perl vcf2bimbam_thread\.pl @ARGV\n\n";

my $path; my @vcfs; my $list; my $out; my $chr; my $ran; my $ioff = 0;
my $thread_in = 1; my $reg;
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-vcf"){
		if (-d $ARGV[$i+1]){
            $path = $ARGV[$i+1];
            if ($path =~ /\/$/){
                $path =~ s/\/$//;
            }
            @vcfs = <$path\/*.vcf.gz>;
            unless (@vcfs){
                @vcfs = <$path\/*.vcf>;
                unless(@vcfs){
                    die "Cannot find vcfs.\n";
                }
            }
		}
        elsif (-e $ARGV[$i+1]){
            @vcfs = $ARGV[$i+1];
            if ($ARGV[$i+1] =~ /\//){
                my @tmps = split(/\//, $ARGV[$i+1]);
                pop(@tmps);
                $path = join("\/", @tmps);
            }
        }
		else {
            die "Cannot find the file\(s\).\n";
		}
	}
	if ($ARGV[$i] eq "\-list"){
        if (-e $ARGV[$i+1]){
            $list = $ARGV[$i+1];
        }
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]\.";
	}
	if ($ARGV[$i] eq "\-t"){ #only extract the target chromosome (for multi-contig vcf)
		$sep = 1;
		if (length($ARGV[$i+1]) > 0){
			$chr = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-ioff"){ #imputation off (don't use imputation values)
		$ioff = 1;
	}
	if ($ARGV[$i] eq "\-n"){ #threads
        $thread_in = "$ARGV[$i+1]";
        if ($thread_in =~ /[^0-9]/){
        	die "-n parameter should be an integer number.\n";
        }
	}
	#for test only
	if ($ARGV[$i] eq "\-reg"){
		$reg = $ARGV[$i+1];
	}
	#end of test
}

#print "debug: @vcfs\n";

unless(@vcfs){
	die "-vcf argument is required.\n";
}
my @chrs = &chr_name($vcfs[0]);
my @chr_lengths = &chr_lengths($vcfs[0]);
my $chr_names = join(" ", @chrs);
if ($chr){
	if ($chr_names !~ /\b$chr\b/){
		$chr = "";
		die "-t: Cannot find the target chromosome\/contig.\n";
	}
}
#print "debug: @chrs\n";

#my $begin_time = time();
my $thread;
foreach my $cnt (0..$#vcfs){
	$thread = $thread_in;
	if ($vcfs[$cnt] !~ /gz$/){
		system("bgzip $vcfs[$cnt]");
		$vcfs[$cnt] = $vcfs[$cnt]."\.gz";
	}
    my $out = $vcfs[$cnt];
    $out =~ s/vcf\.gz$/$ran\Qbimbam\E.gz/;
    if (-e $out){
    	system("rm $out");
    }
    unless (-e "$vcfs[$cnt].tbi"){
    	system("tabix $vcfs[$cnt]");
    }
    if ($sep == 1 && length($chr) > 0){
    	$out =~ s/bimbam.gz$/$chr.$ran\Qbimbam\E.gz/;
    }
    elsif (scalar(@vcfs) == 1){
    	my @check;
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			push(@check, $_);
    		}
    	}
    	if (scalar(@check) > 1){
    		$thread = 1;
    		print "-n option only support the single chromosome\/contig vcf currently.\n";
    		print "Reset thread to 1.\n";
    	}
    	else {
    		$chr = $check[0];
    	}
    }
    else {
    	my @check;
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			$chr = $_;
    			$out =~ s/bimbam.gz$/$chr.$ran\Qbimbam\E.gz/;
    		}
    	}
    	if (scalar(@check) > 1){
    		$thread = 1;
    		print "-n option only support the single chromosome\/contig vcf currently.\n";
    		print "Reset thread to 1.\n";
    	}     	
    }
    #for test only
    if ($reg){
    	$chr_lengths[$cnt] = $reg;
    }
    #end of test
	my $division; my $results;
	#print "debug: $thread\n";
	if ($thread > 1){
		if ($vcfs[$cnt] =~ /\.gz$/){
    		unless (-e "$vcfs[$cnt].tbi"){
    			system("tabix $vcfs[$cnt]");
    		}
		}
		else {
    		system("bgzip -\@ $thread $vcfs[$cnt]");
    		system("tabix $vcfs[$cnt].gz");
    		$vcf = "$vcfs[$cnt].gz";
		}
    	$division = $chr_lengths[$cnt] / $thread;
    	if ($division - int($division) > 0.5){
    		$thread = $thread - 1;
    	}
    	my @thr = ();
    	foreach my $d (1..$thread){
    		my $region_start; my $region_end;
    		if ($d < $thread){
    			if ($d == 1){
    				$region_start = 1;
    			}
    			else {
    				$region_start = int($chr_lengths[$cnt]/$thread*($d-1))+1;
    			}
    			$region_end = int($chr_lengths[$cnt]/$thread*($d));
    		}
    		else {
    			$region_start = int($chr_lengths[$cnt]/$thread*($d-1))+1;
    			$region_end = $chr_lengths[$cnt];
    		}
    		#print "debug: $vcfs[$cnt] $sep $chr $region_start $region_end $list $ioff\n";
    		$thr[$d] = threads->create('v2b', ($vcfs[$cnt],$sep,$chr,$region_start,$region_end,$list,$ioff));
    	}
    	foreach my $d (1..$thread){
    		$results .= $thr[$d]->join();
    	}
    	undef(@thr);
	}
	else {
		my $region_start = 0; my $region_end = 0;
		$results = &v2b($vcfs[$cnt],$sep,$chr,$region_start,$region_end,$list,$ioff);
    }
    open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
    print OUT $results;
    close(OUT);
}
#my $end_time = time();
#printf("%.2f\n", $end_time - $begin_time);
print "Done.\n";

sub v2b {
    my $vcf = $_[0]; my $sep = $_[1]; my $chr = $_[2]; my $region_start = $_[3]; my $region_end = $_[4]; my $list = $_[5]; my $ioff = $_[6];
    my @lists;
    #print "debug: $chr\:$region_start\-$region_end\n";
    if ($list){
        open (LIST, "<$list") || die BOLD "Cannot open $list: $!", RESET, "\n";
        @lists = <LIST>;
        chomp(@lists);
        shift(@lists);
        close(LIST);
        foreach my $j (0..$#lists){
            my @sp_lines = split(/\t+|\s+/, $lists[$j]);
            $lists[$j] = $sp_lines[0];
        }
    }
    if ($vcf =~ /\.gz$/){
    	unless (-e "$vcf.tbi"){
    		system("tabix $vcf");
    	}
    	if ($region_start != 0){
        	open(INPUT, "-|", "tabix -h $vcf $chr\:$region_start\-$region_end") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    	else {
    		open(INPUT, "-|", "bgzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    }
    else {
    	die BOLD "$vcf must be compressed first.", RESET, "\n";
    }
    my $out;
    my @samples; my @index; my $sep_check = 0;
    my @header;
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        if ($line =~ /^\#/){
        	#print "debug: $line\n";
            if ($line =~ /^\#CHROM/){
                @samples = split(/\t/, $line);
                if ($list){
                	foreach (@lists){
                		foreach my $k (0..$#samples){
                			if ($_ eq $samples[$k]){
                				push(@index, $k);
                				push(@header, $samples[$k]);
                			}
                		}
                	}
                }
                splice(@samples, 0, 9);
                unless ($list){
                	@header = @samples;
                }
            }
            next;
        }
        else { #skip lines with un-targeted chromosomes
        	if ($sep == 1 && length($chr) > 0){
        		if ($line !~ /^\b$chr\b/){
        			if ($sep_check == 1){
        				last;
        			}
        			next;
        		}
        		else {
        			$sep_check = 1;
        		}
        	}
        }
        my @eles = split(/\t/, $line);
        my @ele8 = split(/\:/, $eles[8]); #split FORMAT info
        my $ds;
        foreach my $i (1..$#ele8){ #determine the position of DS in FORMAT column
            if ($ele8[$i] eq "DS"){
                $ds = $i;
            }
        }
        unless ($ds){
            $ds = 0;
        }
        #skip only one allele type position
        my @check;
        if ($list){
            @check = @index;
        }
        else {
            @check = @samples;
        }
        foreach my $j (9..$#eles){
            my $sum;
            my @format = split(/\:/, $eles[$j]);
            @allele = split(/\/|\|/, $format[0]);
            if ($ds == 0 || $ioff == 1){ #if there is no "DS" field or turn off imputed data, use sum of "GT" field
                foreach (@allele){
                    if ($_ =~ /[1-9]/){
                    	$_ = 1;
                        $sum += $_;
                    }
                    elsif ($_ == 0){
                        $sum += $_;
                    }
                }
            }
            else {
                if ($format[$ds] ne "."){
                    $sum = $format[$ds];
                }
                else {
                    foreach (@allele){
                        if ($_ =~ /[1-9]/){
                        	$_ = 1;
                            $sum += $_;
                        }
                        elsif ($_ == 0){
                            $sum += $_;
                        }
                    }
                }
            }
            if ($list){
                foreach my $k (0..$#index){
                    if ($j == $index[$k]){
                        $check[$k] = $sum;
                    }
                }
            }
            else {
                $check[$j-9] = $sum;
            }
        }
        @check = do { my %seen; grep { !$seen{$_}++ } @check };
        if (scalar(@check) == 1){
            next;
        }
        #end skipping process
		my $out_line; my @list_array;
		if ($list){
			foreach (@index){
				push(@list_array, 0);
			}
		}
        foreach my $j (9..$#eles){
            my @format = split(/\:/, $eles[$j]);
            my @allele; my $sum; my $miss = 0;
            @allele = split(/\/|\|/, $format[0]);
            if ($ds == 0 || $ioff == 1){ #if there is no "DS" field or turn off imputed data, use sum of "GT" field
                foreach (@allele){
                    if ($_ eq '.'){
                        $miss ++;
                    }
                }
                if ($miss == 2){
                    $sum = "NA";
                }
                else {
                    foreach (@allele){
                        if ($_ =~ /[1-9]/){
                            $_ = 1;
                            $sum += $_;
                        }
                        elsif ($_ == 0){
                            $sum += $_;
                        }
                    }
                    unless ($sum =~ /[0-9]/){
                        $sum = "NA";
                    }
                }
            }
            else { #if there is a "DS" field, but has no value, use sum of "GT" field
                $sum = $format[$ds];
                if ($sum eq "."){
                    $sum = "";
                    foreach (@allele){
                        if ($_ eq '.'){
                            $miss ++;
                        }
                    }
                    if ($miss == 2){
                        $sum = "NA";
                    }
                    else {
                        foreach (@allele){
                            if ($_ =~ /[1-9]/){
                        	    $_ = 1;
                                $sum += $_;
                            }
                            elsif ($_ == 0){
                                $sum += $_;
                            }
                        }
                        unless ($sum >= 0 && $sum <= 2){
                            $sum = "NA";
                        }
                    }
                }
            }
            if ($sum ne "NA" && $sum > 2){
                print "The input vcf is not bi-allele vcf, please filter first.\n";
                close(OUT);
                system("rm $out");
                close(INPUT);
                exit;
            }
			unless ($list){ #no list
                push(@list_array, $sum);
            }
            else { #with list
            	foreach my $l (0..$#index){
            		if ($j == $index[$l]){
            			$list_array[$l] = $sum;
            		}
            	}
            }
        }
        my $n_0 = 0;
        my $n_1 = 0;
        my $n_2 = 0;
        foreach my $i (0..$#list_array){
            if ($list_array[$i] >= 0 && $list_array[$i] <= 0.5){
                $n_0++;
            }
            if ($list_array[$i] > 0.5 && $list_array[$i] < 1.5){
                $n_1++;
            }
            if ($list_array[$i] >= 1.5 && $list_array[$i] <= 2){
                $n_2++;
            }
        }
        my $alleles = "$n_0 $n_1 $n_2";
        my $count = () = $alleles =~ /\b0\b/g;
        if ($count >= 2){
            next;
        }
        $out .= "$eles[2]\:$eles[0]\:$eles[1]\, $eles[3]\, $eles[4]\, ";
        $out .= join(", ", @list_array);
        $out .= "\n";
    }
    close(INPUT);
	return $out;
}

sub usage {
    print "Usage: perl vcf2bimbam_thread.pl -vcf PATH [-list LIST_FILE] [-t CHR] [-ioff] [-sn SN] [-n]\n";
}
