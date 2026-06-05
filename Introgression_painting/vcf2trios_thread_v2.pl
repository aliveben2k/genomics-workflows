#!/usr/bin/perl

use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
use File::Temp qw(tempfile);
use Fcntl qw(:flock);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
our @TMP_CHUNKS;
$SIG{INT} = sub { cleanup_tmp_chunks(); die "Interrupted.\n"; };
$SIG{TERM} = sub { cleanup_tmp_chunks(); die "Terminated.\n"; };
END { cleanup_tmp_chunks(); }
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

print "Input command line:\n";
print "perl vcf2trios_thread_v2\.pl @ARGV\n\n";

my $path; my @vcfs; my $list; my $out; my $chr; my $ran; my $skm = 0;
my $thread_in = 1; my $reg; my $sep = 0; my $path_l;
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
	if ($ARGV[$i] eq "\-list"){ #trio_list
        if (-e $ARGV[$i+1]){
            $list = $ARGV[$i+1];
        }
        if ($list =~ /\//){
			my @tmps = split(/\//, $list);
			pop(@tmps);
			$path_l = join("\/", @tmps);
        }
		else {
			$path_l = ".";
		}
		$path_l = &check_path($path_l);
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]";
	}
	if ($ARGV[$i] eq "\-t"){ #only extract the target chromosome (for multi-contig vcf)
		$sep = 1;
		if (length($ARGV[$i+1]) > 0){
			$chr = $ARGV[$i+1];
		}
	}
	if ($ARGV[$i] eq "\-skm"){ #skip monosite check
		$skm = 1;
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

unless(@vcfs){
	die "-vcf argument is required.\n";
}
unless ($path_l){
	$path_l = $path || ".";
	$path_l = &check_path($path_l);
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
my $info_out;
if ($ran){
	$info_out = "$path_l\/$ran\_genome_info.txt";
}
else {
	$info_out = "$path_l\/genome_info.txt";
}
&update_genome_info($info_out, \@chrs, \@chr_lengths);

my $thread; my @check;
foreach my $cnt (0..$#vcfs){
	@check = ();
	$thread = $thread_in;
	if ($vcfs[$cnt] !~ /gz$/){
		system("bgzip $vcfs[$cnt]");
		$vcfs[$cnt] = $vcfs[$cnt]."\.gz";
	}
    my $out = $vcfs[$cnt];
    $out =~ s/vcf\.gz$/$ran.\Qtrios\E.gz/;
    if (-e $out){
    	system("rm $out");
    }
    unless (-e "$vcfs[$cnt].tbi"){
    	system("tabix $vcfs[$cnt]");
    }
    if ($sep == 1 && length($chr) > 0){
    	$out =~ s/trios.gz$/$chr.$ran.\Qtrios\E.gz/;
    }
    my $out_base = $out;
    if (scalar(@vcfs) == 1){
    	foreach (@chrs){
    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
    		if ($return =~ /[a-z0-9]/i){
    			push(@check, $_);
    		}
    	}
    }
    else {
	    	foreach (@chrs){
	    		my $return = `tabix $vcfs[$cnt] $_ \| head -n 1`;
	    		if ($return =~ /[a-z0-9]/i){
	                push(@check, $_);
	    		}
	    	}   	
	    }
	    if ($sep == 1 && length($chr) > 0){
	    	@check = grep { $_ eq $chr } @check;
	    }
    my $line = `zcat $vcfs[$cnt] \| head -n 10000 \| grep \"\#CHROM"`;
    @vcf_samples = split(/\t/, $line);
    chomp(@vcf_samples);
    my @header_for_print; my @lists;
    if ($list){
    	open(LIST, "<$list") || die "Cannot open $list: $!\n";
    	my @list_tmp = <LIST>;
    	close(LIST);
    	chomp(@list_tmp);
    	foreach (@list_tmp){
    		my @line_eles = split(/\t|\s+/, $_);
    		shift(@line_eles);
    		push(@lists, @line_eles);
    	}
        foreach (@lists){
            foreach my $k (0..$#vcf_samples){
                if ($_ eq $vcf_samples[$k]){
                	push(@header_for_print, $vcf_samples[$k]);
                }
            }
        }
    }
    splice(@vcf_samples, 0, 9);
    unless ($list){
    	foreach (@vcf_samples){
    		push(@header_for_print, $_);
    	}
    }
    unless ($sep == 1 && scalar(@check) > 1){
        open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
        print OUT "Allele\t";
        print OUT join("\t", @header_for_print), "\n";
    }
    foreach my $i (0..$#check){
        if ($sep == 1 && scalar(@check) > 1){
    	    $out = $out_base;
    	    $out =~ s/trios.gz$/$check[$i].$ran\Qtrios\E.gz/;
            open(OUT, "|-", "gzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
            print OUT "Allele\t";
        	print OUT join("\t", @header_for_print), "\n";
        }
        my $chr_length; #get chromosome length
        foreach my $j (0..$#chrs){
            if ($check[$i] eq $chrs[$j]){
                my @vcf_header;
                if ($vcfs[$cnt] =~ /gz$/){
                    @vcf_header = `zcat $vcfs[$cnt] \| head -n 10000 \| grep \"\#\#contig\=\"`;
                }
                else {
                    @vcf_header = `cat $vcfs[$cnt] \| head -n 10000 \| grep \"\#\#contig\=\"`;
                }
                chomp(@vcf_header);
                foreach my $k (0..$#vcf_header){
                    my @tmp = split(/\,/, $vcf_header[$k]);
                    if ($tmp[0] =~ /$check[$i]\b/){
                        my @tmp2 = split(/\=/, $tmp[1]);
                        chomp(@tmp2);
                        $chr_length = $tmp2[1];
                    }
                }
            }
	        }
	        #for test only
	        if ($reg){
	    	    $chr_length = $reg;
	        }
	        #end of test
		    my $division;
		    if ($thread > 1){
		    	my $thread_for_chr = $thread;
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
    	    $division = $chr_length / $thread_for_chr;
	    	    if ($division - int($division) > 0.5){
    		    $thread_for_chr = $thread_for_chr - 1;
    	    }
    	    my @thr; my @chunk_files; my @chunk_regions;
    	    foreach my $d (1..$thread_for_chr){
    		    my $region_start; my $region_end;
    		    if ($d < $thread_for_chr){
    			    if ($d == 1){
    				    $region_start = 1;
    			    }
    			    else {
    				    $region_start = int($chr_length/$thread_for_chr*($d-1))+1;
    			    }
    			    $region_end = int($chr_length/$thread_for_chr*($d));
    		    }
    		    else {
    			    $region_start = int($chr_length/$thread_for_chr*($d-1))+1;
    			    $region_end = 0; # open-ended final chunk; avoids truncating if contig header length is shorter than records
    		    }
    		    my $chunk_chr = $check[$i];
    		    $chunk_chr =~ s/[^A-Za-z0-9_.-]/_/g;
    		    my ($chunk_fh, $chunk_file) = tempfile(
    		    	"vcf2trios.$chunk_chr.$d.XXXXXX",
    		    	DIR => $path_l,
    		    	SUFFIX => ".gz",
    		    	UNLINK => 0
    		    );
    		    close($chunk_fh);
    		    $chunk_files[$d] = $chunk_file;
    		    $chunk_regions[$d] = $region_end ? "$check[$i]:$region_start-$region_end" : "$check[$i]:$region_start-";
    		    push(@TMP_CHUNKS, $chunk_file);
    		    $thr[$d] = threads->create('v2b', ($vcfs[$cnt],$sep,$check[$i],$region_start,$region_end,$list,$skm,$chunk_file,1));
    	    }
    	    foreach my $d (1..$thread_for_chr){
    		    my $written_lines = $thr[$d]->join();
    		    die BOLD "Worker $d did not return a valid status for $check[$i].\n", RESET unless defined $written_lines;
    		    if ($written_lines > 0){
    		    	my $chunk_in = IO::Uncompress::Gunzip->new($chunk_files[$d]) || die BOLD "Cannot read gzip temporary output $chunk_files[$d]: $GunzipError", RESET, "\n";
    		    	while (my $chunk_line = <$chunk_in>){
    		    		print OUT $chunk_line;
    		    	}
    		    	close($chunk_in);
    		    }
    		    elsif ($d == $thread_for_chr){
    		    	warn "Final chunk $chunk_regions[$d] wrote zero lines for $check[$i]. Check whether the VCF has variants in this tail region after filtering.\n";
    		    }
    		    unlink($chunk_files[$d]) || warn "Cannot remove temporary output $chunk_files[$d]: $!\n" if -e $chunk_files[$d];
    		    undef($thr[$d]);
    	    }
    	    @TMP_CHUNKS = grep { -e $_ } @TMP_CHUNKS;
    	    undef(@thr);
    	    undef(@chunk_files);
	    }
		    else {
			    my $region_start = 0; my $region_end = 0;
			    &v2b($vcfs[$cnt],$sep,$check[$i],$region_start,$region_end,$list,$skm,\*OUT);
        }
        if ($sep == 1 && scalar(@check) > 1){
            close(OUT);
        }
    }
    unless ($sep == 1 && scalar(@check) > 1){
        close(OUT);
    }
}
print "Done.\n";

sub v2b {
    my $vcf = $_[0]; my $sep = $_[1]; my $chr = $_[2]; my $region_start = $_[3]; my $region_end = $_[4]; my $list = $_[5]; my $skm = $_[6]; my $sink = $_[7]; my $compress_sink = $_[8];
    my @lists;
    if ($list){
    	open(LIST, "<$list") || die "Cannot open $list: $!\n";
    	my @list_tmp = <LIST>;
    	close(LIST);
    	chomp(@list_tmp);
    	foreach (@list_tmp){
    		my @line_eles = split(/\t|\s+/, $_);
    		shift(@line_eles);
    		push(@lists, @line_eles);
    	}
    }
    if ($vcf =~ /\.gz$/){
    	unless (-e "$vcf.tbi"){
    		system("tabix $vcf");
    	}
    	if ($region_start != 0){
    		my $region;
    		if ($region_end != 0){
    			$region = "$chr\:$region_start\-$region_end";
    		}
    		else {
    			$region = "$chr\:$region_start\-";
    		}
        	open(INPUT, "-|", "tabix", "-h", $vcf, $region) || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    	else {
    		open(INPUT, "-|", "bgzip", "-dc", $vcf) || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    	}
    }
    else {
    	die BOLD "$vcf must be compressed first.", RESET, "\n";
    }
    my $out_fh;
    if (defined $sink){
    	if (ref($sink)){
    		$out_fh = $sink;
    	}
    	elsif ($compress_sink){
    		$out_fh = IO::Compress::Gzip->new($sink) || die BOLD "Cannot write gzip temporary output $sink: $GzipError", RESET, "\n";
    	}
    	else {
    		open($out_fh, ">$sink") || die BOLD "Cannot write temporary output $sink: $!", RESET, "\n";
    	}
    }
    else {
    	open($out_fh, ">", \$sink) || die BOLD "Cannot write output buffer: $!", RESET, "\n";
    }
    my $written_lines = 0;
    my @samples; my @index; my @sample_cols; my $sep_check = 0;
    my @header;
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        if ($line =~ /^\#/){
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
	                	my $last_sample_col = $#samples + 9;
	                	@sample_cols = (9..$last_sample_col);
	                }
	                else {
	                	@sample_cols = @index;
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
        if ($skm == 0){
        my @check = (0) x scalar(@sample_cols);
        foreach my $pos (0..$#sample_cols){
            my $j = $sample_cols[$pos];
            my $sum;
            my @format = split(/\:/, $eles[$j]);
            @allele = split(/\/|\|/, $format[0]);
            foreach (@allele){
                $sum += $_;
            }
            $check[$pos] = $sum;
        }
        @check = do { my %seen; grep { !$seen{$_}++ } @check };
        if (scalar(@check) == 1){
            next;
        }
        }
        #end skipping process
		my $out_line; my @list_array;  my @list_array2;
			if ($list){
				@list_array = (0) x scalar(@index);
				@list_array2 = (0) x scalar(@index);
			}
        foreach my $pos (0..$#sample_cols){
            my $j = $sample_cols[$pos];
            my @format = split(/\:/, $eles[$j]);
            my @allele; my $allele1; my $allele2;
            @allele = split(/\/|\|/, $format[0]);
            foreach my $x (0..$#allele){
                if ($allele[$x] =~ /[1-9]/){
                    $allele[$x] = $allele[$x];
                }
                elsif ($allele[$x] == 0){
                    $allele[$x] = $allele[$x];
                }
                unless ($allele[$x] =~ /[0-9]/){
                    $allele[$x] = "NA";
                }
            }
			unless ($list){ #no list
                push(@list_array, "$allele[0]");
                push(@list_array2, "$allele[1]");
            }
            else { #with list
            	$list_array[$pos] = "$allele[0]";
            	$list_array2[$pos] = "$allele[1]";
            }
        }
        print $out_fh "$eles[1]\_$eles[0]\_a1\t";
        print $out_fh join("\t", @list_array);
        print $out_fh "\n";
        print $out_fh "$eles[1]\_$eles[0]\_a2\t";
        print $out_fh join("\t", @list_array2);
        print $out_fh "\n";
        $written_lines += 2;
    }
    close(INPUT);
    close($out_fh) unless ref($_[7]);
	return $written_lines;
}

sub cleanup_tmp_chunks {
	return if threads->tid();
	foreach my $chunk_file (@TMP_CHUNKS){
		unlink($chunk_file) if defined $chunk_file && -e $chunk_file;
	}
	@TMP_CHUNKS = ();
}

sub update_genome_info {
	my ($info_out, $chr_ref, $len_ref) = @_;
	my $lock_file = "$info_out.lock";
	open(my $lock_fh, ">", $lock_file) || die "Cannot lock genome info $info_out: $!\n";
	flock($lock_fh, LOCK_EX) || die "Cannot lock genome info $info_out: $!\n";

	my %seen;
	my @lines;
	if (-e $info_out){
		open(my $in_fh, "<", $info_out) || die "Cannot read genome info $info_out: $!\n";
		while (my $line = <$in_fh>){
			chomp($line);
			next if $line eq "";
			if ($line =~ /^Chr\tLength$/){
				next;
			}
			my ($chr_name) = split(/\t/, $line);
			next if $seen{$chr_name}++;
			push(@lines, $line);
		}
		close($in_fh);
	}

	foreach my $i (0..$#$chr_ref){
		next if $seen{$chr_ref->[$i]};
		push(@lines, "$chr_ref->[$i]\t$len_ref->[$i]");
		$seen{$chr_ref->[$i]} = 1;
	}

	open(my $out_fh, ">", $info_out) || die "Cannot write the genome info to $info_out: $!\n";
	print $out_fh "Chr\tLength\n";
	foreach my $line (@lines){
		print $out_fh "$line\n";
	}
	close($out_fh);
	flock($lock_fh, LOCK_UN);
	close($lock_fh);
}

sub usage {
    print "Usage: perl vcf2trios_thread_v2.pl -vcf PATH [-list LIST_FILE] [-t CHR] [-skm] [-sn SN] [-n]\n";
}
