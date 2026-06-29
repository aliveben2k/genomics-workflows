#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
my $dir = getcwd;
use Time::HiRes qw(time);

chomp(@ARGV);
if ($#ARGV == -1){
	&usage;
	exit;
}

print "Input command line:\n";
print "perl vcf2cnv\.pl @ARGV\n\n";

my $path; my @vcfs; my $list; my @lists; my $out; my $sep = 0; my $chr; my $ran;
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
            open (LIST, "<$ARGV[$i+1]") || die BOLD "Cannot open $ARGV[$i+1]: $!", RESET, "\n";
            @lists = <LIST>;
            chomp(@lists);
            shift(@lists);
            close(LIST);
            foreach my $j (0..$#lists){
            	my @sp_lines = split(/\t+|\s+/, $lists[$j]);
            	$lists[$j] = $sp_lines[0];
            }
            $list = $ARGV[$i+1];
        }
	}
	if ($ARGV[$i] eq "\-sn"){
        $ran = "$ARGV[$i+1]\.";
	}
	if ($ARGV[$i] eq "\-sep"){
		$sep = 1;
		if (length($ARGV[$i+1]) > 0){
			$chr = $ARGV[$i+1];
		}
	}
}

unless(@vcfs){
	die "-vcf argument is required.\n";
}
my @chrs = &chr_name($vcfs[0]);
my $chr_names = join(" ", @chrs);
if ($chr_names !~ /\b$chr\b/ || length($chr) <= 0){
	$chr = "";
}
my $begin_time = time();
foreach my $vcf (@vcfs){
    if ($vcf =~ /\.gz$/){
        open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    else {
        open (INPUT, "<$vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    $out = $vcf;
    $out =~ s/\.gz$//;
    $out =~ s/vcf$/$ran\Qcnv\E.txt/;
    if (-e $out){
    	system("rm $out");
    }
    if ($sep == 1 && length($chr) > 0){
    	$out =~ s/cnv.txt$/$chr.$ran\Qcnv\E.txt/;
    }
    elsif ($sep == 1){
    	my @check;
    	if ($vcf =~ /gz$/){
    		@check = `zcat $vcf \| head -n 10000`;
    	}
    	else {
    		@check = `cat $vcf \| head -n 10000`;
    	}
        foreach (@check){
            if ($_ !~ /^\#/){
                my @line_eles = split(/\t/, $_);
                $out =~ s/cnv.txt$/$line_eles[0].$ran\Qcnv\E.txt.gz/;
                last;
            }
        }    	
    }
    open(OUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
    print OUT "CNV_ID\tCHR\tPOS\tEND.POS\t";
    my $cnt = 0; my @samples; my @index; my $sep_check = 0;
    while (my $line = <INPUT>){
        chomp($line);
        $line =~ s/[\x0A\x0D]//g;
        if ($line =~ /^\#/){
            if ($line =~ /^\#CHROM/){
                @samples = split(/\t/, $line);
                if ($list){
                	foreach (@lists){ #sample name
                		foreach my $k (0..$#samples){
                			if ($_ eq $samples[$k]){
                				push(@index, $k);
                			}
                		}
                	}
                	print OUT join("\t", @lists), "\n";
                }
                splice(@samples, 0, 9);
                unless ($list){
                	print OUT join("\t", @samples), "\n";
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
        my $cn;
        foreach my $i (1..$#ele8){ #determine the position of DS in FORMAT column
            if ($ele8[$i] eq "CN"){
                $cn = $i;
            }
        }
        unless ($cn){
            die "Cannot find the CN tag.\n";
        }
		my $out_line; my @list_array; my @check;
        @check = @samples;
        foreach my $j (9..$#eles){
            my @format = split(/\:/, $eles[$j]);
            my $allele; my $sum;
            $sum = $format[$cn];
            unless ($sum =~ /[0-9]/){
                $sum = "NA";
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
        $eles[7] =~ s/^END\=//i;
        my @eles7_tmp = split(/\;/, $eles[7]);
        $eles[7] = $eles7_tmp[0];
        undef(@eles7_tmp);
        print OUT "$eles[2]\t$eles[0]\t$eles[1]\t$eles[7]\t";
        print OUT join("\t", @list_array), "\n";
        $cnt++;
        print "\rProcessing $vcf: $cnt CNVs done.";
    }
    print "\n";
    close(INPUT);
    close(OUT);
}
my $end_time = time();
printf("%.2f\n", $end_time - $begin_time);
print "Done.\n";

sub usage {
    print "Usage: perl vcf2cnv.pl -vcf PATH [-list FILE] [-sep CHR]\n";
}
sub chr_name {
	$time = scalar localtime();
	my $file = shift;
	my @content; my @line; my @id;
#	print "chr_name\n";
	if (-e $file){}
	else {
		return "no";
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `gzip \-cd $file \| head \-n 10000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 10000 $file`;
	}
	else {
		return "no";
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
            push(@id, $line[3]);
		}
	}
	return (@id);
} #get chromosome name
