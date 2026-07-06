#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
my $home = (getpwuid $>)[7];
if (-e "$home\/softwares\/qsub_subroutine.pl"){
	require "$home\/softwares\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
else {
	die "Cannot find required subroutine file: qsub_subroutine.pl\n";
}

my $o_path; my $hap = 0; my $in_path; my $list; my $pos; my $out_chr; my $tag = "";
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-i"){ 
		$in_path = $ARGV[$i+1]; #input path+prefix
		if ($in_path =~ /\/$/){
			$in_path =~ s/\/$//;
		}
	}
	if ($ARGV[$i] eq "-list"){ #subset sample list
		$list = $ARGV[$i+1];
		unless (-e $list){
		    die "Cannot find the $list.\n";
		}
	}
	if ($ARGV[$i] eq "\-o"){
		$o_path = $ARGV[$i+1];
		if ($o_path =~ /\/$/){
			$o_path =~ s/\/$//;
		}
		unless (-d $o_path){
			system("mkdir $o_path");
		}
		$o_path = &check_path($o_path);
	}
	if ($ARGV[$i] eq "\-pos"){
        $pos = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-chr"){
        $out_chr = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-tag"){
        $tag = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-hap"){ #if haploid
        $hap = 1;
	}
}

my @sublist;
if ($list){
    open(LIST, "<$list") || die "Cannot open $list: $!\n";
    @sublist = <LIST>;
    chomp(@sublist);
    close(LIST);
}
my @index;
if (-e "$in_path.sample"){
    open(SAM, "<$in_path.sample") || die "Cannot open $in_path.sample: $!\n";
    my @tmplist = <SAM>;
    chomp(@tmplist);
    close(SAM);
    splice(@tmplist, 0, 2);
    foreach my $i (0..$#tmplist){
        if (@sublist){
            foreach my $j (0..$#sublist){
                if ($tmplist[$i] =~ /\b$sublist[$j]\b/){
                    push(@index, $i);
                }
            }
        }
        else {
            push(@index, $i);
        }
    }
}
else {
    die "Cannot find the $in_path.sample.\n";
}

if (-e "$in_path.haps"){
    open(HAPS, "<$in_path.haps") || die "Cannot open $in_path.haps: $!\n";
    my @outdata; my $chr;
    while (my $line = <HAPS>){
        chomp($line);
        my @eles = split(/\t|\s+/, $line);
        if ($eles[2] eq $pos){
            $chr = $eles[0];
            splice(@eles, 0, 5); #remove information part
            foreach my $i (0..$#index){
                if ($hap == 0){
                    $order_idx = $index[$i] * 2;
                    push(@outdata, $eles[$order_idx]);
                    push(@outdata, $eles[$order_idx+1]);
                }
                else {
                    push(@outdata, $eles[$index[$i]]);
                }
            }
            last;
        }
    }
    close(HAPS);
    die "Cannot find position $pos in $in_path.haps.\n" unless @outdata;
    $chr = $out_chr if defined $out_chr && $out_chr ne "";
    my $suffix = "";
    $suffix = "_$tag" if defined $tag && $tag ne "";
    open(OUT, ">$o_path\/chr$chr\_$pos$suffix.txt") || die "Cannot write $o_path\/chr$chr\_$pos$suffix.txt: $!\n";
    print OUT join("\n", @outdata), "\n";
    close(OUT);
    #calculate allele frequency
    my $derive_no = 0;
    foreach my $i (0..$#outdata){
        if ($outdata[$i] != 0){
            $derive_no++;
        }
    }
    my $ratio = sprintf "%.4f", $derive_no / scalar(@outdata);
    open(OUT2, ">$o_path\/chr$chr\_$pos$suffix\_freq.txt") || die "Cannot write $o_path\/chr$chr\_$pos$suffix\_freq.txt: $!\n";
    print OUT2 $ratio;
    close(OUT2);
}


