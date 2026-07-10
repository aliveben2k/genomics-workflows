#!/usr/local/bin/perl
use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
my $home = (getpwuid $>)[7];
if (-e "$home\/software\/qsub_subroutine.pl"){
	require "$home\/software\/qsub_subroutine.pl";
}
elsif (-e "$home\/qsub_subroutine.pl"){
	require "$home\/qsub_subroutine.pl";
}
my $time = scalar localtime();

chomp(@ARGV);
my $input = $ARGV[0];
my @chrs = &chr_name($input);
my $chr_names = join(" ", @chrs);
my @chr_lengths = &chr_lengths($input);
#get available chrs
my @avail_chr; my $ck_chr;
unless ($input =~ /gz$/){
    system("bgzip -\@ 8 $input");
    $input = "$input.gz";
} 
unless (-e "$input.tbi"){
    system("tabix -\@ 8 $input");
}
foreach my $n (0..$#chr_names){
    $ck_chr = `tabix -\@ 8 $input $chr_names[$n] \| head -n 1`;
    if ($ck_chr =~ /\b$chr_names[$n]\b/){
        push(@avail_chr, $chr_names[$n]);
        push(@avail_lengths, $chr_lengths[$n]);
    }
}

my $out = $input;

my @tmp_outs = split(/\//, $out);

my $chr;
if ($tmp_outs[-1] =~ /raw/){
	$tmp_outs[-1] =~ s/raw/star_filtered/;
}
else {
	$tmp_outs[-1] =~ s/\.vcf/\.star_filtered\.vcf/;
}
$out = join("\/", @tmp_outs);

print "\[$time\]\: Filtering start...\n";

open (OUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
my $thread = 8;
foreach my $n (0..$#avail_chr){
    my $division; my $results;
    #define lengths processed by each thread
    $division = $avail_lengths[$n] / $thread;
    if ($division - int($division) > 0.5){
    	$thread = $thread - 1;
    }
    my @thr;
    foreach my $d (1..$thread){
        my $region_start; my $region_end;
        if ($d < $thread){
            if ($d == 1){
                $region_start = 1;
            }
            else {
                $region_start = int($avail_lengths[$n]/$thread*($d-1))+1;
            }
            $region_end = int($avail_lengths[$n]/$thread*($d));
        }
        else {
            $region_start = int($avail_lengths[$n]/$thread*($d-1))+1;
            $region_end = $avail_lengths[$n];
        }
        $thr[$d] = threads->create('v2f', ($input,$avail_chr[$n],$region_start,$region_end));
    }
    foreach my $d (1..$thread){
    	$results .= $thr[$d]->join();
    }
    undef(@thr);    
    print OUT $results;
}
close(OUT);
$time = scalar localtime();
print "\[$time\]\: Filtering is done.\n";

sub v2f {
    my $vcf = $_[0]; my $chr = $_[1]; my $region_start = $_[2]; my $region_end = $_[3];

    if ($region_start != 0){
        open(INPUT, "-|", "tabix -h $vcf $chr\:$region_start\-$region_end") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    else {
    	open(INPUT, "-|", "bgzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
    }
    my $num = 1; 
    my @sample_two_alleles;
    my $no_DP = 0; my $less_min = 0; my $has_aster = 0; my $min_depth = 3;
    my $out;
    while (my $line = <INPUT>){
        $fixed = 0;
	    my @ALT_alleles_noAster; my $which_allele_aster; my $where_in_FORMAT_is_DP;
	    chomp $line;
	    if ($line =~ /^#/){
		    $out .= "$line\n";
	    }
        else {
            my @line_vec = split('\t', $line);
            my @ALT_alleles = split(',', $line_vec[4]); #original ALT column
            for (my $i=0; $i<=$#ALT_alleles; $i++){
                if ($ALT_alleles[$i] !~ /\*/){ #if the ALT allele is not *, put it in @ALT_alleles_noAster
                    push(@ALT_alleles_noAster, $ALT_alleles[$i]);
                }
                else {
                    $which_allele_aster = $i+1; #allele number of '*'
                }
            }
            if (@ALT_alleles_noAster eq () || @ALT_alleles_noAster == ()){ #if no ALT allele after filtering, skip the position
                next;
            }
            my $ALT_join = join(',', @ALT_alleles_noAster); #new ALT column
            if ($line !~ /\*/){
                $which_allele_aster = 500;
            }
            for (my $j=0; $j<=8; $j++){ #print information columns of the position except INFO
                if ($j == 4){
                    $line_vec[$j] = $ALT_join;
                }
                if ($j == 0){
                    $out .= "$line_vec[$j]";
                }
                else {
                    $out .= "\t$line_vec[$j]";
                }
            }
            my @format = split (':', $line_vec[8]);
            for (my $k=0; $k<=$#format; $k++){
                if ($format[$k] =~ /DP/){ #get DP index in the INFO column
                    $where_in_FORMAT_is_DP = $k;
                }
            }
            for (my $l=9; $l<=$#line_vec; $l++){
                my @sample_vec;
                @sample_vec = split(':', $line_vec[$l]);
                @sample_two_alleles = split(/\/|\|/, $sample_vec[0]);
                if ($sample_two_alleles[0] eq "." && $sample_two_alleles[1] eq "."){}
                elsif ($sample_vec[$where_in_FORMAT_is_DP] eq () || $sample_vec[$where_in_FORMAT_is_DP] == ()){
                    $sample_two_alleles[0] = ".";
                    $sample_two_alleles[1] = ".";
                    $fixed = 1;
                    $no_DP += 1;
                }
                elsif ($sample_vec[$where_in_FORMAT_is_DP] eq "." || int($sample_vec[$where_in_FORMAT_is_DP]) < $min_depth){
                    $sample_two_alleles[0] = ".";
                    $sample_two_alleles[1] = ".";
                    $less_min += 1;
                }
                elsif (int($sample_two_alleles[0]) eq $which_allele_aster || int($sample_two_alleles[1]) eq $which_allele_aster){
                    $sample_two_alleles[0] = ".";
                    $sample_two_alleles[1] = ".";
                    $has_aster += 1;
                }
                foreach my $m (0..$#sample_two_alleles){
                    if ($sample_two_alleles[$m] ne "."){
                        if (int($sample_two_alleles[$m]) > $which_allele_aster){
                            $sample_two_alleles[$m] = int($sample_two_alleles[$m]) - 1;
                        }
                    }
                }
                if ($sample_vec[0] =~ /\|/){
                    if ($sample_two_alleles[0] ne "." && $sample_two_alleles[1] ne "."){
                        $sample_vec[0] = join("\|", @sample_two_alleles);
                    }
                    else {
                        $sample_vec[0] = join("\/", @sample_two_alleles);
                    }
                }
                else {
                    $sample_vec[0] = join("\/", @sample_two_alleles);
                }
                if ($sample_vec[0] eq '.'){
                    $sample_vec[0] = './.';
                }
                my $join = join(':', @sample_vec);
                $out .= "\t$join";	
            }
            $out .= "\n";
        }
    }
    close(INPUT);
    return $out;
}
