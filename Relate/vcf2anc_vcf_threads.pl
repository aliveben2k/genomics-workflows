#!/usr/bin/perl

use Cwd qw(getcwd);
use Term::ANSIColor qw(:constants);
use threads;
use threads::shared;
use File::Temp qw(tempfile);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);

our @TMP_CHUNKS;
$SIG{INT} = sub { cleanup_tmp_chunks(); die "Interrupted.\n"; };
$SIG{TERM} = sub { cleanup_tmp_chunks(); die "Terminated.\n"; };
END { cleanup_tmp_chunks(); }

chomp(@ARGV);
print "The script is written by Ben Chien. Jan. 2023.\n";
print "Usage: perl vcf2anc_vcf_threads.pl -vcf TARGET_VCF -aid ANCESTRAL_ID\|A_LIST_FILE [-o OUTPUT_FILE_NAME] [-list A_LIST_FILE] [-bi] [-rchr RENAME_FILE] [-keep] [-hap] [-nm] [-sf] [-thap] [-map GENETIC_MAP] [-shapeit CONDA_ENV] [-n THREADS]\n";
print "-aid: an ancestral id. Multiple samples could be seperated by comma, or listed in a file \(one sample per line\).\n";
print "-o: output file name without extension.\n";
print "-list: only keep the listed samples.\n";
print "-bi: only keep bi-allele SNPs. Default: false\n";
print "-rchr: rename CHR column. Default: false\n";
print "-keep: also keep the ancestral genotype. Default: false\n";
print "-hap: output one haplotype per sample for Relate. With -sf, SHAPEIT5 is run first and the first phased allele is used. Default: false\n";
print "-nm: no missing, the missing haploid will be imputed as the other available haploid or the ancient one.\n";
print "-sf: output Relate haps/sample format. Genotypes are phased/imputed with SHAPEIT5 before conversion.\n";
print "-thap: outputting rehh\:thap format instead of vcf format.\n";
print "-map: genetic map for SHAPEIT5_phase_common when -sf is used.\n";
print "-shapeit: conda environment name for SHAPEIT5. Default: shapeit5\n";
print "-n: number of threads. Default: 1. Multi-thread mode requires a bgzip/tabix indexed VCF.\n";
print "The output file is the same as the input vcf, and the name is PREFIX.anc.vcf.gz\n";
print "This is a local script, not a server script.\n\n";
print "Input command line:\n";
print "perl vcf2anc_vcf_threads\.pl @ARGV\n\n";

my $vcf; my @aids; my $hap = 0; my $kaid = 0; my $list; my @lists; my $bi = 0; my @rm_lists; my $nm = 0; my $sf = 0; my $out; my $thap = 0; my $thread_in = 1; my $path = "."; my $shapeit_map; my $shapeit = "shapeit5";
for (my $i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] eq "\-vcf"){
		if (-e $ARGV[$i+1]){
            $vcf = $ARGV[$i+1];
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
	if ($ARGV[$i] eq "\-aid"){
		if (-e $ARGV[$i+1]){
			open(ALIST, "<$ARGV[$i+1]") || die BOLD "Cannot open aid list, $ARGV[$i+1]: $!", RESET, "\n";
			@aids = <ALIST>;
			chomp(@aids);
			close(ALIST);
            foreach my $j (0..$#aids){
            	my @sp_lines = split(/\t+|\s+/, $aids[$j]);
            	$aids[$j] = $sp_lines[0];
            }	  
		}
		elsif ($ARGV[$i+1] =~ /\,/){
			@aids = split(/\,/, $ARGV[$i+1]);
		}
		else {
			push(@aids, $ARGV[$i+1]);
		}
	}
	if ($ARGV[$i] eq "\-keep"){ #keep ancestral samples
        $kaid = 1;
	}
	if ($ARGV[$i] eq "\-hap"){
        $hap = 1;
	}
	if ($ARGV[$i] eq "\-bi"){
        $bi = 1;
	}
	if ($ARGV[$i] eq "\-nm"){ #no missing
        $nm = 1;
	}
	if ($ARGV[$i] eq "\-rchr"){
		if (-e $ARGV[$i+1]){
			open(RE, "<$ARGV[$i+1]") || die BOLD "Cannot open the rename list, $ARGV[$i+1]: $!", RESET, "\n";
			@rm_lists = <RE>;
			chomp(@rm_lists);
			close(RE);
        }
        else {
        	die BOLD "Cannot find the rename list.", RESET, "\n";
        }
	}
	if ($ARGV[$i] eq "\-list"){
        if (-e $ARGV[$i+1]){
            open (LIST, "<$ARGV[$i+1]") || die BOLD "Cannot open sample list, $ARGV[$i+1]: $!", RESET, "\n";
            @lists = <LIST>;
            chomp(@lists);
            if ($ARGV[$i+1] =~ /poplabels/){
            	shift(@lists);
            }
            close(LIST);
            foreach my $j (0..$#lists){
            	my @sp_lines = split(/\t+|\s+/, $lists[$j]);
            	$lists[$j] = $sp_lines[0];
            }
            $list = $ARGV[$i+1];
        }
	}
	if ($ARGV[$i] eq "\-sf"){ #shapeit format for Relate
        $sf = 1;
	}
	if ($ARGV[$i] eq "\-thap"){ #thap format: for rehh
        $thap = 1;
	}	
	if ($ARGV[$i] eq "\-map"){
		if (-e $ARGV[$i+1]){
			$shapeit_map = $ARGV[$i+1];
		}
		else {
			die BOLD "Cannot find the genetic map for SHAPEIT5: $ARGV[$i+1]", RESET, "\n";
		}
	}
	if ($ARGV[$i] eq "\-shapeit"){
		$shapeit = $ARGV[$i+1];
		unless ($shapeit && $shapeit !~ /^\-/){
			die "-shapeit requires a conda environment name.\n";
		}
	}
	if ($ARGV[$i] eq "\-o"){
        $out = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "\-n"){
        $thread_in = $ARGV[$i+1];
        if ($thread_in =~ /[^0-9]/ || $thread_in < 1){
        	die "-n parameter should be an integer number >= 1.\n";
        }
	}
}

unless (-e $vcf){
	die "Cannot find target vcf: $vcf.\n";
}

unless ($out){
	$out = $vcf;
	$out =~ s/.gz$//;
	if ($list){
		if ($list =~ /\//){
			my @tmp = split(/\//, $list);
			$list = $tmp[-1];
		}
		$list =~ s/\.txt$|\.list$//;
		$out =~ s/\.vcf$/.$list.vcf/;
	}
	$out =~ s/\.vcf$/.anc.vcf.gz/;
	if ($sf == 1){
		$out =~ s/vcf\.gz$/haps/;
	}
	elsif ($thap == 1){
		$out =~ s/vcf\.gz$/thap/;
	}
}
else {
	if ($thap == 1){
		$out .= ".thap";
	}
	elsif ($sf == 0){ #vcf
		$out .= ".vcf.gz";
	}
	else { #haps format
		$out .= ".haps";
	}
}
my $shapeit_haps_out; my $shapeit_sample_out; my $shapeit_unphased_vcf; my $shapeit_tagged_vcf; my $shapeit_phased_bcf; my $run_shapeit = 0; my $output_hap = $hap;
if ($sf == 1){
	die BOLD "-sf requires -map so SHAPEIT5_phase_common can phase before haps/sample conversion.\n", RESET unless $shapeit_map;
	$run_shapeit = 1;
	$shapeit_haps_out = $out;
	$shapeit_sample_out = $out;
	$shapeit_sample_out =~ s/\.haps$/.sample/;
	$shapeit_unphased_vcf = "$out.unphased.vcf.gz";
	$shapeit_tagged_vcf = "$out.shapeit5_input.vcf.gz";
	$shapeit_phased_bcf = "$out.shapeit5_phased.bcf";
	push(@TMP_CHUNKS, $shapeit_unphased_vcf, "$shapeit_unphased_vcf.tbi", $shapeit_tagged_vcf, "$shapeit_tagged_vcf.tbi", $shapeit_phased_bcf, "$shapeit_phased_bcf.csi");
	$out = $shapeit_unphased_vcf;
	$sf = 0;
	$hap = 0;
}
if ($sf == 1){
	my $out2 = $out;
	$out2 =~ s/\.haps$/.sample/;
	open(OUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
	open(SAMPLE, ">$out2") || die BOLD "Cannot write $out2: $!", RESET, "\n";
}
elsif ($thap == 1){
	my $out2 = $out;
	$out2 =~ s/\.thap$/.map/;
	open(OUT, ">$out") || die BOLD "Cannot write $out: $!", RESET, "\n";
	open(SAMPLE, ">$out2") || die BOLD "Cannot write $out2: $!", RESET, "\n";
}
else {
	open(OUT, "|-", "bgzip \> $out") || die BOLD "Cannot write $out: $!", RESET, "\n";
}
print "Start processing vcf...\n";
my ($header_lines, $ids_ref, $contig_lengths_ref) = read_vcf_header($vcf);
my %context = prepare_header_context($ids_ref, \@aids, \@lists, $list, $kaid, $sf, $thap, $hap);
my @aid_indexes = @{$context{aid_indexes}};
my @sample_indexes = @{$context{sample_indexes}};
my @sample_names = @{$context{sample_names}};
if ($sf == 0 && $thap == 0){
	foreach my $header_line (@$header_lines){
		print OUT "$header_line\n";
	}
	print OUT join("\t", @{$context{fixed_fields}}, @sample_names), "\n";
}
elsif ($sf == 1){
	print SAMPLE "ID_1 ID_2 missing\n0 0 0\n";
	foreach my $sample (@sample_names){
		if ($hap == 0){
			print SAMPLE "$sample $sample 0\n";
		}
		else {
			print SAMPLE "$sample NA 0\n";
		}
	}
	close(SAMPLE);
}

if ($thread_in > 1){
	die BOLD "-n > 1 requires bgzip-compressed VCF input ending in .gz.\n", RESET unless $vcf =~ /\.gz$/;
	unless (-e "$vcf.tbi"){
		system("tabix", $vcf) == 0 || die BOLD "Cannot index $vcf with tabix.\n", RESET;
	}
	my @contigs = get_tabix_contigs($vcf);
	my @thr; my @chunk_files; my @map_chunk_files; my @chunk_regions;
	foreach my $contig (@contigs){
		$contig =~ s/^\s+|\s+$//g;
		my $chr_length = $contig_lengths_ref->{$contig};
		my $thread_for_chr = $thread_in;
		if (!$chr_length){
			my @known_contigs = sort keys %$contig_lengths_ref;
			my $known = @known_contigs ? join(",", @known_contigs[0..($#known_contigs < 9 ? $#known_contigs : 9)]) : "none";
			warn "Cannot find contig length for $contig in VCF header. Parsed contigs with lengths: $known. Processing this contig as one whole region instead of splitting across threads.\n";
			$thread_for_chr = 1;
		}
		else {
			$thread_for_chr = $chr_length if $chr_length < $thread_for_chr;
		}
		@thr = (); @chunk_files = (); @map_chunk_files = (); @chunk_regions = ();
		foreach my $d (1..$thread_for_chr){
			my $region_start;
			my $region_end;
			if (!$chr_length){
				$region_start = 1;
				$region_end = 0;
			}
			elsif ($d < $thread_for_chr){
				$region_start = ($d == 1) ? 1 : int($chr_length/$thread_for_chr*($d-1))+1;
				$region_end = int($chr_length/$thread_for_chr*$d);
			}
			else {
				$region_start = int($chr_length/$thread_for_chr*($d-1))+1;
				$region_end = 0;
			}
			my $chunk_chr = $contig;
			$chunk_chr =~ s/[^A-Za-z0-9_.-]/_/g;
			my ($chunk_fh, $chunk_file) = tempfile("vcf2anc.$chunk_chr.$d.XXXXXX", DIR => $path, SUFFIX => ".gz", UNLINK => 0);
			close($chunk_fh);
			push(@TMP_CHUNKS, $chunk_file);
			$chunk_files[$d] = $chunk_file;
			if ($thap == 1){
				my ($map_fh, $map_file) = tempfile("vcf2anc.$chunk_chr.$d.map.XXXXXX", DIR => $path, SUFFIX => ".gz", UNLINK => 0);
				close($map_fh);
				push(@TMP_CHUNKS, $map_file);
				$map_chunk_files[$d] = $map_file;
			}
			$chunk_regions[$d] = $region_end ? "$contig:$region_start-$region_end" : "$contig:$region_start-";
			$thr[$d] = threads->create(
				'process_region',
				$vcf, $contig, $region_start, $region_end, $chunk_file, $map_chunk_files[$d],
				\@aid_indexes, \@sample_indexes, \@rm_lists,
				$bi, $hap, $kaid, $nm, $sf, $thap
			);
		}
		foreach my $d (1..$thread_for_chr){
			my $written_lines = $thr[$d]->join();
			die BOLD "Worker $d did not return a valid status for $contig.\n", RESET unless defined $written_lines;
			if ($written_lines > 0){
				stream_gzip_file($chunk_files[$d], \*OUT);
				stream_gzip_file($map_chunk_files[$d], \*SAMPLE) if $thap == 1;
			}
			elsif ($d == $thread_for_chr){
				warn "Final chunk $chunk_regions[$d] wrote zero lines for $contig after filtering.\n";
			}
			unlink($chunk_files[$d]) || warn "Cannot remove temporary output $chunk_files[$d]: $!\n" if -e $chunk_files[$d];
			unlink($map_chunk_files[$d]) || warn "Cannot remove temporary output $map_chunk_files[$d]: $!\n" if $thap == 1 && -e $map_chunk_files[$d];
			undef($thr[$d]);
		}
		@TMP_CHUNKS = grep { -e $_ } @TMP_CHUNKS;
	}
}
else {
	my $input_fh = open_vcf_stream($vcf);
	while (my $line = <$input_fh>){
		next if $line =~ /^\#/;
		chomp($line);
		my ($out_line, $map_line) = process_record_line($line, \@aid_indexes, \@sample_indexes, \@rm_lists, $bi, $hap, $kaid, $nm, $sf, $thap);
		next unless defined $out_line && length($out_line) > 0;
		print OUT "$out_line\n";
		print SAMPLE "$map_line\n" if $thap == 1 && defined $map_line && length($map_line) > 0;
	}
	close($input_fh);
}
if ($thap == 1){
	close(SAMPLE);
}
close(OUT);
if ($run_shapeit == 1){
	&run_shapeit5_phase_common($shapeit_unphased_vcf, $shapeit_tagged_vcf, $shapeit_phased_bcf, $shapeit_map, $thread_in, $shapeit);
	&convert_phased_bcf_to_haps($shapeit_phased_bcf, $shapeit_haps_out, $shapeit_sample_out, $output_hap);
	foreach my $tmp ($shapeit_unphased_vcf, "$shapeit_unphased_vcf.tbi", $shapeit_tagged_vcf, "$shapeit_tagged_vcf.tbi", $shapeit_phased_bcf, "$shapeit_phased_bcf.csi"){
		unlink($tmp) if defined $tmp && -e $tmp;
	}
	@TMP_CHUNKS = grep { -e $_ } @TMP_CHUNKS;
}
print "Done.\n";

sub open_vcf_stream {
	my $vcf = shift;
	my $fh;
	if ($vcf =~ /\.gz$/){
		open($fh, "-|", "gzip", "-dc", $vcf) || die BOLD "Cannot open vcf $vcf: $!", RESET, "\n";
	}
	else {
		open($fh, "<$vcf") || die BOLD "Cannot open vcf $vcf: $!", RESET, "\n";
	}
	return $fh;
}

sub read_vcf_header {
	my $vcf = shift;
	my $fh;
	if ($vcf =~ /\.gz$/){
		open($fh, "-|", "bcftools", "view", "-h", $vcf) || die BOLD "Cannot read VCF header from $vcf with bcftools view -h: $!", RESET, "\n";
	}
	else {
		$fh = open_vcf_stream($vcf);
	}
	my @header_lines;
	my @ids;
	my %contig_lengths;
	while (my $line = <$fh>){
		$line =~ s/[\x0A\x0D]+$//;
		if ($line =~ /^##contig=/){
			if ($line =~ /ID=([^,\s>]+)/){
				my $contig_id = $1;
				$contig_id =~ s/^\s+|\s+$//g;
				if ($line =~ /length\s*=\s*([0-9]+)/i){
					$contig_lengths{$contig_id} = $1;
				}
			}
		}
		if ($line =~ /^#CHROM/){
			@ids = split(/\t/, $line);
			last;
		}
		push(@header_lines, $line) if $line =~ /^#/;
	}
	close($fh);
	die "Cannot find #CHROM header in $vcf.\n" unless @ids;
	return (\@header_lines, \@ids, \%contig_lengths);
}

sub prepare_header_context {
	my ($ids_ref, $aids_ref, $lists_ref, $list, $kaid, $sf, $thap, $hap) = @_;
	my @ids = @$ids_ref;
	my @aid_names = @$aids_ref;
	my @aid_indexes;
	my %aid_name_seen;
	my %list_seen = map { $_ => 1 } @$lists_ref;
	foreach my $i (9..$#ids){
		foreach my $aid (@aid_names){
			if ($ids[$i] eq $aid){
				push(@aid_indexes, $i);
				$aid_name_seen{$aid} = 1;
			}
		}
	}
	die "Cannot find any ancestral ID in the vcf.\n" unless @aid_indexes;
	my %aid_index_seen = map { $_ => 1 } @aid_indexes;
	my @sample_indexes;
	my @sample_names;
	foreach my $k (9..$#ids){
		next if $kaid == 0 && $aid_index_seen{$k};
		next if $list && !$list_seen{$ids[$k]};
		push(@sample_indexes, $k);
		push(@sample_names, $ids[$k]);
	}
	my @fixed_fields = @ids[0..8];
	return (
		aid_indexes => \@aid_indexes,
		sample_indexes => \@sample_indexes,
		sample_names => \@sample_names,
		fixed_fields => \@fixed_fields,
	);
}

sub get_tabix_contigs {
	my $vcf = shift;
	open(my $fh, "-|", "tabix", "-l", $vcf) || die BOLD "Cannot list contigs from $vcf: $!", RESET, "\n";
	my @contigs = <$fh>;
	chomp(@contigs);
	close($fh);
	die "Cannot find contigs from tabix index for $vcf.\n" unless @contigs;
	return @contigs;
}

sub process_region {
	my ($vcf, $chr, $region_start, $region_end, $chunk_file, $map_chunk_file, $aids_ref, $sample_indexes_ref, $rm_lists_ref, $bi, $hap, $kaid, $nm, $sf, $thap) = @_;
	my $region = $region_end ? "$chr:$region_start-$region_end" : "$chr:$region_start-";
	open(my $input, "-|", "tabix", "-h", $vcf, $region) || die BOLD "Cannot open $vcf region $region: $!", RESET, "\n";
	my $out_fh = IO::Compress::Gzip->new($chunk_file) || die BOLD "Cannot write gzip temporary output $chunk_file: $GzipError", RESET, "\n";
	my $map_fh;
	if ($thap == 1){
		$map_fh = IO::Compress::Gzip->new($map_chunk_file) || die BOLD "Cannot write gzip temporary map $map_chunk_file: $GzipError", RESET, "\n";
	}
	my $written_lines = 0;
	while (my $line = <$input>){
		next if $line =~ /^\#/;
		chomp($line);
		my ($out_line, $map_line) = process_record_line($line, $aids_ref, $sample_indexes_ref, $rm_lists_ref, $bi, $hap, $kaid, $nm, $sf, $thap);
		next unless defined $out_line && length($out_line) > 0;
		print $out_fh "$out_line\n";
		print $map_fh "$map_line\n" if $thap == 1 && defined $map_line && length($map_line) > 0;
		$written_lines++;
	}
	close($input);
	close($out_fh);
	close($map_fh) if $thap == 1;
	return $written_lines;
}

sub stream_gzip_file {
	my ($file, $out_fh) = @_;
	return unless defined $file && -e $file;
	my $in = IO::Uncompress::Gunzip->new($file) || die BOLD "Cannot read gzip temporary output $file: $GunzipError", RESET, "\n";
	while (my $line = <$in>){
		print $out_fh $line;
	}
	close($in);
}

sub run_shapeit5_phase_common {
	my ($unphased_vcf, $tagged_vcf, $phased_bcf, $shapeit_map, $threads, $shapeit) = @_;
	system("tabix", "-f", "-p", "vcf", $unphased_vcf) == 0 || die BOLD "Cannot index $unphased_vcf with tabix.\n", RESET;
	system("bcftools", "+fill-tags", $unphased_vcf, "-Oz", "-o", $tagged_vcf, "--", "-t", "AC,AN") == 0 || die BOLD "Cannot add AC/AN tags to $unphased_vcf with bcftools +fill-tags.\n", RESET;
	system("tabix", "-f", "-p", "vcf", $tagged_vcf) == 0 || die BOLD "Cannot index $tagged_vcf with tabix.\n", RESET;
	my @contigs = get_tabix_contigs($tagged_vcf);
	die BOLD "SHAPEIT5 phasing expects one chromosome per converter run; found: ", join(",", @contigs), "\n", RESET if scalar(@contigs) != 1;
	my $region = $contigs[0];
	system("conda", "run", "-n", $shapeit, "SHAPEIT5_phase_common", "--input", $tagged_vcf, "--map", $shapeit_map, "--region", $region, "--output", $phased_bcf, "--thread", $threads) == 0 || die BOLD "SHAPEIT5_phase_common failed for $tagged_vcf.\n", RESET;
}

sub convert_phased_bcf_to_haps {
	my ($phased_bcf, $haps_out, $sample_out, $output_hap) = @_;
	open(my $input_fh, "-|", "bcftools", "view", $phased_bcf) || die BOLD "Cannot read phased BCF $phased_bcf with bcftools view: $!", RESET, "\n";
	open(my $haps_fh, ">$haps_out") || die BOLD "Cannot write $haps_out: $!", RESET, "\n";
	open(my $sample_fh, ">$sample_out") || die BOLD "Cannot write $sample_out: $!", RESET, "\n";
	my @sample_names;
	while (my $line = <$input_fh>){
		chomp($line);
		if ($line =~ /^#CHROM/){
			my @ids = split(/\t/, $line);
			@sample_names = @ids[9..$#ids];
			print $sample_fh "ID_1 ID_2 missing\n0 0 0\n";
			foreach my $sample (@sample_names){
				if ($output_hap == 1){
					print $sample_fh "$sample NA 0\n";
				}
				else {
					print $sample_fh "$sample $sample 0\n";
				}
			}
			close($sample_fh);
			next;
		}
		next if $line =~ /^\#/;
		die "Cannot find sample header in phased BCF $phased_bcf.\n" unless @sample_names;
		my @eles = split(/\t/, $line);
		my @alt = split(/\,/, $eles[4]);
		my @hap_alleles;
		foreach my $i (9..$#eles){
			my @gt = split(/\:/, $eles[$i]);
			my @alleles = split(/\/|\|/, $gt[0]);
			if (scalar(@alleles) == 1){
				push(@alleles, $alleles[0]);
			}
			my @alleles_to_output = $output_hap == 1 ? ($alleles[0]) : @alleles[0..1];
			foreach my $allele (@alleles_to_output){
				if (!defined $allele || $allele eq "."){
					push(@hap_alleles, ".");
				}
				elsif ($allele =~ /^[0-9]+$/ && $allele > 1){
					push(@hap_alleles, 1);
				}
				else {
					push(@hap_alleles, $allele);
				}
			}
		}
		print $haps_fh join(" ", $eles[0], $eles[2], $eles[1], $eles[3], $alt[0], @hap_alleles), "\n";
	}
	close($input_fh);
	die BOLD "bcftools view failed for phased BCF $phased_bcf.\n", RESET if $?;
	close($haps_fh);
}

sub process_record_line {
	my ($line, $aids_ref, $sample_indexes_ref, $rm_lists_ref, $bi, $hap, $kaid, $nm, $sf, $thap) = @_;
	my @aids = @$aids_ref;
	my @sample_indexes = @$sample_indexes_ref;
	my @rm_lists = @$rm_lists_ref;
	my @eles = split(/\t/, $line);
	my $nucl = &get_ancestral_allele($line, \@aids);
	$nucl = uc($nucl);
	$eles[2] = "$eles[0]\_$eles[1]";
	my $check_existance = $eles[3];
	$check_existance .= ",$eles[4]";
	if ($nucl =~ /\*/){
		if ($check_existance !~ /\*/){
			return;
		}
		$nucl =~ s/\*/B/g;
	}
	if ($check_existance =~ /\*/){
		$check_existance =~ s/\*/B/g;
	}
	if ($check_existance !~ /$nucl/){
		return;
	}
	$nucl =~ s/B/\*/g;
	$check_existance =~ s/B/\*/g;
	my $anc; my @derivs;
	my @sorted_nucls = split(/\,/, $check_existance);
	foreach my $k (0..$#sorted_nucls){
		if ($nucl eq $sorted_nucls[$k]){
			$anc = $k;
		}
		else {
			push(@derivs, $k);
		}
	}
	return unless defined $anc;
	unshift(@derivs, $anc);
	@sorted_nucls = @sorted_nucls[@derivs];
	if ($bi == 1 && $hap == 0){
		if (scalar(@derivs) > 2){
			return;
		}
	}
	my @out_eles;
	foreach my $j (@sample_indexes){
		my @out_alleles;
		my @gt = split(/\:/, $eles[$j]);
		my $gt_sep = ($gt[0] =~ /\|/) ? "\|" : "\/";
		my @alleles = split(/\/|\|/, $gt[0]);
		if (scalar(@alleles) != 2 && $hap == 0){
			if ($nm == 0){
				@alleles = ('.', '.');
				push(@out_eles, join($gt_sep, @alleles));
			}
			if ($nm == 1){
				if (scalar(@alleles) == 1){
					@alleles = ($alleles[0], $alleles[0]);
				}
				else {
					@alleles = (0, 0);
				}
				push(@out_eles, join($gt_sep, @alleles));
			}
		}
		elsif (scalar(@alleles) == 1 && $hap == 1){
			if ($alleles[0] !~ /[^0-9]/){
				foreach my $m (0..$#derivs){
					if ($alleles[0] == $derivs[$m]){
						$out_alleles[0] = $m;
					}
				}
				push(@out_eles, $out_alleles[0]);
			}
			else {
				if ($nm == 0){
					push(@out_eles, '.');
				}
				else {
					push(@out_eles, 0);
				}
			}
		}
		else {
			foreach my $l (0..1){
				if ($alleles[$l] !~ /[^0-9]/){
					foreach my $m (0..$#derivs){
						if ($alleles[$l] == $derivs[$m]){
							push(@out_alleles, $m);
						}
					}
				}
				else {
					if ($nm == 0){
						push(@out_alleles, '.');
					}
					else {
						push(@out_alleles, 0);
					}
				}
			}
			if ($hap == 0){
				push(@out_eles, join($gt_sep, @out_alleles));
			}
			else {
				push(@out_eles, $out_alleles[0]);
			}
		}
	}
	foreach my $q (0..8){
		if ($q == 0 && @rm_lists){
			$eles[$q] = &rename_chr($eles[$q], \@rm_lists);
		}
		if ($q == 3){
			$eles[$q] = $sorted_nucls[0];
			if ($thap == 1){
				$eles[$q] = 0;
			}
			shift(@sorted_nucls);
		}
		if ($q == 4){
			if ($bi == 0){
				if ($thap == 1){
					foreach my $x (1..scalar(@sorted_nucls)){
						$sorted_nucls[$x-1] = $x;
					}
				}
				$eles[$q] = join(",", @sorted_nucls);
			}
			else {
				if ($thap == 1){
					$eles[$q] = 1;
				}
			}
		}
		if ($q == 8){
			$eles[$q] = "GT";
		}
	}
	if ($bi == 1){
		my @unique = do {my %seen; sort {$a <=> $b} grep {$_ ne "." && !$seen{$_}++} map {split(/\/|\|/, $_)} @out_eles};
		return unless scalar(@unique) == 2;
		$eles[4] = $sorted_nucls[$unique[1]-1] if $thap == 0;
	}
	my $joint_out;
	my $out_line;
	my $map_line;
	if ($sf == 1 || $thap == 1){
		if ($sf == 1){
			$out_line = "$eles[0] $eles[2] $eles[1] $eles[3] $eles[4] ";
		}
		elsif ($thap == 1){
			$map_line = "$eles[2] $eles[0] $eles[1] $eles[3] $eles[4]";
			$out_line = "";
		}
		$joint_out = join(" ", @out_eles);
		$joint_out =~ s/\/|\|/ /g;
		if ($bi == 1){
			$joint_out =~ s/[2-9]/1/;
		}
		$out_line .= $joint_out;
	}
	else {
		my @fixed = @eles[0..8];
		$joint_out = join("\t", @out_eles);
		if ($bi == 1){
			$joint_out =~ s/[2-9]/1/;
		}
		$out_line = join("\t", @fixed, $joint_out);
	}
	return ($out_line, $map_line);
}

sub cleanup_tmp_chunks {
	return if threads->tid();
	foreach my $chunk_file (@TMP_CHUNKS){
		unlink($chunk_file) if defined $chunk_file && -e $chunk_file;
	}
	@TMP_CHUNKS = ();
}

sub get_ancestral_allele {
	my $line = shift; my @aids = @{$_[-1]};
	my @values;
	my @eles = split(/\t/, $line);
	my @info = split(/\:/, $eles[8]);
	my @nucls = split(/\,/, $eles[4]); #ALT
	unshift(@nucls, $eles[3]); #REF
	my $ad;
	#my $anc_allele = $eles[$aid];
	my @anc_alleles;
	foreach my $n (0..$#aids){
		if ($aids[$n] !~ /[^0-9]/){
			push(@anc_alleles, $eles[$aids[$n]]);
		}
	}
	foreach my $k (0..$#info){
		if ($info[$k] eq "AD"){
			$ad = $k;
		}
	}
	my @anc_nucl_array;
	foreach my $o (0..$#anc_alleles){
		@anc_info = split(/\:/, $anc_alleles[$o]);
		my $anc_number;
		my @alleles = split(/\/|\|/, $anc_info[0]);
		my $anc_nucl;
		if ($ad =~ /[0-9]/){
			my @depths = split(/\,/, $anc_info[$ad]);
			foreach my $j (0..$#depths){
				if ($depths[$j] eq "."){
					$depths[$j] = 0;
				}
			}
			if (scalar(@alleles) == 2){
				my @sorted_depths = sort {$depths[$b] <=> $depths[$a]} 0..$#depths;
				my @sorted_alleles = @alleles[@sorted_depths];
				if ($sorted_depths[0] == 0){
					$anc_number = "N";
				}
				$anc_number = @sorted_alleles[0];
			}
			else {
				$anc_number = $alleles[0];
			}
			if ($anc_number eq "N"){
				$anc_nucl = $anc_number;
			}
			else {
				$anc_nucl = $nucls[$anc_number];
			}
		}
		else {
			$anc_number = $alleles[0];
			$anc_nucl = $nucls[$anc_number];
		}
		push(@anc_nucl_array, $anc_nucl);
	}
	#count numbers of the same elements in the array
	my %counts = ();
	foreach (@anc_nucl_array){
		$counts{$_}++;
	}
	my @keys = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
	my @vals = @counts{@keys};
	return $keys[0];
}

sub rename_chr {
	my $name = shift; my @lists = @{$_[-1]};
	my $ori_name = $name;
	foreach (@lists){
		my @eles = split(/\t+|\s+/, $_);
		if ($name eq $eles[0]){
			$name = $eles[1];
			last;
		}
	}
	if ($name eq $ori_name){
		print RED "$ori_name is unchanged.\n", RESET;
	}
	return $name;
}

sub modify_chr {
	my $name = shift;
	my $ori_name = $name;
	if ($name !~ /\D/){
		if ($name =~ /^0+/){
			$name =~ s/^0+//;
		}
		return $name;
	}
	else {
		$name =~ s/\D+/ /g;
		if ($name =~ /\s+/){
			my @tmps = split(/\s+/, $name);
			$name = $tmps[-1];
			if ($name =~ /^0+/){
				$name =~ s/^0+//;
			}
			if ($name !~ /\d/){
				$name = $tmps[-2];
					if ($name =~ /^0+/){
						$name =~ s/^0+//;
					}
				if ($name !~ /\d/){
					$name = $ori_name;
				}
			}
		}
		#print "debug: Name: $name\n";
		return $name;
	}
}
