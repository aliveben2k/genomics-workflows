#!/usr/bin/perl
#list fo the function: pbs_setting, status, get_1_eles, chr_name, chr_lengths, sample_name, check_path, rnd_str, name_exists_in_dir...... 

=functions for pbs_setting
[-cj_local] [-cj_env PATH] [-cj_conda ENV_NAME] [-cj_node INT] [-cj_ppn INT] [-cj_mem INT] [-cj_qname JOB_NAME] [-cj_proj PROJECT_ID] [-cj_module MODULE] [-cj_docker PATH] [-cj_qout PATH] [-cj_sn SN] [-cj_exc] [-cj_quiet] [-cj_queue QUEUE_NAME] [-cj_mail EMAIL_ADDRESS] [-cj_server FILE]
-cj_local	Run the script locally.
-cj_env		Environment path that need to be set in $PATH. Can be used for multiple times.
-cj_conda	Conda environment name. Only works for h71 and h81 servers.
-cj_node	Node that will be used in the job. Default: 1.
-cj_ppn		Core that will be used in the job. Default: 1.
-cj_time	walltime of the job. Default: system default.
-cj_mem		Momory that will be used in the job. Default: system default.
-cj_gpu		GPU to use. (NARO server only)
-cj_qname	Name of the job, if defined, the job name will be {SN}_{job_name}. Default: {SN}_cj
-cj_queue	Select queue for the job. Only works for Taiwania 1.
-cj_proj	ID of the project. Only works for Taiwania 1.
-cj_module	A module that need to be loaded. Can be used for multiple times.
-cj_docker  docker image path. Only works for NARO's server. 
-cj_qout	The output path where the job execution info should be stored.
-cj_sn		Serial number {SN} of the job. Default: 4-digit random characters.
-cj_mail	E-mail address. It will send the notice when the job starts and is done. Only works for Taiwania 1.
-cj_exc		Send the job for execution.
-cj_quiet	Supress the message.
-cj_server	Server description file. Format: IP_ADDRESS_OR_PREFIX SERVER_TYPE.
=cut

=functions for rnd_str (for old program)
input: none
main function: generate a four-character serial number
=cut

=functions for name_exists_in_dir
input: path($dir), serial_number($ran)
main function: check if the $ran already exists
=cut

=functions for rnd_str2 (for new program)
input: path($dir)
main function: generate a four-character serial number and check if the serial number exists
=cut

our @QUERY_RETRY_SLEEPS = ([30, 60], [60, 120], [120, 180]);
our @SUBMIT_RETRY_SLEEPS = ([60, 120], [120, 180], [180, 240]);
our @QUEUE_BUSY_SLEEPS = ([120, 180], [180, 300], [300, 420]);
our @ACTIVE_JOB_SLEEPS = ([90, 150], [180, 300], [300, 420]);

=functions for status_slurm
input: user_ID, job_file_name, quiet mode, partition
main function: send and track the jobs if the job number reaches the limitation
=cut

=functions for status_other
input: Serial_number($ran), user_ID, job_file_name, quiet mode, partition
main function: send and track the jobs if the job number reaches the limitation (>200)
=cut

=functions for status
input: Serial_number($ran), user_ID
main function: tracking the job status by the serial number given
=cut

=function for get_1_eles
input: array_of_vcf_files, contig_id
function: giving a series of vcf file paths, and reorder them by contig_id order in the vcf
=cut

=function for chr_name
input: a_vcf_file, prefix_name_of_the_contigs
function: get the contig names in the vcf, if the second argument is given, it will return a list only with the prefixed name
=cut

=function for chr_lengths
input: a_vcf_file, prefix_name_of_the_contigs
function: get the contig lengths in the vcf, if the second argument is given, it will return a list of lengths only with the prefixed name
=cut

=function for sample_name
input: a_vcf_file
function: get the sample IDs in the vcf
=cut

=function for check_path
input: a_path(to a file or a folder)
function: change relative path to absolute path
=cut

use Term::ANSIColor qw(:constants);
use Cwd qw(getcwd);

our $SERVER_DESCRIPTION_FILE;

sub _server_type {
	my $type = lc(shift // "");
	$type =~ s/[\s-]+/_/g;
	my %types = (
		pbs => [1, "pbs"],
		pbs_pro => [2, "pbs_pro"],
		ntu => [2, "pbs_pro"],
		slurm => [3, "slurm"],
		taiwania => [3, "slurm"],
		pbs_special => [4, "pbs_special"],
		naro => [4, "pbs_special"],
		sge => [4, "pbs_special"],
		slurm_scion => [5, "slurm_scion"],
		scion => [5, "slurm_scion"],
	);
	return unless exists $types{$type};
	return @{$types{$type}};
}

sub _local_ip_addresses {
	my @ips;
	if (open(my $fh, "-|", "ip", "route", "get", "1.2.3.4")){
		while (my $line = <$fh>){
			push(@ips, $1) if $line =~ /\bsrc\s+(\S+)/;
		}
		close($fh);
	}
	return @ips;
}

sub detect_server {
	my $description_file = shift;
	my $home = (getpwuid $>)[7];
	$description_file ||= $SERVER_DESCRIPTION_FILE;
	$description_file ||= $ENV{"QSUB_SERVER_DESCRIPTION"};
	$description_file ||= "$home/.qsub_server.conf" if -e "$home/.qsub_server.conf";
	$description_file ||= "$home/software/qsub_server.conf" if -e "$home/software/qsub_server.conf";
	$description_file ||= "qsub_server.conf" if -e "qsub_server.conf";

	my @rules;
	if ($description_file){
		open(my $fh, "<", $description_file) or die "Cannot open server description file $description_file: $!\n";
		my $line_number = 0;
		while (my $line = <$fh>){
			$line_number++;
			chomp($line);
			$line =~ s/\r$//;
			$line =~ s/#.*$//;
			$line =~ s/^\s+|\s+$//g;
			next unless length($line);
			my ($ip_pattern, $type, @extra) = split(/[\t,\s]+/, $line);
			next if $ip_pattern =~ /^IP(?:_ADDRESS|_PREFIX)?$/i && $type =~ /^SERVER_TYPE$/i;
			die "Invalid server description at $description_file line $line_number.\n" if !$type || @extra;
			my ($serv, $canonical_type) = _server_type($type);
			die "Unknown server type '$type' at $description_file line $line_number.\n" unless defined $serv;
			push(@rules, [$ip_pattern, $serv, $canonical_type]);
		}
		close($fh);
		die "No server mappings found in $description_file.\n" unless @rules;
		$SERVER_DESCRIPTION_FILE = $description_file;
	}
	else {
		@rules = (
			["172.28.111", 3, "slurm"],
			["150.26.186", 5, "slurm_scion"],
			["140.112.2", 2, "pbs_pro"],
			["150.26.179", 4, "pbs_special"],
		);
	}

	my @ips = _local_ip_addresses();
	foreach my $rule (@rules){
		my ($ip_pattern, $serv, $canonical_type) = @{$rule};
		if ($ip_pattern eq "*"){
			return ($serv, $canonical_type, $ips[0] // "unknown", $description_file);
		}
		foreach my $ip (@ips){
			if ($ip =~ /^\Q$ip_pattern\E(?:\.|$)/){
				return ($serv, $canonical_type, $ip, $description_file);
			}
		}
	}
	my $source = $description_file ? " in $description_file" : "";
	my $detected = @ips ? join(", ", @ips) : "none";
	die "Cannot match local IP address ($detected) to a server type$source.\n";
}

sub pbs_setting {
my $arg = shift;
my @args = split(/\s/, $arg);
my $nodes = 1; #set nodes used for qsub job
my $ppn = 1; #set ppn used for qsub job
my $mem = 0; #set memory used for qsub job; 0 means no memory defined.
my $gpu = 0; #set the GPU usage;
my $home1 = (getpwuid $>)[7]; #get the $HOME path

my $description_file;
for (my $i=0; $i<=$#args; $i++){
	if ($args[$i] eq "-cj_server"){
		die "-cj_server requires a description file.\n" unless $args[$i+1];
		$description_file = $args[$i+1];
		$args[$i] = "";
		$args[$i+1] = "";
		last;
	}
}
my ($serv, $server_type, $server_ip) = detect_server($description_file);
my $scr_dir = '#PBS'; 
my $par = '-q'; #which job class to use
my $snode; #node to use
my $wall = '-l walltime='; #running time 
my $jenv = '-V'; #copy current environment
my $smail = '-M '; #send e-mail
my $mtype = '-m abe'; #event of notification by e-mail
my $jn = '-N ';  #job name
my $sppn; #total task (CPU) to use
my $pname = '-P '; #account to charge (Taiwania)
my $otype = '-j oe';
my $naro_server = 'hostos_c1';
my $docker = '';

if ($serv == 3){ #Slurm system
			$scr_dir = '#SBATCH';
			$par = '-p';
		$snode = '-N ';
		$wall = '-t ';
		$jenv = '--export=ALL';
		$smail = '--mail-user=';
		$mtype = '-–mail-type=ALL';
		$jn = '-J '; #job name
		$sppn = '--ntasks-per-node=';
			$pname = '-A ';
			$otype = '';
}
elsif ($serv == 5){ #Scion-style Slurm system
			$scr_dir = '#SBATCH';
		$par = '--partition';
		$snode = '--nodes=';
		$wall = '--time=';
		$jenv = '--export=ALL';
		$smail = '--mail-user=';
		$mtype = '-–mail-type=ALL';
		$jn = '--job-name '; #job name
		$sppn = '--ntasks-per-node=';
			#$pname = '-A ';
			$otype = '';
}
elsif ($serv == 4){ #NARO-style PBS_special system
			$scr_dir = '#$';
		$par = '-jc';
		$jenv = '-cwd';
		$wall = '-mods l_hard h_rt ';
			$otype = '-j y';
			$jn = '-N ';
}

if ($#args == -1){
	exit;
}

my $query;

my $exc; my $sn; my @envs; my $ran; my $proj; my $mail; my $user_queue; my $home;
my $qname = "cj"; my $conda; my @module; my $quiet; my $local; my $timel; my $docker;
for (my $i=0;$i<=$#args;$i++){
	if ($args[$i] eq "\-cj_exc"){
		$exc = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_quiet"){
		$quiet = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_local"){
		$local = 1;
		$args[$i] = "";
	}
	if ($args[$i] eq "\-cj_sn"){
		if ($args[$i+1] && $args[$i+1] !~ /^\-/){
			$ran = $args[$i+1];
			$args[$i] = "";
			$args[$i+1] = "";
			$sn = 1;
		}
		else {
			exit;
		}
	}
	if ($args[$i] eq "\-cj_env"){
		my $env = 'export PATH='.$args[$i+1].':$PATH';
		push(@envs, $env);
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_conda"){
		$conda = $args[$i+1];
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_node"){
		if ($args[$i+1] !~ /[^0-9]/){
			$nodes = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_ppn"){
		if ($args[$i+1] !~ /[^0-9]/){
			$ppn = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_time"){
		if ($args[$i+1] !~ /[^0-9\:]/){
			$timel = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_mem"){
		if ($args[$i+1] !~ /[^0-9]/){
			$mem = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_gpu"){
		if ($args[$i+1] !~ /[^0-9]/){
			$gpu = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_qname"){
		if ($args[$i+1] =~ /\w/){
			$qname = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_proj"){
		if ($args[$i+1] !~ /[^a-z0-9]/i){
			$proj = $args[$i+1];
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_module"){
		push(@module, "module load $args[$i+1]");
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_qout"){
		if (-d $args[$i+1]){
			if ($args[$i+1] =~ /\/$/){
				$args[$i+1] =~ s/\/$//;
			}
			$home = $args[$i+1];
			if ($home =~ /\/qsub_files$/){
				$home =~ s/\/qsub_files$//;
			}
			unless ($home){
				$home = '.';
			}
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_mail"){
		if ($args[$i+1] =~ /\@/){
			$mail = "$scr_dir $smail$args[$i+1]\n$scr_dir $mtype\n";
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_queue"){
		if ($args[$i+1] =~ /^trans|^ct|^intr|^single|^normal|^large|^ai/){
			$user_queue = $args[$i+1];
		}
		else {
			print "$user_queue is incorrect, skipped.\n";
			$user_queue = "";
		}
		$args[$i] = "";
		$args[$i+1] = "";
	}
	if ($args[$i] eq "\-cj_docker"){
	    $docker = $args[$i+1];
		$args[$i] = "";
		$args[$i+1] = "";	    
	}
}
unless ($home){
	$home = $home1;
}
my $c_line = join(" ", @args);
$c_line =~ s/\s+/ /g;
$c_line =~ s/^\s+//;

unless ($local == 1){
    unless (-d "$home\/qsub_files"){
        system ("mkdir $home\/qsub_files");
    }
    unless (-d "$home\/qsub_files\/out"){
        system ("mkdir $home\/qsub_files\/out");
    }
    my $path = "$home\/qsub_files";
    unless (defined $ran) {
        $ran = rnd_str2($path);
    }
    #check memory
    my $adjust; my $check_ppn; my $check_nodes;
    if ($mem == 0){
        $mem = "";
    }
    else {
        if ($serv == 1 || $serv == 2){ #every core only has 6 gb memory, so if request large memory, adjust core number based on memory value.
            $check_ppn = $mem / 6;
            $check_ppn++ if ($check_ppn > int($check_ppn));
            if ($ppn < int($check_ppn)){
                $ppn = int($check_ppn);
                $adjust = "Adjusted ncpus\/mpiprocs to $ppn based on memory request.\n";
            }
            if ($serv == 1){
                $mem = "$scr_dir \-l mem\=$mem\gb\n";
            }
            if ($serv == 2){
                $mem = "\:mem\=$mem\gb";
            }
            if ($serv == 3 || $serv == 5){
            }
        }
        if ($serv == 3){
            $check_ppn = $mem / 12;
            $check_ppn++ if ($check_ppn > int($check_ppn));
            if ($ppn < int($check_ppn)){
                $ppn = int($check_ppn);
                $adjust = "Adjusted ncpus\/mpiprocs\/ntasks-per-node to $ppn based on memory request.\n";
            }
            if ($mem == 0){
                $mem = "";
            }
            else {
                $mem = "$scr_dir --mem\=$mem\G\n";
            }
        }
        if ($serv == 5){
            if ($user_queue eq "large"){
                #$check_ppn = $mem / 7;
                if ($mem > 1536){
                    $check_nodes = 2;
                    $adjust = "Adjusted $check_nodes to 2 based on memory request.\n";
                    if ($mem > 3072){
                        $adjust .= "Adjusted $mem to 3072, because the request exceed the resources.\n";
                        $mem = 3072;
                    }
                }
            }
            elsif ($user_queue eq "ai"){
                #$check_ppn = $mem / 10;
                if ($mem > 1024){
                    $mem = 1024;
                    $check_nodes = 1;
                    $adjust = "Adjusted $check_nodes to 1 based on memory request.\n";
                    $adjust .= "Adjusted $mem to 1024, because the request exceed the resources.\n";
                }
            }
            else {
                #$check_ppn = $mem / 3.5;
                $check_nodes = $mem / 768;
                $check_nodes++ if ($check_nodes > int($check_nodes));
                $check_nodes = int($check_nodes);
                if ($check_nodes > 2){
                    $check_nodes = 2;
                    $mem = 1536;
                }
            }
            #$check_ppn++ if ($check_ppn > int($check_ppn));
            #if ($ppn < int($check_ppn)){
            #    $ppn = int($check_ppn);
            #    $adjust = "Adjusted ncpus\/mpiprocs\/ntasks-per-node to $ppn based on memory request.\n";
            #}
            if ($mem == 0){
                $mem = "";
            }
            else {
                $mem = "$scr_dir --mem\=$mem\G\n";
            }
        }
    }
    unless (defined $check_nodes){
        $check_nodes = 1;
    }
    #check ppn
    my $check_nodes_tmp;
    if ($serv == 1){
        $check_nodes = $ppn / 12;
    }
    elsif ($serv == 2){
        $check_nodes = $ppn / 20;
    }
    elsif ($serv == 3){ #Taiwania 3
        my @avail_proj = `get_su_balance`; my $check_proj = 0; my $balance;
        chomp(@avail_proj);
        @avail_proj = grep {$_ ne ""} @avail_proj;
        foreach (@avail_proj){
            if ($_ =~ /$proj/){
                $check_proj = 1;
                $_ =~ s/[\{\}\"]//g;
                my @tmp = split(/\,/, $_); #balance is the last one
                chomp(@tmp);
                my @proj_tmp = split(/\:/, $tmp[0]);
                my @bal_tmp = split(/\:/, $tmp[-1]);
                $balance = $proj_tmp[1];
                $proj = $proj_tmp[1];
            }
        }
        if ($check_proj == 0){
            unless ($quiet == 1){
                print "\nProject is not found\/defined, selecting the best project for the job...";
            }
            foreach (@avail_proj){
                $_ =~ s/[\{\}\"]//g;
                my @tmp = split(/\,/, $_);
                chomp(@tmp);
                my @proj_tmp = split(/\:/, $tmp[0]); #project ID is the first one
                my @bal_tmp = split(/\:/, $tmp[-1]); #balance is the last one
                unless ($proj){
                    $proj = $proj_tmp[1];
                    $balance = $bal_tmp[1];
                }
                else {
                    if ($balance < $bal_tmp[1]){
                        $balance = $bal_tmp[1];
                        $proj = $proj_tmp[1];
                    }
                }
            }
        }
        unless ($quiet == 1){
            print BOLD "\nProject: $proj is used. $balance balance is available.\n", RESET;
        }
        if ($balance <= 0){
            die "There is no balance for the project.\n";
            print "Project condition\(s\):\n";
            system("get_su_balance");
            exit;
        }
        elsif ($balance < 10){
            print "WARNING: The available balance for the project is less than 10.\n";
            print "Project condition\(s\):\n";
            system("get_su_balance");
        }
        elsif ($balance < 5){
            print "WARNING: The available balance for the project is less than 5.\n";
            print "Project condition\(s\):\n";
            system("get_su_balance");
        }
        $proj = "$scr_dir $pname$proj\n";
        if ($user_queue){
            $query = "$scr_dir $par $user_queue\n";
        }
        if ($ppn <= 56){
            $query = "$scr_dir $par ct56\n";
        }
        elsif ($ppn <= 224){
            $query = "$scr_dir $par ct224\n";
        }
        elsif ($ppn <= 560){
            $query = "$scr_dir $par ct560\n";
        }
        elsif ($ppn <= 2240){
            $query = "$scr_dir $par ct2k\n";
        }
        else {
            $query = "$scr_dir $par ct8k\n";
        }
        $check_nodes = $ppn / 56; #56 threads in 1 node
    }
    elsif ($serv == 5){
        if ($user_queue){
            $query = "$scr_dir $par $user_queue\n";
            if ($user_queue eq "ai"){
                $check_nodes_tmp = $ppn / 96;
                if ($check_nodes_tmp > 1){
                    $ppn = 96;
                    $check_nodes_tmp = 1;
                }
            }
            else {
                $check_nodes_tmp = $ppn / 192;
                if ($check_nodes_tmp > 2){
                    $ppn = 384;
                    $check_nodes_tmp = 2;
                }
            }
        }
        else {
            $query = "$scr_dir $par normal\n";
            $check_nodes_tmp = $ppn / 192;
            if ($check_nodes_tmp > 2){
                $ppn = 384;
                $check_nodes_tmp = 2;
            }            
        }
    }
    if (defined $check_nodes_tmp){
        if ($check_nodes_tmp > $check_nodes){
            $check_nodes = $check_nodes_tmp;
        }
    }
    $check_nodes++ if ($check_nodes > int($check_nodes));
    if ($nodes < int($check_nodes)){
        $nodes = int($check_nodes);
    }
    if ($serv == 4){ #NARO server
        my $gpu_core = 1; my $cpu_lv = 1; my $mem_lv = 1;
        my @s_resource = (1,1,1); #cpu,gpu,memory
        if ($gpu == 0 && $docker){
            $gpu = 1;
        }
        if ($gpu == 0){
            if ($ppn == 1){}
            elsif ($ppn == 2){$s_resource[0] = 2;}
            elsif ($ppn > 2 && $ppn <= 4){$s_resource[0] = 4;}
            elsif ($ppn > 4 && $ppn <= 8){$s_resource[0] = 8;}
            elsif ($ppn > 8 && $ppn <= 16){$s_resource[0] = 16;}
            elsif ($ppn > 16 && $ppn <= 32){$s_resource[0] = 32;}
            else {
                $nodes = $ppn / 32;
                $nodes++ if ($nodes > int($nodes));
                $naro_server = "parallel_gic";
                $s_resource[0] = 32;
            }
            if ($mem > 0) {
                if ($mem <=14){}
                elsif ($mem > 14 && $mem <= 28){$s_resource[2] = 2;}
                elsif ($mem > 28 && $mem <= 56){$s_resource[2] = 4;}
                elsif ($mem > 56 && $mem <= 112){$s_resource[2] = 8;}
                elsif ($mem > 112 && $mem <= 225){$s_resource[2] = 16;}
                elsif ($mem > 225 && $mem <= 450){$s_resource[2] = 32;}
                elsif ($mem > 450 && $mem <= 1000){$s_resource[2] = 64;}
            }
            if ($naro_server ne "parallel_gic"){
                my @class = sort {$b <=> $a} @s_resource;
                my $lv = 1;
                if ($class[0] == 64){
                    $lv = "large";
                }
                else {
                    $lv = $class[0];
                }
                $naro_server = "hostos_c$lv";
            }
        }
        else { #using GPU
            if ($ppn <= 8){}
            elsif ($ppn > 8 && $ppn <= 36){$s_resource[0] = 4;}
            elsif ($ppn > 36 && $ppn <= 72){$s_resource[0] = 8;}
            else {
                $nodes = $ppn / 72;
                $nodes++ if ($nodes > int($nodes));
                $naro_server = "parallel";
                $s_resource[0] = $ppn;
                if ($nodes > 2){
                    $nodes = 2;
                }
                if ($s_resource[0] > 144){
                    $s_resource[0] = 144;
                }
            }
            if ($gpu == 1){$s_resource[1] = 1;}
            elsif ($gpu > 1 && $gpu <= 4){$s_resource[1] = 4;}
            elsif ($gpu > 4 && $gpu <= 8){$s_resource[1] = 8;}
            else {
                $nodes = $gpu / 8;
                $nodes++ if ($nodes > int($nodes));
                $naro_server = "parallel";
                if ($nodes > 2){
                    $nodes = 2;
                }
                if ($nodes == 1){
                    $s_resource[1] = 8;
                }
                else {
                    $s_resource[1] = 16;
                }
            }
            if ($mem > 0) {
                if ($mem <= 180){$s_resource[2] = 1;}
                elsif ($mem > 180 && $mem <= 720){$s_resource[2] = 4;}
                elsif ($mem > 720 && $mem <= 1440){$s_resource[2] = 8;}
            }
            my $lv;
            if ($naro_server ne "parallel"){
                my @class = sort {$b <=> $a} @s_resource;
                $lv = $class[0];
                $naro_server = "hostos_g$lv";
            }
            if ($docker){
                if ($naro_server eq "parallel"){
                    $naro_server = "docker_g8";
                }
                else {
                    $naro_server = "docker_g$lv";
                }
            }	
        }
    }
    
    my $m_out = $mem;
    if ($m_out eq ""){
        $m_out = "system default";
    } else {
        $m_out =~ s/\:mem\=|$scr_dir \-l mem\=|$scr_dir --mem\=//g;
    }
    
    my $t_out = $timel;
    if ($timel){
        if ($serv == 4){
            my @tmp = split(/\:/, $timel);
            chomp(@tmp);
            if (length($tmp[1]) == 1){
                $tmp[1] = "0".$tmp[1];
                $tmp[1] =~ s/\s+//g;
            }
            if (length($tmp[2]) == 1){
                $tmp[2] = "0".$tmp[2];
                $tmp[2] =~ s/\s+//g;
            }
            $timel = join("\:", @tmp)
        }
        $t_out = "$scr_dir $wall$timel\n";
    }
    else {
        $timel = "system default";
        if ($serv == 4){
            $timel = '320:0:0'; #just set it as 168 hours
            my @tmp = split(/\:/, $timel);
            chomp(@tmp);
            if (length($tmp[1]) == 1){
                $tmp[1] = "0".$tmp[1];
                $tmp[1] =~ s/\s+//g;
            }
            if (length($tmp[2]) == 1){
                $tmp[2] = "0".$tmp[2];
                $tmp[2] =~ s/\s+//g;
            }
            $timel = join("\:", @tmp);
            $t_out = "$scr_dir $wall$timel\n";
        }
    }
    
    my $q_type = $query;
    $q_type =~ s/$scr_dir|$par|\s//g;
    unless ($q_type){
        $q_type = "system default";
    }
    unless ($quiet == 1){
        print BOLD "\nInfo of the job:\n", RESET;
        print "serial number: $ran\njob name: $qname\n";
    }
    if ($serv == 3){
        my $p_name = $proj;
        $p_name =~ s/$scr_dir|-P|-A|\s//gi;
        unless ($quiet == 1){
            print "project name: $p_name\n";
        }
    }
    unless ($quiet == 1){
        if ($serv == 4){
            print "job class: $naro_server\n";
        }
        print "queue name: $q_type\nnode\(select\): $nodes\nppn\(ncpus\/mpiprocs\): $ppn\nmem: $m_out\nwalltime: $timel\n";
        if ($adjust){
            print $adjust;
        }
    }
    
    open (INPUT, ">$home\/qsub_files\/$ran\_$qname\.q") || die BOLD "Cannot write $home\/qsub_files\/$ran\_$qname\.q: $!", RESET, "\n";
	if ($serv == 1){
    	print INPUT "\#\!\/bin\/bash\n$scr_dir \-l nodes\=$nodes\:ppn\=$ppn\n$mem$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out \-j oe\n";
	}
	if ($serv == 2){
    	print INPUT "\#\!\/bin\/bash\n$scr_dir \-l select\=$nodes\:ncpus\=$ppn\:mpiprocs=$ppn$mem\n$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out \-j oe\n$scr_dir $jenv\n";
	}
	if ($serv == 3 || $serv == 5){
		unless ($local == 1){
			my $j_name = "$scr_dir $jn$ran\_$qname\n";
    		print INPUT "\#\!\/bin\/bash\n$proj$query$j_name$scr_dir $snode$nodes\n$scr_dir $sppn$ppn\n$mem$t_out$mail$scr_dir \-o $home\/qsub_files\/out\/$ran\_$qname\.out\n$scr_dir $jenv\n\n";
   	 	}
	}
	if ($serv == 4){
	    my $parti = "";
	    if ($naro eq "parallel"){
	        $parti = "$scr_dir -par 72\n";
	    }
	    if ($docker){
		    $docker = "$scr_dir -ac d\=$docker\n";
		}
	    my $j_name = "$scr_dir $jn$ran\_$qname\n";
		print INPUT "\#\!\/bin\/bash\n$scr_dir -S /bin/bash\n$docker$t_out$j_name$mail$scr_dir $jenv\n$scr_dir $par $naro_server\n$parti$scr_dir $otype\n$scr_dir -o $home\/qsub_files\/out\/$ran\_$qname\.out\n$scr_dir -V\n\n";
	}
}
if (@module){
    if ($local == 1){
    	foreach (@module){
    		system("$_");
    	}
    }
    else {
		print INPUT join("\n", @module), "\n";
	}
}
unless ($local == 1){
	if ($serv == 1 || $serv == 2){
		print INPUT "\ncd \$PBS_O_WORKDIR\n";
	}
}
if ($serv == 5){
    print INPUT "module load miniconda\n";
    print INPUT "\ncd \$SLURM_SUBMIT_DIR\n";
}
if ($serv == 4){
    print INPUT "source /home/chinc518/.bashrc\n";
    print INPUT "source /home/chinc518/miniforge3/etc/profile.d/conda.sh\n";
    unless ($docker){
        print INPUT "source /etc/profile.d/modules.sh\n";
        my $env_ck = `echo \$gatk4 2\>\&1`;
        #if ($env_ck =~ /gatk-4\.2\.2\.0/ || $c_line =~ /canu/){
        #    print INPUT "export JAVA_HOME\=\$HOME\/softwares\/openjdk_8\n";
        #    print INPUT "export PATH\=\$JAVA_HOME\/bin\:\$PATH\n";
        #}
    }
    print INPUT "export LD_LIBRARY_PATH\=/home/chinc518/miniforge3/envs/libcurl/lib:/home/chinc518/softwares/libs\:\$LD_LIBRARY_PATH\n";
    print INPUT "\ncd \$SGE_O_WORKDIR\n";
    
}
if ($conda){
	unless ($quiet == 1){
		print "conda env: $conda\n";
    }
    if ($serv == 1 || $serv == 2){
    	if ($local == 1){
    		system("conda activate $conda");
    	}
    	else {
    		print INPUT "source activate $conda\n";
   		}
    }
    if ($serv == 3 || $serv == 4 || $serv == 5){
    	if ($local == 1){
    		system("conda activate $conda");
    	}
    	else {
    		print INPUT "conda activate $conda\n";
   		}    
    }
}
if (@envs){
	foreach (@envs){
		unless ($quiet == 1){
			print "env setting: $_\n";
		}
		if ($local == 1){
			system("$_");
		}
		else {
			print INPUT "$_\n";
		}
	}
}

my @c_lines;
if ($c_line =~ /\\n/){
	my @c_lines = split(/\\n/, $c_line);
	foreach (@c_lines){
        $_ =~ s/^\s+|\s+$//;
        if ($local == 1){
        	system("$_");
        }
        else {
        	print INPUT "$_\n";
        }
	}
}
else {
	if ($local == 1){
		system("$c_line");
	}
	else {
		print INPUT "$c_line\n";
	}
}
if ($conda  && ($serv == 1 || $serv == 2 || $serv == 3 || $serv == 4 || $serv == 5)){
	if ($local == 1){
		system("conda deactivate");
	}
	else {
    	print INPUT "conda deactivate\n";
    }
}
unless ($local == 1){
	close(INPUT);
	unless ($quiet == 1){
		print BOLD "\nThe qsub file is: $home\/qsub_files\/$ran\_$qname\.q\n", RESET;
	}
}
if ($exc == 1 && $local != 1){
	my @tmp = split(/\//, $home1);
	my $uid = $tmp[-1];
	if ($serv == 3 || $serv == 5){
		status_slurm($uid, "sbatch $home\/qsub_files\/$ran\_$qname\.q", $quiet, $q_type, $serv);
	}
	else {
		status_other($ran, $uid, "qsub $home\/qsub_files\/$ran\_$qname\.q", $quiet, $q_type, $serv);
	}
}
#always delete qsub files and log files older than 90 days
system("find $home\/qsub_files \! -type d -mtime \+90 -exec rm -f \{\} \+");
system("find $home\/qsub_files\/out \! -type d -mtime \+90 -exec rm -f \{\} \+");
return "qsub $home\/qsub_files\/$ran\_$qname\.q";
}

sub rnd_str {
    my @letters = grep { /^[A-Za-z]$/ } @_;
	join("", $letters[rand @letters], @_[map{rand@_} 1..($_[0] - 1)]);
} #generate serial number (for old programs)

sub name_exists_in_dir {
    my ($dir, $target) = @_;
    opendir my $dh, $dir or die "Cannot open $dir: $!";
    while (my $entry = readdir $dh) {
        next if $entry eq '.' or $entry eq '..';

        if ($entry =~ /\Q$target\E/) {
            closedir $dh;
            return 1;
        }
    }
    closedir $dh;
    return 0;
} #check if the serial number exists

sub rnd_str2 {
    my ($path) = @_;
    my $try = 0;
    my $tmp;
    while ($try++ < 1000) {
        $tmp = rnd_str(4, "A".."Z", 0..9);
        return $tmp unless name_exists_in_dir($path, $tmp);
    }
    die "Failed to generate unique run ID";
} #this is for returning a serial number and check existence. (for new programs)

sub scheduler_backoff {
	my ($mode, $attempt) = @_;
	$attempt = 1 unless defined $attempt && $attempt > 0;
	my $idx = $attempt <= 3 ? 0 : ($attempt <= 10 ? 1 : 2);
	my $ranges_ref;
	if ($mode eq "query_retry"){
		$ranges_ref = \@QUERY_RETRY_SLEEPS;
	}
	elsif ($mode eq "submit_retry"){
		$ranges_ref = \@SUBMIT_RETRY_SLEEPS;
	}
	elsif ($mode eq "queue_busy"){
		$ranges_ref = \@QUEUE_BUSY_SLEEPS;
	}
	else {
		$ranges_ref = \@ACTIVE_JOB_SLEEPS;
	}
	my ($low, $high) = @{$ranges_ref->[$idx]};
	return $low + int(rand($high - $low + 1));
}

sub status_slurm {
	my $time = scalar localtime();
	my $uid = shift; my $job_q = shift; my $quiet = shift; my $q_type = shift; my $serv = shift;
	my $job_count; my $check_sent = 0; my $wait_round = 0;# my $core_count; 
	my $home;
	unless (defined $uid){
	    $home = (getpwuid $>)[7]; #get the $HOME path
	    my @tmp = split(/\//, $home);
	    $uid = $tmp[-1];
	}
	unless ($quiet == 1){
		print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
		print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
		print "\[$time\]\: Looking for $uid job\(s\)...\n";
	}
	do {
		$time = scalar localtime();
		my @stat;
		$job_count = 0;
		do {
			@stat = `squeue -u $uid 2\>\&1`;
		} until ($stat[0] !~ /Socket timed out/i);
		foreach (@stat){
			if ($_ =~ /$uid/i && $_ =~ /$q_type\b/){
				$job_count++; #running and padding jobs are counted together
			}
		}
		if ((($q_type eq "ct224" && $job_count >= 75) || ($q_type eq "ct56" && $job_count >= 80) || ($q_type eq "ct560" && $job_count >= 45) || ($q_type eq "ct2k" && $job_count >= 18) || ($q_type eq "ct8k" && $job_count >= 6) || ($q_type eq "trans" && $job_count >= 30)) && $serv == 3){
			unless ($quiet == 1){
				print "\rYour request is over the limitation in the waiting\/running list. Some jobs will be sent later.";
			}
			$wait_round++;
			sleep(scheduler_backoff("queue_busy", $wait_round));
		}
		elsif ($serv == 5 && $job_count >= 80){
			unless ($quiet == 1){
				print "\rYour request is over the limitation in the waiting\/running list. Some jobs will be sent later.";
			}
			$wait_round++;
			sleep(scheduler_backoff("queue_busy", $wait_round));		    
		}
		else {
			my $repeat = 0; my $violate = 0;
			do {
				$violate = 0;
				my $tmp; my $resend = 0;
				do {
					$tmp = `$job_q 2\>\&1`;
					$resend++;
					sleep(scheduler_backoff("submit_retry", $resend)) if ($tmp =~ /submission failed|error/i && $resend < 50);
				} until ($tmp !~ /submission failed|error/i || $resend == 50);
				print "\n$job_q job is sent.\n";
				my $show_job_count = $job_count + 1;
				if ($tmp =~ /violates/){
					$violate = 1;
				}
				else {
					print "$tmp";
				}
				if ($show_job_count == 1){
					print "$show_job_count job is on the list.\n";
				}
				else {
					print "$show_job_count jobs are on the list.\n";
				}
				$repeat++;
			} until ($repeat == 3 || $violate == 0);
			if ($violate == 1){
				print "Job: $job_q violates queue and\/or server resource limits.\n";
				exit;
			}
			$check_sent = 1;
		}
	} while ($check_sent == 0);
	1;
}
sub status_other {
	my $time = scalar localtime();
	my $ran = shift; my $uid = shift; my $job_q = shift; my $quiet = shift; my $q_type = shift; my $serv = shift;
	my $job_count; my $check_sent = 0; my $wait_round = 0; #my $core_count;
	my $home;
	unless (defined $uid){
	    $home = (getpwuid $>)[7]; #get the $HOME path
	    my @tmp = split(/\//, $home);
	    $uid = $tmp[-1];
	}
	unless ($quiet == 1){
		print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
		print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
		print "\[$time\]\: Looking for $ran job\(s\)...\n";
	}
	do {
		$time = scalar localtime();
		my @stat;
        open(my $fh, "-|", "qstat", "-u", $uid) or die "Failed to run qstat: $!";
        while (my $line = <$fh>) {
            chomp $line;
            push @stat, $line;
        }
        close($fh);
		my @temp;
		$job_count = 0;
		foreach (@stat){
		    my $short = substr($uid, 0, 8);
			if ($_ =~ /$short/i){
				$job_count++;
				@temp = split(/\s+/, $_);
			}
		}
		if ($job_count >= 200){
			unless ($quiet == 1){
				print "\rYour jobs are too many. Some jobs will be sent later.";
			}
			if ($serv == 4){
                serv4chgrp();
		    }
			$wait_round++;
			sleep(scheduler_backoff("queue_busy", $wait_round));
		}
		else {
			my $job_cmd = $job_q;
			my $job_file = $job_q;
			$job_file =~ s/^\s*qsub\s+//;
			my $tmp = do {
                open(my $fh, "-|", "qsub", $job_file) or die "Failed to run qsub: $!";
                local $/;
                <$fh>;
            };
            chomp $tmp;
			print "\n$job_cmd job is sent.\n";
			my $show_job_count = $job_count + 1;
			print "$tmp";
			if ($show_job_count == 1){
				print "$show_job_count job is on the list.\n";
			}
			else {
				print "$show_job_count jobs are on the list.\n";
			}
			if ($serv == 4){
                serv4chgrp();
		    }
			$check_sent = 1;
		}
	} while ($check_sent == 0);
	serv4chgrp();
	1;
}
sub serv4chgrp {
    #my @resp = `ls -l`;
    my @resp;
    open(my $fh, "-|", "ls", "-l") or die "Failed to run ls: $!";
    @resp = <$fh>;
    close($fh);
    my $chg = 0; my $main_grp;
    chomp(@resp);
    foreach (@resp){
        if ($_ =~ /^total/){
            next;
        }
        my @grps = split(/\s+|\t+/, $_);
        if ($grps[3] =~ /chinc/){
            $chg = 1;
            next;
        }
        else {
            $main_grp = $grps[3];
        }
    }
    if ($chg == 1 && $main_grp){
        system("chgrp -R $main_grp \*");
    }
    1;
}
sub status {
	my $time = scalar localtime();
	my $ran = shift; my $uid = shift;
	my $job_count; my $wait_round = 0;
	my $home;
	unless (defined $uid){
	    $home = (getpwuid $>)[7]; #get the $HOME path
	    my @tmp = split(/\//, $home);
	    $uid = $tmp[-1];
	}
		my ($configured_serv) = detect_server();
		my $serv;
		if ($configured_serv == 3 || $configured_serv == 5){
			$serv = 1; #Slurm
		}
		elsif ($configured_serv == 4){
			$serv = 3; #NARO PBS_special
		}
		else {
			$serv = 2; #PBS/PBS Pro
		}
	print "\[$time\]\: Press ctrl \+ c to terminate this script.\n";
	print "\[$time\]\: WARNING: If you terminate this script, the following qsub job\(s\)\/step\(s\) will not be generated and executed. If you run this script in background using \"nohup\", you can ignore this message.\n";
	print "\[$time\]\: Looking for $ran job\(s\)...\n";
	do {
		$time = scalar localtime();
		my @stat; my $resend = 0; $err;
			if ($serv == 1){
				if ($uid){
					do {
						open(my $fh, "-|", "squeue", "-u", $uid) or do {
							$err = 1;
							@stat = ();
							last;
						};
						@stat = <$fh>;
						close($fh);
						$err = 0;
						$resend++;
						foreach (@stat){
							if ($_ =~ /error/i){
								$err = 1;
							}
						}
						sleep(scheduler_backoff("query_retry", $resend)) if ($err == 1);
					} until ($err == 0 || $resend == 50);
				}
				else {
					do {
						open(my $fh, "-|", "squeue") or do {
							$err = 1;
							@stat = ();
							last;
						};
						@stat = <$fh>;
						close($fh);
						$err = 0;
						$resend++;
						foreach (@stat){
							if ($_ =~ /error/i){
								$err = 1;
							}
						}
						sleep(scheduler_backoff("query_retry", $resend)) if ($err == 1);
					} until ($err == 0 || $resend == 50);
				}
			}
			elsif ($serv == 2) {
			@stat = do {
                open(my $fh, "-|", "qstat", "-G") or die $!;
                <$fh>;
            };
		}
		else {
		    #@stat = `qstat`;
		    @stat = do {
                open(my $fh, "-|", "qstat") or die $!;
                <$fh>;
            };
		}
		my @temp;
		my @jobs;
		$job_count = 0;
		foreach (@stat){
			if ($uid){
				my $uid_cut;
				if (length($uid) > 8){
					$uid_cut = substr $uid, 0, 8;
				}
				else {
					$uid_cut = $uid;
				}
				if ($_ =~ /$uid_cut/ && $_ =~ /$ran/){
					$job_count += 1;
					@temp = split(/\s+/, $_);
					push(@jobs, $temp[1]);
				}
			}
			elsif ($_ =~ /$ran/){
				$job_count += 1;
				@temp = split(/\s+/, $_);
				push(@jobs, $temp[1]);
			}
		}
		if ($job_count != 0){
			my $job = join("\t", @jobs);
			if ($job_count == 1){
				print "\r\[$time\]\: $job_count "; 
				print "$ran";
				print " job is still running!          ";
			}
			elsif ($job_count > 1){
				print "\r\[$time\]\: $job_count "; 
				print "$ran";
				print " jobs are still running!      ";
			}
			if ($serv == 3){
                serv4chgrp();
		    }
			$wait_round++;
			sleep(scheduler_backoff("active_jobs", $wait_round));
		}
		else {
			print "\r\[$time\]\: no job is running!                           ";
			if ($serv == 3){
                serv4chgrp();
		    }
		}
		@jobs = ();
	} while ($job_count != 0);
	print "\n";
	1;
} #present status of each step
sub get_1_eles{
     my @vcfs = @{$_[0]}; my @ids = @{$_[1]};
     my @ordered;
     foreach $id (@ids){
        foreach $vcf (@vcfs){
        	if ($vcf =~ /gz$/){
            	open(INPUT, "-|", "gzip -dc $vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
            }
            else {
            	open(INPUT, "<$vcf") || die BOLD "Cannot open $vcf: $!", RESET, "\n";
            }
            while (my $line = <INPUT>){
                chomp($line);
                if ($line =~ /^\#/){
                    next;
                }
                else {
                    my @eles = split(/\t+|\s+/, $line);
                    if ($eles[0] eq $id){
                        push(@ordered, $vcf);
                    }
                    last;
                }
            }
            close(INPUT);
        }
    }
    return @ordered;
} #reorder files as the contig order in the vcf file
sub chr_name {
	$time = scalar localtime();
	my $file = shift; my $pre = shift;
	my @content; my @line; my @id;
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
			if ($pre){
                if ($line[3] =~ /^$pre/i){
                	my @tmp = split(/\s+/, $line[3]);
                	$line[3] = $tmp[0];
                    push(@id, $line[3]);
                }
			}
			else {
                my @tmp = split(/\s+/, $line[3]);
                $line[3] = $tmp[0];
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
		@content = `gzip \-cd $vcf \| head \-n 5000`;
	}
	elsif ($vcf =~ /\.vcf$/){
		@content = `head -n 1000 $vcf`;
	}
	foreach (@content){
		if ($_ =~ /\#\#contig\=/){
			@line = split(/\<|\>|\=|\,/, $_);
			if ($pre){
                if ($line[3] =~ /^$pre/i){
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
sub sample_name {
	$time = scalar localtime();
	my $file = shift;
	my @content; my @line;
	if (-e $file){}
	else {
		return 2;
	}
	if ($file =~ /\.vcf\.gz/){
		@content = `gzip \-cd $file \| head \-n 10000`;
	}
	elsif ($file =~ /\.vcf$/){
		@content = `head -n 10000 $file`;
	}
	else {
		return 2;
	}
	foreach (@content){
		if ($_ =~ /\#CHROM/){
			@line = split(/\t/, $_);
			for (my $i=9; $i<=$#line; $i++){
                push(@samples, $line[$i]);
			}
			last;
		}
	}
	chomp(@samples);
	return (@samples);
} #get sample name
sub check_path {
	my $path = shift;
	my $dir = getcwd;
	if ($path =~ /\// || $path =~ /\.\./){
		my @path_eles = split(/\//, $path);
		my @dir_eles = split(/\//, $dir);
		if ($path_eles[0] eq "."){
			$path =~ s/^.//;
			$path = "$dir$path";
		}
		else {
			my $dot_cnt = -1;
			foreach (@path_eles){
				if ($_ eq ".."){
					$dot_cnt++;
				}
			}
			if ($dot_cnt == -1){
				if (-e "$dir\/$path" || -d "$dir\/$path"){
					$path = "$dir\/$path";
				}
			}
			else {
				for (my $i=0; $i<=$dot_cnt; $i++){
					shift(@path_eles);
					pop(@dir_eles);
				}
				$path = join("\/", @dir_eles)."\/".join("\/", @path_eles);
			}
		}
	}
	else {
		if ($path eq "."){
			$path = "$dir";
		}
		else {
			$path = "$dir\/$path";
		}
	}
	if ($path =~ /\/$/){
        $path =~ s/\/$//;
	}
	return($path);
} #relative path to absolute path

1;
