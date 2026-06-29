#!/usr/bin/perl
#Usage: perl create_job.pl COMMAND_LINE [-cj_help] [-cj_env PATH] [-cj_conda ENV_NAME] [-cj_node INT] [-cj_ppn INT] [-cj_mem INT] [-cj_qname JOB_NAME] [-cj_queue QUEUE_NAME] [-cj_proj PROJECT_ID] [-cj_module MODULE] [-cj_docker PATH] [-cj_mail EMAIL_ADDRESS] [-cj_qout PATH] [-cj_sn SN] [-cj_exc]

=functions for pbs_setting
[-cj_local] [-cj_env PATH] [-cj_conda ENV_NAME] [-cj_node INT] [-cj_ppn INT] [-cj_mem INT] [-cj_qname JOB_NAME] [-cj_proj PROJECT_ID] [-cj_module MODULE] [-cj_qout PATH] [-cj_sn SN] [-cj_exc] [-cj_quiet] [-cj_queue QUEUE_NAME] [-cj_mail EMAIL_ADDRESS]
-cj_local	Run the script locally.
-cj_env		Environment path that need to be set in $PATH. Can be used for multiple times.
-cj_conda	Conda environment name. Only works for h71 and h81 servers.
-cj_node	Node that will be used in the job. Default: 1.
-cj_ppn		Core that will be used in the job. Default: 1.
-cj_time	walltime of the job. Default: system default.
-cj_mem		Momory that will be used in the job. Default: system default.
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
=cut

use Cwd qw(getcwd);
use FindBin;
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


chomp(@ARGV);
print "This script is written by Ben Chien. Apr.2022\n";

my @server = `ip route get 1.2.3.4 \| awk \'\{print \$7\}\'`;
chomp(@server);
my $serv = -1;
foreach (@server){
	if ($_ =~ /140.112.2/){
		$serv = 2;
		print BOLD "NTU server is detected. Set to the server\'s PBS job setting.\n\n", RESET;
	}
    if ($_ =~ /172.28.111/){
		$serv = 3;
		print BOLD "Taiwania server is detected. Set to the server\'s Slurm job setting.\n\n", RESET;
	}
	if ($_ =~ /150.26.179/){
		$serv = 4;
		print BOLD "Shiho server is detected. Set to the server\'s PBS_special job setting.\n\n", RESET;
	}
}
if ($serv == -1){
	print "Cannot find the right server.\n";
	exit;
}
&usage;
if ($#ARGV == -1){
	exit;
}
print "Input command line:\nperl create_job.pl @ARGV\n";

my $exc; my @envs; my $ran; my $proj; my $mail; my $query; my $user_queue;
my $qname; my $conda; my @module; my $qout; my $local;
for (my $i=0;$i<=$#ARGV;$i++){
	if ($ARGV[$i] eq "\-cj_exc"){
		$exc = "-cj_exc ";
		$ARGV[$i] = "";
	}
	if ($ARGV[$i] eq "\-cj_local"){
		$local = "-cj_local ";
		$ARGV[$i] = "";
	}
	if ($ARGV[$i] eq "\-cj_sn"){
		if ($ARGV[$i+1] && $ARGV[$i+1] !~ /^\-/){
			$ran = "-cj_sn $ARGV[$i+1] ";
			$ARGV[$i] = "";
			$ARGV[$i+1] = "";
		}
		else {
			&usage;
			exit;
		}
	}
	if ($ARGV[$i] eq "\-cj_env"){
		push(@envs, "-cj_env $ARGV[$i+1] ");
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
#		splice(@ARGV, $i, 1);
	}
	if ($ARGV[$i] eq "\-cj_help"){
		&usage;
		&help;
		exit;	
	}
	if ($ARGV[$i] eq "\-cj_conda"){
		$conda = "-cj_conda $ARGV[$i+1] ";
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_node"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$nodes = "-cj_node $ARGV[$i+1] ";
		}
		else {
			print "-cj_node value is wrong, use 1 as default.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_ppn"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$ppn = "-cj_ppn $ARGV[$i+1] ";
		}
		else {
			print "-cj_ppn value is wrong, use 1 as default.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_mem"){
		if ($ARGV[$i+1] !~ /[^0-9]/){
			$mem = "-cj_mem $ARGV[$i+1] ";
		}
		else {
			print "-cj_mem value is wrong, use default setting.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_qname"){
		if ($ARGV[$i+1] =~ /\w/){
			$qname = "-cj_qname $ARGV[$i+1] ";
		}
		else {
			print "-cj_qname value is wrong, use \"cj\" as default.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_proj"){
		if ($ARGV[$i+1] !~ /[^a-z0-9]/i){
			$proj = "-cj_proj $ARGV[$i+1] ";
		}
		else {
			print "-cj_proj value is wrong, skipped.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_module"){
		push(@module, "-cj_module $ARGV[$i+1] ");
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_qout"){
		if (-d $ARGV[$i+1]){
			if ($ARGV[$i+1] =~ /\/$/){
				$ARGV[$i+1] =~ s/\/$//;
			}
			$qout = "-cj_qout $ARGV[$i+1] ";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_mail"){
		if ($ARGV[$i+1] =~ /\@/){
			$mail = "-cj_mail $ARGV[$i+1] ";
		}
		else {
			print "E-mail format is incorrect, skipped.\n";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_queue"){
		if ($ARGV[$i+1] =~ /^trans|^ct/){
			$user_queue = "-cj_queue $ARGV[$i+1] ";
		}
		else {
			print "$user_queue is incorrect, skipped.\n";
			$user_queue = "";
		}
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
	if ($ARGV[$i] eq "\-cj_docker"){
		$docker = "-cj_docker $ARGV[$i+1] ";
		$ARGV[$i] = "";
		$ARGV[$i+1] = "";
	}
}

my $c_line = join(" ", @ARGV);
$c_line =~ s/^\s+|\s+$//;

&pbs_setting("$proj$local@envs$conda$docker$ppn$mem$nodes$ran$mail$qout$query$user_queue$qname@module$exc$c_line");
#print "debug1: $proj$local@envs$conda$docker$ppn$mem$nodes$ran$mail$qout$query$user_queue$qname@module$exc$c_line\n";

sub usage {
	print BOLD "Usage: perl create_job.pl COMMAND_LINE [-cj_help] [-cj_env PATH] [-cj_conda ENV_NAME] [-cj_node INT] [-cj_ppn INT] [-cj_mem INT] [-cj_qname JOB_NAME] [-cj_queue QUEUE_NAME] [-cj_proj PROJECT_ID] [-cj_module MODULE] [-cj_qout PATH] [-cj_sn SN] [-cj_mail EMAIL_ADDRESS] [-cj_exc]\n", RESET;
	print "This script can detect which server you are using and use its setting.\n";
	print "Only gc3, h71, h81 and Taiwania 1 servers are supported.\n";
	print "-cj_conda only supports h71 and h81 servers.\n";
	print "-cj_env supports all servers, and conda env path can be set here for gc3 and Taiwania 1.\n";
	print "If -cj_qname is set, job name is random serial number + job name. Default job name is \"cj\".\n";
	print "-cj_queue only works for Taiwania 1.\n";
	print "-cj_proj only works for Taiwania 1.\n";
	print "-cj_module only works for Taiwania 1.\n";
	print "If you have multiple command lines, please use \'\' to quote your command lines, lines are separated by \"\\n\".\n";
	print "For Taiwania 1, the script will select the best queue automatically. Unless you use the option -cj_queue.\n";
	print "If -cj_proj is not set for Taiwania 1, the script will select the best project.\n";
	print "Please use -cj_help to see the details.\n";
	return;
}
sub help {
	print "-cj_env\t\tEnvironment path that need to be set in \$PATH. Can be used for multiple times.\n";
	print "-cj_conda\tConda environment name.\n";
	print "-cj_node\tNode that will be used in the job. Default: 1.\n";
	print "-cj_ppn\t\tCore that will be used in the job. Default: 1.\n";
	print "-cj_mem\t\tMomory that will be used in the job. Default: system default.\n";
	print "-cj_qname\tName of the job. Default: \{SN\}_cj\n";
	print "-cj_queue\tSelect queue for the job.\n";
	print "-cj_proj\tID of the project.\n";
	print "-cj_module\tA module that need to be loaded. Can be used for multiple times.\n";
	print "-cj_docker\tA docker image path. Only works for NARO's server.\n";
	print "-cj_qout\tThe output path where the job execution info should be stored.\n";
	print "-cj_sn\t\tSerial number \{SN\} of the job. Default: 4-digit random characters.\n";
	print "-cj_mail\tE-mail address. It will send the notice when the job starts and is done.\n";
	print "-cj_exc\t\tSend the job for execution.\n";
}