package LogOut;

use File::Spec;
use Projet;

my %config = Projet::load_config();
my $general_log_folder = $config{racine}."/log";

sub printLogByChr {
    $log_message = shift;
    $subroutine_caller = shift;
    $run_name = shift;
    $sample_name = shift;
    $chr = shift;
    $run_log_folder = $general_log_folder."/".$run_name;
    mkdir($run_log_folder);
    
    if ($subroutine_caller =~ /([\w\W_]+)::([\w\W_]+)/){
	$subroutine = $2;
    }
    $current_log_file = $run_log_folder."/".$sample_name."_".$chr."_".$subroutine.".log";

    open (LOG, ">>$current_log_file") || die "cannot open $current_log_file log file:$!\n";
    print LOG $log_message."\n";
    close LOG;
}

sub printLog {
    $log_message = shift;
    $subroutine_caller = shift;
    $run_name = shift;
    $sample_name = shift;
    $run_log_folder = $general_log_folder."/".$run_name;
    mkdir($run_log_folder);
    
    if ($subroutine_caller =~ /([\w\W_]+)::([\w\W_]+)/){
	$subroutine = $2;
    }
    $current_log_file = $run_log_folder."/".$sample_name."_".$subroutine.".log";

    open (LOG, ">>$current_log_file") || die "cannot open $current_log_file log file:$!\n";
    print LOG $log_message."\n";
    close LOG;
}

sub printLogStat {
    $log_message = shift;
    $subroutine_caller = shift;
    $run_name = shift;
    $log_hash = shift;
    $run_log_folder = $general_log_folder."/".$run_name;
    mkdir($run_log_folder);
    
    if ($subroutine_caller =~ /([\w\W_]+)::([\w\W_]+)/){
	$subroutine = $2;
    }
    $current_log_file = $run_log_folder."/".$run_name."_".$subroutine.".log";

    open (LOG, ">$current_log_file") || die "cannot open $current_log_file log file:$!\n";
    print LOG $log_message."\n";
    foreach my $keys (%$log_hash){
	print LOG "***** $keys *****\n$$log_hash{$keys}\n\n";
    }
    close LOG;
}


1;
