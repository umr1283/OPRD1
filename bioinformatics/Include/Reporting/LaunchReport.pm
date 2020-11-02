package LaunchReport;

use Reporting::FastQC;
use Reporting::Stats;
use Reporting::Html;
use Reporting::Tables;
use Reporting::PDF;
use Data::Dumper;

sub launch_all_reports {
    my $dna_list = shift; # array list of samples
    my $run_data = shift; # object containing samples info
    my $fastqc_script = shift;
    my $samtools_script = shift;
    my $result_folder = shift; # results folder in anges
    my $run_disk = shift; # disk containing the run data
    my $run_name = shift;
    my $run_nas = shift;
    my $run_from = shift;
    my $bed_type = shift; # normal bed = 0, alternative bed = 1, whole genome = 2
    my $alt_bed = shift;
    my $pdf = shift;
    my $sexDetscript = shift;
    my ($sample_report_folder);

    my $report_folder = $run_disk."/".$run_name."/Reporting";
    unless (-d $report_folder){system("mkdir $report_folder")};
    
    #### for each sample
    for my $dna (@{$dna_list}){
	
        ## retrieve coverage, duplicates and variations stats
        Stats::recup_all_stats($dna, $$run_data{$dna}, $sample_report_folder, $result_folder, $bed_type, $run_name, $alt_bed);

    }
   
    ## retrieve duplicates stats
    if ($bed_type ne "1"){
	Stats::duplicates_stats($dna_list, $run_data, $report_folder, $run_name, $run_disk, $run_nas);
    }

    ## export stats results into csv and xls
    if($bed_type ne "1"){
	Tables::build_tab_files($dna_list, $run_data, $report_folder, $result_folder, $run_name, $run_nas, $bed_type, $samtools_script, $sexDetscript);
    } else {
	Tables::stat_cov_in_alt_bed($dna_list, $run_data, $report_folder, $result_folder, $run_name, $run_nas, $alt_bed, $bed_type, $samtools_script, $sexDetscript);
    }

    ## when files are generated, load stats into sql table
    if ($bed_type ne "1"){
	Stats::stats_to_sql($run_disk, $run_name);
    }
    ## general html of the run
    # if Hiseq4000 or NextSeq
    #if ($bed_type ne "1"){
	    #Html::run_html_hiseq4000($dna_list, $run_data, $report_folder, $run_disk, $run_name, $result_folder, $run_from);
    #}
    
    ### PDF
    if ($pdf eq "1"){
	PDF::generate_pdf($run_disk, $run_name, $dna_list);
    }
    
}

1;
