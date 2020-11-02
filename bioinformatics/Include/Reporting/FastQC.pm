package FastQC;

sub launch_FastQC {
    my $dna = shift;
    my $run_disk = shift;
    my $run_name = shift;
    my $sample_report_folder = shift;
    my $fastqc_script = shift;
    my $report_folder = shift;
    my $run_nas = shift;
    my $demultiplex_folder;

    #$dna=~s/_/-/g;
    my $FastQC_folder = $sample_report_folder."/FastQC_report";
    unless (-d $FastQC_folder) { system("mkdir $FastQC_folder");}
    
    ## peu importe la capture, on va chercher les fastQ dans le dossier demultiplexage du run
    #if (($run_nas =~/Hiseq3/)||($run_nas =~/NEXTSEQ/)){
    $demultiplex_folder = $run_disk."/".$run_name."/demultiplexage/fastq";
    system ("$fastqc_script $demultiplex_folder/$dna\_R1.fastq.gz --outdir=$FastQC_folder -t 6 > $FastQC_folder/fastqc.log" );
    system ("$fastqc_script $demultiplex_folder/$dna\_R2.fastq.gz --outdir=$FastQC_folder -t 6 > $FastQC_folder/fastqc.log" );
    #} else {
#	$demultiplex_folder = $run_disk."/".$run_name."/demultiplexage/Project_FC/Sample_".$dna;
#	system ("$fastqc_script $demultiplex_folder/\*$dna\*.gz --outdir=$FastQC_folder --casava -t 6 > $FastQC_folder/fastqc.log" );
#    }
    
    
}

1;
