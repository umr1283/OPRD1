package Link;

use Data::Dumper;

sub link_fastq {
        my($run,$disque,$zip) = @_;

	#Creation du dossier run
	$target = "/path/$run";
	unless ( -d $target) { system ("mkdir $target"); }

	#Creation du lien symbolique
	if ($zip eq "1"){ # run_presta, lien vers zip
	    $fastq_source = "$disque/$run/demultiplexage/";
	} else { # run labo, lien vers fastq
	    $fastq_source = "$disque/$run/demultiplexage/fastq/";
	}

	$fastq_link = "/path/$run/fastq";
	system ("ln -s $fastq_source $fastq_link");

	#Statistique
	$stat_link = "/path/$run/demultiplexage";
	$stat_source = "$disque/$run/demultiplexage/Reports/html/";

	system ("ln -s $stat_source $stat_link");
}

sub link_bam {
        my($run,$disque) = @_;

        #Creation du dossier run
        $target = "/path/$run";
        unless ( -d $target) { system ("mkdir $target"); }

        #Creation du lien symbolique
        $bam_source = "$disque/$run/alignement/gatk";
        $bam_link = "/path/$run/bam";
        system ("ln -s $bam_source $bam_link");
}

sub link_vcf {
        my($run,$disque) = @_;

        #Creation du dossier run
        $target = "/path/$run";
        unless ( -d $target) { system ("mkdir $target"); }

        #Creation du lien symbolique
        $vcf_source = "$disque/$run/variation/vcf_gatk";
        $vcf_link = "/path/$run/vcf";
        system ("ln -s $vcf_source $vcf_link");
}

sub link_stats {
        my($run,$disque) = @_;

        #Creation du dossier run
        $target = "/path/$run";
        unless ( -d $target) { system ("mkdir $target"); }

        #Creation du lien symbolique
        $report_source = "$disque/$run/Reporting";
        $report_link = "/path/$run/reporting";
        system ("ln -s $report_source $report_link");
}

sub zip_fastq {
    my ($run,$disque) = @_;

    $fastq_folder = "$disque/$run/demultiplexage";
    $out_zip = "$disque/$run/demultiplexage/$run.tar.gz";
    system("tar -czvf $out_zip $fastq_folder/*.fastq.gz");
}

sub link_annot {
    my $run = shift;
    my $list = shift;
    my $racine = shift;

    #print Dumper(\$list);
    my $target = "/path/$run";
    my $snp_target = $target."/SNPs";
    my $indel_target = $target."/INDELs";

    my ($snp_source,$indel_source,$snp_link,$indel_link);
    
    unless ( -d $target) { system ("mkdir $target"); }
    unless ( -d $snp_target) { system ("mkdir $snp_target"); }
    unless ( -d $indel_target) { system ("mkdir $indel_target"); }
    
    foreach my $l (@$list){
#	$snp_source = $racine."/results/SNP_filter/F_basecount/F_basecount_$l";
	$snp_source = $racine."/results/SNP_filter/F_disease_inher/F_disease_inher_$l";
	if (-f $snp_source."_hg19.txt"){
	    $snp_link = $snp_target."/F_basecount_$l\_hg19.txt";
	    system("ln -s $snp_source\_hg19.txt $snp_link");
	} elsif (-f $snp_source."_hg19.txt.gz"){
	    $snp_link = $snp_target."/F_basecount_$l\_hg19.txt.gz";
	    system("ln -s $snp_source\_hg19.txt.gz $snp_link");
	}
	
	if (-f $snp_source."_hg38.txt"){
#	    $snp_link = $snp_target."/F_basecount_$l\_hg38.txt";
	    $snp_link = $snp_target."/F_disease_inher_$l\_hg38.txt";
	    system("ln -s $snp_source\_hg38.txt $snp_link");
	} elsif (-f $snp_source."_hg38.txt.gz"){
#	    $snp_link = $snp_target."/F_basecount_$l\_hg38.txt.gz";
	    $snp_link = $snp_target."/F_disease_inher_$l\_hg38.txt.gz";
	    system("ln -s $snp_source\_hg38.txt.gz $snp_link");
	} 

	
#	$indel_source = $racine."/results/INDEL_filter/F_basecount/F_basecount_$l";
	$indel_source = $racine."/results/INDEL_filter/F_disease_inher/F_disease_inher_$l";
	if (-f $indel_source."_hg19.txt"){
	    $indel_link = $indel_target."/F_basecount_$l\_hg19.txt";
	    system("ln -s $indel_source\_hg19.txt $indel_link");
	}
	
	if (-f $indel_source."_hg38.txt"){
#	    $indel_link = $indel_target."/F_basecount_$l\_hg38.txt";
	    $indel_link = $indel_target."/F_disease_inher_$l\_hg38.txt";
	    system("ln -s $indel_source\_hg38.txt $indel_link");
	}
	
    }
}

sub link_cnv {
    my($run,$disque) = @_;
    
    #Creation du dossier run
    $target = "/path/$run";
    unless ( -d $target) { system ("mkdir $target"); }
    
    #Creation du lien symbolique
    $cnv_source = "$disque/$run/cnv";
    $cnv_link = "/path/$run/cnv";
    system ("ln -s $cnv_source $cnv_link");
}



1;
