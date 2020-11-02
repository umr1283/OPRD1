package Rnaseq;

use Data::Dumper;

sub launch_rnaseq {
    my $run_disk = shift;
    my $run_name = shift;
    my $run = shift;
    my $config = shift;
    my $dna_list = shift;
    
    my ($sample_folder, $tophat_out, $cufflinks_out, $hisat_out, $host, $dataset);
    
    my $run_folder = "$run_disk/$run_name";
    my $genome = $$run{reference};
    my $ref_genome = $$run{genome_RNA};
    my $ref_transcriptome = $$run{transcriptome};
    my $mrna_gtf = $$run{mrna};
    my $rrna_gtf = $$run{rrna};
    my $tophat_path = $$config{tophat2};
    my $cufflinks_path = $$config{cufflink2};
    my $trim_galore = $$config{trim_galore};
    my $cutadapt = $$config{cutadapt};
    my $fastqc_script = $$config{fastqc};
    my $bowtie2_folder = $$config{bowtie2};
    my $bowtie2 = $bowtie2_folder."/bowtie2";
    my $RSeQC = $$config{rseqc};
    my $samtools = $$config{samtools};
    my $hisat2 = $$config{hisat2};
    my $ref_hisat = $$run{genome_hisat};
    my $ref_STAR = $$run{genome_STAR};
    my $mrna_STAR = $$run{mrna_STAR};
    my $STAR = $$config{STAR};
    my $ref_RSEM = $$run{reference_RSEM};
    my $typerun = $$run{typerun};
    my $RSEM_addgenename = $$config{racine}."/Include/rnaseq_add_genename.R";

    my $ngsbioinfo = "/path";
    my $link_folder = $ngsbioinfo."/$run_name";
    unless (-d $link_folder) {system("mkdir $link_folder");}
    unless (-d "$link_folder/analysis") {system("mkdir $link_folder/analysis");}

    my $fastq_folder = "$run_disk/$run_name/demultiplexage/fastq";
    system("ln -s $fastq_folder $link_folder/fastq");
    system("ln -s $run_folder/demultiplexage/Reports/html $link_folder/demultiplexage");
    
    my $rnaseq_folder;
    if ($$run{trim} eq "1"){
	$rnaseq_folder = "$run_disk/$run_name/rnaseq_analysis_trim";
    } else {
	$rnaseq_folder = "$run_disk/$run_name/rnaseq_analysis_noTrim";
    }
    unless (-d $rnaseq_folder) {system("mkdir $rnaseq_folder");}
    my $bam_folder = "$rnaseq_folder/analysis";
    unless (-d $bam_folder) {system("mkdir $bam_folder");}
    
    for my $adn (@$dna_list){
	print "RNAseq analysis of $adn...\n";
	$sample_folder = $rnaseq_folder."/Sample_$adn";
	unless (-d $sample_folder) {system("mkdir $sample_folder");}
	
	## unzip fastq
	system("cp $fastq_folder/$adn\_*.fastq.gz $sample_folder");
	

### with trim
	if ($$run{trim} eq "1"){
## trim adaptors
	    system("gunzip $sample_folder/*.gz");
	    if ($typerun eq "PE"){
		system("$trim_galore --illumina -o $sample_folder --clip_R1 12 --clip_R2 12 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --paired --length 10 --path_to_cutadapt $cutadapt $sample_folder/*R1.fastq $sample_folder/*R2.fastq");
	    } elsif ($typerun eq "SR"){
		system("$trim_galore --illumina -o $sample_folder --clip_R1 12 --clip_R2 12 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --paired --length 10 --path_to_cutadapt $cutadapt $sample_folder/*R1.fastq");
	    }
	    # launch tophat
	    if ($$run{mapping} eq "tophat2"){
		print "tophat mapping\n";
		$tophat_out = "$sample_folder/tophat";
		system("mkdir $tophat_out");
		if ($typerun eq "PE"){
		    system("$tophat_path -p 8 --fusion-search -o $tophat_out $ref_genome/$genome $sample_folder/*R1_val_1.fq $sample_folder/*R2_val_2.fq > $tophat_out/tophat.log\n");
		    # -G transcriptome
		} elsif ($typerun eq "SR"){
		    system("$tophat_path -p 8 --fusion-search -o $tophat_out $ref_genome/$genome $sample_folder/*R1_val_1.fq > $tophat_out/tophat.log\n");
		}
	    }
	    
            # launch Hisat2
	    elsif ($$run{mapping} eq "hisat2"){
		print "hisat mapping\n";
		$hisat_out = "$sample_folder/hisat";
		system("mkdir $hisat_out");
		if ($typerun eq "PE"){
		    system("$hisat2 -q -t -p 8 --reorder -x $ref_hisat -1 $sample_folder/*R1_val_1.fq -2 $sample_folder/*R2_val_2.fq -S $hisat_out/$adn.sam > $hisat_out/hisat.log 2>&1");
		} elsif ($typerun eq "SR"){
		    system("$hisat2 -q -t -p 8 --reorder -x $ref_hisat -U $sample_folder/*R1_val_1.fq -S $hisat_out/$adn.sam > $hisat_out/hisat.log 2>&1");
		}
		system("$samtools view -Sb $hisat_out/$adn.sam > $hisat_out/$adn.bam");
		system("rm $hisat_out/$adn.sam");
		system("$samtools sort $hisat_out/$adn.bam $hisat_out/$adn.sorted");
	    }

	    # launch STAR
	    elsif ($$run{mapping} eq "STAR"){
		print "STAR mapping\n";
		my $star_out = "$sample_folder/STAR";
		system("mkdir $star_out");
		chdir($star_out);
		if ($typerun eq "PE"){
		    system("$STAR/STAR --runMode alignReads --runThreadN 25 --genomeDir $ref_STAR --readFilesIn $sample_folder/*R1_val_1.fq $sample_folder/*R2_val_2.fq --sjdbGTFfile $mrna_STAR --outFileNamePrefix $adn --readFilesCommand zcat --outTmpDir tmp --outReadsUnmapped Fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic");
		} elsif ($typerun eq "SR"){
		    system("$STAR/STAR --runMode alignReads --runThreadN 25 --genomeDir $ref_STAR --readFilesIn $sample_folder/*R1_val_1.fq --sjdbGTFfile $mrna_STAR --outFileNamePrefix $adn --readFilesCommand zcat --outTmpDir tmp --outReadsUnmapped Fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic");
		}
	    }
	    
	    # launch RSEM/STAR
	    elsif ($$run{mapping} eq "RSEM"){
		print "RSEM mapping\n";
		my $rsem_out = "$sample_folder/RSEM";
		system("mkdir $rsem_out");
		if ($typerun eq "PE"){
		    system("$$config{RSEM}/rsem-calculate-expression -p 20 --star --star-path $STAR --star-gzipped-read-file --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end $sample_folder/*R1_val_1.fq $sample_folder/*R2_val_2.fq $ref_RSEM $rsem_out/$adn");
		} elsif ($typerun eq "SR"){
		    system("$$config{RSEM}/rsem-calculate-expression -p 20 --star --star-path $STAR --star-gzipped-read-file --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end $sample_folder/*R1_val_1.fq $ref_RSEM $rsem_out/$adn");
		}
		if ($genome eq "hg38"){
		    $host = "www.ensembl.org";
		    $dataset = "hsapiens_gene_ensembl";
		} elsif ($genome eq "hg19"){
		    $host = "feb2014.archive.ensembl.org";
		    $dataset = "hsapiens_gene_ensembl";
		} elsif ($genome eq "mm10"){
		    $host = "may2017.archive.ensembl.org";
		    $dataset = "mmusculus_gene_ensembl";
		} elsif ($genome eq "Sscrofa11"){
		    $host = "www.ensembl.org";
		    $dataset = "sscrofa_gene_ensembl";
		} else {
		    print "!!!Warning : cannot add gene names : reference unknown (different from hg19, hg38, mm10 or Sscrofa11)\n";
		    exit;
		}
		system("Rscript --vanilla $RSEM_addgenename $rsem_out/$adn $host $dataset");

		system("ln -s $rsem_out/$adn.genes_genename.results $link_folder/analysis/$adn.genes_genename.results");
		system("ln -s $rsem_out/$adn.isoforms_genename.results $link_folder/analysis/$adn.isoforms_genename.results");
	    }

            # launch FastQC
	    my $fastQC_folder = $sample_folder."/FastQC_report";
	    system("mkdir $fastQC_folder");
	    system("$fastqc_script $sample_folder/*R1_val_1.fq $sample_folder/*R2_val_2.fq --outdir=$fastQC_folder -t 30 > $fastQC_folder.log");
	    
            # launch mapping VS ncRNA
	    system("$bowtie2 -x $$run{ncrna} -U $sample_folder/$adn\_R1_val_2.fq -S $sample_folder/$adn\_R1.NCRNA.sam --un $sample_folder/$adn\_R1.minusNCRNA.fastq > $sample_folder/$adn\_R1.align_summary.txt 2>&1");
	    system("rm $sample_folder/$adn\_R1.NCRNA.sam");
	    if ($typerun eq "PE"){
		system("$bowtie2 -x $$run{ncrna} -U $sample_folder/$adn\_R2_val_2.fq -S $sample_folder/$adn\_R2.NCRNA.sam --un $sample_folder/$adn\_R2.minusNCRNA.fastq > $sample_folder/$adn\_R2.align_summary.txt 2>&1");
		system("rm $sample_folder/$adn\_R2.NCRNA.sam");
	    }
	    ###system("$$config{samtools} view -Sb $sample_folder/$adn\_R1.NCRNA.sam > $sample_folder/$adn\_R1.NCRNA.bam");
	    ###system("$$config{samtools} view -Sb $sample_folder/$adn\_R2.NCRNA.sam > $sample_folder/$adn\_R2.NCRNA.bam");
	    system("gzip $sample_folder/*.fastq");

	    #### create Rscript file and launch Rsubread
	    
	} else {
#### no trim
	    # launch tophat
	    if ($$run{mapping} eq "tophat2"){
		print "tophat mapping\n";
		$tophat_out = "$sample_folder/tophat";
		system("mkdir $tophat_out");
		if ($typerun eq "PE"){
		    system("$tophat_path -p 8 --fusion-search -o $tophat_out $ref_genome/$genome $fastq_folder/$adn\_*R1.fastq.gz $fastq_folder/$adn\_*R2.fastq.gz > $tophat_out/tophat.log\n");
		} elsif ($typerun eq "SR"){
		    system("$tophat_path -p 8 --fusion-search -o $tophat_out $ref_genome/$genome $fastq_folder/$adn\_*R1.fastq.gz > $tophat_out/tophat.log\n");
		}
	    }

	    # launch Hisat2
	    elsif ($$run{mapping} eq "hisat2"){
		print "hisat mapping\n";
		$hisat_out = "$sample_folder/hisat";
		system("mkdir $hisat_out");
		if ($typerun eq "PE"){
		    system("$hisat2 -q -t -p 8 --reorder -x $ref_hisat -1 $fastq_folder/$adn\_*R1.fastq.gz -2 $fastq_folder/$adn\_*R2.fastq.gz -S $hisat_out/$adn.sam > $hisat_out/hisat.log 2>&1");
		} elsif ($typerun eq "SR"){
		    system("$hisat2 -q -t -p 8 --reorder -x $ref_hisat -U $fastq_folder/$adn\_*R1.fastq.gz -S $hisat_out/$adn.sam > $hisat_out/hisat.log 2>&1");
		}
		system("$samtools view -Sb $hisat_out/$adn.sam > $hisat_out/$adn.bam");
		system("rm $hisat_out/$adn.sam");
		system("$samtools sort $hisat_out/$adn.bam $hisat_out/$adn.sorted");
            }
	    
	    # launch STAR
	    elsif ($$run{mapping} eq "STAR"){
		print "STAR mapping\n";
		my $star_out = "$sample_folder/STAR";
		system("mkdir $star_out");
		chdir($star_out);
		# paired-end
		if ($typerun eq "PE"){
		    system("$STAR/STAR --runMode alignReads --runThreadN 25 --genomeDir $ref_STAR --readFilesIn $fastq_folder/$adn\_R1.fastq.gz $fastq_folder/$adn\_R2.fastq.gz --sjdbGTFfile $mrna_STAR --outFileNamePrefix $adn --readFilesCommand zcat --outTmpDir tmp --outReadsUnmapped Fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic");
		} elsif ($typerun eq "SR"){
		# single end
		    system("$STAR/STAR --runMode alignReads --runThreadN 25 --genomeDir $ref_STAR --readFilesIn $fastq_folder/$adn\_R1.fastq.gz --sjdbGTFfile $mrna_STAR --outFileNamePrefix $adn --readFilesCommand zcat --outTmpDir tmp --outReadsUnmapped Fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic");
		}
	    }

            # launch RSEM/STAR
	    elsif ($$run{mapping} eq "RSEM"){
		print "RSEM mapping\n";
		my $rsem_out = "$sample_folder/RSEM";
		system("mkdir $rsem_out");
		# paired end
		if ($typerun eq "PE"){
		    system("$$config{RSEM}/rsem-calculate-expression -p 10 --star --star-path $STAR --star-gzipped-read-file --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end $fastq_folder/$adn\_R1.fastq.gz $fastq_folder/$adn\_R2.fastq.gz $ref_RSEM $rsem_out/$adn");
		} elsif ($typerun eq "SR"){
		    # single end
		    system("$$config{RSEM}/rsem-calculate-expression -p 10 --star --star-path $STAR --star-gzipped-read-file --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files $fastq_folder/$adn\_R1.fastq.gz $ref_RSEM $rsem_out/$adn");
		}
		if ($genome eq "hg38"){
		    $host = "www.ensembl.org";
		    $dataset = "hsapiens_gene_ensembl";
		} elsif ($genome eq "hg19"){
		    $host = "feb2014.archive.ensembl.org";
		    $dataset = "hsapiens_gene_ensembl";
		} elsif ($genome eq "mm10"){
		    $host = "may2017.archive.ensembl.org";
		    $dataset = "mmusculus_gene_ensembl";
		} elsif ($genome eq "Sscrofa11"){
		    $host = "www.ensembl.org";
		    $dataset = "sscrofa_gene_ensembl";
		} else {
		    print "!!!Warning : cannot add gene names : reference unknown (different from hg19, hg38, mm10 or Sscrofa11)\n";
		    exit;
		}
		system("Rscript --vanilla $RSEM_addgenename $rsem_out/$adn $host $dataset");

		system("ln -s $rsem_out/$adn.genes_genename.results $link_folder/analysis/$adn.genes_genename.results");
		system("ln -s $rsem_out/$adn.isoforms_genename.results $link_folder/analysis/$adn.isoforms_genename.results");
		### add links to results files
	    }

            # launch fastQC
	    my $fastQC_folder = $sample_folder."/FastQC_report";
	    system("mkdir $fastQC_folder");
	    system("$fastqc_script $sample_folder/*.fastq.gz --outdir=$fastQC_folder -t 10 > $fastQC_folder.log");
            
	    # launch mapping VS ncRNA
	    system("$bowtie2 -x $$run{ncrna} -U $sample_folder/$adn\_R1.fastq.gz -S $sample_folder/$adn\_R1.NCRNA.sam --un $sample_folder/$adn\_R1.minusNCRNA.fastq  > $sample_folder/$adn\_R1.align_summary.txt 2>&1");
	    my $fastq_R2 = "$sample_folder/$adn\_R2.fastq.gz";
	    if(-f $fastq_R2){
		system("$bowtie2 -x $$run{ncrna} -U $fastq_R2 -S $sample_folder/$adn\_R2.NCRNA.sam --un $sample_folder/$adn\_R2.minusNCRNA.fastq > $sample_folder/$adn\_R2.align_summary.txt 2>&1");
	    }
	    system("rm $sample_folder/$adn\_R1.NCRNA.sam");
	    my $R2_sam = "$sample_folder/$adn\_R2.NCRNA.sam";
	    if (-f $R2_sam){
		system("rm $R2_sam");
	    }
	    system("gzip $sample_folder/*.fastq");

	    ## create Rscript file and launch Rsubread
	    
	}

	
        ## regroup BAMs from tophat for analysis
	if ($$run{mapping} eq "tophat2"){
	    system("ln -s $sample_folder/tophat/accepted_hits.bam $bam_folder/$adn\_tophat.bam");
	    system("$samtools index $bam_folder/$adn\_tophat.bam");
	} elsif ($$run{mapping} eq "hisat2"){
	    system("ln -s $sample_folder/hisat/$adn.bam $bam_folder/$adn\_hisat.bam");
	    system("$samtools index $bam_folder/$adn\_hisat.bam");
	} elsif ($$run{mapping} eq "STAR"){
	    system("ln -s $sample_folder/STAR/$adn\Aligned.sortedByCoord.out.bam $bam_folder/$adn\_STAR.bam");
	    system("$samtools index $bam_folder/$adn\_STAR.bam");
	} elsif ($$run{mapping} eq "RSEM"){
	    system("ln -s $sample_folder/RSEM/$adn.genome.sorted.bam $bam_folder/$adn\_RSEM.bam");
	    system("$samtools index $bam_folder/$adn\_RSEM.bam");
	}


    }
    print "RNAseq analysis of $adn END...\n";

}

sub rna_stats {
    my $run_options = shift;
    my $list_adn = shift;
    
    my $typerun = $$run_options{typerun};
    my $run_folder = $$run_options{disque}."/".$$run_options{name};

    my $ngsbioinfo = "/path";
    my $link_folder = $ngsbioinfo."/".$$run_options{name};
    unless (-d $link_folder) {system("mkdir $link_folder");}
    
    my $rna_folder;
    
    print "### launch STATS RNA...\n";
    
    if ($$run_options{trim} eq "1"){
	$rna_folder = $run_folder."/rnaseq_analysis_trim";
    } else {
	$rna_folder = $run_folder."/rnaseq_analysis_noTrim";
    }
    my $stats_folder = $run_folder."/Reporting";
    system("mkdir $stats_folder");
    my $bam_folder = $rna_folder."/analysis";
    my $mapper = $$run_options{mapping};

    my (%R1_dup_stat, %R2_dup_stat, %tot_reads_nb, %ncrna1_reads, %ncrna2_reads, %mapped_reads);


    foreach my $adn (@$list_adn){
	my $sample_folder = $rna_folder."/Sample_$adn";
	my $adn_folder = "$stats_folder/Sample_$adn";
	unless (-d $adn_folder){system("mkdir $adn_folder");}
	# FastQC
	my $FastQC_folder = "$sample_folder/FastQC_report";
	system("mv $FastQC_folder $adn_folder");
	system("unzip $adn_folder/FastQC_report/$adn\_R1_fastqc.zip -d $adn_folder/FastQC_report");
	if ($typerun eq "PE"){
	    system("unzip $adn_folder/FastQC_report/$adn\_R2_fastqc.zip -d $adn_folder/FastQC_report");
	}
	# recup duplicates
	print "##### launch DUP stats\n";
	duplicates_stats($stats_folder, $adn, \%R1_dup_stat, \%R2_dup_stat);
	# recup nb reads total/ncRNA mapped/RNA mapped
	print "##### launch MAPPED stats\n";
	mapped_stats($stats_folder, $adn, $sample_folder, $run_folder, \%tot_reads_nb, \%ncrna1_reads, \%ncrna2_reads, \%mapped_reads, $mapper);
    
    }

    
    # recup stats_demultiplexing
    
    print "#### print STAT FILE... \n";
    my $stat_out_file = $stats_folder."/".$$run_options{name}."_stats.csv";
    open (OUT,">$stat_out_file") || die "Cannot open Stat file $stat_out_file :$!\n";
    print OUT "sample\t\% R1_duplicates\t\% R2_duplicates\t\# total reads\t\% R1 mapping ncRNA\t\% R2 mapping ncRNA\t\% total mapped reads\n";
    foreach my $adn (@$list_adn){
	my $dup_R1 = sprintf("%.2f",$R1_dup_stat{$adn});
	print OUT "$adn\t$dup_R1\t";
	if ($R2_dup_stat{$adn}){
	    my $dup_R2 = sprintf("%.2f",$R2_dup_stat{$adn});
	    print OUT "$dup_R2\t";
	} else {
	     print OUT "NA\t";
	}
	print OUT "$tot_reads_nb{$adn}\t$ncrna1_reads{$adn}\t";
	if ($ncrna2_reads{$adn}){
	    print OUT "$ncrna2_reads{$adn}\t";
	} else {
	     print OUT "NA\t";
	}
	print OUT "$mapped_reads{$adn}\n";
    }
    close OUT;

    system("ln -s $stats_folder $link_folder/Reporting");

}

sub filt_transcript{
    my $folder_in = shift;
    my $sample = shift;
    my $fileIn = "$folder_in/transcripts.gtf";
    my $fileOut = "$folder_in/$sample.transcripts.gtf";
    my ($FPKM, $cov, $FPKM_nb, $cov_nb, $total_filt_nb, $total_nb);


    open (IN, "$fileIn") || die "Cannot open $fileIn $!\n";
    open (OUT, ">$fileOut") || die "Cannot open $fileOut $!\n";
    while (my $l = <IN>){
	chomp $l;
	$total_nb++;
	if ($l=~/^.*FPKM \"([\d.]+)\";.*cov \"([\d.]+)\";.*/){
		$FPKM = $1; $cov = $2;
		# print "$FPKM , $cov \n";
		if ($FPKM !~ /^([0\.]+)$/){
			#print "++ $FPKM\n";
			$FPKM_nb++;
		}
		if ($cov !~ /^([0\.]+)$/){
		        #print "++ $cov\n";
			$cov_nb++;
		}
		if (($FPKM !~/^([0\.]+)$/) && ($cov !~/^([0\.]+)$/)) {
			$total_filt_nb++;
			print OUT "$l\n";
		}
	} else {
		die "line read error";
	}
    }
    close IN;
    close OUT;
    print "************************** Sample $sample ************************\n
* total lines : $total_nb\n
* FPKM = 0 : ".($total_nb-$FPKM_nb)."
* filt FPKM : $FPKM_nb\n
* cov = 0 : ".($total_nb-$cov_nb)."
* filt cov : $cov_nb\n
* filt total : $total_filt_nb\n
************************************************************
";

}

sub duplicates_stats {
    my $stats_folder = shift;
    my $dna = shift;
    my $R1_stat = shift;
    my $R2_stat = shift;
    
    my ($R1_dup, $R2_dup);
    
    my $R1_file = "$stats_folder/Sample_$dna/FastQC_report/$dna\_R1_fastqc/fastqc_data.txt";
    my $R2_file = "$stats_folder/Sample_$dna/FastQC_report/$dna\_R2_fastqc/fastqc_data.txt";
    
    $R1_dup = count_dup($R1_file);
    $$R1_stat{$dna} = $R1_dup;
    if (-f $R2_file){
	$R2_dup = count_dup($R2_file);
	$$R2_stat{$dna} = $R2_dup;
    } else {
	$R2_dup = "NA";
	$$R2_stat{$dna} = $R2_dup;
    }

}

sub count_dup {
    my $input_file = shift;
    my ($dup,$dedup);
    open (IN, $input_file) || die "Cannot open file dup_stats $input_file : $! \n";
    while (my $l=<IN>){
	#if ($l=~/^#Total Duplicate Percentage\s+([\d.]+)/){
	if ($l=~/^#Total Deduplicated Percentage\s+([\d.]+)/){
	    $dedup = $1;
	    $dup = 100 - $dedup;
	}	
    }    
    close IN;
    return $dup;
}

sub mapped_stats {
    my $stats_folder = shift;
    my $dna = shift;
    my $sample_folder = shift;
    my $run_folder = shift;
    my $tot_reads_nb = shift;
    my $ncrna1_reads = shift;
    my $ncrna2_reads = shift;
    my $mapped_reads = shift;
    my $mapper = shift;
    
    ### total reads number
    my $demultiplexing_folder = "$run_folder/demultiplexage/fastq";
    $$tot_reads_nb{$dna} = `zcat $demultiplexing_folder/$dna\_R*.fastq.gz | echo \$((\`wc -l\`/4))`;
    $$tot_reads_nb{$dna} =~ s/\012//;
    $$tot_reads_nb{$dna} =~ s/\015//;

    ### ncrna mapped reads
    my $ncrna1_file = "$sample_folder/$dna\_R1.align_summary.txt";
    my $ncrna2_file = "$sample_folder/$dna\_R2.align_summary.txt";
    $$ncrna1_reads{$dna} = count_ncrna($ncrna1_file);
    if (-f $ncrna2_file){
	$$ncrna2_reads{$dna} = count_ncrna($ncrna2_file);
    } else {
	$$ncrna2_reads{$dna} = "NA";
    }
        
    ### mapped reads
    if ($mapper eq "tophat2"){
	my $log_file = "";
    } elsif ($mapper eq "hisat2"){
	my $log_file = "";
    } elsif ($mapper eq "STAR"){
	my $log_file = $sample_folder."/STAR/".$dna."Log.final.out";
	$$mapped_reads{$dna} = star_log($log_file);
    } elsif ($mapper eq "RSEM"){
	my $log_file = $sample_folder."/RSEM/".$dna.".temp/".$dna."Log.final.out";
	$$mapped_reads{$dna} = star_log($log_file);
    }
    
}

sub count_ncrna {
    my $input_file = shift;
    my $count;
    open (IN, $input_file) || die "Cannot open NCrna_stats file $input_file : $! \n";
    while (my $l=<IN>){
	if ($l=~/^([\d\.]+)\% overall alignment rate/){
	    $count = $1;
	}	
    }    
    close IN;
    return $count;
}

sub star_log{
    my $input_file = shift;
    my ($total, $count1, $count2, $count3);
    open (IN, $input_file) || die "Cannot open STAR log file $input_file : $! \n";
    while (my $l=<IN>){
	if ($l=~/Uniquely mapped reads \% \|\s+([\d\.]+)\%/){
	    $count1 = $1;
	} elsif($l=~/\% of reads mapped to multiple loci \|\s+([\d\.]+)\%/) {
	    $count2 = $1;
	} elsif($l=~/\% of reads mapped to too many loci \|\s+([\d\.]+)\%/) {
	    $count3 = $1;
	}
    }    
    close IN;
    $total = $count1 + $count2 + $count3;
    return $total;
    
}

sub link_rna {
    
    
}

1;
