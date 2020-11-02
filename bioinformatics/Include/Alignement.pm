package Alignement;
use Data::Dumper;
use Sys::Hostname;

sub casava {
	my ($disque,$name,$gerald) = @_;
	system("mkdir $disque/$name/alignement");
	system("mkdir $disque/$name/alignement/casava");
	system("$gerald ./Config/gerald_config_PE_Hg19.txt --EXPT_DIR $disque/$name/demultiplexage --OUT_DIR $disque/$name/alignement/casava --make 1>$disque/$name/log/gerald.txt 2>$disque/$name/log/gerald.txt");
	system("cd $disque/$name/alignement/casava/;make -j80 all 1>>$disque/$name/log/gerald.txt 2>>$disque/$name/log/gerald.txt;cd -");
	#system("chmod -R 775 $disque/$name/alignement");
}

sub align_new_bwa {
	my ($bwa,$disque,$name,$adn,$genome_bwa_07,$samtools) = @_;
        unless (-d "$disque/$name/alignement") { system("mkdir $disque/$name/alignement"); }
        unless (-d "$disque/$name/alignement/bwa") { system("mkdir $disque/$name/alignement/bwa"); }
	#$adn =~ s/_/-/g;
	# Alignement
	#my $pro = `nproc`;
	#chomp $pro;
	my $pro = 20;
        system("$bwa mem -M -t $pro $genome_bwa_07 $disque/$name/demultiplexage/fastq/".$adn."_R1.fastq.gz $disque/$name/demultiplexage/fastq/".$adn."_R2.fastq.gz > $disque/$name/alignement/bwa/$adn.sam");
	# Sam2bam
	system("$samtools view -bS -q1 -@ 20 -o $disque/$name/alignement/bwa/$adn.bam $disque/$name/alignement/bwa/$adn.sam");
        # sort .bam
        system("$samtools sort $disque/$name/alignement/bwa/$adn.bam -@ 20 -o $disque/$name/alignement/bwa/$adn.sorted.bam");
        # index .bam
        system("$samtools index -b $disque/$name/alignement/bwa/$adn.sorted.bam");
}

sub gatk_realign {
    my $disque = shift;
    my $run_name = shift;
    my $adn = shift ;
    my $config = shift;
    my $run_data = shift;
    my $ref_genome = shift;
    my $align_tool = shift;

    my ($sample_data, $bam_folder, $bam_file, $bam_RG_file, $bam_RG_ordered_file, $realigned_file, $dedup_file,$recal_csv, $recal_file, $output_dir, $bedfile, $capture, $cov_file, $bed_option);

    my $host = hostname;
    my $gatk_launch;


    $gatk_launch = "gatk3";


    my $alignement_dir = "$disque/$run_name/alignement";
    my $gatk_dir = $alignement_dir."/gatk";
    unless (-d "$alignement_dir") { system("mkdir $alignement_dir"); }
    unless (-d "$gatk_dir") { system("mkdir $gatk_dir"); }
    my $log_folder = "$disque/$run_name/log";
    my $tmp_folder = $log_folder."/tmp";
    unless (-d "$tmp_folder") { system("mkdir $tmp_folder"); }

    foreach my $dna (@{ $adn }){

	#$dna =~ s/_/-/g;

	$sample_data = $$run_data->{$dna}; #####

	$output_dir = "$gatk_dir/Sample_$dna";
	unless (-d "$output_dir") { system("mkdir $output_dir"); }

	$capture = $sample_data->{capture};
	print "Capture utilisee : \n\n$capture\n\n\n !!!";
	# en fonction capture /bed, lancer pour exome, haloplex ou raindance (fichier de depart change)
	if ($capture =~ /raindance/){
	    $bam_folder = "$disque/$run_name/alignement/bwa";
	    $bam_file = "$bam_folder/Sample_$dna/$dna.bam";
	} elsif ($capture =~ /Haloplex/){
	    $bam_folder = "$disque/$run_name/alignement/bwa";
	    $bam_file = "$bam_folder/Sample_$dna/$dna.bam";
	} elsif ($capture eq "Agilent_84MB_Met"){
	    print "capture Methylation... please launch methylation analysis\n";
	} else {
	    $bam_folder = "$disque/$run_name/variation/casava";
	    $bam_file = "$bam_folder/Sample_$dna/genome/bam/sorted.bam";
	}

	# si exome, en fonction casava ou bwa, changer bam_folder et bam_file
	if ($align_tool eq "casava"){
	    $bam_folder = "$disque/$run_name/variation/casava";
	    $bam_file = "$bam_folder/Sample_$dna/genome/bam/sorted.bam";
	} elsif ($align_tool eq "bwa"){
	    $bam_folder = "$disque/$run_name/alignement/bwa";
	    $bam_file = "$bam_folder/$dna.sorted.bam";
	}

	if (!(-f $bam_file)){
	    print "BAM file $bam_file does not exist!!!\n";
	    next;
	}
	$bedfile = $$config{$capture};
	if ($bedfile eq "/media/Data/bed/capture/Twist/Exome_V_1_3_0_RefSeq_Gencode_Probe_Targets_hg38.bed"){
		$bedfile = "/media/Data/bed/capture/Twist/Exome+50bp/Exome_V_1_3_0_RefSeq_Gencode_Probe_Targets_hg38_+50bp.bed";
	}
	###
 print "BED utilise : \n\n$bedfile\n\n\n !!!";
	if ($bedfile eq "WGS"){
	    $bed_option = " ";
	} else {
	    $bed_option = "-L $bedfile";
	}
	$bam_RG_file = "$output_dir/$dna\_RG.bam";
	$bam_RG_ordered_file = "$output_dir/$dna\_RG.dedup.sorted.bam";
	$realigned_file = "$output_dir/$dna\_RG.dedup.sorted.realigned.bam";
	$dedup_file = "$output_dir/$dna\_RG.dedup.bam";
	$recal_csv = "$output_dir/$dna\_RG.dedup.sorted.realigned.recal.csv";
	$second_recal_csv = "$output_dir/$dna\_RG.dedup.sorted.realigned.second_recal.csv";
	$recal_file = "$output_dir/$dna\_RG.dedup.sorted.realigned.recal.bam";
	$cov_file = "$output_dir/$dna.cov";

	### Picard : AddReadGroups
	system ("java -Xms12g -Xmx24g -XX:ParallelGCThreads=20 -jar $$config{picard} AddOrReplaceReadGroups INPUT=$bam_file OUTPUT=$bam_RG_file RGID=1 RGLB=library RGPL=illumina RGPU=runBarcode RGSM=$dna CREATE_INDEX=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_folder > $log_folder/$dna\_AddOrReplaceReadGroups.log");

	### Picard : MarkDuplicates (enlever les duplicates)
	system ("java -XX:ParallelGCThreads=20 -Xms12g -Xmx24g -jar $$config{picard} MarkDuplicates INPUT=$bam_RG_file REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=$log_folder/$dna\_markDuplicate.dups OUTPUT=$dedup_file TMP_DIR=$tmp_folder");

	### Samtools : index (reindex unduplicated file)
	system ("$$config{samtools} index $dedup_file");

	### Picard : ReorderSam (reordonner les bam et generer index .bai)
	system ("java -XX:ParallelGCThreads=20 -Xms12g -Xmx24g -jar $$config{picard} ReorderSam INPUT=$dedup_file OUTPUT=$bam_RG_ordered_file REFERENCE=$ref_genome VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$tmp_folder > $log_folder/$dna\_ReorderSam.log");

        ### GATK : RealignTargetCreator (creation fichier intervals pour le realignement local)
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T RealignerTargetCreator -nt 20 -R $ref_genome -I $bam_RG_ordered_file $bed_option -known $$config{hg38_Mills_indels_file} -known $$config{hg38_kg_indels_file} --disable_auto_index_creation_and_locking_when_reading_rods -log $log_folder/$dna\_intervals.log -o $bam_RG_ordered_file.intervals");


	### GATK : IndelRealigner (realignement local)
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T IndelRealigner -allowPotentiallyMisencodedQuals -R $ref_genome -I $bam_RG_ordered_file -targetIntervals $bam_RG_ordered_file.intervals -known $$config{hg38_Mills_indels_file} -known $$config{hg38_kg_indels_file} --disable_auto_index_creation_and_locking_when_reading_rods -log $log_folder/$dna\_realign.log -o $realigned_file");

	### GATK : BaseRecalibrator (créer table de recalibration des bases) puis PrintReads
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T BaseRecalibrator -nct 20 -allowPotentiallyMisencodedQuals -I $realigned_file -R $ref_genome -knownSites $$config{hg38_dbsnp_file} -knownSites $$config{hg38_Mills_indels_file} -knownSites $$config{hg38_kg_indels_file} --disable_auto_index_creation_and_locking_when_reading_rods -o $recal_csv");
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T BaseRecalibrator -nct 20 -allowPotentiallyMisencodedQuals -I $realigned_file -R $ref_genome -knownSites $$config{hg38_dbsnp_file} -knownSites $$config{hg38_Mills_indels_file} -knownSites $$config{hg38_kg_indels_file} -BQSR $recal_csv --disable_auto_index_creation_and_locking_when_reading_rods -o $second_recal_csv");
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T PrintReads -nct 20 -allowPotentiallyMisencodedQuals -R $ref_genome -I $realigned_file -BQSR $second_recal_csv --disable_auto_index_creation_and_locking_when_reading_rods -o $recal_file");

	### Samtools : index (reindex recalibrated file)
	system ("$$config{samtools} index $recal_file");

        #retrieve coverage from GATK bam
	system("$$config{samtools} depth $recal_file > $cov_file");
	system("gzip $cov_file");
    }

}


sub bwa {
    my ($config,$disque,$name,$genome_bwa,$adn,$bed,$run_seq) = @_;
    my $alignement_dir = "$disque/$name/alignement";
    my $bwa_dir = $alignement_dir."/bwa";
    unless (-d "$alignement_dir") { system("mkdir $alignement_dir"); }
    unless (-d "$bwa_dir") { system("mkdir $bwa_dir"); }

    foreach my $dna (@{ $adn }){
        print $dna."\n";
        my $sample_fastq_dir;
	if (($run_seq=~/Hiseq3/) || ($run{nas}=~/NEXTSEQ/)){
	    $sample_fastq_dir = "$disque/$name/demultiplexage/fastq";
	} else {
	    $sample_fastq_dir = "$disque/$name/demultiplexage/Project_FC/Sample_$dna";
	}
        # create sample result folder
        system("mkdir $bwa_dir/Sample_$dna");

	# if haloplex -> cutadapt
	if ($bed=~ /haloplex/) {
	    print "HALOPLEX align\n";
	    # launch cutadapt on each R1/R2
	    system("cat $sample_fastq_dir/*$dna*R1*.fastq.gz > $sample_fastq_dir/$dna\_all_R1.fastq.gz");
	    system("cat $sample_fastq_dir/*$dna*R2*.fastq.gz > $sample_fastq_dir/$dna\_all_R2.fastq.gz");


	    opendir(DIR, $sample_fastq_dir) || die "Cannot open $sample_fastq_dir :$!\n";
	    my @files = grep {/\.fastq\.gz/} readdir(DIR);
	    closedir(DIR);
	    foreach my $file (@files){
		if ($file=~ /all_R1/){
		    system("$$config{cutadapt} -a $$config{R1} $sample_fastq_dir/$file -o $bwa_dir/Sample_$dna/$dna.R1.fastq -m 5"); #adapteur R1
		}elsif ($file=~ /all_R2/){
		    system("$$config{cutadapt} -a $$config{R2} $sample_fastq_dir/$file -o $bwa_dir/Sample_$dna/$dna.R2.fastq -m 5"); #adapteur R2
		}
	    }

	    # cat R1 - R2 fastq
	    system("cat $bwa_dir/Sample_$dna/$dna.R1.fastq $bwa_dir/Sample_$dna/$dna.R2.fastq > $bwa_dir/Sample_$dna/$dna.all.fastq");

	} elsif ($bed=~ /raindance/){
	    # if raindance -> no cutadapt
	    $sample_fastq_dir = "$disque/$name/demultiplexage/Project_FC/Sample_$dna";
	    print "RAINDANCE align\n";
	    system("zcat $sample_fastq_dir/*$dna*fastq.gz > $bwa_dir/Sample_$dna/$dna.all.fastq");
	}

        # bwa
        system("$$config{bwa} bwasw -t 10 $genome_bwa $bwa_dir/Sample_$dna/$dna.all.fastq > $bwa_dir/Sample_$dna/$dna.sam");
        # .sai to .bam
        system("$$config{samtools} view -bS -q1 -o $bwa_dir/Sample_$dna/$dna.bam $bwa_dir/Sample_$dna/$dna.sam");
        # sort .bam
        system("$$config{samtools} sort $bwa_dir/Sample_$dna/$dna.bam  $bwa_dir/Sample_$dna/$dna.sorted");
        # index .bam
        system("$$config{samtools} index $bwa_dir/Sample_$dna/$dna.sorted.bam");
    }
}



1;
