package Stats;

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use Projet;
use DBconnect;
use Detection;

sub recup_all_stats {
    my $sample_name = shift;
    my $sample_data = shift;
    my $sample_report_folder = shift;
    my $result_folder = shift;
    my $bed_type = shift; # normal bed = 0, alternative bed = 1, whole genome = 2
    my $run_name = shift;
    my $alt_bed = shift;
    my @arr = split("/",$alt_bed);
    my $bedname = $arr[(scalar(@arr))-1];
    $bedname =~ s/.bed//;
    my  $sample_cover_in_capture;
    if (($bed_type eq "0") || ($bed_type eq "2")){
	$sample_cover_in_capture = $result_folder."/SEQ_cover/InCapture/".$sample_data->{result_dir};
    } elsif ($bed_type eq "1"){
	$sample_cover_in_capture = $result_folder."/SEQ_cover/InBed/".$run_name."/".$bedname."/".$sample_data->{result_dir};
#    } elsif ($bed_type eq "2"){
#	print "wgs\n";
#	$sample_cover_in_capture = $result_folder."/SEQ_cover/all/".$sample_data->{result_dir};
    }

    # compute_mean_cover
    compute_mean_cover($sample_data, $sample_cover_in_capture);

    # compute_deep_nbp
    compute_deep_nbp($sample_data, $sample_cover_in_capture);

    # compute_deep_nbp_cumul
    compute_deep_nbp_cumul($sample_data, $sample_cover_in_capture);

    # compute_stats_regions
    compute_stat_regions($sample_data, $sample_cover_in_capture);

}

sub compute_mean_cover {
    my $sample_data = shift;
    my $sample_cover_in_capture = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)  
    
    my ($log, $line, $pos, $cov, %n_in, %n_out, %sum, $mean, $chrom, @chrom, $i, @tabcov, $long, $Q1, $Q3, $med, $valmed, $valmed1, $valmed2, $valQ1, $valQ3, $espinterQ,@tabcov1,@tabcovall,@tabcovall1,$longall, $Q1all, $Q3all, $medall, $valmedall, $valmedall1, $valmedall2, $valQ1all, $valQ3all, $espinterQall,$j,$k);
    my $sex = $sample_data->{sex};
    
#    if ($sex eq 'F'){
#	@chrom = (1 .. 22, X);
#    } elsif ($sex eq 'M'){
#	@chrom = (1 .. 22, X, Y);
#    } else {
	@chrom = (1 .. 22, X, Y);
#    }
    
    my $fileOut = $sample_cover_in_capture."/meanX.txt";
#    if (-f $fileOut) {
#	$log = scalar(localtime())."$fileOut already exists\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
#	return;
#   }

    $log = "compute mean cover for base $base meth $meth dna $dna => $fileOut\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open(FHO, "> $fileOut") or die("cannot create $fileOut:$!\n");
    
    $k = 0;
    my %cov_count;
    my $total;

    for $chrom (@chrom) {
	@tabcov = ();
	@tabcov1 = ();
	$j = 0;
	my $fileInCover = $sample_cover_in_capture."/chr$chrom.cov";

	if (! -f $fileInCover) {
	    $log = " !!!!!  skip cover file $fileInCover\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    next;
	}

	#si fichier est vide
	if (-z $fileInCover) {
	    $log = " !!!!! skip cover file $fileInCover because empty\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    next;
	}

	$log = "open cover file $fileInCover\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	open(FHI, $fileInCover) || die("cannot open $fileInCover:$!\n");

	while (my $line = <FHI>) {
	    next if ($line =~ /^\#/);
	    $n_in{$chrom}++;
	    $n_in{'all'}++;
	    chomp $line;
	    ($pos, $cov) = split(/\s+/, $line);
	    next if ($cov == 0);
	    $j++;
	    #$k++;
	    $n_out{'all'}++;
	    $sum{'all'} += $cov;
	    $n_out{$chrom}++;
	    $sum{$chrom} += $cov;

 	    $tabcov1[$j] = $cov;
 	    #$tabcovall1[$k] = $cov;
	    
	    $cov_count{$cov}++;
	    $total++;
	}

	close(FHI);

	if ($n_out{$chrom} == 0) {
	    $mean = 0;
	} else {
	    $mean = sprintf("%.2f", $sum{$chrom} / $n_out{$chrom});
	}

	@tabcov = sort {$a <=> $b} @tabcov1;

	$long = @tabcov;
	$Q1 = $long/4;
	$med = ($long+1)/2;
	$Q3 = $Q1*3;

 	if (int($med) == $med){
 	    $valmed = $tabcov[$med];
 	} else{
	    $valmed1 = $tabcov[int($med)];
	    $valmed2 = $tabcov[int($med)+1];
	    $valmed = ($valmed1+$valmed2)/2;
	}

	if (int($Q1) == $Q1){
	    $valQ1 = $tabcov[$Q1];
	    $valQ3 = $tabcov[$Q3];
	} elsif(($Q1 =~ /\.5/)||($Q1 =~ /\.75/)){
	    $valQ1 = $tabcov[int($Q1)+1];
	    $valQ3 = $tabcov[int($Q3)];
	} else {
	    $valQ1 = $tabcov[int($Q1)];
	    $valQ3 = $tabcov[int($Q3)+1];
	}

	print FHO "$chrom\t$mean\t$valmed\t$valQ1\t$valQ3\n";
	$log = "$chrom\t$mean\t$valmed\t$valQ1\t$valQ3\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    }

    if ($n_out{'all'} == 0) {
	$mean = 0;
    } else {
	$mean = sprintf("%.2f", $sum{'all'} / $n_out{'all'});
    }

    #@tabcovall = sort {$a <=> $b} @tabcovall1;

    #$longall = @tabcovall;
    #$Q1all = $longall/4;
    #$medall = ($longall+1)/2;
    #$Q3all = $Q1all*3;

    #if (int($medall) == $medall){
	#$valmedall = $tabcovall[$medall];
    #} else {
	#$valmedall1 = $tabcovall[int($medall)];
	#$valmedall2 = $tabcovall[int($medall)+1];
	#$valmedall = ($valmedall1+$valmedall2)/2;
    #}

     #if (int($Q1all) == $Q1all){
	#$valQ1all = $tabcovall[$Q1all];
	#$valQ3all = $tabcovall[$Q3all];
    #} elsif(($Q1all =~ /\.5/)||($Q1all=~ /\.75/)){
	#$valQ1all = $tabcovall[int($Q1all)+1];
	#$valQ3all = $tabcovall[int($Q3all)];
    #} else {
	#$valQ1all = $tabcovall[int($Q1all)];
	#$valQ3all = $tabcovall[int($Q3all)+1];
    #}
    
    my $t = my $cumul = 0;
    my $q1 = my $med = my $q3 = 0;
    for my $i (sort {$a<=>$b} keys %cov_count){
	$t += $cov_count{$i};
	$cumul = ($t/$total)*100;
	if ($cumul > 25){
	    if ($q1 eq "0"){
		$valQ1all = $i;
		#print "$cumul - $i\n";
	    }
	    $q1 = 1;
	}
	if ($cumul > 50){
	    if ($med eq "0"){
		$valmedall = $i;
		#print "$cumul - $i\n";
	    }
	    $med = 1;
	}
	if ($cumul > 75){
	    if ($q3 eq "0"){
		$valQ3all = $i;
		#print "$cumul - $i\n";
	    }
	    $q3 = 1;
	}
    }
    
    
    
    
    print FHO "all\t$mean\t$valmedall\t$valQ1all\t$valQ3all\n";
    print "all\t$mean\t$valmedall\t$valQ1all\t$valQ3all\n";
    
    close FHO;
    
#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

sub compute_deep_nbp {
    my $sample_data = shift;
    my $sample_cover_in_capture = shift;
    my $sex = $sample_data->{sex};
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)  

    my ($log, @chrom, $chrom, $fileIn, $line, $pos, $cov, %cptr, $fileOut, $n);

#    if ($sex eq 'F'){
#	@chrom = (1 .. 22, X);
#    } elsif ($sex eq 'M'){
#	@chrom = (1 .. 22, X, Y);
#    } else {
	@chrom = (1 .. 22, X, Y);
#    }
    
    $fileOut = $sample_cover_in_capture."/deep_nbp";
#    if (-f $fileOut) {
#	$log = "$fileOut already exists\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
#   	return;
#    }
    
    for $chrom (@chrom) {
	$fileIn = $sample_cover_in_capture."/chr$chrom.cov";
	if (! -f $fileIn) {
	    $log = " !!!!! skip cover file $fileIn\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    next;
	}
	$log = "open cover file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	$n = 0;
	open(FHI, $fileIn) or die("cannot open $fileIn:$!\n");
	while ($line = <FHI>) {
	    $n++;
	    next if ($line =~ /^\#/);
	    chomp $line;
	    ($pos, $cov) = split(/\s+/, $line);
	    $cptr{$cov}++;
	    $log.= $sample_data->{dna}." $chrom pos $n\n" if (($n % 1000000) == 0);
	}
	close(FHI);
    }

    open(FHO, "> $fileOut") or die("cannot create $fileOut:$!\n");
    print FHO "cov\tn_bp\n";
    for $cov (sort {$a <=> $b} (keys(%cptr))) {
	print FHO "$cov\t$cptr{$cov}\n";
    }
    close(FHO);
    
#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    
}

sub compute_deep_nbp_cumul {
    my $sample_data = shift;
    my $sample_cover_in_capture = shift;
    my $fileIn = $sample_cover_in_capture."/deep_nbp";
    my $fileOut = $sample_cover_in_capture."/deep_nbp_cumul";
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)  

    my ($log, $line, $cov, $nbp, @thr, $thr, %cptr, $all, $pct);

    if (! -f $fileIn) {
	$log = "no cover results for base $base meth $meth dna $dna ($fileIn)\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }

#    if (-f $fileOut) {
#	$log = "$fileOut already exists\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
#  	return,
#    }

    $log = "compute deep nbp cumul for base $base meth $meth dna $dna => $fileOut\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    my @thr = (100, 50, 40, 30, 20, 15, 10, 8, 5, 1);

    open(FHI, $fileIn) or die("cannot open $fileIn:$!\n");
    <FHI>;
    while ($line = <FHI>) {
	chomp $line;
	($cov, $nbp) = split(/\t/, $line);
	$all += $nbp;
	for $thr (@thr) {
	    if ($cov >= $thr) {
		$cptr{$thr} += $nbp;
	    }
	}
    }
    close(FHI);

    open(FHO, "> $fileOut") or die("cannot create $fileOut:$!\n");
    print FHO "thr\tn_bp\tpct\n";
    for $thr (@thr) {
	if (!($cptr{$thr})){
	    $nbp = 0;
	} else {
	    $nbp = $cptr{$thr};
	}
	if ($all == 0) {
	    $pct = 0;
	} else {
	    $pct = $cptr{$thr} / $all * 100;
	}
	print FHO "$thr\t$nbp\t$pct\n";
    }
    close(FHO);

#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

sub compute_stat_regions {
    my $sample_data = shift;
    my $sample_cover_in_capture = shift;
    my $sex = $sample_data->{sex};
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)  

    my (@chrom, $chrom, $fileIn, $fileOut, $line, $start, $end, $size, $mi_depth, $max_depth, $mean_depth);
    my (@thr, $thr,$log);
    my (%region, %bp, $n_region, $n_bp, $pct_region, $pct_bp);
    
    $fileOut = $sample_cover_in_capture."/stat_regions";
#    if (-f $fileOut) {
#	$log = "$fileOut already exists\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
#  	return;
#    }
    $log = "compute stat_region for base $base meth $meth dna $dna => $fileOut\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    @thr = (100, 50, 40, 30, 20, 15, 10, 8, 5, 1);
#    if ($sex eq 'F'){
#	@chrom = (1 .. 22, X);
#    } elsif ($sex eq 'M'){
#	@chrom = (1 .. 22, X, Y);
#    } else {
	@chrom = (1 .. 22, X, Y);
#    }
    
    $region{'all'}= 0;
    for $chrom (@chrom) {
	$fileIn =  $sample_cover_in_capture."/chr$chrom.stat";
	if (! -f $fileIn) {
	    $log = "skip cover file $fileIn\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    next;
	}
	$log = "open cover file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	open(FHI, $fileIn) or die("cannot open $fileIn:$!\n");
	<FHI>;
	
	while ($line = <FHI>) {
	    chomp $line;
	    ($chrom, $start, $end, $size, $mi_depth, $max_depth, $mean_depth) = split(/\t/, $line);
	    $mean_depth =~ s/O/0/g;
	    $region{'all'}++;
	    $bp{'all'} += $size;
	    for $thr (@thr) {
		if ($mean_depth >= $thr) {
		    $region{$thr} ++;
		    $bp{$thr} += $size;
		}
	    }
	}
	close(FHI);
    }
   

    open(FHO, "> $fileOut") or die("cannot create $fileOut:$!\n");
    print FHO "#thr\tn_region\tpct_region\tn_bp\tpct_bp\n";
    for $thr (@thr) {
	if(!($region{$thr})){
	    $n_region = 0;
	    $pct_region = 0;
	} else {
	    $n_region = $region{$thr};
	    $pct_region = $region{$thr} / $region{'all'} * 100;
	}
	if(!($bp{$thr})){	
	    $n_bp = 0;
	    $pct_bp = 0;
	} else {
	    $n_bp = $bp{$thr};
	    $pct_bp = $bp{$thr} / $bp{'all'} * 100;
	}
	print FHO "$thr\t$n_region\t$pct_region\t$n_bp\t$pct_bp\n";
    }
    print FHO "all\t$region{'all'}\n";
    close(FHO);

#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

sub reads_stats {
    my $dna = shift;
    my $run_data = shift; # object with samples info
    my $report_folder = shift; # run folder Reporting
    my $run_name = shift;
    my $run_disk = shift;
    my $run_nas = shift; # hiseq/miseq
    my $samtools_script = shift;
    
    print "$dna\n$run_name\n$run_disk\n$run_nas\n$samtools_script\n$report_folder\n";
    
    my ($demultiplexing_folder, $alignment_folder, $gatk_folder, $dna_fastq, $bed, $sample_bam, $sample_dedup_bam, $sample_bam_out, $sample_dedup_bam_out);
    my (%tot_reads_nb, %mapped_reads_nb, %inbed_reads_nb, %mapped_dedup_reads_nb, %inbed_dedup_reads_nb, %gatk_reads_nb, %gatk_inbed_reads_nb);

       
    my $ontarget_bam_folder = $run_disk."/".$run_name."/alignement/ontarget";
    unless (-d $ontarget_bam_folder) { system("mkdir $ontarget_bam_folder"); }
    

    $bed = $run_data->{$dna}{bed_file};
    print "BED file : ".$bed."\n";
    #if (-d $run_disk."/".$run_name."/alignement/gatk/Sample_$dna"){
    $sample_bam = "$run_disk/$run_name/alignement/bwa/$dna.sorted.bam";
    $sample_dedup_bam = "$run_disk/$run_name/alignement/gatk/Sample_$dna/$dna\_RG.dedup.sorted.realigned.recal.bam";
    #} elsif (-d $run_disk."/".$run_name."/alignement/bwa") {
	#my $dnabis = $dna;
	#$dnabis =~ s/_/-/g;
	#print $dnabis."\n";
	#if (-f "$run_disk/$run_name/alignement/bwa/$dna.sorted.bam"){
	#    $sample_bam = "$run_disk/$run_name/alignement/bwa/$dna.sorted.bam";
	#} elsif (-f "$run_disk/$run_name/alignement/bwa/$dnabis.sorted.bam"){
	#    $sample_bam = "$run_disk/$run_name/alignement/bwa/$dnabis.sorted.bam";
	#} else {
	#    $sample_bam = "$run_disk/$run_name/alignement/bwa/Sample_$dna/$dna.sorted.bam";
	#}
    #} elsif (-d $run_disk."/".$run_name."/alignement/casava") {
	#$sample_bam = "$run_disk/$run_name/variation/casava/Sample_$dna/genome/bam/sorted.bam";
    #}
    
    ###$sample_gatk_bam = "$run_disk/$run_name/alignement/gatk/Sample_$dna/$dna\_RG.sorted.realigned.dedup.recal.bam";
    $sample_bam_out = "$ontarget_bam_folder/$dna\_inbed.bam";
    $sample_dedup_bam_out = "$ontarget_bam_folder/$dna\_dedup_inbed.bam";
    
    ## writeRegion : extract inbed reads, mapped and dedup inbed
    system("/media/Script/prog/bamUtil/bin/bam writeRegion --in $sample_bam --bed $bed --out $sample_bam_out");
    system("/media/Script/prog/bamUtil/bin/bam writeRegion --in $sample_dedup_bam --bed $bed --out $sample_dedup_bam_out");
    
    ## flagstat : aligned, dedup, inbed, dedup inbed
    system("$samtools_script flagstat $sample_bam > $ontarget_bam_folder/$dna\_flagstat.txt");
    system("$samtools_script flagstat $sample_dedup_bam > $ontarget_bam_folder/$dna\_dedup_flagstat.txt");
    system("$samtools_script flagstat $sample_bam_out > $ontarget_bam_folder/$dna\_inbed_flagstat.txt");
    system("$samtools_script flagstat $sample_dedup_bam_out > $ontarget_bam_folder/$dna\_inbed_dedup_flagstat.txt");
    

    ## total reads number
    $demultiplexing_folder = $run_disk."/".$run_name."/demultiplexage/fastq";
    $tot_reads_nb{$dna} = `zcat $demultiplexing_folder/$dna\_R*.fastq.gz | echo \$((\`wc -l\`/4))`;
    $tot_reads_nb{$dna} =~ s/\012//;
    $tot_reads_nb{$dna} =~ s/\015//;


    ## write csv : read flagstats + total reads number
    my $flag_file_bam = "$ontarget_bam_folder/$dna\_flagstat.txt";
    my $flag_dedup_bam = "$ontarget_bam_folder/$dna\_dedup_flagstat.txt";
    my $flag_inbed_file_bam = "$ontarget_bam_folder/$dna\_inbed_flagstat.txt";
    my $flag_inbed_dedup_bam = "$ontarget_bam_folder/$dna\_inbed_dedup_flagstat.txt";
    
    $mapped_reads_nb{$dna} = recup_flagstat($flag_file_bam);
    $inbed_reads_nb{$dna} = recup_flagstat($flag_inbed_file_bam);
    $mapped_dedup_reads_nb{$dna} = recup_flagstat($flag_dedup_bam);
    $inbed_dedup_reads_nb{$dna} = recup_flagstat($flag_inbed_dedup_bam);

    my $outfile = "$run_disk/$run_name/Reporting/stat_reads_$run_name.txt";
    if (-f $outfile){
	open (OUT,">>$outfile") || die "Cannot open file $outfile : $! \n";
	my $percent_mapped = sprintf("%.2f",(($mapped_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_dedup_mapped = sprintf("%.2f",(($mapped_dedup_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_ontarget = sprintf("%.2f",(($inbed_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_dedup_ontarget = sprintf("%.2f",(($inbed_dedup_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	print OUT "$dna\t$tot_reads_nb{$dna}\t$mapped_reads_nb{$dna}\t$percent_mapped\t$inbed_reads_nb{$dna}\t$percent_ontarget\t$mapped_dedup_reads_nb{$dna}\t$percent_dedup_mapped\t$inbed_dedup_reads_nb{$dna}\t$percent_dedup_ontarget\n";
	close OUT;
    } else {
	open (OUT,">$outfile") || die "Cannot open file $outfile : $! \n";
	print OUT "Sample\tnb reads total\tnb total reads mapped\t\% reads mapped\tnb ontarget reads\t\% ontarget reads\tnb uniq reads mapped\t\% uniq reads mapped\tnb uniq ontarget reads\t\% uniq ontarget reads\n";  
	my $percent_mapped = sprintf("%.2f",(($mapped_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_dedup_mapped = sprintf("%.2f",(($mapped_dedup_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_ontarget = sprintf("%.2f",(($inbed_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	my $percent_dedup_ontarget = sprintf("%.2f",(($inbed_dedup_reads_nb{$dna} * 100)/$tot_reads_nb{$dna}));
	print OUT "$dna\t$tot_reads_nb{$dna}\t$mapped_reads_nb{$dna}\t$percent_mapped\t$inbed_reads_nb{$dna}\t$percent_ontarget\t$mapped_dedup_reads_nb{$dna}\t$percent_dedup_mapped\t$inbed_dedup_reads_nb{$dna}\t$percent_dedup_ontarget\n";
	close OUT;
    }
}

sub recup_flagstat{

#Correction Mehdi 160719 car stats reads_mapped supperieurs à 100%.

    #my $infile = shift;
    #my $read_count;

    #open(IN, $infile) || die "Cannot open file $infile : $!\n";
    #my $first_line = <IN>;
    #if ($first_line =~ /^([\d]+)\s\+\s([\d]+)\sin total/){
	#$read_count = $1;
	#print "flagstat from ".$infile." : ".$read_count."\n";
   # }
   # close IN;
    #return $read_count;


    my $infile = shift;

    my $read_count;
    my $total_count;
    my $secondary;

    open(IN, $infile) || die "Cannot open file $infile : $!\n";
    while (my $line = <IN>) {
   
    if ($line =~ /^([\d]+)\s\+\s([\d]+)\sin total/){
	$total_count = $1;
    }
   if ($line =~ /^([\d]+)\s\+\s([\d]+)\ssecondary/){
    $secondary = $1;
    }
     
  }
   close IN;

   $read_count = $total_count - $secondary;
   return $read_count;
}

sub duplicates_stats {
    my $dna_list = shift;
    my $run_data = shift; # object with samples info
    my $report_folder = shift; # run folder Reporting
    my $run_name = shift;
    my $run_disk = shift;
    my $run_nas = shift; # hiseq/miseq
    
    my (%R1_stat, %R2_stat, $R1_dup, $R2_dup);
    
    foreach $dna (@{$dna_list}) {
	# fastqc last version : add unzip files
	system("unzip $report_folder/Sample_$dna/FastQC_report/$dna\_R1_fastqc.zip -d $report_folder/Sample_$dna/FastQC_report");
	system("unzip $report_folder/Sample_$dna/FastQC_report/$dna\_R2_fastqc.zip -d $report_folder/Sample_$dna/FastQC_report");
	my $R1_file = "$report_folder/Sample_$dna/FastQC_report/$dna\_R1_fastqc/fastqc_data.txt";
	my $R2_file = "$report_folder/Sample_$dna/FastQC_report/$dna\_R2_fastqc/fastqc_data.txt";
	
	$R1_dup = count_dup($R1_file);
	$R2_dup = count_dup($R2_file);
	
	$R1_stat{$dna} = $R1_dup;
	$R2_stat{$dna} = $R2_dup;
	
    }
    
    my $outfile = "$report_folder/dup_rate_$run_name.csv";
    open (OUT,">$outfile") || die "Cannot open file $outfile : $! \n";
    print OUT "#dna\tR1_duplicates_percent\tR2_duplicates_percent\n";
    foreach $dna (@{$dna_list}){
	my $dup_R1 = sprintf("%.2f",$R1_stat{$dna});
	my $dup_R2 = sprintf("%.2f",$R2_stat{$dna});
	print OUT "$dna\t$dup_R1\t$dup_R2\n";
    }
    close OUT;
}

sub count_dup {
    my $input_file = shift;
    my ($dedup,$dup);
    open (IN, $input_file) || die "Cannot open file $input_file : $! \n";
    while (my $l=<IN>){
	##if ($l=~/^#Total Duplicate Percentage\s+([\d.]+)/){ #old version fastqc
	if ($l=~/^#Total Deduplicated Percentage\s+([\d.]+)/){
	    $dedup = $1;
	    $dup = 100 - $dedup;
	}
    }  
    close IN;
    return $dup;
}

sub stats_to_sql {
    my $run_disk = shift;
    my $run_name = shift;
    
    my ($stat_cov, $stat_mut, $dup_rate, $stat_reads, $cov, $mut, $dup, $reads, $samplesheet, @samples, $sample, $stat_sample);
    my $dbh = DBconnect::connect_db('RUNS');
    $samplesheet = "$run_disk/$run_name/samplesheet.csv";
    $stat_cov = "$run_disk/$run_name/Reporting/stat_cov_in_capture_$run_name.txt";
    $stat_mut = "$run_disk/$run_name/Reporting/stat_mut_$run_name.txt";
    $dup_rate = "$run_disk/$run_name/Reporting/dup_rate_$run_name.csv";
    $stat_reads = "$run_disk/$run_name/Reporting/stat_reads_$run_name.txt";

    print "***************************\n";
    @samples = Detection::read_sample($run_disk,$run_name);
    $cov = read_cov_file($samples,$stat_cov);
    $mut = read_mut_file($samples,$stat_mut);
    $dup = read_dup_file($samples,$dup_rate);
    $reads = read_stat_reads_file($samples,$stat_reads);

    #print Dumper(\@samples);
    foreach $sample (@samples){
	add_stats($dbh, $sample, $$cov{$sample},$$mut{$sample},$$dup{$sample},$$reads{$sample});
    }
    
    $dbh->disconnect();
}


sub add_stats() {
    my ($dbh, $name_sample, $cov, $mut, $dup, $reads) = @_;
    
    #Insertion dans la table
    #~ my $requete_sql_run = "INSERT IGNORE INTO samples_stats (name_sample,nb_total_reads,nb_mapped_reads,nb_reads_ontarget,meanx,chrX_meanX,chrY_meanX,aut_meanX,mediane,quartile1,quartile3,mbp_total,mbp8,pmbp8,mbp20,pmbp20,mbp50,pmbp50,mbp100,pmbp100,regions_total,regions8,pregions8,regions20,pregions20,regions50,pregions50,regions100,pregions100,nb_snps,nb_indels,r1_dup,r2_dup) VALUES (\'$name_sample\',$reads,$cov,$mut,$dup) ON DUPLICATE KEY UPDATE \'name_sample\' = \'$name_sample\' ;";
	#my $requete_sql_run = "INSERT IGNORE INTO samples_stats (name_sample,nb_total_reads,nb_mapped_reads,nb_reads_ontarget,meanx,chrX_meanX,chrY_meanX,aut_meanX,mediane,quartile1,quartile3,mbp_total,mbp8,pmbp8,mbp20,pmbp20,mbp50,pmbp50,mbp100,pmbp100,regions_total,regions8,pregions8,regions20,pregions20,regions50,pregions50,regions100,pregions100,nb_snps,nb_indels,r1_dup,r2_dup) VALUES (\'$name_sample\',$reads,$cov,$mut,$dup) ON DUPLICATE KEY UPDATE nb_total_reads = VALUES(nb_total_reads),nb_mapped_reads=VALUES(nb_mapped_reads),nb_reads_ontarget=VALUES(nb_reads_ontarget),meanx=VALUES(meanx),chrX_meanX= VALUES(chrX_meanX),chrY_meanX = VALUES(chrY_meanX),aut_meanX=VALUES(aut_meanX),mediane=VALUES(mediane),quartile1=VALUES(quartile1),quartile3=VALUES(quartile3),mbp_total=VALUES(mbp_total),mbp8=VALUES(mbp8),pmbp8=VALUES(pmbp8),mbp20=VALUES(mbp20),pmbp20=VALUES(pmbp20),mbp50 = VALUES(mbp50),pmbp50=VALUES(pmbp50),mbp100=VALUES(mbp100),pmbp100=VALUES(pmbp100),regions_total=VALUES(regions_total),regions8=VALUES(regions8),pregions8=VALUES(pregions8),regions20=VALUES(regions20),pregions20=VALUES(pregions20),regions50=VALUES(regions50),pregions50=VALUES(pregions50),regions100=VALUES(regions100),pregions100=VALUES(pregions100),nb_snps=VALUES(nb_snps),nb_indels=VALUES(nb_indels),r1_dup=VALUES(r1_dup),r2_dup=VALUES(r2_dup) ;";
	my $requete_sql_run = "REPLACE INTO samples_stats (name_sample,nb_total_reads,nb_mapped_reads,nb_reads_ontarget,nb_dedup_mapped_reads,nb_dedup_reads_ontarget,meanx,chrX_meanX,chrY_meanX,aut_meanX,mediane,quartile1,quartile3,mbp_total,mbp8,pmbp8,mbp20,pmbp20,mbp50,pmbp50,mbp100,pmbp100,regions_total,regions8,pregions8,regions20,pregions20,regions50,pregions50,regions100,pregions100,nb_snps,nb_indels,r1_dup,r2_dup) VALUES (\'$name_sample\',$reads,$cov,$mut,$dup);";

    my $sth_run = $dbh->prepare($requete_sql_run) or die $dbh->errstr;
    
    $sth_run->execute() or die "Echec Requete $requete_sql_run : $DBI::errstr";
    $sth_run->finish();
    print "$requete_sql_run a ete executee correctement. \n";
}


sub read_cov_file() {
    my $samples = shift;
    my ($input) = shift;
    my (%result, $sample, @tab, $stat);
    my $firstLine = 1;

    if (-f $input){
	open(FH,"$input")|| die "Cannot open $input : $!\n";
	while (my $l=<FH>){
	    #skip first line
	    next if $. < 2;
	    #each line on @tab
	    # last cov files
	    if ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'$10\',\'$11\',\'$12\',\'$13\',\'$14\',\'$15\',\'$16\',\'$17\',\'$18\',\'$19\',\'$20\',\'$21\',\'$22\',\'$23\',\'$24\',\'$25\',\'$26\',\'$27\',\'$28\'";
		$result{$sample} = $stat;
	    # old cov files with chrX and chrY NA
	    } elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d.NA]+)\s+([\d.NA]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'$5\',\'$6\',\'NA\',\'$7\',\'$8\',\'$9\',\'$10\',\'$11\',\'$12\',\'$13\',\'$14\',\'$15\',\'$16\',\'$17\',\'$18\',\'$19\',\'$20\',\'$21\',\'$22\',\'$23\',\'$24\',\'$25\',\'$26\',\'$27\'";
		$result{$sample} = $stat;
	    # old cov files
	    } elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'NA\',\'NA\',\'NA\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'$10\',\'$11\',\'$12\',\'$13\',\'$14\',\'$15\',\'$16\',\'$17\',\'$18\',\'$19\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		$result{$sample} = $stat;
	    } elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'NA\',\'NA\',\'NA\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'10\',\'$11\',\'12\',\'$13\',\'14\',\'$15\',\'16\',\'$17\',\'$18\',\'19\',\'$20\',\'21\',\'$22\',\'23\',\'NA\',\'NA\'";
		$result{$sample} = $stat;
	    } elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d\w]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'NA\',\'NA\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'$10\',\'$11\',\'$12\',\'$13\',\'$14\',\'$15\',\'$16\',\'$17\',\'$18\',\'$19\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		$result{$sample} = $stat;
	    } elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'NA\',\'NA\',\'NA\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'$10\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'$11\',\'$12\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		$result{$sample} = $stat;
	    }
	    elsif ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d.]+)\s+([\d]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'NA\',\'NA\',\'NA\',\'$5\',\'$6\',\'$7\',\'$8\',\'$9\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		$result{$sample} = $stat;
	    }
	}

	close FH;
    } else {
		foreach my $samplelist (@$samples){
			$result{$samplelist} = "\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		}
    }
    return \%result;
}

sub read_mut_file() {
    my $samples = shift;
    my ($input) = shift;
    my (%result, $sample, @tab, $stat);
    my $firstLine = 1;

    if (-f $input){
	open(FH,"$input")|| die "Cannot open $input : $!\n";
	while (my $l=<FH>){
	    #skip first line
	    next if $. < 2;
	    #each line on @tab
	    if ($l=~/^([\d\w-]+)\s+([\w]+)\s+([\d])\s+([\d\w]+)\s+([\d\w]+)\s+([\d\w.]+)\s+([\d\w.]+)\s+([\d\w.]+)\s+([\d\w.]+)\s+([\d\w.]+)\s+([\d\w.]+)$/){
		$sample = $1;
		$stat = "\'$4\',\'$5\'";
		$result{$sample} = $stat;
	    }
	}
	close FH;
    } else {
		foreach my $samplelist (@$samples){
			$result{$samplelist} = "\'NA\',\'NA\'";
		}
    }
	return \%result;
}

sub read_dup_file() {
    my $samples = shift;
    my ($input) = shift;
    my (%result, $sample, @tab, $stat);
    my $firstLine = 1;

    if (-f $input){
	open(FH,"$input")|| die "Cannot open $input : $!\n";
	while (my $l=<FH>){
	    #skip first line
	    next if $. < 2;
	    #each line on @tab
	    if ($l=~/^([\d\w-]+)\s+([\d.]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$2\',\'$3\'";
		$result{$sample} = $stat;
	    }
	}
	close FH;
    } else {
		foreach my $samplelist (@$samples){
			$result{$samplelist} = "\'NA\',\'NA\'";
		}
    }
    return \%result;
}

sub read_stat_reads_file() {
    my $samples = shift;
    my ($input) = shift;
    my (%result, $sample, @tab, $stat);
    my $firstLine = 1;

    if (-f $input){
	open(FH,"$input")|| die "Cannot open $input : $!\n";
	while (my $l=<FH>){
	    #skip first line
	    next if $. < 2;
	    #each line on @tab
	    if ($l=~/^([\d\w-_]+)\s+([\d]+)\s+([\d]+)\s+([\d.]+)\s+([\d]+)\s+([\d.]+)\s+([\d]+)\s+([\d.]+)\s+([\d]+)\s+([\d.]+)$/){
		$sample = $1;
		$stat = "\'$2\',\'$3\',\'$5\',\'$7\',\'$9\'";
		#$stat = "\'$2\',\'$3\',\'$5\'";
		$result{$sample} = $stat;
	    }
	}
	close FH;
    } else {
		foreach my $samplelist (@$samples){
		    $result{$samplelist} = "\'NA\',\'NA\',\'NA\',\'NA\',\'NA\'";
		    #$result{$samplelist} = "\'NA\',\'NA\',\'NA\'";
		}
    }
    return \%result;
}



1;
