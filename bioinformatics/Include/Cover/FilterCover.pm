package FilterCover;


use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use LogOut;
use Spreadsheet::WriteExcel;


sub filter_cover_in_capture {
    my $sample_data = shift;
    my $result_folder = shift;
    my $config = shift;
    
    my $seq_cover_folder = $result_folder."/SEQ_cover/InCapture/".$sample_data->{result_dir};
    my $seq_cover_all_folder = $result_folder."/SEQ_cover/all/".$sample_data->{result_dir};
    my (@chrs, @PID); 
    my $bedFile;
    
    ## en fonction de la capture
    if ($sample_data->{capture} =~ /raindance/){
 	$bedFile = $sample_data->{bed_file};
        print "*** filter sample raindance\n";
        #if ($sample_data->{sex} eq 'F'){
        #    @chrs = (1 .. 22, 'X');
        #} elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
        #}
        print "start PIDs fork filter $sample_data->{dna}\n";
        foreach my $chr (@chrs){
            my $pid = fork();
            if($pid){
                push (@PID, $pid);
            } elsif ($pid == 0){
                filter_cover($sample_data, $seq_cover_folder, $seq_cover_all_folder, $chr, $bedFile);
                exit 0;
            } else {
                die "unable to fork: $!";
            }
        }
  
    } elsif ($sample_data->{capture} =~ /Agilent_Haloplex/){
	$bedFile = $sample_data->{bed_file};
        print "*** filter sample Haloplex\n";
        #if ($sample_data->{sex} eq 'F'){
        #    @chrs = (1 .. 22, 'X');
        #} elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
        #}
        print "start PIDs fork filter $sample_data->{dna}\n";
        foreach my $chr (@chrs){
            my $pid = fork();
            if($pid){
                push (@PID, $pid);
            } elsif ($pid == 0){
                filter_cover($sample_data, $seq_cover_folder, $seq_cover_all_folder, $chr, $bedFile);
                exit 0;
            } else {
                die "unable to fork: $!";
            }
        }
    } elsif ($sample_data->{capture} =~ /WGS/){
	print "*** WGS : no filtering";
	system("mkdir $seq_cover_folder") unless (-d "$seq_cover_folder");
	#if ($sample_data->{sex} eq 'F'){
	#    @chrs = (1 .. 22, 'X');
	#} elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
	#} else {
	#    @chrs = (1 .. 22, 'X', 'Y');
	#}
	print "start PIDs fork filter $sample_data->{dna}\n";
	foreach my $chr (@chrs){
	    my $pid = fork();
	    if($pid){
		push (@PID, $pid);
	    } elsif ($pid == 0){
		system("cp $seq_cover_all_folder/chr$chr.cov $seq_cover_folder/chr$chr.cov");
		exit 0;
	    } else {
		die "unable to fork: $!";
	    }	    
	}
	
    } else {
	$bedFile = $$config{$sample_data->{capture}};
	my @recup_exome_jobs;   # pour lancer en parallèle
	#if ($sample_data->{sex} eq 'F'){
	#    @chrs = (1 .. 22, 'X');
	#} elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
	#} else {
	#    @chrs = (1 .. 22, 'X', 'Y');
	#}
	print "start PIDs fork filter $sample_data->{dna}\n";
	foreach my $chr (@chrs){
	    my $pid = fork();
	    if($pid){
		push (@PID, $pid);
	    } elsif ($pid == 0){
		filter_cover($sample_data, $seq_cover_folder, $seq_cover_all_folder, $chr, $bedFile);
		exit 0;
	    } else {
		die "unable to fork: $!";
	    }	    
	}
		
    }
    for my $pid (@PID) {
	#print "waitpid $pid\n";
	my $tmp = waitpid($pid, 0);
	print "done with pid $tmp\n";
    }
    print "END of fork PIDs filter\n";
}

sub filter_cover {
    my $sample_data = shift;
    my $seq_cover_in_capture_folder = shift;
    my $seq_cover_all_folder = shift;
    my $chr = shift;
    my $bedFile = shift;
    
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log
        
    my ($log);
    my ($line_B, $line_C);
    my ($bed_chrom, $bed_start, $bed_end);
    my ($i, $cover_pos);
    my ($in_bed_region);
    my ($cov_depth, $min, $max, $n, $sum, $mean);
    my ($bed_n);
    
    #my $seq_cover_all_folder = $seq_cover_folder."/all/".$sample_data->{result_dir};
    #my $seq_cover_in_capture_folder = $seq_cover_folder."/InCapture/".$sample_data->{result_dir};
    mkdir($seq_cover_in_capture_folder);
    
    my $bed_chrom_file = $bedFile."_chr".$chr;
    my $seq_cover_all_file = $seq_cover_all_folder."/chr".$chr.".cov";
    my $seq_cover_in_capture_file = $seq_cover_in_capture_folder."/chr".$chr.".cov";
    my $seq_stat_in_capture_file = $seq_cover_in_capture_folder."/chr".$chr.".stat";
    
#    if (-f $seq_cover_in_capture_file){
#	$log = scalar(localtime())." : ".$seq_cover_in_capture_file." already exists!\n";
#	LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
#	return;
#    }

    $log = "[dna] = $sample_data->{dna}, [chr] = $chr\n";
    $log.= "will do $seq_cover_in_capture_file\n";
    $log.= "open bed file $fileInBed\n";
    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );   
    open(FHIB, $bed_chrom_file) || die("cannot open $bed_chrom_file:$!\n");

    $log = "open cover file $seq_cover_all_file\n";
    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
    open(FHIC, $seq_cover_all_file) || die("cannot open $seq_cover_all_file:$!\n");
    <FHIC>;
    <FHIC>;

    $log = "open out.cov  file $seq_cover_in_capture_file\n";
    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
    open(FHOC, ">$seq_cover_in_capture_file") || die("cannot create $seq_cover_in_capture_file:$!\n");

    $log = "open out.stat file $seq_stat_in_capture_file\n";
    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
    open(FHOS, ">$seq_stat_in_capture_file") || die("cannot create $seq_stat_in_capture_file:$!\n");
    print FHOS "chrom\tstart\tend\tsize_bp\tmin_depth\tmax_depth\tmean_depth\n";

    $in_bed_region = 0;

    # 1ère region
    $line_B = <FHIB>;
    $bed_n++;
    $line_B =~ s/\012//;
    $line_B =~ s/\015//;
    ($bed_chrom, $bed_start, $bed_end) = split(/\t/, $line_B);
    $bed_chrom =~ s/chr//;
    $log.= "first bed $bed_chrom $bed_start $bed_end\n";

    while ($line_C = <FHIC>) {
	chomp $line_C;
	$i = index($line_C, " ");
	$cover_pos = substr($line_C, 0, $i);
	$cov_depth = substr($line_C, $i + 1);
	$cov_depth =~ s/ //g;

 	# print "pos : $cover_pos || bed pos : $bed_start\n";
	# print "\t* pos $cover_pos cov_depth $cov_depth\n";

	if ($cover_pos < $bed_start) {
	    # avant la 1ère region
	    # print "\t- cover_pos $cover_pos < bed_start $bed_start\n";
	    next;
	}
	if (($cover_pos >= $bed_start) and ($cover_pos <= $bed_end)) {
	    # dans une region
	    if ($in_bed_region == 0) {
		# debut de la region
		# print "\tenter bed region $bed_chrom $bed_start $bed_end\n";
		$in_bed_region = 1;
		$n = 0;
		$min = 0;
		$max = 0;
		$sum = 0;
	    }
	    # print "\t+ cover_pos $cover_pos [ bed_start $bed_start and bed_end $bed_end]\n";
	    print FHOC "$cover_pos\t$cov_depth\n";
	    $n++;
	    $sum += $cov_depth;
	    $min = $cov_depth if ($cov_depth < $min);
	    $max = $cov_depth if ($cov_depth > $max);
	    next;
	}
	if ($cover_pos > $bed_end) {
	    # après la region bed
	    if ($in_bed_region == 1) {
		# region en cours
		# print "\tleave bed region $bed_chrom $bed_start $bed_end\n";
		$mean = sprintf("%.1f", $sum / $n);
		print FHOS "$bed_chrom\t$bed_start\t$bed_end\t$n\t$min\t$max\t$mean\n";
		$in_bed_region = 0;
	    } else {
		# pas de region en cours
		# print "\tno cover at all for bed region $bed_chrom $bed_start $bed_end\n";
		print FHOS "$bed_chrom\t$bed_start\t$bed_end\t0\t0\t0\t0\n";
		# passe au bed suivant
	    }
	    $line_B = <FHIB>;
	    goto FIN unless ($line_B) ;
	    $bed_n++;
	    $line_B =~ s/\012//;
	    $line_B =~ s/\015//;
	    ($bed_chrom, $bed_start, $bed_end) = split(/\t/, $line_B);
	    $bed_chrom =~ s/chr//;
	    # print "next bed $bed_chrom $bed_start $bed_end\n";
	    next;
	}
    }
  FIN:
    if ($in_bed_region == 1) {
	# region en cours
	$log = "\tleave last bed region $bed_chrom $bed_start $bed_end\n";
	LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
	$mean = sprintf("%.1f", $sum / $n);
	print FHOS "$bed_chrom\t$bed_start\t$bed_end\t$n\t$min\t$max\t$mean\n";
	$log = "next bed $bed_chrom $bed_start $bed_end\n";
	LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
    } else {
	# pas de region en cours
	$log = "\tno cover at all for last bed region $bed_chrom $bed_start $bed_end\n";
	LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
	# passe au bed suivant
    }
    close(FHIB);
    close(FHIC);
    close(FHOC);
    close(FHOS);
    
#   LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );    
}

sub filter_cover_in_bed {
    my $sample_data = shift;
    my $result_folder = shift;
    my $bedFile = shift;

    my $selection_name = $bedFile;
    $selection_name =~ s/\/media\///;
    $selection_name =~ s/Data\/bed\/capture\///;
    $selection_name =~ s/Data\/bed\///;
    $selection_name =~ s/Data\/ngs\/data\/raindance\///;
    $selection_name =~ s/_Covered\.bed//;
    $selection_name =~ s/\.bed//;
    $selection_name =~ s/_sorted//;
    $selection_name =~ s/_uniq//;
    $selection_name =~ s/_gene_ids//;
    $selection_name =~ s/\//_/;

    my $seq_cover_folder = $result_folder."/SEQ_cover";
    my $depth_thr = 8;

    my ($chrom, $line_B, $bed_n, $bed_chrom, $bed_start, $bed_end, $bed_gene, $gene, $n_exon);
    my ($origDir, $dir1, $dir2, $dir3, $dir4, $dir5);
    my ($fileInCover, $fileOutCov, $fileOutStat);

    print "**** $selection_name\n";

    my $dna = $sample_data->{result_dir};
    $origDir = "$seq_cover_folder/all/$dna";
    $dir2 = "$seq_cover_folder/InBed";
    system("mkdir $dir2") unless (-d $dir2);
    $dir3 = "$dir2/$selection_name";
    system("mkdir $dir3") unless (-d $dir3);
    $dir4 = "$dir3/$dna";
    system("mkdir $dir4") unless (-d $dir4);

    print "open bed file $bedFile\n";
    open(FHIB, $bedFile) or die("cannot open $bedFile:$!\n");

    ## ! Choice of output dir !
    $fileOutStat = "$dir4/exon_positions_${dna}_thr_${depth_thr}.stat";
    print "open out.stat file $fileOutStat\n";
    open(FHOS, ">$fileOutStat") || die("cannot create $fileOutStat");
    print FHOS "chrom\tstart\tend\tgene\texon_size\tcovered_size_bp\tmin_depth\tmax_depth\tmean_depth\tmedian_depth\tQ1_depth\tQ3_depth\tmissing_bp\tpercent_missing_bp\tR20\tR50\tR100\n";

    $bed_n = 0;
    $gene = "init";
    $n_exon = 1;

    while ($line_B = <FHIB>) {
	next if $line_B =~ /#/;
	$bed_n++;
	$line_B =~ s/\012//;
	$line_B =~ s/\015//;
	($bed_chrom, $bed_start, $bed_end, $bed_gene) = split(/\t/, $line_B);
#	($bed_chrom, $bed_start, $bed_end, $bed_gene, $bed_type) = split(/\t/, $line_B);
	$bed_chrom =~ s/chr//;
	$bed_start =~ s/ //g;
	$bed_end =~ s/ //g;
	if ($bed_start > $bed_end) {
	    my $bed_temp = $bed_end;
	    $bed_end = $bed_start;
	    $bed_start = $bed_temp;
	}
	print " bed $bed_n : $bed_chrom $bed_start $bed_end\n";

	unless (($bed_gene =~ /-exon/) || ($bed_gene =~ /_cds/)) {
	    if ($bed_gene eq $gene) {
		$n_exon ++;
	    } else {
		$gene = $bed_gene;
		$n_exon = 1;
	    }
	    $bed_gene = $bed_gene . "-exon". $n_exon;
	}

	$fileInCover = "$origDir/chr$bed_chrom.cov";
	## ! Choice of output dir !
	$fileOutCov = "$dir4/chr$bed_chrom\_$bed_start\_$bed_end.cov";

#	coverExonsPositions_region($bed_chrom, $bed_start, $bed_end, $fileInCover, $fileOutCov);
	coverExonsPositions_region($bed_chrom, $bed_start, $bed_end, $bed_gene, $fileInCover, $fileOutCov, $depth_thr);
#	coverExonsPositions_region($bed_chrom, $bed_start, $bed_end, $bed_gene, $bed_type, $fileInCover, $fileOutCov, $depth_thr);
    }
    close(FHIB);
    close(FHOS);

}

sub coverExonsPositions_region {
#   my ($bed_chrom, $bed_start, $bed_end, $fileInCover, $fileOutCov) = @_;
    my ($bed_chrom, $bed_start, $bed_end, $bed_gene, $fileInCover, $fileOutCov, $depth_thr) = @_;
#   my ($bed_chrom, $bed_start, $bed_end, $bed_gene, $bed_type, $fileInCover, $fileOutCov, $depth_thr) = @_;
    my ($line_C);
    my ($i, $cover_pos);
#   my ($in_bed_region, $missing_bp, $last_cover_pos_inbed);
    my ($in_bed_region, $exon_size, $missing_bp_count, $missing_bp, $percent_missing_bp, $last_cover_pos_inbed);
    my ($cov_depth, $min, $max, $n, %n_R, $sum, $mean, $read_number, @depths, @depths_sorted, $long, $med, $Q1, $Q3, $valmed, $valmed1, $valmed2, $valQ1, $valQ3);
    my ($bed_n);

#    if (-f $fileOutCov) {
#	print "$fileOutCov already exists\n";
#	return;
#    }

#    print "will do $fileOutCov\n";

    print "open cover file $fileInCover\n";
    open(FHIC, $fileInCover) or die("cannot open $fileInCover:$!\n");

    print "open out.cov  file $fileOutCov\n";
    open(FHOC, ">$fileOutCov") || die("cannot creat $fileOutCov");
	
    $in_bed_region = 0;
    $missing_bp_count = 0;
    $percent_missing_bp = 0;
    @depths = ();
    @depths_sorted = ();

    while ($line_C = <FHIC>) {
	next if ($line_C =~ '\#');
	chomp $line_C;
	$i = index($line_C, " ");
	$cover_pos = substr($line_C, 0, $i);
	$cov_depth = substr($line_C, $i + 1);
	$cov_depth =~ s/ //g;
	if ($cover_pos < $bed_start) { # avant la region
	    next;
	}
	if (($cover_pos >= $bed_start) and ($cover_pos <= $bed_end)) { # dans la region
	    if ($in_bed_region == 0) { # debut de la region
		$in_bed_region = 1;
		$n = 0;
		$n_R{'R20'} = 0;
		$n_R{'R50'} = 0;
		$n_R{'R100'} = 0;
		$min = 1000;
		$max = 0;
		$sum = 0;
		if ($cover_pos > $bed_start) {
		    $missing_bp_count = $cover_pos - $bed_start;
#		    print "$missing_bp_count bp not covered at beginning of region\n";
		}
	    }

	    if ($cov_depth < $depth_thr) {
#		print "depth below threshold at position $cover_pos\n";
		$missing_bp_count ++;
	    } else {
		print FHOC "$cover_pos\t$cov_depth\n";
		$n++;
		if ($cov_depth >= 20) { # nb bp with depth > 20
		    $n_R{'R20'}++;
		}
		if ($cov_depth >= 50) { # nb bp with depth > 50
		    $n_R{'R50'}++;
		}
		if ($cov_depth >= 100) { # nb bp with depth > 100
		    $n_R{'R100'}++;
		}
		$sum += $cov_depth;
		$min = $cov_depth if ($cov_depth < $min);
		$max = $cov_depth if ($cov_depth > $max);
		push(@depths, $cov_depth);
	    }

	    if ($cover_pos == $bed_end) {
#		print "end of bed reached\n";
	    }

	    $last_cover_pos_inbed = $cover_pos;
	    next;
	}
	if ($cover_pos > $bed_end) { # après la region bed
		last;
	}
    }

    if ($last_cover_pos_inbed) {
	if ($last_cover_pos_inbed < $bed_end) {
	    my $end_missing_bp = $bed_end - $last_cover_pos_inbed;
	    $missing_bp_count = $missing_bp_count + $end_missing_bp;
#	    print "$end_missing_bp bp not covered at end of region\n";
	}
    } else {
	$missing_bp_count = $bed_end - $bed_start + 1;
#	print "region not covered at all\n";
	$n = 0;
	$min = 0;
	$max = 0;
	$n_R{'R20'} = 0;
	$n_R{'R50'} = 0;
	$n_R{'R100'} = 0;
	$sum = 0;
    }

    if ($n == 0) {
	$min = 0;
    }

    unless ($sum == 0) {
	$mean = sprintf("%.1f", $sum / $n);
    } else {
	$mean = 0;
    }

    $exon_size = $bed_end - $bed_start + 1;
    $missing_bp = $exon_size - $n;
    $percent_missing_bp = sprintf("%.1f", ($missing_bp / $exon_size) * 100);

    unless ($bed_gene) {
	$bed_gene = "NA";
    }

    @depths_sorted = sort {$a <=> $b} @depths;

    $long = @depths;
    if ($long == 0) {
	$Q1 = 0;
	$med = 0;
	$Q3 = 0;
    } else {
	$Q1 = $long/4;
	$med = ($long+1)/2;
	$Q3 = $Q1*3;
    }

    if ($med == 0) {
	$valmed = 0;
    } else {
	if (int($med) == $med){
	    $valmed = $depths_sorted[$med];
	} else{
	    $valmed1 = $depths_sorted[int($med)];
	    $valmed2 = $depths_sorted[int($med)+1];
	    $valmed = ($valmed1+$valmed2)/2;
	}
    }

    if ($Q1 == 0) {
	$valQ1 = 0;
	$valQ3 = 0;
    } else {
	if (int($Q1) == $Q1){
	    $valQ1 = $depths_sorted[$Q1];
	    $valQ3 = $depths_sorted[$Q3];
	} elsif(($Q1 =~ /\.5/)||($Q1 =~ /\.75/)){
	    $valQ1 = $depths_sorted[int($Q1)+1];
	    $valQ3 = $depths_sorted[int($Q3)];
	} else {
	    $valQ1 = $depths_sorted[int($Q1)];
	    $valQ3 = $depths_sorted[int($Q3)+1];
	}
    }

    print FHOS "$bed_chrom\t$bed_start\t$bed_end\t$bed_gene\t$exon_size\t$n\t$min\t$max\t$mean\t$valmed\t$valQ1\t$valQ3\t$missing_bp\t$percent_missing_bp\t$n_R{'R20'}\t$n_R{'R50'}\t$n_R{'R100'}\n";

    close(FHIC);
    close(FHOC);

}

sub exon_covers_xls_summary {
    my (@DNA, %EXON, %PERCENT, %MISSING, %MEAN, %MIN, %MAX, %MEDIAN, %Q1, %Q3);
    my $adn = shift;
    my $result_folder = shift;
    my $run = shift;

    my $selection = $$run{"bed"};
    $selection =~ s/\/media\///;
    $selection =~ s/Data\/bed\/capture\///;
    $selection =~ s/Data\/bed\///;
    $selection =~ s/Data\/ngs\/data\/raindance\///;
    $selection =~ s/_Covered\.bed//;
    $selection =~ s/\.bed//;
    $selection =~ s/_sorted//;
    $selection =~ s/_uniq//;
    $selection =~ s/_gene_ids//;
    $selection =~ s/\//_/;
    my $depth_thr = 8;

    @DNA = @{$adn};

    #splice @DNA, 44, 4;

    print "dna=[@DNA]\n";

    my $selection_dir = "$result_folder/SEQ_cover/InBed/$selection";
    system("mkdir $selection_dir") unless (-d $selection_dir);
    my $name_run = $$run{"name"};
    my $result_dir = "$selection_dir/$name_run";
    system("mkdir $result_dir") unless (-d $result_dir);

    my $dna_workbook;
    my $dna_header_format;
    my $color1;
    my $color2;
    my $color3;
    my $color4;
    my $thr1 = 15;
    my $thr2 = 30;

    foreach my $dna(@DNA) {

	my $cover_file = "$result_folder/SEQ_cover/InBed/$selection/${dna}_hg38/exon_positions_${dna}_hg38_thr_${depth_thr}.stat";

	my ($gene, @split_gene, $transcript, $exon, $max_exon, %GENE, %TRANSCRIPT);
	$max_exon = 0;

	next unless (-e $cover_file); # in case you want to generate xls files for a partial run

	open(FILE, $cover_file) or  die ("Cannot open file $cover_file") ;
	while(my $line = <FILE>) {
	    chomp($line);
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    next if ($line =~ "chrom");
	    my @tab = split(/\t/, $line);
	    my $chrom = $tab[0];
	    my $start = $tab[1];
	    my $end = $tab[2];
	    $gene = $tab[3];
#	    my $exon_size = $tab[4];
#	    my $covered_size = $tab[5];
	    my $min_depth = $tab[6];
	    my $max_depth = $tab[7];
	    my $mean_depth = $tab[8];
	    my $median_depth = $tab[9];
	    my $Q1_depth = $tab[10];
	    my $Q3_depth = $tab[11];
	    my $missing_bp = $tab[12];
	    my $percent_missing_bp = $tab[13];
	    my $key = $chrom."-".$start."-".$end."-".$gene;

	    if ($gene =~ /-exon/) {
		@split_gene = split (/-exon/, $gene);
	    } elsif ($gene =~ /_cds_/) {
		@split_gene = split (/_cds_/, $gene);
	    }# else {
#		$split_gene[0] = $gene;
#	    }
	    $transcript = $split_gene[0];
	    $exon = $split_gene[1];
	    push @{$GENE{$transcript}}, $exon;
	    if ($#{$GENE{$transcript}} > $max_exon) {
		$max_exon = $#{$GENE{$transcript}};
	    }
	    my $coord = $chrom."-".$start."-".$end;
	    $TRANSCRIPT{$transcript, $exon} = $coord;

	    $EXON{$key} = 1;
	    $PERCENT{$dna, $key} = $percent_missing_bp;
	    $MISSING{$dna, $key} = $missing_bp;
	    $MEAN{$dna, $key} = $mean_depth;
	    $MIN{$dna, $key} = $min_depth;
	    $MAX{$dna, $key} = $max_depth;
	    $MEDIAN{$dna, $key} = $median_depth;
	    $Q1{$dna, $key} = $Q1_depth;
	    $Q3{$dna, $key} = $Q3_depth;
	}

	$dna_workbook = Spreadsheet::WriteExcel->new("$result_folder/SEQ_cover/InBed/$selection/$name_run/cover_exons_${dna}_hg38_thr_${depth_thr}.xls");
	die "Problems creating new Excel file: $!" unless defined $dna_workbook;

	$dna_header_format = $dna_workbook->add_format();
	$dna_header_format->set_bold();

	# Green
#	$color1 = $dna_workbook->add_format();
#	$color1->set_bg_color(17);
	$color1 = $dna_workbook->add_format(fg_color=>17);
	# Yellow
#	$color2 = $dna_workbook->add_format();
#	$color2->set_bg_color(13);
	$color2 = $dna_workbook->add_format(fg_color=>13);
	# Orange
#	$color3 = $dna_workbook->add_format();
#	$color3->set_bg_color(53);
	$color3 = $dna_workbook->add_format(fg_color=>53);
	# Red
#	$color4 = $dna_workbook->add_format();
#	$color4->set_bg_color(10);
	$color4 = $dna_workbook->add_format(fg_color=>10);

	# Header
	my $dna_worksheet = $dna_workbook->add_worksheet('exon_covers');
#	my $dna_worksheet = $dna_workbook->add_worksheet('amplicon_covers');
	$dna_worksheet->write(0, 0, 'Gene/Exon', $dna_header_format);
#	$dna_worksheet->write(0, 0, 'Gene/Amplicon', $dna_header_format);
	for (my $x = 1; $x<=$max_exon+1; $x++) {
	    $dna_worksheet->write(0, $x, $x, $dna_header_format);
	}

	# Values
	my $l = 1;
	my @sorted_transcripts = sort keys %GENE;
	foreach my $transcr(@sorted_transcripts) {
	    my $c = 1;
	    $dna_worksheet->write($l, 0, $transcr);
	    foreach my $x(@{$GENE{$transcr}}) {
		# Choice according to gene column in bed data
		my $key;
		if ($gene =~ /-exon/) {
		    $key = $TRANSCRIPT{$transcr, $x}."-".$transcr."-exon".$x;
		} elsif ($gene =~ /_cds_/) {
		    $key = $TRANSCRIPT{$transcr, $x}."-".$transcr."_cds_".$x; # works for data from UCSC table
		}
#		if ($MEDIAN{$dna,$key} >= 100) {
		if ($PERCENT{$dna,$key} == 0) {
#		    $dna_worksheet->write($l, $x, $TRANSCRIPT{$transcr, $x}, $color1);
		    $dna_worksheet->write($l, $c, $TRANSCRIPT{$transcr, $x}, $color1);
#		} elsif (($MEDIAN{$dna,$key} >= 8) && ($MEDIAN{$dna,$key} < 100)) {
		} elsif (($PERCENT{$dna,$key} > 0) && ($PERCENT{$dna,$key} <= $thr1)) {
#		    $dna_worksheet->write($l, $x, $TRANSCRIPT{$transcr, $x}, $color2);
		    $dna_worksheet->write($l, $c, $TRANSCRIPT{$transcr, $x}, $color2);
#		} elsif (($MEDIAN{$dna,$key} >= 8) && ($MEDIAN{$dna,$key} < 100)) {
		} elsif (($PERCENT{$dna,$key} > $thr1) && ($PERCENT{$dna,$key} <= $thr2)) {
#		    $dna_worksheet->write($l, $x, $TRANSCRIPT{$transcr, $x}, $color3);
		    $dna_worksheet->write($l, $c, $TRANSCRIPT{$transcr, $x}, $color3);
#		} elsif ($MEDIAN{$dna,$key} < 8) {
		} elsif ($PERCENT{$dna,$key} > $thr2) {
#		    $dna_worksheet->write($l, $x, $TRANSCRIPT{$transcr, $x}, $color4);
		    $dna_worksheet->write($l, $c, $TRANSCRIPT{$transcr, $x}, $color4);
		}
		$c++;
	    }
	    $l++;
	}
    }

    my $workbook = Spreadsheet::WriteExcel->new("$result_folder/SEQ_cover/InBed/$selection/$name_run/cover_exons_${selection}_thr_${depth_thr}.xls");
    die "Problems creating new Excel file: $!" unless defined $workbook;

    my $header_format = $workbook->add_format();
    $header_format->set_bold();

    # Percent missing bp
    # Header
    my $worksheet = $workbook->add_worksheet('percent_missing_bp');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    my $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    my $l = 1;
    my @sorted_exons = sort keys %EXON;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $PERCENT{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Missing bp
    # Header
    $worksheet = $workbook->add_worksheet('missing_bp');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $MISSING{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Mean depth
    # Header
    $worksheet = $workbook->add_worksheet('mean_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $MEAN{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Min depth
    # Header
    $worksheet = $workbook->add_worksheet('min_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $MIN{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Max depth
    # Header
    $worksheet = $workbook->add_worksheet('max_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $MAX{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Median depth
    # Header
    $worksheet = $workbook->add_worksheet('median_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $MEDIAN{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Q1 depth
    # Header
    $worksheet = $workbook->add_worksheet('Q1_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $Q1{$dna, $exon});
	    $c++;
	}
	$l++;
    }

    # Q3 depth
    # Header
    $worksheet = $workbook->add_worksheet('Q3_depth');
    $worksheet->write(0, 0, 'Exon/DNA', $header_format);
    $c = 1;
    foreach my $dna(@DNA) {
	$worksheet->write(0, $c, $dna, $header_format);
	$c++;
    }

    # Values
    $l = 1;
    foreach my $exon(@sorted_exons) {
	$c = 1;
	$worksheet->write($l, 0, $exon);
	foreach my $dna(@DNA) {
	    $worksheet->write($l, $c, $Q3{$dna, $exon});
	    $c++;
	}
	$l++;
    }
}


sub filter_cover_in_alternative_bed {
    my $sample_data = shift;
    my $result_folder = shift;
    my $alt_bed = shift;
    my $run_name = shift;
    my $config = shift;
    my @arr = split("/",$alt_bed);
    my $bedname = $arr[(scalar(@arr))-1];
    $bedname =~ s/.bed//;

    my $seq_cover_folder = $result_folder."/SEQ_cover";
    my $seq_cover_all_folder = $seq_cover_folder."/all/".$sample_data->{result_dir};
    my $seq_cover_in_alt_bed = $seq_cover_folder."/InBed/".$run_name;
    my $seq_cover_in_alt_bed_folder = $seq_cover_in_alt_bed."/".$bedname;
    my $seq_cover_in_alt_bed_folder_sample = $seq_cover_in_alt_bed_folder."/".$sample_data->{result_dir};
    mkdir($seq_cover_in_alt_bed);
    mkdir($seq_cover_in_alt_bed_folder);
    
    my (@chrs, @PID); 
    
    
    my @recup_exome_jobs;   # pour lancer en parallèle
    if ($sample_data->{sex} eq 'F'){
	@chrs = (1 .. 22, 'X');
    } elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
    } else {
	@chrs = (1 .. 22, 'X', 'Y');
    }
    print "start PIDs fork filter $sample_data->{dna} in alternative bed\n";
    foreach my $chr (@chrs){
	my $pid = fork();
	if($pid){
	    push (@PID, $pid);
	} elsif ($pid == 0){
	    filter_cover($sample_data, $seq_cover_in_alt_bed_folder_sample, $seq_cover_all_folder, $chr, $alt_bed);
	    exit 0;
	    } else {
		die "unable to fork: $!";
	    }	    
    }
    

    for my $pid (@PID) {
	#print "waitpid $pid\n";
	my $tmp = waitpid($pid, 0);
	print "done with pid $tmp\n";
    }
    print "END of fork PIDs filter\n";
}

1;
