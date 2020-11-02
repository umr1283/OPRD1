package Tables;


use Data::Dumper;
use Spreadsheet::WriteExcel::Big;


sub build_tab_files{
    my $dna_list = shift;
    my $run_data = shift;
    my $report_folder = shift;
    my $result_folder = shift;
    my $run_name = shift;
    my $run_nas = shift;
    my $bed_type = shift;
    my $samtools = shift;
    my $sexDet = shift;
    my $thr = 8;



    $stat_cov_in_capture_file = $report_folder."/stat_cov_in_capture_".$run_name.".txt";
    launch_stat_cov_in_capture($dna_list, $stat_cov_in_capture_file, $thr, $run_data, $result_folder, $run_name, $bed_type, $samtools, $sexDet);
    # stat_mut : nb snps, nb indels, comp MA
    my $stat_mut_file = $report_folder."/stat_mut_".$run_name.".txt";
    launch_stat_mut($dna_list, $stat_mut_file, $thr, $run_data, $result_folder, $run_name, $run_nas);

    # stat_ontarget : nb total reads, nb reads mapped, nb reads on target
    my $stat_ontarget_file = $report_folder."/stat_reads_".$run_name.".txt";

    # stats to xls
    my $stats_to_xls_file =  $report_folder."/stat_QC_".$run_name.".xls";
    convert_to_xls($stat_cov_in_capture_file, $stat_mut_file, $stat_ontarget_file, $stats_to_xls_file);

}

sub launch_stat_cov_in_capture {
    my $dna_list = shift;
    my $stat_cov_in_capture_file = shift;
    my $thr = shift;
    my $run_data = shift;
    my $dir = shift;
    my $run_name = shift;
    my $bed_type = shift;
    my $samtools = shift;
    my $sexDet = shift;

    my ($log, $base, $dna, %dna, %data, $k, @k, %title, %log, $bp, $region, %bp_total);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my $bedname = "";
    my $fileOut = $stat_cov_in_capture_file;
    $log .= "will create file $fileOut\n";
#   LogOut::printLogStat($log, $current_sub, $run_name, %log);

    my @thr = (8, 20, 50, 100);

    for $dna (@{ $dna_list }) {
 	($data{$dna, 'meanX'}, $data{$dna, 'chrX_meanX'},$data{$dna, 'chrY_meanX'}, $data{$dna, 'aut_meanX'}, $data{$dna, 'mediane'}, $data{$dna, 'Q1'}, $data{$dna, 'Q3'}, $log{$dna}) = get_mean_median_cover($stat_cover_in_capture, $$run_data{$dna}->{result_dir}, \@thr, $dir, $bed_type, $run_name, $bedname, $samtools, $$run_data{$dna}->{bed_file}, $dna, $sexDet);
	($bp_total{$dna}, $bp{$dna}, $log{$dna}) = get_bp($stat_cover_in_capture, $$run_data{$dna}->{result_dir}, \@thr, $dir, $bed_type, $run_name, $bedname);
 	($region{$dna}, $log{$dna}) = get_region($stat_cov_in_capture, $$run_data{$dna}->{result_dir}, \@thr, $dir, $bed_type, $run_name, $bedname);
    }
    #print Dumper(\%data);
    #print Dumper(\%region);
    #print Dumper(\$bp_total);
    #print Dumper(\%bp);

### Ajout Mehdi column sex 250219

    $log .= "write $fileOut\n";
    #LogOut::printLogStat($log, $current_sub, $run_name, %log);
    open(FHO, ">$fileOut") or die("cannot create $fileOut:$!");
    print FHO "dna\tmethod\tlane\tsex\tmeanX\tchrX_meanX\tSRY_meanX\tsex_ratio\tflag_SRY\taut_meanX\tmediane\tQ1\tQ3\tMbp total";

    for my $i (@thr) {
	print FHO "\tMbp>=$i\t\%bp>=$i";
    }
    print FHO "\tnb_region_total";
    for my $i (@thr) {
	print FHO "\tnb_regions_meanX >=$i\t\%regions_meanX >=$i";
    }
    print FHO "\n";
    for $dna (@{$dna_list}) {
	print FHO "$dna\t$$run_data{$dna}->{method}\t$$run_data{$dna}->{lane}\t$$run_data{$dna}->{sex}";
	print FHO "\t$data{$dna,'meanX'}\t$data{$dna,'chrX_meanX'}\t$data{$dna,'chrY_meanX'}";

  my $ratio_Y = ($data{$dna,'chrY_meanX'} / $data{$dna,'meanX'});
  my $flag_Y = "NA";
  if ((($ratio_Y > 1) && ($$run_data{$dna}->{sex} eq "M")) || (($ratio_Y < 1) && ($$run_data{$dna}->{sex} eq "F"))){
    $flag_Y = "OK";
  }
  else {
    $flag_Y = "WARNING";
  }
  print FHO "\t$ratio_Y\t$flag_Y";

  print FHO "\t$data{$dna,'aut_meanX'}\t$data{$dna,'mediane'}\t$data{$dna,'Q1'}\t$data{$dna,'Q3'}";
	print FHO "\t$bp_total{$dna}";
	for my $i (@thr){
	    my $bp_Mb = sprintf("%.2f", ($bp{$dna}{$i}[1]) / 1000000);
	    my $bp_pct = sprintf("%.2f", ($bp{$dna}{$i}[2]));
	    print FHO "\t$bp_Mb\t$bp_pct";
	}
	print FHO "\t$region{$dna}{'all'}[1]";
	for my $i (@thr){
	    my $reg_pct = sprintf("%.2f", ($region{$dna}{$i}[2]));
	    print FHO "\t$region{$dna}{$i}[1]\t$reg_pct";
	}
	print FHO "\n";

    }


#    LogOut::printLogStat($log, $current_sub, $run_name, %log);

}

sub get_region {
    my ($base, $dna, $thr, $dir, $bed_type, $run_name, $bedname) = @_;
    my ($file, $line, @tab, $region, $region_pct, $log);
    my %regions;

    if ($bed_type eq "0"){
	$file = "$dir/SEQ_cover/InCapture/$dna/stat_regions";
    } elsif ($bed_type eq "1") {
	$file = "$dir/SEQ_cover/InBed/$run_name/$bedname/$dna/stat_regions";
    } elsif ($bed_type eq "2") {
	$file = "$dir/SEQ_cover/all/$dna/stat_regions";
    }
    $log.= $file."\n";

    if (-f $file) {
	open(FHI, $file) or die("cannot open $file:$!");
	<FHI>;
	while($line = <FHI>) {
	    undef @tab;
	    chomp $line;
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    @tab = split(/\t/, $line);
	    @{$regions{$tab[0]}} = @tab;
	}
	close(FHI);
	#print Dumper(\%regions);
    } else {
	$log.= "no file $file\n";
    }

    return (\%regions, $log);
}


sub get_bp {
    my ($base, $dna, $thr, $dir, $bed_type, $run_name, $bedname) = @_;
    my ($file, $line, $deep, $tmp1, $tmp2, $bp_total, $bp_deep, $bp_deep_pct,$log, @tab);
    my %bp;

    if ($bed_type eq "0"){
	$file = "$dir/SEQ_cover/InCapture/$dna/deep_nbp_cumul";
    } elsif ($bed_type eq "1"){
	$file = "$dir/SEQ_cover/InBed/$run_name/$bedname/$dna/deep_nbp_cumul";
    } elsif ($bed_type eq "2") {
	$file = "$dir/SEQ_cover/all/$dna/deep_nbp_cumul";
    }
    $bp_total = "NA";

    if (-f $file) {
	open(FHI, $file) or die("cannot open $file:$!");
	<FHI>;
	while($line = <FHI>) {
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    @tab = split(/\t/, $line);
	    if ($tab[0] == 1) {
		$bp_total = $tab[1];
	    }
	    @{$bp{$tab[0]}} = @tab;
	}
	close(FHI);
    } else {
	$bp_total = 0;
	$log.= "no file $file\n";
    }
    unless ($bp_total eq "NA") {
	$bp_total = sprintf("%.2f", $bp_total / 1000000);
    }
  #  unless ($bp_deep8 eq "NA") {
#	$bp_deep = sprintf("%.2f", $bp_deep / 1000000);
 #   }
    print $bp_total."\n";
    return ($bp_total, \%bp, $log);

}


sub get_mean_median_cover {
    my ($base, $dna, $thr, $dir, $bed_type, $run_name, $bedname, $samtools, $bed_file, $dna_name, $sexDet) = @_;
    my ($file, %data, $meanX, $median, $Q1, $Q3, $log,$chrX_meanX,$chrY_meanX,$aut_meanX,$chr);

    if ($bed_type eq "0"){
	$file = "$dir/SEQ_cover/InCapture/$dna/meanX.txt";
    } elsif ($bed_type eq "1"){
	$file = "$dir/SEQ_cover/InBed/$run_name/$bedname/$dna/meanX.txt";
    } elsif ($bed_type eq "2") {
	$file = "$dir/SEQ_cover/all/$dna/meanX.txt";
    }

   if (! -f $file) {
       $meanX = "NA";
       $chrX_meanX = "NA";
       $chrY_meanX = "NA";
       $aut_meanX = "NA";
       $median = "NA";
       $Q1 = "NA";
       $Q3 = "NA";
    } else {
    	$log.= "open $file\n";

	open(FILE, "$file") or  die ("Erreur d'ouverture du fichier : cannot open file $file $!\n") ;
	my (%data, $nb);

	#Recuperation des stats pour cet adn
	while(my $line = <FILE>) {
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    chomp($line);
	    if ($line=~/^all/){
		($nb, $data{'meanX'} , $data{'mediane'} , $data{'Q1'}, $data{'Q3'}) = split ("\t", $line);
	    } else {
		next;
	    }
	}
	close FILE;

	#chrX, chrY and aut meanX calculation
	#my $in_file = "/media/Run/$run_name/alignement/gatk/Sample_$dna_name/$dna_name\_RG.dedup.sorted.realigned.recal.bam";
	#my @result = `$samtools depth -q30 -Q37 -a -b $bed_file $in_file | awk -f $sexDet`;
	#foreach my $i (@result){
	 #   print "result : ".$i."\n";
	 #   chomp $i;
	 #   @res2 = split("\t",$i);
	 #   foreach my $j (@res2){
		#print "second split : ".$j."\n";
		#$data{$res2[0]} = $res2[1];
	    #}
	   #
	#}


	if ($bed_type eq "0"){
		$fichier = "$dir/SEQ_cover/InCapture/$dna/CoverX_Y.txt";
    } elsif ($bed_type eq "1"){
		$fichier = "$dir/SEQ_cover/InBed/$run_name/$bedname/$dna/CoverX_Y.txt";
    } elsif ($bed_type eq "2") {
		$fichier = "$dir/SEQ_cover/all/$dna/CoverX_Y.txt";
    }

	open (FICHIER, $fichier);
	while($line = <FICHIER>) {
		#print "result : ".$line."\n";
		chomp $line;
		@donnees = split("\t",$line);
		foreach my $line2 (@donnees){
			#print "second split : ".$line2."\n";
			$data{$donnees[0]} = $donnees[1];
	   }
	}
	close FICHIER;

	if (defined($data{'meanX'})) {
	    $meanX = $data{'meanX'};
	} else {
	    $meanX = "NA";
	}

	if (defined($data{'xCoverage'})) {
	    $chrX_meanX = $data{'xCoverage'};
	} else {
	    $chrX_meanX = "NA";
	}

	if (defined($data{'yCoverage'})) {
	    $chrY_meanX = $data{'yCoverage'};
	} else {
	    $chrY_meanX = "NA";
	}

	if (defined($data{'autCoverage'})) {
	    $aut_meanX = $data{'autCoverage'};
	} else {
	    $aut_meanX = "NA";
	}

	if (defined($data{'mediane'})) {
	    $median = $data{'mediane'};
	} else {
	    $median = "NA";
	}

	if (defined($data{'Q1'})) {
	    $Q1 = $data{'Q1'};
	} else {
	    $Q1 = "NA";
	}

	if (defined($data{'Q3'})) {
	    $Q3 = $data{'Q3'};
	} else {
	    $Q3 = "NA";
	}
    }
    #print "$meanX, $chrX_meanX, $chrY_meanX, $median, $Q1, $Q3, $log";
    return ($meanX, $chrX_meanX, $chrY_meanX, $aut_meanX, $median, $Q1, $Q3, $log);

}

sub launch_stat_mut {
    my $dna_list = shift;
    my $fileOut = shift; # $stat_mut_file
    my $thr = shift;
    my $run_data = shift;
    my $dir = shift; # $result_folder
    my $run_name = shift;
    my $run_nas = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($base, $mutMeth, $dna, %dna, %data, $k, @k, %title, $log, %log);

    $log .= "will create file $fileOut\n";
#   LogOut::printLogStat($log, $current_sub, $run_name, %log);

    for $dna (@{$dna_list}) {
    	my $capture = $$run_data{$dna}->{capture};
	if($capture=~  /raindance/) {
	    $mutMeth = 'targeted';
	} elsif($capture eq 'Agilent_Haloplex') {
	    $mutMeth = 'targeted';
	} else {
	    #if (($run_nas eq 'Hiseq3') || ($run_nas eq 'NextSeq')){
		$mutMeth = 'gatk';
	    #} else {
		#$mutMeth = 'casava';
	    #}
	}
	$data{$dna, 'indels'} = get_indels($mutMeth, $$run_data{$dna}, $dir);
	$data{$dna, 'rs'} = get_snps($mutMeth, $$run_data{$dna}, $dir);

	next if ($$run_data{$dna}->{MA} eq '0');
	($data{$dna, 'MA_found_pct'},
	 $data{$dna, 'MA_ok_pct'},
	 $data{$dna, 'MA_homo_found_pct'},
	 $data{$dna, 'MA_homo_ok_pct'},
	 $data{$dna, 'MA_het_found_pct'},
	 $data{$dna, 'MA_het_ok_pct'}) = get_MA_cmp($dna, $dir);
    }

    @k = ('rs', 'indels',
	  'MA_found_pct', 'MA_ok_pct',
	  'MA_homo_found_pct', 'MA_homo_ok_pct',
	  'MA_het_found_pct', 'MA_het_ok_pct'
	  );

    %title = (
	       'rs' => "nb de snps detectes",
	       'indels' => "nb d'insertions/deletions detectes",
	       'MA_found_pct' => "\% snps MA trouve",
	       'MA_ok_pct' => "\% snps MA ok",
	       'MA_homo_found_pct' => "\% homo snps MA trouve",
	       'MA_homo_ok_pct' => "\% homo snps MA ok",
	       'MA_het_found_pct' => "\% het snps MA trouve",
	       'MA_het_ok_pct' => "\% het snps MA ok"
	  );

    $log .= "write $fileOut\n";
#   LogOut::printLogStat($log, $current_sub, $run_name, %log);
    open(FHO, ">$fileOut") or die("cannot create $fileOut:$!");
    print FHO "dna\tmethod\tlane";
    for $k (@k) {
	print FHO "\t$title{$k}";
    }
    print FHO "\n";
    for $dna (@{$dna_list}) {
	print FHO "$dna\t$$run_data{$dna}->{method}\t$$run_data{$dna}->{lane}";
	for $k (@k) {
	    if (defined($data{$dna, $k})) {
		print FHO "\t$data{$dna, $k}";
	    } else {
		print FHO "\tNA";
	    }
	}
	print FHO "\n";
    }

    #LogOut::printLogStat($log, $current_sub, $run_name, %log);

    close(FHO);
}

sub get_indels {
    my ($mutMeth, $dna_data, $dir) = @_;
    my ($indels);
    my ($file, $line, %data, @tab);

    if($mutMeth eq 'targeted') {
	$file = "$dir/SEQ_indel/$mutMeth/$dna_data->{result_dir}/$dna_data->{dna}_indels.vcf";
    } elsif ($mutMeth = 'gatk'){
	$file = "$dir/SEQ_indel/$mutMeth/gatk_indels_$dna_data->{dna}_hg38.txt";
    } else {
	$file = "$dir/SEQ_indel/$mutMeth/$dna_data->{result_dir}/score3_indels_tot_sorted.txt";
    }


    if (! -f $file) {
	$indels = "NA";
    } else {
	open(FHI, $file) or die("cannot open $file:$!");
	while($line = <FHI>) {
	    next if ($line =~ /^\#/);
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    $indels++;
	}
	close(FHI);
    }
    return $indels;
}

sub get_snps {
    my ($mutMeth, $dna_data, $dir) = @_;
    my ($snps);

    my ($file, $line, %data, @tab);

	if($mutMeth eq 'targeted') {
		$file = "$dir/SEQ_snp/$mutMeth/$dna_data->{result_dir}/$dna_data->{dna}_snps.vcf";
	} elsif ($mutMeth = 'gatk') {
	    $file = "$dir/SEQ_snp/$mutMeth/gatk_snps_$dna_data->{dna}_hg38.txt";
	} else {
		if($dna =~ 'sida') {
		    $file = "$dir/SEQ_snp/$mutMeth/$dna_data->{result_dir}/score3_snp_tot_sorted_ori.txt";
		} else {
		    $file = "$dir/SEQ_snp/$mutMeth/$dna_data->{result_dir}/score3_snp_tot_sorted.txt";
		}
   	}

    if (! -f $file) {
	$snps = "NA";
    } else {
	open(FHI, $file) or die("cannot open $file:$!");
	while($line = <FHI>) {
	    next if ($line =~ /^\#/);
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    $snps++;
	}
	close(FHI);
    }
    return $snps;

}

sub get_MA_cmp {
    my ($dna, $dir) = @_;
    my ($file, $line, %data, @tab);

    $file = "$dir/MA_stat/$dna";

    if (! -f $file) {
	print "missing $file\n";
	@tab = ("NA", "NA", "NA", "NA", "NA", "NA");
    } else {
	print "open $file\n";
	open(FHI, $file) or die("cannot open $file:$!");
	while($line = <FHI>) {
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    @tab = split(/\t/, $line);
	    $data{$tab[0]} = $tab[1];
	}
	close(FHI);
	if (defined($data{'pct_geno_found'})) {
	    @tab = (
		    $data{'rs_seq'},
		    $data{'pct_geno_found'},
		    $data{'pct_geno_ok'},
		    $data{'pct_homo_found'},
		    $data{'pct_homo_ok'},
		    $data{'pct_het_found'},
		    $data{'pct_het_ok'})
	} else {
	    @tab = ("NA", "NA", "NA", "NA", "NA", "NA");
	}
    }
    return @tab;
}

sub convert_to_xls {
    my $stat_cov_in_capture_file = shift;
    my $stat_mut_file = shift;
    my $stat_ontarget_file = shift;
    my $stats_to_xls_file = shift;

    my ($workbook, $fileName, $sheetNamelong, $sheetName, $worksheet);
    my (@files);
    my ($row, $col, $line, $el, @tab, $i);

    @files = ($stat_cov_in_capture_file , $stat_mut_file, $stat_ontarget_file);

    $workbook = Spreadsheet::WriteExcel::Big->new("$stats_to_xls_file");

    foreach $fileName (@files) {
	my (@tmp) = split(/\//, $fileName);
	$sheetNamelong = pop(@tmp);
	if ($sheetNamelong =~ /:/) {
	    $sheetNamelong =~ s/^.*://;
	    $fileNamelong =~ s/:.*$//;
	    print "filename $fileName => sheetname $sheetNamelong\n";
	}

	#$sheetName =~ s/\.txt//;
	#if (length($sheetName) > 31) {
	#    print "!!!! sheet=$sheetName > 31 caracters\n";
	#}

	if ($sheetNamelong =~/^([\w-_]+)_run_([\w\d-_]+)/ ){
	    $sheetName = $1;
	}

	print "work on file $fileName, sheet=$sheetName\n";

	$worksheet = $workbook->add_worksheet($sheetName);

	open (FH, $fileName) || die ("cannot open $fileName\n");
	$row = 0;
	$col = 0;
	while ($line = <FH>) {
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    @tab = split(/\t/, $line);

	    if ($row == 0 ) {
		foreach $el (@tab) {
		    $worksheet->write($row, $col, $el);
		    $col++;
		}
		$row++;
	    } else {
		$col = 0;
		$el = shift @tab;
		$worksheet->write($row, $col, $el);
		$col++;
		foreach $el (@tab) {
		    $worksheet->write($row, $col, $el);
		    $col++;
		}
		$row++;
	    }
	}
	close(FH);
    }
    $workbook->close() or die "Error closing fileName: $!";
}

sub stat_cov_in_alt_bed {
    my $dna_list = shift;
    my $run_data = shift;
    my $report_folder = shift;
    my $result_folder = shift;
    my $run_name = shift;
    my $run_nas = shift;
    my $alt_bed = shift;
    my $bed_type = shift;
    my $samtools = shift;
    my $sexDet = shift;

    my $dir = $result_folder;
    my @thr = (8, 20, 50, 100);
    my @arr = split("/",$alt_bed);
    my $bedname = $arr[(scalar(@arr))-1];
    $bedname =~ s/.bed//;
    my $result_bed_folder = $report_folder."/".$bedname;
    #print "Report folder : ".$result_bed_folder."\n";
    my $stat_cov_in_capture_file = "$result_bed_folder/stat_cov_in_bed.txt";
    mkdir($result_bed_folder);

    my ($log, $base, $dna, %dna, %data, $k, @k, %title, %log, $bp_total);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my $fileOut = $stat_cov_in_capture_file;
    $log .= "will create file $fileOut\n";
#   LogOut::printLogStat($log, $current_sub, $run_name, %log);

    for $dna (@{ $dna_list }) {
 	($data{$dna, 'meanX'}, $data{$dna, 'chrX_meanX'}, $data{$dna, 'chrY_meanX'}, $data{$dna, 'mediane'}, $data{$dna, 'Q1'}, $data{$dna, 'Q3'}, $log{$dna}) = get_mean_median_cover($stat_cover_in_capture, $$run_data{$dna}->{result_dir}, $thr, $dir, $bed_type, $run_name, $bedname, $samtools, $$run_data{$dna}->{bed_file}, $dna, $sexDet);
	($bp_total{$dna}, $bp{$dna}, $log{$dna}) = get_bp($stat_cover_in_capture, $$run_data{$dna}->{result_dir}, \@thr, $dir, $bed_type, $run_name, $bedname);
 	($region{$dna}, $log{$dna}) = get_region($stat_cov_in_capture, $$run_data{$dna}->{result_dir}, \@thr, $dir, $bed_type, $run_name, $bedname);
    }

    $log .= "write $fileOut\n";
    LogOut::printLogStat($log, $current_sub, $run_name, %log);
    open(FHO, ">$fileOut") or die("cannot create $fileOut:$!");
    print FHO "dna\tmethod\tlane\tmeanX\tchrX_meanX\tchrY_meanX\taut_meanX\tmediane\tQ1\tQ3\tMbp total";

    for my $i (@thr) {
	print FHO "\tMbp>=$i\t\%bp>=$i";
    }
    for my $i (@thr) {
	print FHO "\tnb_regions_meanX >=$i\t\%regions_meanX >=$i";
    }
    print FHO "\n";
    for $dna (@{$dna_list}) {
	print FHO "$dna\t$$run_data{$dna}->{method}\t$$run_data{$dna}->{lane}";
	print FHO "\t$data{$dna,'meanX'}\t$data{$dna,'chrX_meanX'}\t$data{$dna,'chrY_meanX'}\t$data{$dna,'aut_meanX'}\t$data{$dna,'mediane'}\t$data{$dna,'Q1'}\t$data{$dna,'Q3'}";
	print FHO "\t$bp_total{$dna}";
	for my $i (@thr){
	    my $bp_Mb = sprintf("%.2f", ($bp{$dna}{$i}[1]) / 1000000);
	    my $bp_pct = sprintf("%.2f", ($bp{$dna}{$i}[2]));
	    print FHO "\t$bp_Mb\t$bp_pct";
	}
	for my $i (@thr){
	    my $reg_pct = sprintf("%.2f", ($bp{$dna}{$i}[2]));
	    print FHO "\t$region{$dna}{$i}[1]\t$reg_pct";
	}
	print FHO "\n";

    }

  #  LogOut::printLogStat($log, $current_sub, $run_name, %log);

}

1;
