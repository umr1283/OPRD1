package FilterSNP;

use Data::Dumper;
use strict;
no strict 'refs';
use DBI;

##################################################
## Filter SNP lists based on score
sub score_filter_snp {
    my $sample_data = shift;
    my $result_folder = shift;
    my $capture;
## en fonction de la capture
    if ($sample_data->{capture} =~ /raindance/) {
        print "filter SNP raindance\n";
        convert_vcf_snp_to_fscore($sample_data,$result_folder);
    } elsif ($sample_data->{capture} =~ /Haloplex/) {

        convert_vcf_snp_to_fscore($sample_data,$result_folder);

    } elsif ($sample_data->{capture} =~ /Met/){
        print "Methylation\n";                       ### rajouter sub methylation

    } else {   # exome
        score_filter_snp_exome($sample_data,$result_folder);
    }

}

##################################################
## Filter SNP lists based on score for exome
sub score_filter_snp_exome {
    my $sample_data = shift;
    my $result_folder = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

### $threshold ??

    my ($log, $score, $depth, $chrom, $pos, $n, $type, $AR, $tmp1, $tmp2, $A1, $A2, $ref_allele, $hh);
    my $thr_depth = 8;
    my @tab;
    my $F_score_folder = $result_folder."/SNP_filter/F_score";
    my $fileOut = $F_score_folder."/F_score_".$sample_data->{result_dir}.".tmp.txt";
    my $fileIn = $result_folder."/SEQ_snp/casava/".$sample_data->{result_dir}."/score3_snp_tot_sorted.txt";

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }
     if (-f $fileOut) {
	 $log = scalar(localtime())."  $fileOut already exists\n";
	 LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	 return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHO, ">$fileOut") or die("cannot create $fileOut");
    print FHO "#chrom\tpos\tAR\tA1\tA2\thh\tdepth\tscore\n";

    $log = "";
    while (my $line = <FHI>) {
	$line =~ s/\012//;
	$line =~ s/\015//;
	next if ($line =~ /^\#/);
	@tab = split(/\t/, $line);
	$score = $tab[9];
	$score =~ s/:.*//;
	$depth = $tab[8];
	next if ($depth < $thr_depth);
	$n++;

	$chrom = $tab[0];
	$pos = $tab[1];

	$type = $tab[11];
	$ref_allele = $tab[10];

	if ($type eq "SNP_diff") {
	    $AR = $ref_allele;
	    ($tmp1, $tmp2) = split(//, $tab[6]); # from casava 1.8 (both alleles in input table)
	    $A1 = $tmp1; # from casava 1.8 (both alleles in input table)#
	    $A2 = $tmp2; # from casava 1.8 (both alleles in input table)
	    $hh = 'hom';
	} elsif ($type =~ /SNP_het1/) {
	    ($tmp1, $tmp2) = split(//, $tab[6]);
	    $AR = $ref_allele;
	    $A1 = $ref_allele;
	    $A2 = $tmp2;
	    $hh = 'het';
	} elsif ($type =~ /SNP_het2/) {
	    ($tmp1, $tmp2) = split(//, $tab[6]);
	    $AR = $ref_allele;
	    $A1 = $ref_allele;
	    $A2 = $tmp1;
	    $hh = 'het';
	} elsif ($type =~ /SNP_other_het/) {
	    ($tmp1, $tmp2) = split(//, $tab[6]);
	    $AR = $ref_allele;
	    $A1 = $tmp1;
	    $A2 = $tmp2;
	    $hh = 'het2';
	} else {
	    $log.= "!! last field=$type [$line]\n";
	    next;
	}
	print FHO "$chrom\t$pos\t$AR\t$A1\t$A2\t$hh\t$depth\t$score\n";
    }

    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);

    $log = " $n writen on $fileOut\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}


##################################################
## Filter SNP lists based on score for haloplex or raindance
sub convert_vcf_snp_to_fscore {
    my $sample_data = shift;
    my $result_folder = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my $log;
    my @chrs = (1..22, 'X', 'Y');

    my $vcf_file = $result_folder."/SEQ_snp/targeted/".$sample_data->{result_dir}."/".$sample_data->{dna}."_snps.vcf";
    my $F_score_folder = $result_folder."/SNP_filter/F_score";
    my $fileOut = $F_score_folder."/F_score_".$sample_data->{result_dir}.".txt";

    if (-f $fileOut) {
        $log= scalar(localtime())."  $fileOut already exists\n";
         LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
        return;
    } else {
        $log = scalar(localtime())." creating $fileOut\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

        open (FH, $vcf_file) || die("cannot open $vcf_file");
        $log = "open $vcf_file dna=[".$sample_data->{dna}."]\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

        open (FHS, ">$fileOut") || die("cannot create $fileOut");
        $log = "open $fileOut dna=[".$sample_data->{dna}."]\n";

        print FHS "#chrom\tpos\tAR\tA1\tA2\thh\tdepth\tscore\n";

        while (my $line = <FH>) {
            $line =~ s/\012//;
            $line =~ s/\015//;
            next if ($line =~ /^\#/);

            my ($chr, $pos, $CIGAR, $ref, $alt1, $alt2, $hh, $depth, $score, $htype, $cov, $qual);

	    if ($line=~/^chr([\d\w]+)\s+(\d+)\s+(.+)\s+([A-Z-]+)\s+([A-Z,-]+)\s+([\d\.]+)\s+(.+)\s+(.*)DP=(\d+);(.*)\s+([A-Z:]+)\s+([0-9\/]+):.*$/){
		$chr = $1;
		$pos = $2;
		$ref = $4;
		$alt2 = $5;
		$htype = $12;
		if ($htype eq "0/0"){
		    $alt1 = $alt2;
		    $hh = "hom";
		} elsif (($htype eq "0/1") || ($htype eq "1/0")) {
		    $alt1 = $ref;
		    $hh = "het";
		} elsif ($htype eq "1/1"){
		    if ($alt2=~/(\w),(\w)/){
			$alt1 = $1;
			$alt2 = $2;
			$hh = "het2";
		    } else {
			$alt1 = $alt2;
			$hh = "hom";
		    }
		} elsif ($htype eq "1/2"){
		    if ($alt2=~/(\w),(\w)/){
			$alt1 = $1;
			$alt2 = $2;
			$hh = "het2";
		    }
		}

		$cov = $9;
		$qual = $6;

	    }

            if($cov >= 8) {
                print FHS "$chr\t$pos\t$ref\t$alt1\t$alt2\t$hh\t$cov\t$qual\n";
            } else {
		next;
	    }
        }
        close(FH);
        close(FHS);
    }

}



##################################################
## Merge casava and GATK snps
sub merge_fscore {
    my $disk_name = shift;
    my $run_name = shift;
    my $dna = shift;
    my $result_dir = shift;
    my $sample_data = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, my $abs_path);
    my @files;
    my %union;

    my $F_score_folder = $result_dir."/SNP_filter/F_score"; # idem F_score casava folder but files are tmp
    my $F_score_gatk_folder = $result_dir."/SEQ_snp/gatk";

    # F_scores files from casava (old or new samples)
    $log = "********** SNP **********\nSample $dna F_score files : \n  F_score : \n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    opendir (FSCORE, $F_score_folder) || die "Cannot open F_score $F_score_folder folder : $!\n";
    while (my $f = readdir(FSCORE)){
	if (($f =~/_$dna\_hg19(.*).txt$/)){    # F_score_$dna(bis?)_hg19(.tmp)?.txt
	    $log = "\t - ".$f."\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    $abs_path = $F_score_folder."/".$f;
	    push (@files, $abs_path);
	}
    }
    closedir FSCORE;

    # F_score from GATK
    $log = "  F_score GATK : \n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    opendir (GATK, $F_score_gatk_folder) || die "Cannot open F_score_GATK $F_score_gatk_folder folder : $!\n";
    while (my $f = readdir(GATK)){
	if ($f=~/gatk_snps_$dna\_hg19(.*).txt$/){ # F_score_$dna(bis?)_hg19(.tmp)?.txt
	    $log = "\t - ".$f."\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    $abs_path = $F_score_gatk_folder."/".$f;
	    push(@files, $abs_path);
	}
    }
    closedir GATK;

    foreach my $file (@files){
	open(FILE, $file) || die "Cannot open $file : $!";
	while (<FILE>){
	    next if ($_ =~ /^#/);
	    next if ($_ =~ /^M/);
	    my ($chrom, $pos, $ref, $alt1, $alt2, $rest) = split(/\t/, $_);
	    $chrom =~ s/chr//;
	    chomp($alt2);
	    unless ($union{$chrom, $pos, $ref, $alt2}) {
		$union{$chrom, $pos, $ref, $alt2} = $_;
	    }
	}
	close FILE;
    }

    # write output
    my $fileout = "$F_score_folder/F_score_".$dna."_hg19.txt";
    if (-f $fileout){
	$log = "$fileout already exists\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    } else {
	$log = "output written at $fileout\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	open(OUT, ">$fileout") || die "Cannot open $fileout : $!";
	print OUT "#chrom\tpos\tAR\tA1\tA2\thh\tdepth\tscore\n";
	foreach my $key (sort keys %union) {
	    print OUT $union{$key};
	}
	close(OUT);
    }

    # move tmp file (F_score from casava)
#   $log = "mv $F_score_folder/F_score*$dna*.tmp.txt";
    my $tmp_folder = $F_score_folder."/tmp_casava/";
    mkdir($tmp_folder) unless (-d $tmp_folder);
    system("mv $F_score_folder/F_score*$dna*.tmp.txt $tmp_folder");
#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

##################################################
## Filter SNP lists based on consequences
sub coding_filter_snp {
    my $sample_data = shift;
    my $ensembl_version = shift;
    my $result_folder = shift;
    my $filter = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $score, $depth, $chrom, $pos, $AR, $A1, $A2, $hh, $thr_consequence, $n_in, $n_out, @tmp);
    my $annot_table = '37_annot_ensembl_'.$ensembl_version.'_full_descr';
    my $F_coding_folder = $result_folder."/SNP_filter/F_coding";
    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/F_coding") unless (-d "$result_folder/SNP_filter/F_coding");
    my $fileOut = $F_coding_folder."/F_coding_".$sample_data->{result_dir}.".txt";
    my $fileIn = $result_folder."/SNP_filter/F_score/F_score_".$sample_data->{result_dir}.".txt";

    my $dbh = DBconnect::connect_db("SNP_annot");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }
     if (-f $fileOut) {
	 $log = scalar(localtime())."  $fileOut already exists\n";
	 LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	 return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    init_best_consequence($ensembl_version);
    if ($ensembl_version < 68) {
	$thr_consequence = $::sort_consequence{"SYNONYMOUS_CODING"};
    } else {
	$thr_consequence = $::sort_consequence{"synonymous_variant"};
    }

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\tENST\tENSG\tgene_symbol\tconsequence_name\tSwissprot_ID\tAA_change\tBiotype\tGene_description\tRefSeq_mRNA\tRefSeq_peptide\n";
	    next;
	}

	($chrom, $pos, $AR, $A1, $A2, $hh, $depth, $score) = split(/\t/, $line);

	if ($hh eq 'het2') {
	    @tmp = get_consequence($dbh, $chrom, $pos, $AR, $A1, $thr_consequence, $filter, $annot_table);
	    if ($tmp[0] ne 'NA') {
		print FHO "$line\t" . join("\t", @tmp) . "\n";
		$n_out++;
	    }
	    @tmp = get_consequence($dbh, $chrom, $pos, $AR, $A2, $thr_consequence, $filter, $annot_table);
	    if ($tmp[0] ne 'NA') {
		print FHO "$line\t" . join("\t", @tmp) . "\n";
		$n_out++;
	    }

	} else {
	    @tmp = get_consequence($dbh, $chrom, $pos, $AR, $A2, $thr_consequence, $filter, $annot_table);
	    if($filter == 0) {
		if (defined(scalar(@tmp))) {
		    print FHO "$line\t" . join("\t", @tmp) . "\n";
		    $n_out++;
		} else {
		    print FHO "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
		    $n_out++;
		}
	    } else {
	    	 if ($tmp[0] ne 'NA') {
		     print FHO "$line\t" . join("\t", @tmp) . "\n";
		     $n_out++;
	   	 }
	     }
	}
    }

    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);
    $dbh->disconnect();
    $log = "$n_out writen on $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

##################################################
## Filter SNP lists based on consequences
sub consequence_filter_all_transcripts_snp {
    my $sample_data = shift;
    my $ensembl_version = shift;
    my $result_folder = shift;
    my $filter = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $score, $depth, $chrom, $pos, $AR, $A1, $A2, $hh, $thr_consequence, $n_in, $n_out, %tmp);
    my $annot_table = '37_annot_ensembl_'.$ensembl_version.'_full_descr';
    my $F_coding_folder = $result_folder."/SNP_filter/F_coding";
    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/F_coding") unless (-d "$result_folder/SNP_filter/F_coding");
    my $fileOut = $F_coding_folder."/F_coding_".$sample_data->{result_dir}.".txt";
    my $fileIn = $result_folder."/SNP_filter/F_score/F_score_".$sample_data->{result_dir}.".txt";

    my $dbh = DBconnect::connect_db("SNP_annot");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }
     if (-f $fileOut) {
	 $log = scalar(localtime())."  $fileOut already exists\n";
	 LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	 return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    init_best_consequence_all_regions($ensembl_version);
    if ($ensembl_version < 68) {
	$thr_consequence = $::sort_consequence{"SYNONYMOUS_CODING"};
    } else {
	$thr_consequence = $::sort_consequence{"synonymous_variant"};
    }

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\tENST\tENSG\tgene_symbol\tconsequence_name\tSwissprot_ID\tAA_change\tBiotype\tGene_description\tRefSeq_mRNA\tRefSeq_peptide\n";
	    next;
	}

	($chrom, $pos, $AR, $A1, $A2, $hh, $depth, $score) = split(/\t/, $line);

	if ($hh eq 'het2') {
	    %tmp = %{get_consequence_all_transcripts($dbh, $chrom, $pos, $AR, $A1, $thr_consequence, $filter, $annot_table)};
	    foreach my $key (sort { $tmp{$a}[1] cmp $tmp{$b}[1] or $a cmp $b} keys %tmp) { # sort on ENSG then on ENST
		if ($tmp{$key}[3] ne 'NA') {
		    print FHO "$line\t" . join("\t", @{$tmp{$key}}) . "\n";
		    $n_out++;
		}
	    }
	    %tmp = %{get_consequence_all_transcripts($dbh, $chrom, $pos, $AR, $A2, $thr_consequence, $filter, $annot_table)};
	    foreach my $key (sort { $tmp{$a}[1] cmp $tmp{$b}[1] or $a cmp $b} keys %tmp) {
		if ($tmp{$key}[3] ne 'NA') {
		    print FHO "$line\t" . join("\t", @{$tmp{$key}}) . "\n";
		    $n_out++;
		}
	    }
	} else {
	    %tmp = %{get_consequence_all_transcripts($dbh, $chrom, $pos, $AR, $A2, $thr_consequence, $filter, $annot_table)};
	    if($filter == 0) {
		if (%tmp) {
		    foreach my $key (sort { $tmp{$a}[1] cmp $tmp{$b}[1] or $a cmp $b} keys %tmp) {
			print FHO "$line\t" . join("\t", @{$tmp{$key}}) . "\n";
			$n_out++;
		    }
		} else {
		    print FHO "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
		    $n_out++;
		}
	    } else {
		foreach my $key (sort { $tmp{$a}[1] cmp $tmp{$b}[1] or $a cmp $b} keys %tmp) {
		    if ($tmp{$key}[3] ne 'NA') {
			print FHO "$line\t" . join("\t", @{$tmp{$key}}) . "\n";
			$n_out++;
		    }
		}
	    }
	}
    }

    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);
    $dbh->disconnect();
    $log = "$n_out writen on $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

##################################################
## Get SNP consequence from local DB SNP_annot; keep one (the worse) for all transcripts
sub get_consequence {
    my ($dbh, $chrom, $pos, $A1, $A2, $thr_consequence, $filter, $annot_table) = @_;
    my ($log, $req, $sth, @row);
    my ($consequence_name, $consequence_n, $new_consequence_name, $new_consequence_n);
    my ($ENST, $ENSG, $gene_symbol, $new_ENST, $new_ENSG, $new_gene_symbol);
    my ($Swissprot_ID, $AA_change, $Gene_description, $Biotype, $refseq_mrna, $refseq_peptide);
    my ($new_Swissprot_ID, $new_AA_change, $new_Gene_description, $new_Biotype, $new_refseq_mrna, $new_refseq_peptide);
    my ($n);

    $req = "select Consequence, Transcript, Swissprot_ID, AA_change, Gene_Id, Gene_name, Biotype, Gene_description, RefSeq_mRNA, RefSeq_peptide from $annot_table where " .
	   " Chrom = \"$chrom\" " .
	   " and Start = \"$pos\" " .
	   " and Local_alleles = \"$A1/$A2\" ";
	#  " and Biotype NOT LIKE  \"%pseudogene\" " .
    $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
    $sth->execute() ||  print ("\nExecution requete: $req echec\n");

    $n = $sth->rows;

    if ($n > 0) {
	# il y a au moins une ligne en réponse
	$consequence_n = 0;
	#si on filtre
	if ($filter == 1) {
	    while(@row = $sth->fetchrow_array()) {
		$new_ENST = $row[1];
		$new_ENSG = $row[4];
		$new_Swissprot_ID = $row[2];
		$new_gene_symbol = $row[5];
		$new_AA_change = $row[3];
		$new_Biotype = $row[6];
		next if ($new_Biotype =~ /pseudogene/);
		$new_Gene_description = $row[7];
		$new_refseq_mrna = $row[8];
		$new_refseq_peptide = $row[9];
		$new_consequence_name = $row[0];
		if ($new_consequence_name =~ /TF_binging_site_variant/) { ## quick fix when several consequences separated by ','; only two cases found so far
		    $new_consequence_name = 'TF_binding_site_variant';
		}
		if ($new_consequence_name =~ /regulatory_region_variant/) {
		    $new_consequence_name = 'regulatory_region_variant';
		}
		$new_consequence_n = $::sort_consequence{$new_consequence_name};
		if ($new_consequence_n > $thr_consequence) {
		    if ($new_consequence_n > $consequence_n) {
			$consequence_n = $new_consequence_n;
			$consequence_name = $new_consequence_name;
			$ENST = $new_ENST;
			$ENSG = $new_ENSG;
			$Swissprot_ID = $new_Swissprot_ID;
			$AA_change = $new_AA_change;
			$Biotype = $new_Biotype;
			$Gene_description = $new_Gene_description;
			$gene_symbol = $new_gene_symbol;
			$refseq_mrna = $new_refseq_mrna;
			$refseq_peptide = $new_refseq_peptide;
		    }
		}
	    }
	} else {
	    while(@row = $sth->fetchrow_array()) {
		$new_ENST = $row[1];
		$new_ENSG = $row[4];
		$new_Swissprot_ID = $row[2];
		$new_gene_symbol = $row[5];
		$new_AA_change = $row[3];
		$new_Biotype = $row[6];
#  		next if ($new_Biotype =~ /pseudogene/);
		$new_Gene_description = $row[7];
		$new_refseq_mrna = $row[8];
		$new_refseq_peptide = $row[9];
		$new_consequence_name = $row[0];
		if ($new_consequence_name =~ /TF_binging_site_variant/) { ## quick fix when several consequences separated by ','; only two cases found so far
		    $new_consequence_name = 'TF_binding_site_variant';
		}
		if ($new_consequence_name =~ /regulatory_region_variant/) {
		    $new_consequence_name = 'regulatory_region_variant';
		}
		$new_consequence_n = $::sort_consequence{$new_consequence_name};
# 		if ($new_consequence_n > $thr_consequence) {
# 		    print "new consequence :$new_consequence_n\n";
# 		    print "consequence $consequence_n \n";
		if ($new_consequence_n > $consequence_n) {
		    $consequence_n = $new_consequence_n;
		    $consequence_name = $new_consequence_name;
		    $ENST = $new_ENST;
		    $ENSG = $new_ENSG;
		    $Swissprot_ID = $new_Swissprot_ID;
		    $AA_change = $new_AA_change;
		    $Biotype = $new_Biotype;
		    $Gene_description = $new_Gene_description;
		    $gene_symbol = $new_gene_symbol;
		    $refseq_mrna = $new_refseq_mrna;
		    $refseq_peptide = $new_refseq_peptide;
		}
	    }
	}

	if ($filter == 1) {
	    if ($consequence_n > 0) {
		return ($ENST, $ENSG, $gene_symbol, $consequence_name, $Swissprot_ID, $AA_change, $Biotype, $Gene_description, $refseq_mrna, $refseq_peptide);
	    } else {
		return ('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
	    }
	} else {
	    if ($consequence_n > 0) {
		return ($ENST, $ENSG, $gene_symbol, $consequence_name, $Swissprot_ID, $AA_change, $Biotype, $Gene_description, $refseq_mrna, $refseq_peptide);
	    } else {
		return ($new_ENST, $new_ENSG, $new_gene_symbol, $new_consequence_name, $new_Swissprot_ID, $new_AA_change, $new_Biotype, $new_Gene_description, $new_refseq_mrna, $new_refseq_peptide);
	    }
	}

    # !!! Le bloc ci-dessous est important: il permet de savoir si l'étape d'annotation a abouti !!!
    } else {
	# pas d'annot : ne devrait pas arriver
	print "!!! erreur : missing SNP_annot for chrom=$chrom, pos=$pos A1=$A1, A2=$A2\n";
#	$log = "!!! erreur : missing SNP_annot for chrom=$chrom, pos=$pos A1=$A1, A2=$A2\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return ('NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA');
    }
}

##################################################
## Get SNP consequence from local DB SNP_annot; keep one (the worst) for each transcript
sub get_consequence_all_transcripts {
    my ($dbh, $chrom, $pos, $A1, $A2, $thr_consequence, $filter, $annot_table) = @_;
    my ($log, $req, $sth, @row);
    my ($consequence_name, $consequence_n, $new_consequence_name, $new_consequence_n);
    my ($ENST, $ENSG, $gene_symbol, $new_ENST, $new_ENSG, $new_gene_symbol);
    my ($Swissprot_ID, $AA_change, $Gene_description, $Biotype, $refseq_mrna, $refseq_peptide);
    my ($new_Swissprot_ID, $new_AA_change, $new_Gene_description, $new_Biotype, $new_refseq_mrna, $new_refseq_peptide);
    my ($n);

    $req = "select Consequence, Transcript, Swissprot_ID, AA_change, Gene_Id, Gene_name, Biotype, Gene_description, RefSeq_mRNA, RefSeq_peptide from $annot_table where " .
	   " Chrom = \"$chrom\" " .
	   " and Start = \"$pos\" " .
	   " and Local_alleles = \"$A1/$A2\" ";
	#  " and Biotype NOT LIKE  \"%pseudogene\" " .
    $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
    $sth->execute() ||  print ("\nExecution requete: $req echec\n");

    $n = $sth->rows;

    my %transcript_consequence;

    if ($n > 0) {     # il y a au moins une ligne en réponse
	my %transcript_consequences;

	while(@row = $sth->fetchrow_array()) {
	    push @{$transcript_consequences{$row[1]}}, [@row];
	}

	#si on filtre
	if ($filter == 1) {
	    foreach my $key (keys %transcript_consequences) {
		$consequence_n = 0;
		$consequence_name = "";
		$ENST = "";
		$ENSG = "";
		$Swissprot_ID = "";
		$AA_change = "";
		$Biotype = "";
		$Gene_description = "";
		$gene_symbol = "";
		$refseq_mrna = "";
		$refseq_peptide = "";
		foreach my $tcelement_ref(@{$transcript_consequences{$key}}) {
		    my @tcelement = @$tcelement_ref;
		    $new_ENST = $tcelement[1];
		    $new_ENSG = $tcelement[4];
		    $new_Swissprot_ID = $tcelement[2];
		    $new_gene_symbol = $tcelement[5];
		    $new_AA_change = $tcelement[3];
		    $new_Biotype = $tcelement[6];
#		    next if ($new_Biotype =~ /pseudogene/);
		    $new_Gene_description = $tcelement[7];
		    $new_refseq_mrna = $tcelement[8];
		    $new_refseq_peptide = $tcelement[9];
		    $new_consequence_name = $tcelement[0];
		    if ($new_consequence_name =~ /TF_binging_site_variant/) { ## quick fix when several consequences separated by ','; only two cases found so far
			$new_consequence_name = 'TF_binding_site_variant';
		    }
		    if ($new_consequence_name =~ /regulatory_region_variant/) {
			$new_consequence_name = 'regulatory_region_variant';
		    }
		    $new_consequence_n = $::sort_consequence{$new_consequence_name};
		    if ($new_consequence_n > $thr_consequence) {
			if ($new_consequence_n >= $consequence_n) {
			    $consequence_n = $new_consequence_n;
			    $consequence_name = $new_consequence_name;
			    $ENST = $new_ENST;
			    $ENSG = $new_ENSG;
			    $Swissprot_ID = $new_Swissprot_ID;
			    $AA_change = $new_AA_change;
			    $Biotype = $new_Biotype;
			    $Gene_description = $new_Gene_description;
			    $gene_symbol = $new_gene_symbol;
			    $refseq_mrna = $new_refseq_mrna;
			    $refseq_peptide = $new_refseq_peptide;
			}
		    }
		}
		if ($consequence_n > 0) {
		    $transcript_consequence{$key} = [$ENST, $ENSG, $gene_symbol, $consequence_name, $Swissprot_ID, $AA_change, $Biotype, $Gene_description, $refseq_mrna, $refseq_peptide];
		} else {
		    $transcript_consequence{'nada'} = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'];
		}
	    }
	    return \%transcript_consequence;
	} else {
	    foreach my $key (keys %transcript_consequences) {
		$consequence_n = 0;
		$consequence_name = "";
		$ENST = "";
		$ENSG = "";
		$Swissprot_ID = "";
		$AA_change = "";
		$Biotype = "";
		$Gene_description = "";
		$gene_symbol = "";
		$refseq_mrna = "";
		$refseq_peptide = "";
		foreach my $tcelement_ref(@{$transcript_consequences{$key}}) {
		    my @tcelement = @$tcelement_ref;
		    $new_ENST = $tcelement[1];
		    $new_ENSG = $tcelement[4];
		    $new_Swissprot_ID = $tcelement[2];
		    $new_gene_symbol = $tcelement[5];
		    $new_AA_change = $tcelement[3];
		    $new_Biotype = $tcelement[6];
#		    next if ($new_Biotype =~ /pseudogene/);
		    $new_Gene_description = $tcelement[7];
		    $new_refseq_mrna = $tcelement[8];
		    $new_refseq_peptide = $tcelement[9];
		    $new_consequence_name = $tcelement[0];
		    if ($new_consequence_name =~ /TF_binging_site_variant/) { ## quick fix when several consequences separated by ','; only two cases found so far
			$new_consequence_name = 'TF_binding_site_variant';
		    }
		    if ($new_consequence_name =~ /regulatory_region_variant/) {
			$new_consequence_name = 'regulatory_region_variant';
		    }
		    $new_consequence_n = $::sort_consequence{$new_consequence_name};
#		    if ($new_consequence_n > $thr_consequence) {
# 		    print "new consequence :$new_consequence_n\n";
# 		    print "consequence $consequence_n \n";
		    if ($new_consequence_n > $consequence_n) {
			$consequence_n = $new_consequence_n;
			$consequence_name = $new_consequence_name;
			$ENST = $new_ENST;
			$ENSG = $new_ENSG;
			$Swissprot_ID = $new_Swissprot_ID;
			$AA_change = $new_AA_change;
			$Biotype = $new_Biotype;
			$Gene_description = $new_Gene_description;
			$gene_symbol = $new_gene_symbol;
			$refseq_mrna = $new_refseq_mrna;
			$refseq_peptide = $new_refseq_peptide;
		    }
		}
		if ($consequence_n > 0) {
		    $transcript_consequence{$key} = [$ENST, $ENSG, $gene_symbol, $consequence_name, $Swissprot_ID, $AA_change, $Biotype, $Gene_description, $refseq_mrna, $refseq_peptide];
		} else {
		    $transcript_consequence{$key} = [$new_ENST, $new_ENSG, $new_gene_symbol, $new_consequence_name, $new_Swissprot_ID, $new_AA_change, $new_Biotype, $new_Gene_description, $refseq_mrna, $refseq_peptide];
		}
	    }
	    return \%transcript_consequence;
	}

    # !!! Le bloc ci-dessous est important: il permet de savoir si l'étape d'annotation a abouti !!!
    } else {
	# pas d'annot : ne devrait pas arriver
	print "!!! erreur : missing SNP_annot for chrom=$chrom, pos=$pos A1=$A1, A2=$A2\n";
#	$log = "!!! erreur : missing SNP_annot for chrom=$chrom, pos=$pos A1=$A1, A2=$A2\n";
#	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	$transcript_consequence{'nada'} = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'];
	return \%transcript_consequence;
    }
}

##################################################
## Sort consequences coding regions
sub init_best_consequence {
    my $ensembl_version = shift;
    my ($consequence, $n);

    if ($ensembl_version < 68) {
	$n = 20;
	@::sort_consequence = (
	    "FRAMESHIFT_CODING",
	    "COMPLEX_INDEL",
	    "NON_SYNONYMOUS_CODING",
	    "STOP_GAINED",
	    "STOP_LOST",
	    "ESSENTIAL_SPLICE_SITE",
	    "SPLICE_SITE",
	    "WITHIN_MATURE_miRNA",
	    "REGULATORY_REGION",
	    "CODING_UNKNOWN",
	    "SYNONYMOUS_CODING",
	    "NMD_TRANSCRIPT",              # Nonsense-Mediated Decay
	    "PARTIAL_CODON",
	    "WITHIN_NON_CODING_GENE",
	    "3PRIME_UTR",
	    "5PRIME_UTR",
	    "UPSTREAM",
	    "DOWNSTREAM",
	    "INTRONIC",
	    "INTERGENIC"
	    );
    } else {
	$n = 34;
	@::sort_consequence = (
	    "frameshift_variant",
	    "missense_variant",
	    "inframe_insertion",
	    "inframe_deletion",
	    "initiator_codon_variant",
	    "stop_gained",
	    "stop_lost",
	    "splice_donor_variant",
	    "splice_acceptor_variant",
	    "splice_region_variant",
	    "mature_miRNA_variant",
	    "TF_binding_site_variant",
	    "regulatory_region_variant",
	    "TFBS_ablation",
	    "TFBS_amplification",
	    "regulatory_region_ablation",
	    "regulatory_region_amplification",
	    "coding_sequence_variant",
	    "synonymous_variant",
	    "stop_retained_variant",
	    "NMD_transcript_variant",
	    "incomplete_terminal_codon_variant",
	    "non_coding_exon_variant",
	    "nc_transcript_variant",
	    "5_prime_UTR_variant",
	    "3_prime_UTR_variant",
	    "upstream_gene_variant",
	    "downstream_gene_variant",
	    "intron_variant",
	    "transcript_ablation",
	    "transcript_amplification",
	    "feature_elongation",
	    "feature_truncation",
	    "intergenic_variant"
	    );
    }

    for $consequence (@::sort_consequence) {
	$::sort_consequence{$consequence} = $n;
	$n--;
    }
}

##################################################
## Sort consequences all regions
sub init_best_consequence_all_regions {
    my $ensembl_version = shift;
    my ($consequence, $n);

    if ($ensembl_version < 68) {
	$n = 20;
	@::sort_consequence = (
	    "FRAMESHIFT_CODING",
	    "COMPLEX_INDEL",
	    "NON_SYNONYMOUS_CODING",
	    "STOP_GAINED",
	    "STOP_LOST",
	    "ESSENTIAL_SPLICE_SITE",
	    "SPLICE_SITE",
	    "WITHIN_MATURE_miRNA",
	    "REGULATORY_REGION",
	    "CODING_UNKNOWN",
	    "NMD_TRANSCRIPT",              # Nonsense-Mediated Decay
	    "PARTIAL_CODON",
	    "WITHIN_NON_CODING_GENE",
	    "3PRIME_UTR",
	    "5PRIME_UTR",
	    "UPSTREAM",
	    "DOWNSTREAM",
	    "INTRONIC",
	    "INTERGENIC",
	    "SYNONYMOUS_CODING"
	    );
    } else {
	$n = 34;
	@::sort_consequence = (
	    "frameshift_variant",
	    "missense_variant",
	    "inframe_insertion",
	    "inframe_deletion",
	    "initiator_codon_variant",
	    "stop_gained",
	    "stop_lost",
	    "splice_donor_variant",
	    "splice_acceptor_variant",
	    "splice_region_variant",
	    "mature_miRNA_variant",
	    "TF_binding_site_variant",
	    "regulatory_region_variant",
	    "TFBS_ablation",
	    "TFBS_amplification",
	    "regulatory_region_ablation",
	    "regulatory_region_amplification",
	    "coding_sequence_variant",
	    "stop_retained_variant",
	    "NMD_transcript_variant",
	    "incomplete_terminal_codon_variant",
	    "non_coding_exon_variant",
	    "nc_transcript_variant",
	    "5_prime_UTR_variant",
	    "3_prime_UTR_variant",
	    "upstream_gene_variant",
	    "downstream_gene_variant",
	    "intron_variant",
	    "transcript_ablation",
	    "transcript_amplification",
	    "feature_elongation",
	    "feature_truncation",
	    "intergenic_variant",
	    "synonymous_variant"
	    );
    }

    for $consequence (@::sort_consequence) {
	$::sort_consequence{$consequence} = $n;
	$n--;
    }
}

##################################################
## Prepare for cosegregation step
sub launch_cosegregation_snp {
    my ($run_name, $analysis_data, $combinaison, $recess, $compos, $gene, $result_folder, $from) = @_;
    my ($log_folder, $log_file);
     if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $analysis_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $coseg_sub;
    my $dest;
    my $suffix;

    my $from_short = $from;
    $from_short =~ s/F_//;

    if($recess == 1) {
	$coseg_sub = "cosegregation_recess_filter_snp";
#	$dest = "cosegregation_recess";
	$dest = "cosegregation";
	$suffix = "_recess";
    } elsif($compos == 1) {
	$coseg_sub = "cosegregation_composites_filter_snp";
#	$dest = "cosegregation_compos";
	$dest = "cosegregation";
	$suffix = "_compos";
    } elsif($gene == 1) {
	$coseg_sub = "cosegregation_gene_filter_snp";
	$dest = "cosegregation";
	$suffix = "_gene";
    } else {
	$coseg_sub ="cosegregation_filter_snp";
	$dest = "cosegregation";
	$suffix = "";
    }

    #On parcourt tous les samples
    my @dna_cases;
    my @dna_controls;
    my @parent_dna_controls;
    my $family;
    my $sickness;
    foreach my $k (keys %$analysis_data) {
	if ($$analysis_data{$k}{'atteint'} == 1) {
	    push @dna_cases, $$analysis_data{$k}{'result_dir'};
	}
	if ($$analysis_data{$k}{'atteint'} == 0) {
	    if (($$analysis_data{$k}{'parent'} eq 'M') || ($$analysis_data{$k}{'parent'} eq 'F')) {
		push @parent_dna_controls, $$analysis_data{$k}{'result_dir'};
	    } else {
		push @dna_controls, $$analysis_data{$k}{'result_dir'};
	    }
	}
	$family = $$analysis_data{$k}{'pop'};
	$sickness = $$analysis_data{$k}{'maladie'};
    }
    my $cases = join('_',  sort @dna_cases);
    my ($controls);

    #Cas normal : cas versus controles : ecrit dans "cosegregation"
    if (@parent_dna_controls && @dna_controls) {
	$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @dna_controls);
	if ($recess == 1) {
	    ${coseg_sub}->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}", $from, $dest, $result_folder, \@dna_cases, \@parent_dna_controls, \@dna_controls);     # calling the sub using a variable
	} else {
	    my @all_dna_controls = (@parent_dna_controls,@dna_controls);
	    ${coseg_sub}->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}", $from, $dest, $result_folder, \@dna_cases, \@all_dna_controls);     # calling the sub using a variable
	}

    } elsif (@parent_dna_controls && !@dna_controls) {
	$controls = join('_', sort @parent_dna_controls);
	if ($recess == 1) {
	    ${coseg_sub}->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}", $from, $dest, $result_folder, \@dna_cases,\@parent_dna_controls, \@dna_controls);     # calling the sub using a variable
	} else {
	    ${coseg_sub}->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}", $from, $dest, $result_folder, \@dna_cases,\@parent_dna_controls);     # calling the sub using a variable
	}
    } elsif (!@parent_dna_controls && @dna_controls) {
	$controls = join('_', sort @dna_controls);
#	if ($recess == 1) {
#	    die "What the fuck : recess analysis without parents !!!";
#	} else {
	    ${coseg_sub}->($log_folder, $log_file,"from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}", $from, $dest, $result_folder, \@dna_cases,\@dna_controls);     # calling the sub using a variable
#	}
    } elsif (!@parent_dna_controls && !@dna_controls) {
#	if ($recess == 1) {
#	    die "What the fuck : recess analysis without controls !!!";
#	} else {
	    ${coseg_sub}->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${cases}${suffix}", $from, $dest, $result_folder, \@dna_cases);
#	}
    }



    #Autres cas qui teste toutes les autres combinaisons
    if($combinaison == 1) {
    	my @bitMask = ();
    	#On cherche tous les subsets possibles pour la liste d'adn
    	while (generate_mask(\@bitMask, \@dna_cases)) {
    	    my $d_cases = get_one_subset(\@bitMask, \@dna_cases);

    	    my @d_controls;
    	    my $bool = 0;

    	    #On cherche la liste complementaire
    	    foreach my $set (@dna_cases) {
    		foreach my $sub (@$d_cases) {
    		    if($sub eq $set) {
    			$bool = 1;
    		    }
    		}
    		if($bool == 0) {
    		    push @d_controls, $set;
    		} else {
    		    $bool = 0;
    		}
    	    }

    	    #si controls n'est pas vide = cas ou le subset trouve est l'ensemble de tous les adn
    	    if(@d_controls) {
    		my @cases_are_controls = @d_controls;
    		#On ajoute nos "vrais" controles
    		# a cette liste de controles
		# les parents controles
    		if(@parent_dna_controls) {
    		    push @d_controls, @parent_dna_controls;
    		}
		# et les autres controles
		if(@dna_controls) {
    		    push @d_controls, @dna_controls;
    		}
    		my $ca = join('_', sort @$d_cases);
    		my $co = join('_', sort @d_controls);
#    		print "perl $coseg_sub $name_bf_file"."from_coding_in_family_${family}_adn_${ca}_vs_${co} F_coding ${dest}_combi -cases @$d_cases -controls @d_controls > $log_dir/family_${family}_adn_${ca}_vs_${co}\n";
#    		$coseg_sub("${suffix}from_coding_in_family_${family}_adn_${ca}_vs_${co}", $from, "${dest}_combi", $result_folder, @$d_cases,@d_controls,$log_dir/family_${family}_adn_${ca}_vs_${co});
#    		$coseg_sub->("${suffix}from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}", $from, "${dest}_combi", $result_folder, $d_cases, \@d_controls);
    		$coseg_sub->($log_folder, $log_file, "from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}", $from, "${dest}_combi", $result_folder, $d_cases, \@d_controls);
    	    }
    	}
    }
}

##################################################
## Filter SNP lists based on disease status - dominant case
sub cosegregation_filter_snp {
#   my ($sample_data, $name, $src, $dest, $result_folder, $ref_cases, $ref_controls) = @_;
    my ($logFolder, $logFile, $name, $src, $dest, $result_folder, $ref_cases, $ref_controls) = @_;
    my ($log, @cases, @controls, $n_cases, %controls);
    my ($dna, $fileIn, $fileOut, $n_in, $n_out);
    my ($big_line, $no_depth, $no_score, $line, $chrom, $pos, $A1, $A2, $hh, $depth, $score, $enst, $ensg, $conseq);
    my ($key, $tmp, %snp, %snp_hh, %snp_data, %snp_ind);
    my (@tab, $title);

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    #print "launch Cosegregation_filter_snp... \n";

    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/$dest") unless (-d "$result_folder/SNP_filter/$dest");
    @cases = @{$ref_cases};
    @controls = @{$ref_controls};
    $n_cases = $#cases + 1;

    $fileOut = "$result_folder/SNP_filter/$dest/$dest" . "_$name.txt";
     if (-f $fileOut) {
	 $log= scalar(localtime())."  $fileOut already exists\n";
	 LogOut::printLog($log, $current_sub,$logFolder, $logFile);
	 return;
    }
    $n_out = 0;

    for $dna (@controls) {
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";
	if (! -f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn  [dna]=$dna\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    return;
	}
	$log="open control $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	open (FHI, $fileIn) or die("cannot open $fileIn");
	while ($line = <FHI>) {
	    next if ($line =~ /^\#/);
	    @tab = split(/\t/, $line);
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $key = $chrom . "\t" . $pos;
	    $controls{$key}++;
	}
	close(FHI);
    }

    for $dna (@cases) {
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";
	if (! -f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn  [dna]=$dna\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    return;
	}
	$log="open case $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	open (FHI, $fileIn) or die("cannot open $fileIn");
	$n_in = 0;
	while ($line = <FHI>) {
	    $n_in++;
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    chomp($line);

	    @tab = split(/\t/, $line);
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $A1 = $tab[3];
	    $A2 = $tab[4];
	    $hh = $tab[5];

	    $depth = $tab[6];
	    if($depth eq '') {
		$depth = 'NA';
	    }
	    $score = $tab[7];
	    $enst = $tab[8];
	    $ensg = $tab[9];
	    $conseq = $tab[10];
	    splice(@tab, 3, 5);

	    if ($line =~ /^\#/) {
		$title = join("\t", @tab);
	    } else {
		$key = $chrom . "\t" . $pos;
		next if(defined($controls{$key}));
		my @previous_snp_data;
		if ($snp_data{$key}) {
		    @previous_snp_data = split(/\t/,$snp_data{$key});
		} else {
		    @previous_snp_data = ();
		}
		if (($snp_ind{$key} =~ m/($dna)/) && ($chrom eq $previous_snp_data[0]) && ($pos == $previous_snp_data[1])) {
		    $snp_data{$key} = $chrom."\t".$pos."\t".$tab[2]."\t".$tab[3]."\t".$tab[4]."\t".$tab[5]."\t".$previous_snp_data[6].",".$depth."\t".$score."\t".$previous_snp_data[8].",".$enst."\t".$ensg."\t".$conseq;
		} else {
		    $snp{$key}++;
		    if ($hh eq "het2"){
			$snp_hh{$key, "het"}++;
		    } else {
			$snp_hh{$key, $hh}++;
		    }
		    $snp_data{$key} = join ("\t", @tab);
 		    $tmp = $dna."[".$A1."/".$A2.", ".$depth.", ".$score."]";
		    if (defined($snp_ind{$key})) {
			$snp_ind{$key} .= " $tmp";
		    } else {
			$snp_ind{$key} = $tmp;
		    }
		}
	    }
	}
	close(FHI);
    }

    $log= scalar(localtime())."    create $fileOut cases=[@cases] controls=[@controls]\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    print FHO "$title\tnb_ind\tnb_hom\tnb_het\tind_list[alleles, depth, score]\n";

    my $test = scalar(keys(%snp));
    $log = scalar(localtime())."    $test\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

    for $key (sort(keys(%snp))) {
    my @tab2;
	next if ($snp{$key} < $n_cases);
	$snp_hh{$key, 'hom'} = 0 unless (defined($snp_hh{$key, 'hom'}));
	$snp_hh{$key, 'het'} = 0 unless (defined($snp_hh{$key, 'het'}));
	print FHO "$snp_data{$key}";

 @tab2 = split(/\t/, $snp_data{$key});
  if (($tab2[11] eq '') && ($tab2[12] eq '')){
    print FHO "\tNA\tNA";
    }
  if (($tab2[11] ne '') && ($tab2[12] eq '')){
    print FHO "\tNA";
    }
  #
  print FHO "\t$snp{$key}";
	print FHO "\t$snp_hh{$key, 'hom'}";
	print FHO "\t$snp_hh{$key, 'het'}";
	print FHO "\t$snp_ind{$key}\n";
	$n_out++;
    }
    close(FHO);
    $log= scalar(localtime())."    $n_out writen on $fileOut, $n_in read from $fileIn\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
}

##################################################
## Filter SNP lists based on disease status - recessive case
sub cosegregation_recess_filter_snp {
    my ($logFolder, $logFile, $name, $src, $dest, $result_folder, $ref_cases, $ref_parent_controls, $ref_other_controls) = @_;
    my ($log, @cases, @parent_controls, @other_controls, $n_cases, $n_parent_controls, $n_other_controls, %p_controls, %np_controls);
    my ($dna, $fileIn, $fileOut, $n_in, $n_out);
    my ($big_line, $no_depth, $no_score, $line, $chrom, $pos, $A1, $A2, $hh, $depth, $score, $enst, $ensg, $conseq);
    my ($key, $tmp, %snp, %snp_hh, %snp_data, %snp_ind);
    my (@tab, $title);

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    #print "launch Cosegregation_recess_filter_snp... \n";

    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/$dest") unless (-d "$result_folder/SNP_filter/$dest");
    @cases = @{$ref_cases};
    @parent_controls = @{$ref_parent_controls};
    @other_controls = @{$ref_other_controls};
    $n_cases = $#cases + 1;
    $n_parent_controls = $#parent_controls + 1;
    $n_other_controls = $#other_controls + 1;

    $fileOut = "$result_folder/SNP_filter/${dest}/${dest}_${name}.txt";
     if (-f $fileOut) {
	 $log= scalar(localtime())."  $fileOut already exists\n";
#	 print $log;
	 LogOut::printLog($log, $current_sub, $logFolder, $logFile);
#	 LogOut::printLog($log, $current_osub, $sample_data->{run}, $sample_data->{dna});
	 return;
    }

    $n_out = 0;

    # Parents
    for (my $i = 0; $i < $n_parent_controls; $i++) {
	$dna = $parent_controls[$i];
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";

	if (!-f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
#       LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    print $log;
	    return;
	}
	$log="open control $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	#print "open control $fileIn dna=$dna\n";
	open (FHI, $fileIn) or die("cannot open $fileIn");
	while ($line = <FHI>) {
	    next if ($line =~ /^\#/);
	    @tab = split(/\t/, $line);
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $hh = $tab[5];
	    next if ($hh eq 'hom'); # cannot be as hom in parents
	    $key = $chrom . "\t" . $pos;
#	    $key = $chrom . "\t" . $pos . "\t" . $hh;
	    $p_controls{$key}++;
	}
	close(FHI);
    }

    # Other controls
    if ($n_other_controls > 0) { # il y a des controles en plus des parents
	for (my $i = 0; $i < $n_other_controls; $i++) {
	    $dna = $other_controls[$i];
	    $fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";

	    if (!-f $fileIn) {
		$log= "!!missing file $fileIn\n";
		LogOut::printLog($log, $current_sub, $logFolder, $logFile);
		return;
	    }
	    $log =  "open control $fileIn dna=$dna\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    open (FHI, $fileIn) or die("cannot open $fileIn");
	    while ($line = <FHI>) {
		next if ($line =~ /^\#/);
		@tab = split(/\t/, $line);
		$chrom = $tab[0];
		$pos = $tab[1];
		$hh = $tab[5];
		next if ($hh eq 'het'); # can be as het in non parent controls
		$key = $chrom . "\t" . $pos;
#	        $key = $chrom . "\t" . $pos . "\t" . $hh;
		$np_controls{$key}++;
	    }
	    close(FHI);
	}
    }

    for $dna (@cases) {
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";
	if (! -f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn\n";
	    #print $log;
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    return;
	}
	$log="open case $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	open (FHI, $fileIn) or die("cannot open $fileIn");
	$n_in = 0;
	while ($line = <FHI>) {
	    $n_in++;
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    chomp($line);

	    @tab = split(/\t/, $line);
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $A1 = $tab[3];
	    $A2 = $tab[4];
	    $hh = $tab[5];
	    next if ($hh =~ /het/);

	    $depth = $tab[6];
	    if($depth eq '') {
		$depth = 'NA';
	    }
	    $score = $tab[7];
	    $enst = $tab[8];
	    $ensg = $tab[9];
	    $conseq = $tab[10];
	    splice(@tab, 3, 5);

	    if ($line =~ /^\#/) {
		$title = join("\t", @tab);
	    } else {
		$key = $chrom . "\t" . $pos;
		unless ($n_parent_controls == 0) {
		    next if(!defined($p_controls{$key})); # mutation has to be in parents (as het)
		    next if($p_controls{$key} < $n_parent_controls); # mutation has to be in all sequenced parents (as het)
		}
		unless ($n_other_controls == 0) {
		    next if ((defined($np_controls{$key})) && ($np_controls{$key} > 0)); # mutation cannot be as hom in other controls
		}
		$snp{$key}++;
		$snp_hh{$key, $hh}++;
		$snp_data{$key} = join ("\t", @tab);
		$tmp = $dna."[".$A1."/".$A2.", ".$depth.", ".$score."]";
		if (defined($snp_ind{$key})) {
		    $snp_ind{$key} .= " $tmp";
		} else {
		    $snp_ind{$key} = $tmp;
		}
	    }
	}
	close(FHI);
    }

    $log= scalar(localtime())."    create $fileOut cases=[@cases] parents=[@parent_controls] other controls=[@other_controls]\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    print FHO "$title\tnb_ind\tnb_hom\tnb_het\tind_list[alleles, depth, score]\n";

#   my $test = scalar(keys(%snp));
#   $log .= scalar(localtime())."    $test\n";

    for $key (sort(keys(%snp))) {
      my @tab2;
	next if ($snp{$key} < $n_cases);
	$snp_hh{$key, 'hom'} = 0 unless (defined($snp_hh{$key, 'hom'}));
	$snp_hh{$key, 'het'} = 0 unless (defined($snp_hh{$key, 'het'}));
	print FHO "$snp_data{$key}";

  @tab2 = split(/\t/, $snp_data{$key});
   if (($tab2[11] eq '') && ($tab2[12] eq '')){
     print FHO "\tNA\tNA";
     }
   if (($tab2[11] ne '') && ($tab2[12] eq '')){
     print FHO "\tNA";
     }
   #

	print FHO "\t$snp{$key}";
	print FHO "\t$snp_hh{$key, 'hom'}";
	print FHO "\t$snp_hh{$key, 'het'}";
	print FHO "\t$snp_ind{$key}\n";
	$n_out++;
    }
    close(FHO);
    $log= scalar(localtime())."    $n_out writen on $fileOut, $n_in read from $fileIn\n";
    #print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
}

##################################################
## Filter SNP lists based on disease status - dominant case - gene centric
sub cosegregation_gene_filter_snp {
    my ($logFolder, $logFile, $name, $src, $dest, $result_folder, $ref_cases, $ref_controls) = @_;
    my ($log, @cases, @controls, $n_cases, %controls_gene);
    my ($dna, $fileIn, $fileOut, $n_in, $n_out);
    my ($line, $chrom, $pos, $A1, $A2, $depth, $score, $ensg, $conseq, $swiss_id, $aa_change);
    my (%snp, %snp_hh, $snp_data, %snp_ind);
    my ($key_gene, %gene, %gene_data, %gene_ind);
    my (@tab, $title);

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    #print "launch Cosegregation_gene_filter_snp... \n";

    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/$dest") unless (-d "$result_folder/SNP_filter/$dest");
    @cases = @{$ref_cases};
    @controls = @{$ref_controls};
    $n_cases = $#cases + 1;

    $fileOut = "$result_folder/SNP_filter/$dest/$dest" . "_$name.txt";
     if (-f $fileOut) {
	 $log= scalar(localtime())."  $fileOut already exists\n";
	 LogOut::printLog($log, $current_sub,$logFolder, $logFile);
	 return;
    }
    $n_out = 0;

    for $dna (@controls) {
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";
	if (! -f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    return;
	}
	$log="open control $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	open (FHI, $fileIn) or die("cannot open $fileIn");
	while ($line = <FHI>) {
	    next if ($line =~ /^\#/);
	    @tab = split(/\t/, $line);
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $key_gene = $tab[9];
	    $controls_gene{$key_gene}++;
	}
	close(FHI);
    }

    for $dna (@cases) {
	$fileIn = "$result_folder/SNP_filter/$src/$src" . "_$dna.txt";
	if (! -f $fileIn) {
	    $log= scalar(localtime())."  !!!missing file $fileIn\n";
	    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	    return;
	}
	$log="open case $fileIn  [dna]=$dna\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	open (FHI, $fileIn) or die("cannot open $fileIn");
	$n_in = 0;
	while ($line = <FHI>) {
	    $n_in++;
	    $line =~ s/\012//;
	    $line =~ s/\015//;
	    chomp($line);

	    @tab = split(/\t/, $line);
	    unless ($line =~/#/) {
		next if ($tab[28] > 0.01); # mutations rares seulement
	    }
	    $chrom = $tab[0];
	    $pos = $tab[1];
	    $A1 = $tab[3];
	    $A2 = $tab[4];
	    $depth = $tab[6];
	    if($depth eq '') {
		$depth = 'NA';
	    }
	    $score = $tab[7];
	    $ensg = $tab[9];
	    $conseq = $tab[11];
	    $swiss_id = $tab[12];
	    $aa_change = $tab[13];
	    splice(@tab, 1, 8);
	    splice(@tab, 3, 3);
	    splice(@tab, 5);

	    if ($line =~ /^\#/) {
		$title = join("\t", @tab);
	    } else {
		$key_gene = $ensg;
		next if(defined($controls_gene{$key_gene}));
		$gene{$key_gene}++;
		unless ($dna ~~ @{$gene_ind{$key_gene}}) {
		    push(@{$gene_ind{$key_gene}}, $dna);
		}
		my $snp_data = "[".$pos.", ".$A1."/".$A2.", ".$depth.", ".$score.", ".$conseq.", ".$swiss_id.", ".$aa_change."]";
		$snp{$dna,$key_gene} .= "$snp_data";
		$gene_data{$key_gene} = join ("\t", @tab);
	    }
	}
	close(FHI);
    }

    $log= scalar(localtime())."    create $fileOut cases=[@cases] controls=[@controls]\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

   print FHO "$title\tind_list(snp_list[pos, alleles, depth, score, consequence, swissprot_id, aa_change])\n";

    foreach my $k (sort(keys(%gene))) {
	next if ($#{$gene_ind{$k}} + 1 != $n_cases);
	print FHO "$gene_data{$k}\t";
	my $i = 0;
	foreach my $dude(@cases) {
	    $i++;
	    unless ($i == $n_cases) {
		print FHO "$dude($snp{$dude,$k}),";
	    } else {
		print FHO "$dude($snp{$dude,$k})";
	    }
	}
	print FHO "\n";
	$n_out++;
    }
    close(FHO);
    $log= scalar(localtime())."    $n_out writen on $fileOut, $n_in read from $fileIn\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
}

##################################################
## Get one subset for combinations
sub get_one_subset {
    my ($bitMask, $set) = @_;
    my @oneset;
    for (0 .. @$bitMask-1) {
	if ($bitMask->[$_] == 1) {
	    push (@oneset, $set->[$_]);
	}
    }
    return \@oneset;
}

##################################################
## Generate mask for combinations
sub generate_mask {
    my ($bitMask, $set) = @_;

    my $i;
    for ($i = 0; $i < @$set && $bitMask->[$i]; $i++) {
	$bitMask->[$i] = 0;
    }

    if ($i < @$set) {
	$bitMask->[$i] = 1;
	return 1;
    }
    return 0;
}

############ cp fscore from GATK if Hiseq4000 (no merge with casava)
sub cp_fscore {
    my $disk_name = shift;
    my $run_name = shift;
    my $dna = shift;
    my $result_dir = shift;
    my $sample_data = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, my $abs_path);


    my $F_score_folder = $result_dir."/SNP_filter/F_score";
    my $F_score_gatk_folder = $result_dir."/SEQ_snp/gatk";

    my $F_score_file = $F_score_folder."/F_score_".$dna."_hg38.txt";
    my $F_score_gatk_file = $F_score_gatk_folder."/gatk_snps_".$dna."_hg38.txt";

    system("cp $F_score_gatk_file $F_score_file");

    $log = "********** SNP **********\nSample $dna \n copy from F_score gatk file : $F_score_gatk_file to F_score : $F_score_file\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}


1;
