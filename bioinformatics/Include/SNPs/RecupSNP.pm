package RecupSNP;


use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use LogOut;
use Projet;
use List::Util qw/first/;


sub recup_all_snp {
    my $sample_data = shift;
    my $result_folder = shift;
    
    ## en fonction de la capture
    if ($sample_data->{capture} =~ /raindance/){
 	print "recup SNPs raindance\n";                         
	recup_all_snp_targeted($sample_data, $result_folder);

    } elsif ($sample_data->{capture} =~ /Haloplex/){
	print "recup SNPs Haloplex\n";                          
	recup_all_snp_targeted($sample_data, $result_folder);

    } elsif ($sample_data->{capture} =~ 'Met'){

	print "Methylation\n";                       ### rajouter sub methylation

    } else {
	print "recup SNPs Exome\n";
	recup_all_snp_exome($sample_data, $result_folder);

    }    
}


sub recup_all_snp_exome {
    my $sample_data = shift;
    my $result_folder = shift;
    my ($log);
    my (@tab, @chrs);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)    
    
    my $result_dir = $result_folder."/SEQ_snp/casava/".$sample_data->{result_dir};
    mkdir($result_dir);

    if ($sample_data->{sex} eq 'F'){
	@chrs = (1 .. 22, 'X');
    } elsif ($sample_data->{sex} eq 'M'){
	@chrs = (1 .. 22, 'X', 'Y');
    } else {
	@chrs = (1 .. 22, 'X', 'Y');
    }
    
    $log.=scalar(localtime())."  [dna]= $sample_data->{dna}, [sex]= $sample_data->{sex}, [run]=$sample_data->{run}\n";
   
    my $sample_folder = $sample_data->{disk}."/".$sample_data->{run}."/variation/casava/".$sample_data->{sample_dir}; 
    
    ### récupérer parsed_dir
    my @dir = Projet::listDir($sample_folder);
    my $parsed_dir = first { $_ =~ /^Parsed/ } @dir;
    
    my $fileSNP =  $result_dir."/score3_snp.txt";
    my $fileSNP_removed =  $result_dir."/score3_snp_removed.txt";
    my $fileSNP_tot =  $result_dir."/score3_snp_tot.txt";
    my $fileSNP_tot_sorted =  $result_dir."/score3_snp_tot_sorted.txt";
    
    if (-f $fileSNP) {
	$log.= scalar(localtime())."  $fileSNP already exists\n";
	 LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    } else {
	$log.="$0 $sample_data->{dna} version_casava $sample_data->{casava_version}\n";
	$log.="create $fileSNP\n";
	open(FHO, ">$fileSNP") or die("cannot create $fileSNP");
	print FHO "#chrom\tposition\tA\tC\tG\tT\tmodified_call\ttotal\tused\tscore\treference\ttype\n";
	for my $chr (@chrs){
	    my $fileSeq = $sample_folder."/".$parsed_dir."/chr".$chr.".fa/snps.txt";
	    #print "$fileSeq\n";
	    $log.=  "read $fileSeq\n";
	    open(FHI, $fileSeq) or die("cannot open $fileSeq :$!\n");
	    <FHI>;
	    while (my $line = <FHI>) {
		next if (($sample_data->{sex} eq 'M') && (($chr eq 'X') || ($chr eq 'Y')) && ($line =~ 'het'));
		next if $line =~ /^\#/;
		chomp $line;
		$line =~ s/\012//;
		$line =~ s/\015//;
		@tab = split(/\t/, $line);
		my $type;
		
		if ((substr($tab[6],0,1) eq $tab[4]) && (substr($tab[6],1,1) eq $tab[4])) {
		    next;
		} elsif ((substr($tab[6],0,1) eq $tab[4]) && (substr($tab[6],1,1) ne $tab[4])) {
		    $type = 'SNP_het1';
		} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) eq $tab[4])) {
		    $type = 'SNP_het2';
		} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) ne $tab[4]) && (substr($tab[6],0,1) eq substr($tab[6],1,1))) {
		    $type = 'SNP_diff';
		} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) ne $tab[4]) && (substr($tab[6],0,1) ne substr($tab[6],1,1))) {
		    $type = 'SNP_other_het';
		}
		
		my $totalbcalls = $tab[2] + $tab[3];
		print FHO "$chr\t$tab[1]\t$tab[10]\t$tab[11]\t$tab[12]\t$tab[13]\t$tab[6]\t$totalbcalls\t$tab[2]\t$tab[5]\t$tab[4]\t$type\n";
				
	    }
	    close FHI;
	}
	close FHO;
	$log.= scalar(localtime())."  $fileSNP created\n\n";
    }

    if (-f $fileSNP_removed) {
	$log.= scalar(localtime())."  $fileSNP_removed already exists\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    } else {
	$log.= "will create $fileSNP_removed\n";
	$log.= "$0 $sample_data->{dna} version_casava $sample_data->{casava_version}\n";
	$log.= "create $fileSNP_removed\n";
	open(FHO, ">$fileSNP_removed") or die("cannot create $fileSNP_removed :$!\n");
	print FHO "#chrom\tposition\tA\tC\tG\tT\tmodified_call\ttotal\tused\tscore\treference\ttype\n";
	for my $chr (@chrs) {
	    my $origDir = $sample_folder."/".$parsed_dir."/chr".$chr.".fa";
	    #print $origDir."\n";
	    opendir(FHD, $origDir) or die("cannot open origDir $origDir:$!");
	    while (my $tmp = readdir(FHD)) {
		next if ( $tmp eq "." );
		next if ( $tmp eq ".." );
		next unless ($tmp =~ /^\d+$/);
		my $fileSeq_removed = "$origDir/$tmp/snps.removed.txt";
		if (! -f $fileSeq_removed) {
		    $log.= scalar(localtime())."  $0 : missing $fileSeq_removed\n";
		    next;
		} else {
		    $log.= " read $fileSeq_removed\n";
		    open(FHI, $fileSeq_removed) or die("cannot open $fileSeq_removed");
		    <FHI>;
		    while (my $line = <FHI>) {
			next if (($sample_data->{sex} eq 'M') && (($chr eq 'X') || ($chr eq 'Y')) && ($line =~ 'het'));
			next if $line =~ /^\#/;
			chomp $line;
			$line =~ s/\012//;
			$line =~ s/\015//;
			@tab = split(/\t/, $line);
			my $type;
			
			if ((substr($tab[6],0,1) eq $tab[4]) && (substr($tab[6],1,1) eq $tab[4])) {
			    next;
			} elsif ((substr($tab[6],0,1) eq $tab[4]) && (substr($tab[6],1,1) ne $tab[4])) {
			    $type = 'SNP_het1';
			} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) eq $tab[4])) {
			    $type = 'SNP_het2';
			} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) ne $tab[4]) && (substr($tab[6],0,1) eq substr($tab[6],1,1))) {
			    $type = 'SNP_diff';
			} elsif ((substr($tab[6],0,1) ne $tab[4]) && (substr($tab[6],1,1) ne $tab[4]) && (substr($tab[6],0,1) ne substr($tab[6],1,1))) {
			    $type = 'SNP_other_het';
			}
			
			my $totalbcalls = $tab[2] + $tab[3];
			
			print FHO "$chr\t$tab[1]\t$tab[10]\t$tab[11]\t$tab[12]\t$tab[13]\t$tab[6]\t$totalbcalls\t$tab[2]\t$tab[5]\t$tab[4]\t$type\n";
			
		    }
		    close FHI;
		}
	    }
	}
	close FHO;
	$log.= scalar(localtime())."  $fileSNP_removed created\n\n";
    }

    system("cat $fileSNP $fileSNP_removed > $fileSNP_tot");
    open(FHOT, ">$fileSNP_tot_sorted") or die("cannot create $fileSNP_tot_sorted: $!\n");
    print FHOT "#chrom\tposition\tA\tC\tG\tT\tmodified_call\ttotal\tused\tscore\treference\ttype\n";
    close FHOT;
    for my $i (1..22,"X","Y") {
	system("grep -P '^$i\t' $fileSNP_tot | sort -k 2 -n | uniq >> $fileSNP_tot_sorted");
    }
    
    $log.= scalar(localtime())."  $fileSNP_tot and \n $fileSNP_tot_sorted created!!!\n\n";
    ### printLog
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
}

#################################################################
# convert gatk vcf file to F_score and store into SEQ_snp folder
sub recup_gatk {
    my $run_disk = shift;
    my $run_name = shift;
    my $adn = shift;
    my $result_folder = shift;
    my $sample_data = shift;
    my ($log);
    my $dna = $adn;
    #$dna =~s/_/-/g;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    
    my @chrs = (1..22, 'X', 'Y');

    #my $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/".$sample_data->{sample_dir}."/".$sample_data->{dna}."_snps.vcf";
    my $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/Sample_".$dna."/".$dna."_snps.vcf";
    my $gatk_snp_folder = $result_folder."/SEQ_snp/gatk";
    unless (-d "$gatk_snp_folder") { system("mkdir $gatk_snp_folder"); }
    my $fileOut = $gatk_snp_folder."/gatk_snps_".$sample_data->{result_dir}.".txt";

    if (-f $fileOut) {
        $log= scalar(localtime())."  $fileOut already exists\n";
         LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
        return;
    } else {
        $log = scalar(localtime())." creating $fileOut\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

        open (FH, $vcf_file) || die("cannot open vcf file $vcf_file \n");
        $log = "open $vcf_file dna=[".$sample_data->{dna}."]\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

        open (FHS, ">$fileOut") || die("cannot create $fileOut \n");
        $log = "open $fileOut dna=[".$sample_data->{dna}."]\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

        print FHS "#chrom\tpos\tAR\tA1\tA2\thh\tdepth\tscore\n";

        while (my $l = <FH>) {
	    chomp $l;
	    next if ($l =~ /^#/);
	    if ($l=~/^chr([\d\w]+)\s+(\d+)\s+(.+)\s+([A-Z-]+)\s+([A-Z,-]+)\s+([\d\.]+)\s+(.+)\s+(.*)DP=(\d+);(.*)\s+([A-Z:]+)\s+([0-9\/]+):.*$/){ 
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
		## enlever chromosomes mitochondriaux ##
		if ($chr ne "M"){
		    print FHS "$chr\t$pos\t$ref\t$alt1\t$alt2\t$hh\t$cov\t$qual\n";
		}
	    } else {
		$log= "!!! check vcf file line :\n\"   $l  \"\n";
		LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	    }
        }
        close(FH);
        close(FHS);
    }
}

###############################################################
## recup vcf file from haloplex/raindance
sub recup_all_snp_targeted {
    my $sample_data = shift;
    my $result_folder = shift;
    my ($log);
    my (@tab, @chrs);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my $orig_dir =  $sample_data->{disk}."/".$sample_data->{run}."/variation/bwa/Sample_".$sample_data->{dna};

    my $result_dir = $result_folder."/SEQ_snp/targeted/";
    my $sample_result_dir = $result_dir."/".$sample_data->{result_dir};

    system("mkdir $result_dir") unless (-d "$result_dir");
    system("mkdir $sample_result_dir") unless (-d "$sample_result_dir");

    # Split snp/indels vcf file
    my $file_snps = $orig_dir."/".$sample_data->{dna}."_snps.vcf";
    $log = "SNPs from $orig_dir/".$sample_data->{dna}.".raw.vcf to $file_snps";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    if(!(-f $file_snps)) {
        system('sed \'/INDEL/d\''." $orig_dir/".$sample_data->{dna}.".raw.vcf > $file_snps");
	system("cp $file_snps $sample_result_dir");
    } else {
        $log = "The file $file_snps already exists.\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	system("cp $file_snps $sample_result_dir");
        return;
    }
}



1;

