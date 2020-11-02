package RecupINDEL;


use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use LogOut;
use Projet;
use List::Util qw/first/;
use List::Util qw (max);

sub recup_all_indel {
    my $sample_data = shift;
    my $result_folder = shift;
    
    ## en fonction de la capture
    if ($sample_data->{capture} =~ /raindance/){
 	print "raindance\n";                        
	recup_all_indel_targeted($sample_data, $result_folder);

    } elsif ($sample_data->{capture} =~ /Agilent_Haloplex/){
	print "Haloplex\n";                          
	recup_all_indel_targeted($sample_data, $result_folder);

    } elsif ($sample_data->{capture} =~ 'Met'){

	print "Methylation\n";                       ### rajouter sub methylation

    } else {
	
	recup_all_indel_exome($sample_data, $result_folder);

    }
}


sub recup_all_indel_exome {
    my $sample_data = shift;
    my $result_folder = shift;
    my ($log);
    my (@tab, @chrs);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)    
    
    my $result_dir = $result_folder."/SEQ_indel/casava/".$sample_data->{result_dir};
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
    
    my $fileINDEL =  $result_dir."/score3_indels.txt";
    my $fileINDEL_removed =  $result_dir."/score3_indels_removed.txt";
    my $fileINDEL_tot =  $result_dir."/score3_indels_tot.txt";
    my $fileINDEL_tot_sorted =  $result_dir."/score3_indels_tot_sorted.txt";
    
    if (-f $fileINDEL) {
	$log.= scalar(localtime())."  $fileINDEL already exists\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    } else {
	$log.="$0 $sample_data->{dna} version_casava $sample_data->{casava_version}\n";
	$log.="create $fileINDEL\n";
	open(FHO, ">$fileINDEL") or die("cannot create $fileINDEL");
	print FHO "#\$ COLUMNS seq_name pos type ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) depth alt_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count\n";
	for $chr (@chrs) {
	    my $fileSeq = $sample_folder."/".$parsed_dir."/chr".$chr.".fa/indels.txt";
	    $log.=  "read $fileSeq\n";
	    open(FHI, $fileSeq) or die("cannot open $fileSeq");
	    while (my $line = <FHI>) {
		next if ($line =~ /^\#/);
		next if (($sample_data->{sex} eq 'M') && (($chr eq 'X') || ($chr eq 'Y')) && ($line =~ 'het'));
		$line =~ s/chr//;
		$line =~ s/\.fa//;
		print FHO "$line";
	    }
	    close FHI;
	}
	close FHO;
	$log.= scalar(localtime())."  $fileINDEL created\n\n";
    }
    
    if (-f $fileINDEL_removed) {
	$log.= "$fileINDEL_removed already exists\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    } else {
	$log.= "will create $fileINDEL_removed\n";
	$log.= "$0 $sample_data->{dna} version_casava $sample_data->{casava_version}\n";
	$log.= "create $fileINDEL_removed\n";
	open(FHO, ">$fileINDEL_removed") or die("cannot create $fileINDEL_removed");
	print FHO "#\$ COLUMNS seq_name pos type ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) depth alt_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count\n";
	for $chr (@chrs) {
	    my $origDir = $sample_folder."/".$parsed_dir."/chr".$chr.".fa";
	    opendir(FHD, $origDir) or die("cannot open origDir $origDir:$!");
	    while (my $tmp = readdir(FHD)) {
		next if ( $tmp eq "." );
		next if ( $tmp eq ".." );
		next unless ($tmp =~ /^\d+$/);
		my $fileSeq_removed = "$origDir/$tmp/indels.removed.txt";
		if (! -f $fileSeq_removed) {
		    $log.= scalar(localtime())."  $0 : missing $fileSeq_removed\n";
		    next;
		} else {
		    $log.= " read $fileSeq_removed\n";
		    open(FHI, $fileSeq_removed) or die("cannot open $fileSeq_removed : $!\n");
		    while (my $line = <FHI>) {
			next if ($line =~ /^\#/);
			next if (($sample_data->{sex} eq 'M') && (($chr eq 'X') || ($chr eq 'Y')) && ($line =~ 'het'));
			$line =~ s/chr//;
			$line =~ s/\.fa//;
			$line =~ s/\012//;
			$line =~ s/\015//;
			print FHO "$line\n";
		    }
		    close(FHI);
		}
	    }
	}
	close FHO;
	$log.= scalar(localtime())."  $fileINDEL_removed created\n\n";
    }

    system("cat $fileINDEL $fileINDEL_removed > $fileINDEL_tot"); # sorting foireux; a ameliorer
    open(FHOT, ">$fileINDEL_tot_sorted") or die("cannot create $fileINDEL_tot_sorted :$!\n");
    print FHOT "#\$ COLUMNS seq_name pos type ref_upstream ref/indel ref_downstream Q(indel) max_gtype Q(max_gtype) depth alt_reads indel_reads other_reads repeat_unit ref_repeat_count indel_repeat_count\n";
    close(FHOT);

    for my $i (1..22,"X","Y") {
	system("grep -P '^$i\t' $fileINDEL_tot | sort -k 2 -n | uniq >> $fileINDEL_tot_sorted");
    }
    $log.= scalar(localtime())."  $fileINDEL_tot and \n $fileINDEL_tot_sorted created!!!\n\n";
    ### printLog
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}

################################################################
## recup indels from haloplex 
sub recup_all_indel_targeted {
    my $sample_data = shift;
    my $result_folder = shift;
    my ($log);
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my $orig_dir =  $sample_data->{disk}."/".$sample_data->{run}."/variation/bwa/Sample_".$sample_data->{dna};

    my $result_dir = $result_folder."/SEQ_indel/targeted/";
    my $sample_result_dir = $result_dir."/".$sample_data->{result_dir};

    system("mkdir $result_dir") unless (-d "$result_dir");
    system("mkdir $sample_result_dir") unless (-d "$sample_result_dir");

    my $file_indels = $orig_dir."/".$sample_data->{dna}."_indels.vcf";
    if(!(-f $file_indels)) {
        $log = "SNPs from $orig_dir/".$sample_data->{dna}.".raw.vcf to $file_indels";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
        system('grep -e \'^#\|INDEL\'' . " $orig_dir/".$sample_data->{dna}.".raw.vcf > $file_indels");
	system("cp $file_indels $sample_result_dir");
    } else {
        $log = "le fichier $file_indels existe deja.\n";
        LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	system("cp $file_indels $sample_result_dir");
        return;
    }
}


#################################################################
# convert gatk vcf file to F_score and store into SEQ_indel folder
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
    
    #my $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/".$sample_data->{sample_dir}."/".$sample_data->{dna}."_indels.vcf";
    my $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/Sample_".$dna."/".$dna."_indels.vcf";

    my $gatk_indel_folder = $result_folder."/SEQ_indel/gatk";
    unless (-d "$gatk_indel_folder") { system("mkdir $gatk_indel_folder"); }
    my $fileOut = $gatk_indel_folder."/gatk_indels_".$sample_data->{result_dir}.".txt";

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

        print FHS "#chrom\tpos\tCIGAR\tAR\tA1\tA2\thh\tdepth\tscore\n";

        while (my $l = <FH>) {
	    chomp $l;
	    next if ($l =~ /^#/);
	    if ($l=~/^chr([\d\w]+)\s+(\d+)\s+(.+)\s+([A-Za-z-]+)\s+([A-Za-z,-]+)\s+([\d.]+)\s+(.+)\s+(.*)DP=(\d+);(.*)\s+([A-Z:]+)\s+([0-9\/]+):.*/){ 
		$n_line++;
		$chr = $1;
		$pos = $2;
		$pos++;
		$ref = $4;
		@ref = split(//,$ref);
		$alt2 = $5;
		$htype = $12;
		if ($htype eq "0/0"){
		    $alt1 = $alt2;
		    $hh = "hom";
		} elsif (($htype eq "0/1") || ($htype eq "1/0")) {
		    $alt1 = $ref;
		    $hh = "het";
		} elsif ($htype eq "1/1"){
		    if ($alt2=~/(\w+),(\w+)/){
			$alt1 = $1;
			$alt2 = $2;
			$hh = "het2";
			#print $1."\t".$2."\n";
		    } else {
			$alt1 = $alt2;
			$hh = "hom";
		    }
		} elsif ($htype eq "1/2"){
		    if ($alt2=~/(\w+),(\w+)/){
			$alt1 = $1;
			$alt2 = $2;
			$hh = "het2";
			#print $1."\t".$2."\n";
		    }
		}
		@alt1 = split(//, $alt1);
		@alt2 = split(//, $alt2);
		shift @alt1; shift @alt2; shift @ref;
		$max = max(scalar(@ref), scalar(@alt1), scalar(@alt2));
		
		if (scalar(@ref) == $max){ 
		    $i1 = $max - scalar(@alt1); # nombre de caracteres manquantes par rapport a la sequence la plus longue
		    for (0..$i1-1){
			$r1 .="-";
		    }
		    push(@alt1,$r1);
		    $i2 = $max - scalar(@alt2);
		    for (0..$i2-1){
			$r2 .="-";
		    }
		    push(@alt2,$r2);
		    $i = max($i1,$i2);
		    $cigar = $i."D";
		} else {
		    if (scalar(@alt1)==$max ){
			$i1 = $max - scalar(@ref);
			for (0..$i1-1){
			    $r1 .="-"; 
			}
			push(@ref,$r1);
			$i2 = $max - scalar(@alt2);
			for (0..$i2-1){
			    $r2 .="-";
			}
			push(@alt2,$r2);
			$cigar = $i1."I";
		    } 
		    elsif (scalar(@alt2)==$max ){
			$i1 = $max - scalar(@ref);
			for (0..$i1-1){
			    $r1 .="-"; #ajouter "-" pour ref
			}
			push(@ref,$r1);
			$i2 = $max - scalar(@alt1);
			for (0..$i2-1){
			    $r2 .="-";
			}
			push(@alt1,$r2);
			
			$cigar = $i1."I";
		    }
		}
		$cov = $9;
		$qual = $6;
		$r = $r1 = $r2 = "";
		if ($chr ne "M"){
		    print FHS "$chr\t$pos\t$cigar\t".join("",@ref)."\t".join("",@alt1)."\t".join("",@alt2)."\t$hh\t$cov\t$qual\n";
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

1;
