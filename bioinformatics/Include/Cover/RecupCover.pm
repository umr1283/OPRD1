package RecupCover;


use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use Projet;
use LogOut;
use List::Util qw/first/;


sub recup_all_cover {
    my $sample_data = shift;
    my $result_folder = shift;
    my $sequencer = shift;
    my $seq_cover_all_folder = $result_folder."/SEQ_cover/all/".$sample_data->{result_dir};
    mkdir($seq_cover_all_folder);
    print $sample_data->{capture}."\n";
    my (@chrs, @PID); 
    
    ## en fonction de la capture
    if ($sample_data->{capture} =~ /raindance/){
 	print "*** sample raindance\n";  ### raindance / haloplex => pas besoin parallelisation
        recup_all_cover_targeted($sample_data, $seq_cover_all_folder);
    } elsif ($sample_data->{capture} =~ /Haloplex/){
	print "*** sample Haloplex\n";
	recup_all_cover_targeted($sample_data, $seq_cover_all_folder);
    } else {
	# si hiseq4000
	#if (($sequencer =~/Hiseq3/)||($sequencer =~/NEXTSEQ/)){
	    print "recup cover from gatk with $sequencer...\n";
	    recup_cover_from_gatk($sample_data, $seq_cover_all_folder);
    }
    
    
}


sub recup_all_cover_exome {
    my $sample_data = shift;
    my $seq_cover_all_folder = shift;
    my $chr = shift;
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    
    my $log;
    my ($tmp, $fileIn, $cov, $pos);
    my (@tab);
    
    $log = "[dna] = $sample_data->{dna}, [chr] = $chr\n";
    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );

    my $sample_folder = $sample_data->{disk}."/".$sample_data->{run}."/variation/casava/".$sample_data->{sample_dir};
    
    ### récupérer parsed_dir
    my @dir = Projet::listDir($sample_folder);
    my $parsed_dir = first { $_ =~ /^Parsed/ } @dir;
    my $input_folder = $sample_folder."/".$parsed_dir."/chr".$chr.".fa";
    
    my $file_out = $seq_cover_all_folder."/chr".$chr.".cov";
    if (-f $file_out){
	 $log = scalar(localtime())." : $file_out already exists\n";
	 LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );
	 return;
    }

    open(FHO, ">$file_out") or die("cannot open $file_out:$!");
    print FHO "\# ILL coverage file for chr$chr\n";
    print FHO "\# reference-pos coverage\n";

    opendir(FHD, $input_folder) or die("cannot open origDir $input_folder:$!");
    while ($tmp = readdir(FHD)) {
	next if ( $tmp eq "." );
	next if ( $tmp eq ".." );
	next unless ($tmp =~ /^\d+$/);
	system("gunzip $input_folder/$tmp/sites.txt.gz");
	$fileIn = "$input_folder/$tmp/sites.txt";

	if (-f $fileIn) {
	    $log = scalar(localtime())." : reading $fileIn\n";
	    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );  
	    open(FHI, $fileIn) or die("cannot open $fileIn:$!");

	    if ($sample_data->{casava_version} =~ /^1.8/) {
		# sauter 8 lignes avec la version 1.8
		for $tmp ( 1 .. 8) {
		    <FHI>;
		}
	    }
    
	    while (my $line = <FHI>) {
		@tab = split(/\t/,$line);
		$cov = $tab[2];
		$cov =~ s/\s+//g;
		if ($cov != 0) {
		    $pos = $tab[1];
		    $pos =~ s/\s+//g;
		    print FHO "$pos $cov\n";
		}
	    }
	    close(FHI);
	} else {
	    $log = scalar(localtime())." !!! missing $fileIn\n";
	    LogOut::printLogByChr( $log, $current_sub, $sample_data->{run}, $sample_data->{dna}, $chr );  
	}
	
	system("gzip $input_folder/$tmp/sites.txt");
    }
    closedir FHD;
    close FHO;

}

sub recup_all_cover_targeted {
    my $sample_data = shift;
    my $seq_cover_all_folder = shift;
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my $cov_file = $sample_data->{disk}."/".$sample_data->{run}."/variation/bwa/".$sample_data->{sample_dir}."/".$sample_data->{dna}.".cov";
    my @chrs = (1..22, 'X', 'Y');
    my $n_out = 0;
    my ($log, $last_chrom);

    $log = scalar(localtime())." open file $cov_file\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    open(FHI, $cov_file) or die("cannot open file $cov_file :$!\n");
    while (my $line = <FHI>) {
        $line =~ s/\012//;
        $line =~ s/\015//;
        chomp $line;

        next if ($line =~ /^#/);
        my ($chrom, $pos, $cov) = split(/\t/, $line);
        my $chr = substr($chrom, 3);
        next if ("@chrs" !~ /$chr/);

        if ($chrom ne $last_chrom) {
            my $fileOut = "$seq_cover_all_folder/$chrom.cov";
            if ($last_chrom ne "") {
                close(FHO);
                $log = "close $fileOut ($n_out lines)\n";
                LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
            }

            if (-f $fileOut) {
                $log = "$fileOut already exists\n";
                LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
                return;
            }

            $log = "open file $fileOut\n";
            LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

            open(FHO, ">$fileOut") || die("cannot create $fileOut");
            $last_chrom = $chrom;
        }
        print FHO "$pos $cov\n";
        $n_out++;
    }
    close(FHI);

}

sub recup_cover_from_gatk {
    my $sample_data = shift;
    my $seq_cover_all_folder = shift;
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    
    my $dna = $sample_data->{dna};
    #$dna =~ s/_/-/g;

    my $cov_file = $sample_data->{disk}."/".$sample_data->{run}."/alignement/gatk/Sample_".$dna."/".$dna.".cov";
    my @chrs = (1..22, 'X', 'Y');
    my $n_out = 0;
    my ($log, $last_chrom);

    $log = scalar(localtime())." open file $cov_file\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $dna);

    if (-f "$cov_file.gz"){
	system("gunzip $cov_file.gz");
    }

    open(FHI, $cov_file) or die("cannot open file $cov_file :$!\n");
    while (my $line = <FHI>) {
        $line =~ s/\012//;
        $line =~ s/\015//;
        chomp $line;

        next if ($line =~ /^#/);
        my ($chrom, $pos, $cov) = split(/\t/, $line);
        my $chr = substr($chrom, 3);
        next if ("@chrs" !~ /$chr/);

        if ($chrom ne $last_chrom) {
            my $fileOut = "$seq_cover_all_folder/$chrom.cov";
            if ($last_chrom ne "") {
                close(FHO);
                $log = "close $fileOut ($n_out lines)\n";
                LogOut::printLog($log, $current_sub, $sample_data->{run}, $dna);
            }

            if (-f $fileOut) {
                $log = "$fileOut already exists\n";
                LogOut::printLog($log, $current_sub, $sample_data->{run}, $dna);
                return;
            }

            $log = "open file $fileOut\n";
            LogOut::printLog($log, $current_sub, $sample_data->{run},$dna);

            open(FHO, ">$fileOut") || die("cannot create $fileOut");
            $last_chrom = $chrom;
        }
        print FHO "$pos $cov\n";
        $n_out++;
    }
    close(FHI);
    
    if (-f "$cov_file.gz"){
	system("rm $cov_file.gz");
    }
    system("gzip $cov_file");
}



1;
