package Demultiplexage;

use Detection;
use Data::Dumper;

sub demultiplex_samplesheet {
	my ($disque,$name,$samplesheet,$demultiplexeur) = @_;
	system("mkdir $disque/$name/demultiplexage");
	system("mkdir $disque/$name/log");
	system("mv $samplesheet $disque/$name/samplesheet.csv");
	system("$demultiplexeur --input-dir $disque/$name/raw/Data/Intensities/BaseCalls/ --output-dir $disque/$name/demultiplexage --sample-sheet $disque/$name/samplesheet.csv --force 1>$disque/$name/log/demultiplexage1.txt 2>$disque/$name/log/demultiplexage1.txt");
	system("cd $disque/$name/demultiplexage/;make -j80 all 1>$disque/$name/log/demultiplexage2.txt 2>$disque/$name/log/demultiplexage2.txt;cd -");
	#system("chmod -R 775 $disque/$name/demultiplexage");
}

sub demultiplex_targeted {
        my ($disque,$name,$samplesheet,$demultiplexeur) = @_;
        system("mkdir $disque/$name/demultiplexage");
        system("mkdir $disque/$name/log");
        system("mv $samplesheet $disque/$name/samplesheet.csv");
        system("$demultiplexeur --input-dir $disque/$name/raw/Data/Intensities/BaseCalls/ --output-dir $disque/$name/demultiplexage --sample-sheet $disque/$name/samplesheet.csv --use-bases-mask Y*n,I*,Y*n --force 1>$disque/$name/log/demultiplexage1.txt 2>$disque/$name/log/demultiplexage1.txt");  ### modification --use-bases-mask
        system("cd $disque/$name/demultiplexage/;make -j80 all 1>$disque/$name/log/demultiplexage2.txt 2>$disque/$name/log/demultiplexage2.txt;cd -");
        #system("chmod -R 775 $disque/$name/demultiplexage");
}

sub demultiplex_new_bcl {
        my ($disque,$name,$samplesheet,$bcl2fastq2) = @_;
        system("mkdir $disque/$name/demultiplexage");
        system("mkdir $disque/$name/log");
	if (-e $samplesheet){
        system("mv $samplesheet $disque/$name/samplesheet.csv");}
	my $mask = check_mask($disque,$name);
	
	system("$bcl2fastq2 --input-dir $disque/$name/raw/Data/Intensities/BaseCalls/ --ignore-missing-bcls --runfolder-dir $disque/$name/raw --output-dir $disque/$name/demultiplexage $mask 1>$disque/$name/log/demultiplexage.txt 2>$disque/$name/log/demultiplexage.txt\n");
}

sub group_fastq {
	my ($disque,$name,$adn) = @_;
	unless ( -d "$disque/$name/demultiplexage/fastq" ) { system("mkdir $disque/$name/demultiplexage/fastq"); }
	system("cat $disque/$name/demultiplexage/".$adn."_*_R1_*.fastq.gz > $disque/$name/demultiplexage/fastq/".$adn."_R1.fastq.gz");
	system("cat $disque/$name/demultiplexage/".$adn."_*_R2_*.fastq.gz > $disque/$name/demultiplexage/fastq/".$adn."_R2.fastq.gz");
}

sub check_samplesheet {
    my ($disque,$name,$sequencer) = @_;
    my $run_samplesheet = "$disque/$name/raw/SampleSheet.csv";
    # vérifier si samplesheet dans raw existe
    if (!(-f $run_samplesheet)){
	die "Cannot launch demultiplexing : SampleSheet in raw files doesn't exist.\n";
    }
    
    # vérifier si nom samples idems entre samplesheet formulaire et samplesheet run
    my @samples_form = Detection::read_sample($disque,$name);
    my $samples_run = read_samplesheet($run_samplesheet,$sequencer);

    my ($item, $item2);
    my @isect = ();
    my @diff  = ();
    my %count = ();
    my @not_form = ();
    my @not_run = ();
    foreach $item (@samples_form){
	if (!($item ~~ @$samples_run)){ 
	    push @not_run, $item;
	}
    }
    foreach $item (@$samples_run){
	if (!($item ~~ @samples_form)){
	    push @not_form, $item;
	}
    }
    
    if (@not_form){
	print "Pas dans le formulaire...\n";
	print Dumper(\@not_form);
    }
    if (@not_run){
	print "Pas dans la samplesheet du run...\n";
	print Dumper(\@not_run);
    }

    if (@not_form || @not_run){
	die "Sample names in Database and in Demultiplexing Samplesheet are not the same...\n";
    } else {
	print "Roger, can proceed!! \n";
    }

}

sub read_samplesheet {
    my $samplesheet = shift;
    my $sequencer = shift;
    
    my (@line, @samples, $sample,@adn);
    open (IN, "$samplesheet") || die "Cannot open Samplesheet$!\n";
    my $i = 0;
    LINE : while (my $l =<IN>){
	chomp $l;
	if ($l=~/^Sample/){
	    $i = 1;
	    next;
	}
	if ($i == 1){
	    @line = split(",", $l);
	    $sample = $line[1];
	    push(@samples, $sample);
	} 
	if ($l=~/^Lane/){
	    $i = 2;
	    next;
	}
	if ($i == 2){
	    @line = split(",", $l);
	    $sample = $line[2];
	    push(@samples, $sample);
	}
    }
    close IN;
    
    foreach my $objet (@samples){
	push(@adn, $objet) unless $deja_vu{$objet}++;
    }
    return (\@adn);
}

sub check_mask {
    my ($disque,$name) = @_;
    my $run_info = "$disque/$name/raw/RunInfo.xml";
    my $samplesheet = "$disque/$name/samplesheet.csv";
    my $nb = 0;
    my $nb_index;
    my $len_index;
    my $mask;
    my $isindex;
    
    open (IN, "$run_info") || die "Cannot open RunInfo file$!\n";
    while (my $l = <IN>){
	if ($l=~ /Read Number/){
	    $nb++;
	    if ($l=~ /Read Number="2" NumCycles="([\d]+)"/){
		$nb_index = $1;
	    }
	    if ($l=~ /Read Number="3" NumCycles="([\d]+)" IsIndexedRead="([\w])"/){
                $isindex = $2;
            }
	}
    }
    close IN;

    my (@line, $index);
    open (IN2, "$samplesheet") || die "Cannot open samplesheet$!\n";
    while (my $l = <IN2>){
	next if ($l=~/^FCID/);
	chomp $l;
	@line = split(/\,/, $l);
	$index = $line[4];
	$len_index = length($index);
    }
    
    close IN2;
    
    if ($nb_index eq $len_index){ # pas de base-mask
	$mask = "";
    } elsif ($nb_index gt $len_index){
	if ($nb eq "2"){ # single end
	    $mask = "--use-bases-mask Y*,I6N*";
	} elsif ($nb eq "3"){ # paired end
	    if ($isindex eq "Y"){
		$mask = "--use-bases-mask Y*,I*,I*";
	    } elsif ($isindex eq "N"){
		$mask = "--use-bases-mask Y*,I6N*,Y*";
	    }
	} elsif ($nb eq "4"){ # paired end and dual index
	    $mask = "--use-bases-mask Y*,I*,I*,Y*";
	}
    }
	return $mask;
}


1;
