package DBconnect;

use DBI;
use Data::Dumper;
use Redis;


## récupération des infos d'un run dans la DB runs
sub recup_all_data_run {
    my $name_run = shift;
    my $config = shift;
    my $run = shift;
    my $dbh = connect_db('RUNS');

print "lancement de recup_all_data_run : \n\n";

    ## recup alignements faits sur ce run
    my $version_ali = recup_ali($name_run, $dbh);

    ## recup run infos
    my @run_infos = recup_info_run($name_run, $dbh);

    ## recup_samples et créer objet du run
    my %all_samples = recup_info_sample_run($name_run, $dbh, @run_infos, $config, $run);

    # Deconnexion de la base de donnees
    $dbh->disconnect();

    return \%all_samples;
    #return \$run_obj;

}

### récupère la version du genome utilisée pour l'alignement des samples
sub recup_ali {
    my ($name_run, $dbh) = @_;
    my $prep_ali = $dbh->prepare("SHOW COLUMNS FROM ALIGNMENT") or die "$dbh -> errstr";
    $prep_ali->execute() or die "Echec de la requete \n";
    my $i = 0;
    my @ali;

    #Suppression des deux premieres colonnes (id_ali et run_name)
    while (my @alignments = $prep_ali->fetchrow_array) {
	if( $i > 1) {
	    push (@ali, $alignments[0]);
	}
	$i ++;
    }

    $prep_ali->finish();

    my $alis = join(',',@ali);
    my $prep_ali = $dbh->prepare("SELECT $alis FROM ALIGNMENT WHERE name_run = ?;") or die "$dbh -> errstr";
    $prep_ali->execute($name_run) or die "Echec de la requete \n";

    my $alignments = $prep_ali->fetchrow_hashref;

    $prep_ali->finish();

    return \$alignments;
}

### récupère les infos sur le run (dossiers d'alignement, date du run, nas, ...)
sub recup_info_run {
    my ($name_run, $dbh) = @_;
    my $prep_run = $dbh->prepare('SELECT id_run, name_run as run, read_type as method, sequencer, date_run_start, date_run_end, nas_file_name as NAS, casava_version, ali_dir FROM RUN WHERE name_run = ?;') or die "$dbh -> errstr";
    $prep_run->execute($name_run) or die "Echec de la requete \n";
    my @run = $prep_run->fetchrow_array;
    #print Dumper(\@run);
    $prep_run->finish();
    return \@run;
}

sub recup_info_sample_run {
    my ($name_run, $dbh, $run_infos, $config, $run) = @_;
    my %all_samples;
    my $prep_samples = $dbh->prepare('SELECT family as pop, name_sample as dna, sex, line_sample as lane, microarray as MA, prestation, patient_control as atteint, capture, illness as maladie, organism, control as controle FROM SAMPLE WHERE id_run = ?;') or die "$dbh -> errstr";
    $prep_samples->execute($$run_infos[0]) or die "Echec de la requete \n";

    ## ajouter chaque sample en tant qu'objet dans le hash %all_samples avec clé valeur "nom sample" => "object sample"
    while ( my @sample = $prep_samples->fetchrow_array ) {
	$sample_obj = create_obj_run(\@$run_infos, \@sample, $config, $run);
	## si sample existe dejà, c'est qu'il est sur plusieurs lignes => incrementer nombre de lignes
	if (exists($all_samples{$sample[1]})){
	    $all_samples{$sample[1]}->{'lane'}++;
	} else {
	    $all_samples{$sample[1]} = $sample_obj;
	}
     }
    return %all_samples;

}

## création d'un objet run contenant toutes les infos nécessaires pour le run
sub create_obj_run {
    my ($run_infos, $sample, $config, $run) = @_;

    if ($$sample[6] == 0) {
	$$sample[6] = 1;
    } elsif ($$sample[6] == 1) {
	$$sample[6] = 0;
    }

    # renvoie l'objet du run
    my $this = {
	'pop' => $$sample[0],
	'dna' => $$sample[1],
	'result_dir' => $$sample[1]."_hg38", # recup ali
	'sex' => $$sample[2],
	'lane' => 1,
	'MA' => $$sample[4],
	'prestation' => $$sample[5],
	'atteint' => $$sample[6],
	'capture' => $$sample[7],
	'bed_file' => $$config{$$sample[7]},
	'maladie' => $$sample[8],
	'organism' => $$sample[9],
	'controle' => $$sample[10],
	'disk' => $$run{disque},
	'run' => $$run_infos[1],
	'method' => $$run_infos[2],
	'sequencer' => $$run_infos[3],
	'date_debut_run' => $$run_infos[4],
	'date_fin_run' => $$run_infos[5],
	'NAS' => $$run_infos[6],
	'casava_version' => $$run_infos[7],
	'orga_ref' => "hg38", ## recup ali
	'ali_dir' => $$run_infos[8], ## Gerald ou Gerald_hg38 ??
	'sample_dir' => "Sample_".$$sample[1],
    };

    bless ($this);
    return $this;

}

##################################################################################
## récupération des infos sur les samples d'une analyse dans la DB runs
sub recup_all_data_analysis {
    my $samples = shift;
    my $dbh = connect_db('RUNS');

    ## recup_samples et créer objet du run
    my %all_samples = recup_info_sample_analysis($dbh, $samples);

    # Deconnexion de la base de donnees
    $dbh->disconnect();

    return \%all_samples;

}

##################################################################################
sub recup_info_sample_analysis {
    my ($dbh, $samples) = @_;
    my %all_samples;
    my $prep_samples = $dbh->prepare('SELECT family as pop, name_sample as dna, sex, line_sample as lane, microarray as MA, prestation, patient_control as atteint, capture, illness as maladie, organism, control as controle FROM SAMPLE WHERE name_sample = ?;') or die "$dbh -> errstr";
    foreach my $k (keys %$samples) {
	$prep_samples->execute($k) or die "Echec de la requete \n";

	## ajouter chaque sample en tant qu'objet dans le hash %all_samples avec clé valeur "nom sample" => "object sample"
	while ( my @sample = $prep_samples->fetchrow_array ) {
	    $sample_obj = create_obj_analysis(\@sample);
	    ## si sample existe dejà, c'est qu'il est sur plusieurs lignes => incrementer nombre de lignes
	    if (exists($all_samples{$sample[1]})){
		$all_samples{$sample[1]}->{'lane'}++;
	    } else {
		$all_samples{$sample[1]} = $sample_obj;
	    }
	}
    }
    return %all_samples;

}

##################################################################################
## création d'un objet analyse contenant toutes les infos nécessaires pour une analyse
sub create_obj_analysis {
    my ($sample) = @_;

    if ($$sample[6] == 0) {
	$$sample[6] = 1;
    } elsif ($$sample[6] == 1) {
	$$sample[6] = 0;
    }

    # renvoie l'objet de l'analyse
    my $this = {
	'pop' => $$sample[0],
	'dna' => $$sample[1],
	'result_dir' => $$sample[1]."_hg38", # recup ali
	'sex' => $$sample[2],
	'lane' => 1,
	'MA' => $$sample[4],
	'prestation' => $$sample[5],
	'atteint' => $$sample[6],
	'capture' => $$sample[7],
	'maladie' => $$sample[8],
	'organism' => $$sample[9],
	'controle' => $$sample[10],
	'orga_ref' => "hg38", ## recup ali
	'sample_dir' => "Sample_".$$sample[1],
    };

    bless ($this);
    return $this;

}

##################################################################################
## Insert new mutations in DB (table mut_indiv)
sub mut_indiv_insert {
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $req, $sth, $n, @tab);

    my $fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

    my $dbh = DBconnect::connect_db("SNP_annot");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHI, $fileIn) or die("cannot open $fileIn");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	next if ($line =~ /^\#/);
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	@tab = split(/\t/, $line);

	# Check for mutation in mut_indiv table
	$req = "select * from mut_indiv_hg38 where " .
	    " Individu = \"$sample_data->{result_dir}\" " .
	    " and Chrom = \"$tab[0]\" " .
	    " and Pos = \"$tab[1]\" " .
	    " and ARA1 = \"$tab[2]/$tab[3]\" " .
	    " and ARA2 = \"$tab[2]/$tab[4]\" ";

	$sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec: ".$dbh->errstr."\n");
	$sth->execute()||  print ("\nExecution requete: $req echec : ".$sth->errstr."\n");

	$n = $sth->rows;

	if ($n == 0) {
	    $n_out++;
	    $req = "insert into mut_indiv_hg38 set " .
		" Individu = \"$sample_data->{result_dir}\", " .
		" Chrom = \"$tab[0]\", " .
		" Pos = \"$tab[1]\", " .
		" ARA1 = \"$tab[2]/$tab[3]\", " .
		" ARA2 = \"$tab[2]/$tab[4]\", " .
		" hh = \"$tab[5]\", " .
		" Score = \"$tab[6]\" ";

	    $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	    $sth->execute()||  print ("\nExecution requete: $req echec\n");
	}
	$sth->finish();
    }

    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);

    $dbh->disconnect();

    $log  = "$n_out line(s) added to table, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

}


##################################################################################
## Insert new mutations/indivs in DB Redis
sub mut_indiv_insert_redis {
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $n, @tab, $status);

    my $redis = Redis->new(
	server => 'XXXX:XXX',
	);

    my $ind_ref = $sample_data->{result_dir};
    my $ind = $ind_ref;
    $ind =~s/_hg38//;

    # Chargement du status dans redis
    my $dbh = connect_db("RUNS");
    my $req = "select name_sample,patient_control from SAMPLE where name_sample = \"$ind\"";
    my $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
    $sth->execute()||  print ("\nExecution requete: $req echec\n");
    while (my @ligne = $sth->fetchrow_array) {
	$hash_name  = 'Status_table';
	$redis->hset( $hash_name, $ligne[0] => $ligne[1]);
    }

    my $fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open (FHI, $fileIn) or die("cannot open $fileIn");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	next if ($line =~ /^\#/);
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	@tab = split(/\t/, $line);

	next if ($ind_ref !~ /hg38/);
	next if ($ind =~ /Sample/); # sauf exception
	next if ($ind =~ /_S/);
	next if ($ind =~ /test/);
	next if ($ind =~ /uniq/);
	next if ($ind =~ /Next/);

        $status = $redis->hget( "Status_table", $ind );
        if ($status eq "0"){
		       $list_name = "$tab[0]"."-"."$tab[1]"."-cas-hg38";
                       $redis->sadd( $list_name => "$ind_ref" );
		       $n_out++;
        }
        if ($status eq "1"){
		       $list_name = "$tab[0]"."-"."$tab[1]"."-control-hg38";
                       $redis->sadd( $list_name => "$ind_ref" );
		       	$n_out++;
        }
    }

    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);

    $log  = "$n_out line(s) added to table, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
}


##################################################################################
## connexion à une DB mysql du labo
sub connect_db {
    # Parametres de connexion a la base de donnees
    my ($db) = @_;
    my $host = 'XXX';
    my $login = 'XXX';
    my $pwd  = 'XXX';

    # Connexion a la base de donnees MySQL
    my $dbh = DBI->connect( "dbi:mysql:dbname=$db;host=$host;", $login, $pwd ) or die "Couldn't connect to database $db : ". DBI->errstr;
    return $dbh;
}


1;
