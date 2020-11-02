#! /usr/bin/perl -I Include

################################
# Informations to know :
# - Al packages are in include folder.
# - All variables defined with .cfg extention
#   in Config folder are automatiquely loaded
#
#	@ The bioinfo team
################################


### Loading of project packages
use Projet;
use Folder;
use Copie;
use Demultiplexage;
use Alignement;
use Detection;
use DBconnect;
use Link;
use Pirna;
use Mirna;
use Rnaseq;
use Cnv;
use Cover::RecupCover;
use Cover::FilterCover;
use SNPs::RecupSNP;
use SNPs::FilterSNP;
use SNPs::AnnotSNP;
use INDELs::RecupINDEL;
use INDELs::FilterINDEL;
use INDELs::AnnotINDEL;
use Reporting::LaunchReport;
use Data::Dumper;
use Spreadsheet::WriteExcel;

### Chargement des fichiers de configurations du projet
%config = Projet::load_config();

### Vérification de la variable $ARGV[0]
if ($ARGV[0] eq ""){ print "Ne pas oublier le fichier de configuration du run ! \n";
    exit;}

### Chargement du fichier de configurations du run
%run = Projet::load_run($ARGV[0]);

### Vérification des variables indispensable au run
$erreur = Projet::check_config(%run);
if ($erreur ne ""){ print "Out\n$erreur";
    exit;}

### Affiche un warning si le run existe
$warning = Folder::check_run("$run{disque}/$run{name}");
if ($warning ne ""){ print "$warning";}

### On récupère la liste des serveurs à utiliser
$serveurs = $run{server};
@serveur = split(",",$serveurs);

### Chargement du programme de contrôle du processus
$pid = $$;
system("nohup perl shaker.pl $pid $run{name} $run{email} $run{user} 1>shaker_log 2>shaker_log &");

### Logs
$runlog = $run{name};
$runlog =~ s/\-/\_/g;
($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=multiserveur,status=start,run=$runlog";
system ("echo $add >> $run{logfile}");
$timeall1=time();

### Création du dossier du run
Folder::create_run("$run{disque}/$run{name}");

### Copie du run
if ($run{copie} eq "1"){
	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=copie,status=start,run=$runlog";
	system ("echo $add >> $run{logfile}");

	$timecopie1=time();
	print "##### ".localtime()." : launch run copy...\n";
	Copie::run($run{nas},$run{from},$run{disque},$run{name});
	$timecopie2=time();
	$duration=$timecopie2-$timecopie1;

	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=copie,duration=$duration,run=$runlog";
	system ("echo $add >> $run{logfile}");

	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=copie,status=stop,run=$runlog";
	system ("echo $add >> $run{logfile}");
}

### Copie de la SampleSheet
if (-e $run{samplesheet}){
system("mv $run{samplesheet} $run{disque}/$run{name}/samplesheet.csv");}

### verifier Samplesheets
#Demultiplexage::check_samplesheet($run{disque},$run{name}, $run{nas});

### Demultiplexage avec le nouveau bcl2fastq
if ($run{demultiplex_new_bcl} eq "1"){
	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
       	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=demultiplexage,status=start,run=$runlog";
       	system ("echo $add >> $run{logfile}");
       	$time1=time();
	#if ($run{nas}=~/Novaseq/){
	    $bcl2fastq2 = $config{bcl2fastq219_1};
	#} else {
	#    $bcl2fastq2 = $config{bcl2fastq218};
	#}
	print "##### ".localtime()." : launch demultiplexing...\n";
	Demultiplexage::demultiplex_new_bcl($run{disque},$run{name},$run{samplesheet},$bcl2fastq2);

	$time2=time();
       	$duration=$time2-$time1;

       	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
        $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=demultiplexage,duration=$duration,run=$runlog";
       	system ("echo $add >> $run{logfile}");

        ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
       	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=demultiplexage,status=stop,run=$runlog";
       	system ("echo $add >> $run{logfile}");
}

### Lecture des samples

@adn = Detection::read_sample($run{disque},$run{name});

###
### Regroupement des différents fichiers fastq pour l'analyse
if ($run{group_fastq} eq "1"){
    print "##### ".localtime()." : launch group_fastq...\n";
    print @adn;
    for $adn(@adn){
    print "echantillon en cours : $adn\n";
	$adnlog = $adn;
	$adnlog =~ s/\-/\_/g;
	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
	$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=group_fastq,status=start,run=$runlog,sample=$adnlog";
        system ("echo $add >> $run{logfile}");
	$time1=time();
	Demultiplexage::group_fastq($run{disque},$run{name},$adn);
 print "##### ".localtime()." : finished sample $adn...\n";
	$time2=time();
        $duration=$time2-$time1;
	($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
        $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=group_fastq,duration=$duration,run=$runlog,sample=$adnlog";
        system ("echo $add >> $run{logfile}");
        ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
        $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=group_fastq,status=stop,run=$runlog,sample=$adnlog";
        system ("echo $add >> $run{logfile}");
    print "fin de l'echantillon en cours : $adn\n";
    }

 print "##### ".localtime()." : finished group_fastq...\n";
}

if ($run{demultiplex_only} eq "1"){
    system("mv $config{racine}/$ARGV[0] $config{racine}/finished/runs/$ARGV[0]");
    system("touch $config{racine}/$run{name}_finish");
    exit;
}

### On passe les échantillons déjà analysé
$log_dir = $run{disque}."\/".$run{name}."\/log/finished";
@skip_ok = listFile($log_dir);
for $adn (@adn){
	$vue = "0";
        for $file (@skip_ok){
		$match = $adn.".ok";
		$vue = "1" if ($file eq $match);
		print "Analyse de $adn déjà en cache \n" if ($file eq $match);
        }
	if ($vue eq "0"){
		push @adn2, $adn;
	}
}
@adn = @adn2;

### Analyse pour 1 ADN
$work_dir = $config{racine};
while(@adn){
	$adn = shift(@adn);
		#Lancement de l'analyse
		$date = localtime();
		print "nom de cluster  : $run{cluster} \n";
		if ($run{cluster} eq "mesos") {
		print "Analyse de $adn sur $serveur le $date - Il reste ".scalar(grep {defined $_} @adn)." sample à lancer \n";
		system("perl add_job.pl $run{name} $adn");
		}
		if ($run{cluster} eq "ocean") {
		print "Analyse de $adn sur $serveur le $date - Il reste ".scalar(grep {defined $_} @adn)." sample à lancer \n";
		system("perl ocean-add-job.pl $run{name} $adn");
		}
		if ($run{cluster} eq "slurm") {
		print "Analyse de $adn sur $serveur le $date - Il reste ".scalar(grep {defined $_} @adn)." sample à lancer \n";
		system("srun sh /media/Script/anges/ocean-job.sh $run{name} $adn");
		}
		sleep (3);
}

### Lecture des samples une nouvelle fois
@adn = Detection::read_sample($run{disque},$run{name});

### Verification que tous les ADNs sont analysé
$log_dir = $run{disque}."\/".$run{name}."\/log/finished";
while(@adn){
	@check_file = listFile($log_dir);
	$adn = shift(@adn);
	$ok = "0";
	for $file (@check_file){
		$match = $adn.".ok";
		$ok = "1" if ($file eq $match);
	}
	if ($ok eq "0"){
	unshift @adn, $adn;
	sleep (300);
	}
}

# update aligner and reference on database
if ($run{align_new_bwa} eq "1"){
	system("perl $config{racine}/Include/Insert_db.pl -n $run{name} -c 0.7 -r bwa -a hg38 12>$run{disque}/$run{name}/log/insert_db.txt 2>$run{disque}/$run{name}/log/insert_db.txt");}

### Lecture des samples une nouvelle fois
@list_adn = Detection::read_sample($run{disque},$run{name});

### se connecter à BDD runs et récupérer toutes les infos sur le run et les échantillons
$run_data = DBconnect::recup_all_data_run($run{name}, \%config, \%run);

#### CNV analysis
@adn = Detection::read_sample($run{disque},$run{name});
if ($run{cnv_analysis} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=cnvanalysis,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." Launch CNV analysis...\n";
    Cnv::launch_cnv($run{disque},$run{name},\@adn,\%config,\%run, $run_data);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=cnvanalysis,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=cnvanalysis,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

### Lancer calcul des stats, FastQC et rapport html/xls
# 10e paramètre = calculer stats à partir de incapture (0), de inbed (bed alternatif; 1) ou genome complet (2)
$result_folder = $config{results};
if($run{report} eq "1"){
    print "Launch reporting stats...\n";
    if($run{wgs} eq "1"){
	LaunchReport::launch_all_reports(\@list_adn, $run_data, $config{fastqc}, $config{samtools}, $result_folder, $run{disque}, $run{name}, $run{nas}, $run{from}, 2, $run{bed}, $run{pdf}, $config{sexDet});
    } elsif ($run{rna} eq "1"){
	Rnaseq::rna_stats(\%run, \@list_adn);
    } elsif ($run{methylation}){
	LaunchReport::launch_all_reports(\@list_adn, $run_data, $config{fastqc}, $config{samtools}, $result_folder, $run{disque}, $run{name}, $run{nas}, $run{from}, 0, $run{bed}, $run{pdf}, $config{sexDet});
    } else {
	if($run{alternative_filter_cover} eq "0"){
	    LaunchReport::launch_all_reports(\@list_adn, $run_data, $config{fastqc}, $config{samtools}, $result_folder, $run{disque}, $run{name}, $run{nas}, $run{from}, 0, $run{bed}, $run{pdf}, $config{sexDet});
	} elsif($run{alternative_filter_cover} eq "1"){
	    LaunchReport::launch_all_reports(\@list_adn, $run_data, $config{fastqc}, $config{samtools}, $result_folder, $run{disque}, $run{name}, $run{nas}, $run{from}, 1, $run{bed}, $run{pdf}, $config{sexDet});
	} else {
	    LaunchReport::launch_all_reports(\@list_adn, $run_data, $config{fastqc}, $config{samtools}, $result_folder, $run{disque}, $run{name}, $run{nas}, $run{from}, 0, $run{bed}, $run{pdf}, $config{sexDet});
	}
    }
}

if ($run{cover_xls_summary} eq "1"){FilterCover::exon_covers_xls_summary(\@list_adn, $result_folder, \%run);}

### Liens symboliques
@adn_list = Detection::read_sample($run{disque},$run{name});
if ($run{link_fastq} eq "1"){Link::link_fastq($run{name},$run{disque},$run{zip_fastq});}
if ($run{link_bam} eq "1"){Link::link_bam($run{name},$run{disque});}
if ($run{link_vcf} eq "1"){Link::link_vcf($run{name},$run{disque});}
if ($run{link_stats} eq "1"){Link::link_stats($run{name},$run{disque});}
if ($run{link_annot} eq "1"){Link::link_annot($run{name},\@adn_list,$config{racine});}
if ($run{link_cnv} eq "1"){Link::link_cnv($run{name},$run{disque});}

### zip fastq pour presta ou collab
if ($run{zip_fastq} eq "1"){Link::zip_fastq($run{name},$run{disque});}

### Logs
$timeall2=time();
$duration=$timeall2-$timeall1;

($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=multiserveur,duration=$duration,run=$runlog";
system ("echo $add >> $run{logfile}");

($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=multiserveur,status=stop,run=$runlog";
system ("echo $add >> $run{logfile}");

### Fichier de fin
if ($ARGV[0] !~ "run_test") {
system("mv $config{racine}/$ARGV[0] $config{racine}/finished/runs/$ARGV[0]");
}

system("touch $config{racine}/$run{name}_finish");

###########
# Routine #
###########

sub listFile {
    my($dir) = @_;
    my($fhdir) = 'FHDIR';
    my(@liFi);
    my($fic);

    opendir($fhdir, $dir);
    while ($fic = readdir($fhdir)) {
       next if ( $fic eq "." );
       next if ( $fic eq ".." );
       next unless -f "$dir/$fic";
       push @liFi, $fic;
    }
    closedir($fhdir);
    return @liFi;
}

sub loginfo {
	my ($sec,$min,$hour,$day,$month,$year)=(localtime)[0,1,2,3,4,5];
        my ($host);
        $month = sprintf("%02d",($month+1));
        $day = sprintf("%02d",($day));
        $year = sprintf("%04d",($year+1900));
        $host = `hostname`;
        $host =~ s/\R//g;

        return $day,$month,$year,$hour,$min,$sec,$host;
}
