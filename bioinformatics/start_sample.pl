#! /usr/bin/perl -I Include

################################
# Informations to know :
# - Al packages are in include folder.
# - All variables defined with .cfg extention
#   in Config folder are automatiquely loaded
#
#   @ The bioinfo team
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
use Rnaseq;
use Cover::RecupCover;
use Cover::FilterCover;
use SNPs::RecupSNP;
use SNPs::FilterSNP;
use SNPs::AnnotSNP;
use INDELs::RecupINDEL;
use INDELs::FilterINDEL;
use INDELs::AnnotINDEL;
use Reporting::LaunchReport;
use Reporting::Stats;
use Reporting::FastQC;
use Data::Dumper;
use Spreadsheet::WriteExcel;


### Chargement des fichiers de configurations du projet
%config = Projet::load_config();

### Verification de la variable $ARGV[0]
if ($ARGV[0] eq ""){
    print "Ne pas oublier le fichier de configuraiton du run ! \n";
    exit;
}

### Chargement du fichier de configurations du run
$run_config = "$config{racine}/$ARGV[0]";
%run = Projet::load_run($run_config);

### Verification des variables indispensable au run
$erreur = Projet::check_config(%run);
if ($erreur ne ""){
    print "Out\n$erreur";
    exit;
}

### sample a lancer
$adn = $ARGV[1];
push(@adn, $adn);

$runlog = $run{name};
$runlog =~ s/\-/\_/g;
$adnlog = $adn;
$adnlog =~ s/\-/\_/g;

($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=start_sample,status=start,run=$runlog,sample=$adnlog";
system ("echo $add >> $run{logfile}");
$timesample1=time();

$log_file = $config{racine}."/log/".$run{name}."/".$adn.".log";
open (LOG,">>$log_file");
print LOG "analysis sample $adn\n";
### Alignement bwa 0.7

if ($run{align_new_bwa} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=bwa,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch mapping with bwa0.7...\n";
    Alignement::align_new_bwa($config{bwa0715},$run{disque},$run{name},$adn,$run{genome_bwa},$config{samtools});
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=bwa,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=bwa,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
### se connecter a BDD runs et recuperer toutes les infos sur le run et les echantillons
$run_data = DBconnect::recup_all_data_run($run{name}, \%config, \%run);

### Realignment GATK from bwa
if ($run{gatk_realign} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_realign,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch realignment with GATK...\n";
    Alignement::gatk_realign($run{disque} ,$run{name}, \@adn, \%config, \$run_data, $run{genome_gatk}, "bwa");
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_realign,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_realign,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}


### mutations from GATK
if ($run{gatk_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_calling,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch calling with GATK...\n";
    Detection::gatk_mutation($run{disque} ,$run{name}, \@adn, \%config, \$run_data, $run{genome_gatk});
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_calling,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=gatk_calling,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

##### dossier results
$result_folder = $config{results}; ### quand créé une fois, plus besoin de refaire
mkdir($result_folder) unless (-d $result_folder);

### creer dossier log du run
$log_folder = $config{racine}."/log/".$run{name};
mkdir($log_folder) unless (-d $log_folder);

##### dossier data
$data_folder = $config{data}; ### contient les donnees expression, captures, etc...

$run_data = DBconnect::recup_all_data_run($run{name}, \%config, \%run);

##### Recuperation des couvertures
if ($run{recup_cover} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=recup_cover,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch recup_cover...\n";
    RecupCover::recup_all_cover($$run_data{$adn}, $result_folder, $run{nas});
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=recup_cover,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=recup_cover,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
##### Filtre des couvertures
if ($run{filter_cover} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=filter_cover,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch filter_cover...\n";
    FilterCover::filter_cover_in_capture($$run_data{$adn}, $result_folder, \%config);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=filter_cover,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=filter_cover,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

if ($run{alternative_filter_cover} eq "1"){
    print LOG "##### ".localtime()." launch filter_cover_in_alternative_bed...\n";
    FilterCover::filter_cover_in_alternative_bed($$run_data{$adn}, $result_folder, $run{bed}, $run{name}, \%config);}

##### Filtre et stats de couvertures dans un bed
if ($run{bed_stat_cover} eq "1"){
    print LOG "##### launch bed_stat_cover...\n";
    FilterCover::filter_cover_in_bed($$run_data{$adn}, $result_folder, $run{bed});
}

##### Recuperation SNPs/INDELs from GATK
if ($run{recup_gatk} eq "1"){
    print LOG "##### ".localtime()." launch recup_GATK SNP and INDEL...\n";
    RecupSNP::recup_gatk($run{disque} ,$run{name}, $adn, $result_folder,$$run_data{$adn});
    RecupINDEL::recup_gatk($run{disque} ,$run{name}, $adn, $result_folder,$$run_data{$adn});
}

##### Merge F_score (casava/GATK)
if ($run{recup_gatk} eq "1"){
    FilterSNP::cp_fscore($run{disque} ,$run{name},$adn, $result_folder, $$run_data{$adn});
    FilterINDEL::cp_fscore($run{disque} ,$run{name},$adn, $result_folder, $$run_data{$adn});
}

##### Annotation des SNPs et INDELs par Ensembl et stockage dans table mysql; dernier argument = serveur local (1) ou pas (0)
if ($run{ensembl_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch ensembl_annot snps...\n";
    AnnotSNP::ensembl_annot_snp($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{ensembl_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch ensembl_annot indels...\n";
    AnnotINDEL::ensembl_annot_indel($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=ensembl_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

##### Annot et Filtre (dernier arg = 1) ou pas (dernier arg = 0) des SNPs et INDELs sur base des consequences
if ($run{coding_filter_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch coding_filter snps...\n";
    if($run{wgs} eq "1"){
	FilterSNP::coding_filter_snp($$run_data{$adn}, $config{ensembl_version}, $result_folder, 0);
    } else {
	FilterSNP::coding_filter_snp($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    }
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{coding_filter_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch coding_filter indels...\n";
    if($run{wgs} eq "1"){
	FilterINDEL::coding_filter_indel($$run_data{$adn}, $config{ensembl_version}, $result_folder, 0);
    } else {
	FilterINDEL::coding_filter_indel($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    }
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=coding_filter_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

if ($run{consequence_filter_all_transcripts_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch consequence_filter_all_transcripts snps...\n";
    FilterSNP::consequence_filter_all_transcripts_snp($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{consequence_filter_all_transcripts_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch consequence_filter_all_transcripts indels...\n";
    FilterINDEL::consequence_filter_all_transcripts_indel($$run_data{$adn}, $config{ensembl_version}, $result_folder, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=consequence_filter_all_transcripts_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

##### Annotations supplementaires des SNPs et INDELs (last 5 arguments: from folder, family, filter, combinaison, recess)
if ($run{hg19_pos_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch hg19_pos_annot snps...\n";
    AnnotSNP::launch_hg19_pos_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_coding", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{hg19_pos_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch hg19_pos_annot indels...\n";
    AnnotINDEL::launch_hg19_pos_annot_indel($run{name}, $$run_data{$adn}, $result_folder, "F_coding", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=hg19_pos_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{dbsnp_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch dbSNP_annot snps...\n";
    AnnotSNP::launch_dbSNP_annot_snp($run{name}, $$run_data{$adn}, $config{dbsnp}, $result_folder, "F_hg19_pos", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{dbsnp_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch dbSNP_annot indels...\n";
    AnnotINDEL::launch_dbSNP_annot_indel($run{name}, $$run_data{$adn}, $config{dbsnp}, $result_folder, "F_hg19_pos", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbsnp_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{CG_54genomes_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch 54genomes_annot snps...\n";
    AnnotSNP::launch_CG_54genomes_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_dbSNP", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{CG_54genomes_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch 54genomes_annot indels...\n";
    AnnotINDEL::launch_CG_54genomes_annot_indel($run{name}, $$run_data{$adn}, $result_folder, "F_dbSNP", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=CG_54genomes_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{expr_beta_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch expr_beta_annot snps...\n";
    AnnotSNP::launch_expr_beta_annot_snp($run{name}, $$run_data{$adn}, $result_folder, $data_folder, "F_54genomes", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{expr_beta_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch expr_beta_annot indels...\n";
    AnnotINDEL::launch_expr_beta_annot_indel($run{name}, $$run_data{$adn}, $result_folder, $data_folder, "F_54genomes", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=expr_beta_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{good_annot_snp} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch good_annot snps...\n";
    AnnotSNP::launch_good_annot_snp($run{name}, $$run_data{$adn}, $result_folder, $data_folder, "F_expr", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if ($run{good_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch good_annot indels...\n";
    AnnotINDEL::launch_good_annot_indel($run{name}, $$run_data{$adn}, $result_folder, $data_folder, "F_expr", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=good_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
##### Annotation supplementaire des INDELs (last 6 arguments: from folder, family, filter, combinaison, recess, local server)
if ($run{maf_annot_indel} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=maf_annot_indel,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch maf_annot indels...\n";
    AnnotINDEL::launch_maf_annot_indel($run{name}, $$run_data{$adn}, $config{ensembl_version}, $result_folder, "F_good", 0, 0, 0, 0, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=maf_annot_indel,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=maf_annot_indel,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
# Inhouse analysis (dbNSFP from F_good, mut_indiv lists)
if (($run{dbnsfp_annot_snp} eq "1") && ($run{good_annot_snp} eq "1")){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch dbnsfp_annot...\n";
    AnnotSNP::launch_dbnsfp_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_good", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
### extra last argument = count des indiv (1) ou liste des indiv (0)
if (($run{indiv_annot_snp} eq "1") && ($run{good_annot_snp} eq "1")){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
     print LOG "##### ".localtime()." launch indiv_annot...\n";
    AnnotSNP::launch_indiv_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_dbNSFP", 0, 0, 0, 0, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
# Collab or presta (dbNSFP from 54genomes, mut_indiv counts)
if (($run{dbnsfp_annot_snp} eq "1") && ($run{good_annot_snp} eq "0")){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch dbnsfp_annot...\n";
    AnnotSNP::launch_dbnsfp_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_54genomes", 0, 0, 0, 0);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=dbnsfp_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}
if (($run{indiv_annot_snp} eq "1") && ($run{good_annot_snp} eq "0")){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch indiv_annot...\n";
    AnnotSNP::launch_indiv_annot_snp($run{name}, $$run_data{$adn}, $result_folder, "F_dbNSFP", 0, 0, 0, 0, 1);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=indiv_annot_snp,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

##### Annotation genotype
if ($run{base_count_snp} eq "1"){
    print LOG "##### ".localtime()." launch basecount_SNP...\n";
    AnnotSNP::launch_base_count_snp($run{name}, $run{disque}, $$run_data{$adn}, $result_folder, "F_mut_indiv");
}
if ($run{base_count_indel} eq "1"){
    print LOG "##### ".localtime()." launch basecount INDEL...\n";
    if ($run{good_annot_indel} eq "1"){
	AnnotINDEL::launch_base_count_indel($run{name}, $run{disque}, $$run_data{$adn}, $result_folder, "F_maf");
    } elsif ($run{good_annot_indel} eq "0"){
	AnnotINDEL::launch_base_count_indel($run{name}, $run{disque}, $$run_data{$adn}, $result_folder, "F_maf");
    }
}

##### Annotation disease et inheritance
if ($run{disease_inher_snp} eq "1"){
    print LOG "##### ".localtime()." launch disease inheritance SNP...\n";
    AnnotSNP::launch_disease_inher_snp($run{name}, $run{disque}, $$run_data{$adn}, $result_folder, "F_basecount");
}
if ($run{disease_inher_indel} eq "1"){
    print LOG "##### ".localtime()." launch disease inheritance INDEL...\n";
    AnnotINDEL::launch_disease_inher_indel($run{name}, $run{disque}, $$run_data{$adn}, $result_folder, "F_basecount");
}

### Insertion des nouvelles mutations dans notre DB
if($run{mut_indiv_insert} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=mut_indiv_insert,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." launch mut_indiv_insert...\n";
    DBconnect::mut_indiv_insert($$run_data{$adn}, $result_folder, "F_score");
    DBconnect::mut_indiv_insert_redis($$run_data{$adn}, $result_folder, "F_score");
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=mut_indiv_insert,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=mut_indiv_insert,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

#### RNAseq
if ($run{rna} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=rnaseq,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." Launch RNAseq analysis...\n";
    Rnaseq::launch_rnaseq($run{disque},$run{name},\%run,\%config,\@adn);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=rnaseq,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=rnaseq,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}

#### Methylation
$run_data = DBconnect::recup_all_data_run($run{name}, \%config, \%run);
$sample_data = $run_data->{$adn};
if ($run{methylation} eq "1"){
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=methylseq,status=start,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    $time1=time();
    print LOG "##### ".localtime()." Launch Methylation analysis...\n";
    Methylation::launch_methylation($run{disque},$run{name},\%run,\%config,\@adn, $sample_data);
    $time2=time();
    $duration=$time2-$time1;
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=methylseq,duration=$duration,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
    ($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
    $add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=methylseq,status=stop,run=$runlog,sample=$adnlog";
    system ("echo $add >> $run{logfile}");
}


### Launch stats and reporting

# Lecture des samples
@adn = Detection::read_sample($run{disque},$run{name});
$run_data = DBconnect::recup_all_data_run($run{name}, \%config, \%run);



if ($run{report} eq "1"){
    if ($run{rna} ne "1"){
	print LOG "##### ".localtime()." launch Reporting...\n";
	$report_folder = $run{disque}."/".$run{name}."/Reporting";
	unless (-d $report_folder){system("mkdir $report_folder")};
	$sample_report_folder = $report_folder."/Sample_".$adn;
	unless (-d $sample_report_folder){system("mkdir $sample_report_folder")};
	### lancement FastQC
	print LOG "##### ".localtime()." launch FastQC...\n";
	FastQC::launch_FastQC($adn, $run{disque}, $run{name}, $sample_report_folder, $config{fastqc}, $report_folder, $run{nas});
	### lancement calcul ontarget
	if ($run{methylation} ne "1"){
	    print LOG "##### ".localtime()." launch on target for ".$adn."\n";
	    Stats::reads_stats($adn, $run_data, $report_folder, $run{name}, $run{disque}, $run{nas}, $config{samtools});
	}

	#chrX, chrY and aut meanX calculation
	print LOG "##### ".localtime()." launch meanX_Y for ".$adn."\n";
	meanX_Y($adn, $run_data, \%run, \%config);

	print LOG "##### ".localtime()." END of stats for ".$adn."\n";
    }
}

$timesample2=time();
$duration=$timesample2-$timesample1;
($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=start_sample,duration=$duration,run=$runlog,sample=$adnlog";
system ("echo $add >> $run{logfile}");
($day,$month,$year,$hour,$min,$sec,$host) = loginfo();
$add = "$day-".$month."-".$year." ".$hour.":".$min.":".$sec." [".$host."] "."app=start_sample,status=stop,run=$runlog,sample=$adnlog";
system ("echo $add >> $run{logfile}");

### signal fin d'analyse...
$finished_folder = "$run{disque}/$run{name}/log/finished";
unless (-d "$finished_folder") { system("mkdir $finished_folder"); }
system("touch $finished_folder/$adn.ok");

print LOG "##### ".localtime()." ... end\n";
close LOG;

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

sub meanX_Y {
	my $adn = shift;
	my $run_data = shift;
	my $run = shift;
	my $config = shift;

	my ($bedtype, @res2);
	if($$run{alternative_filter_cover} eq "0"){
	    $bedtype = "InCapture";
	} elsif($$run{alternative_filter_cover} eq "1"){
	    $bedtype = "InBed";
	} else {
	    $bedtype = "InCapture";
	}
	my $capture = $$run_data{$adn}->{capture};
	my $folder = $$run_data{$adn}->{result_dir};
	open (COVER,">$$config{results}/SEQ_cover/$bedtype/$folder/CoverX_Y_old.txt");
	my $in_file = "$$run{disque}/$$run{name}/alignement/gatk/Sample_$adn/$adn\_RG.dedup.sorted.realigned.recal.bam";
  #print "$$config{samtools} depth -q30 -Q37 -a -b $$run{bed} $in_file | awk -f $$config{sexDet}\n";
	my @result = `$$config{samtools} depth -q23 -Q37 -a -b $$config{$capture} $in_file | awk -f $$config{sexDet}`;
### Ajout Mehdi 19/01/19 :
  my $result2 = `$$config{samtools} view -b $in_file chrY:2786855-2787699 | $$config{bedtools}/bin/bedtools coverage -abam /media/Data/sry.bed -b stdin -counts | cut -f4`;
  chomp $result2;
  print $result2;
  open (COVER2,">$$config{results}/SEQ_cover/$bedtype/$folder/CoverX_Y.txt");
###
  foreach my $i (@result){
	    #print "result : ".$i."\n";
	    chomp $i;
	    @res2 = split("\t",$i);
	    #foreach my $j (@res2){
		#print "second split : ".$j."\n";
		print COVER "".$res2[0]."\t".$res2[1]."\n";
	    #}
	}
  close COVER;
  open (COVER3,"$$config{results}/SEQ_cover/$bedtype/$folder/CoverX_Y_old.txt");
  while (my $line = <COVER3>) {
    chomp $line;
    if ( $line =~ "yCoverage" ) {    # changer "=~ " en "eq" si sa marche pas
        $line = "yCoverage\t"."$result2";
      print COVER2 "$line"."\n";
    }
    else {
      print COVER2 "$line"."\n";
    }

    }

	close COVER2;
  close COVER3;
}
