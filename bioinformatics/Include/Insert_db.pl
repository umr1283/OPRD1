#!/usr/bin/perl
use warnings;
use strict;
use DBI;
use Getopt::Long;
use DBconnect;

main();

sub main {	
    my ($name_run, $casava, $dir, @ali_arg, $ali);

    GetOptions( 
	'n=s' => \$name_run,     
	'c=s'  => \$casava, 
	'r=s'   => \$dir,  
	'a=s' => \@ali_arg,		
	);

    if(!(defined $name_run && defined $casava && defined $dir && @ali_arg)) {
	die "Les arguments doivent etre de la forme : perl add_in_bd_after_ali.pl -n run_name -c casava_version -r ali_dir -a hg38";
    }

    my $dbh = DBconnect::connect_db("RUNS");

    #Recuperation des noms des colonnes de la table alignment
    my $alignment = list_ali_column($dbh);

    #recup des ali donnes en argument
    $ali = join(' ', @ali_arg);		
    my $exist = 0;

    #Verification que ces ref existent
    for my $a_argv (@ali_arg) {
	for $a (@$alignment) {
	    if($a_argv eq $a) {
		$exist = 1;	
	    }
	}		
	if($exist == 0) {
	    die "Votre alignement de reference n'existe pas";	
	}
	$exist = 0;
    }

    add_ali($dbh, $name_run, @ali_arg);
    add_casava_and_dir($dbh, $name_run, $casava, $dir);

    # Deconnexion de la base de donnees
    $dbh->disconnect();
}

sub list_ali_column {
    my ($dbh) = @_;

    my $prep_ali = $dbh->prepare("SHOW COLUMNS FROM ALIGNMENT") or die "$dbh -> errstr";
    $prep_ali->execute() or die "Echec de la requete \n";

    my $i = 0;
    my @ali;

    #Suppression des deux premiers colonnes (id_ali et run_name)
    while (my @alignments = $prep_ali->fetchrow_array) {
	if( $i > 1) {
	    push (@ali, $alignments[0]);	
	}
	$i ++;
    }	
    
    $prep_ali->finish();
    return \@ali;
}

#Ajout des alignements
sub add_ali() {
    my ($dbh, $name_run, @alignments) = @_;	
    my $requete_sql_ali;
    
    my $prep = $dbh->prepare('SELECT name_run FROM ALIGNMENT WHERE name_run = ?') or die "$dbh -> errstr";			
    $prep->execute($name_run) or die "Echec de la requete \n";
    my @values = $prep -> fetchrow_array;
    $prep -> finish(); 

    if(scalar(@values) != 0) { 
	my (@sql, $sql);
	if($#alignments == 0) {
	    $sql[0] = $alignments[0].'=1';			
	} else {
	    for (my $i = 0; $i < $#alignments; $i ++) {		
		$sql[$i] = "$alignments[$i]=1";				
	    }			
	}
	$sql = join(",", @sql);		
	$requete_sql_ali = "UPDATE ALIGNMENT set $sql WHERE name_run='".$name_run."'";
    } else {
	my $sql = '';
	my $res = '';		
	#creation de la requete sql
	if($#alignments == 0) {
	    $sql = $alignments[0];
	    $res = '1';
	} else {
	    for (my $i = 0; $i < $#alignments; $i ++) {		
		$sql = $sql . " $alignments[$i],"; 
		$res = $res . '1,';
	    }
	    $sql = $sql . " $alignments[$#alignments]"; 
	    $res = $res . '1';
	}
	$requete_sql_ali = "INSERT INTO ALIGNMENT(name_run, $sql) VALUES ('".$name_run."', $res)";
    }	
    
    my $sth_ali = $dbh -> prepare($requete_sql_ali) or die $dbh->errstr;

    $sth_ali->execute() or die "Echec Requete $requete_sql_ali : $DBI::errstr";	
    $sth_ali->finish();
    print "$requete_sql_ali a ete executee correctement.\n";
}
    
    #Ici update de la table RUN
    my $requete_sql_run = "UPDATE RUN SET casava_version = '".$casava."', ali_dir = '".$dir."' WHERE name_run = '".$name_run."'";
    my $sth_run = $dbh->prepare($requete_sql_run) or die $dbh->errstr;
    
    $sth_run->execute() or die "Echec Requete $requete_sql_run : $DBI::errstr";
    $sth_run->finish();
    print "$requete_sql_run a ete executee correctement. \n";	
}
