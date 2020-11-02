package AnnotSNP;

##################################################################################
# From a tab delimited input file containing ID, chromosome, start, end and allele,
# this module retrieves from Ensembl Local ID, Annotated ID, Local alleles,
# Annotated alleles, Chromosome, Start, End, Consequence, Transcript, AA change,
# Gene ID, and Gene name and writes it to a file
##################################################################################

use strict;
use warnings;
use DBI;
use DBconnect;
use Data::Dumper;
use HTTP::Tiny;
use JSON;

use Bio::Perl;
use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);

use Redis;
use MongoDB;

##################################################################################
###
sub ensembl_annot_snp {
    my $sample_data = shift;
    my $dbversion = shift;
    my $result_folder = shift;
    my $local_server = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my $log;

    my $port;
    if (($dbversion < 48) || ($local_server == 1)) {
	$port = '3306';
    } else {
	$port = '5306';
    }


    my $organism = 'homo_sapiens';

    my $F_score_folder = $result_folder."/SNP_filter/F_score";
    my $allele_file = $F_score_folder."/F_score_".$sample_data->{result_dir}.".txt";

    if (! -f $allele_file) {
	$log = scalar(localtime())."  !!!missing file $allele_file\n";
	LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	return;
    }

    my $dbh = DBconnect::connect_db("SNP_annot");

    my $table = '37_annot_ensembl_'.$dbversion.'_full_descr';

#    my $org = $sample_data->{organism};
#    if (defined $org) {
#	$organism = lc($org);
#    }

    # Get registry
    my $registry = 'Bio::EnsEMBL::Registry';

    if ($local_server == 1) {
	$registry->load_registry_from_db(
	    -host => 'XXX',
	    -user => 'XXX',
	    -pass => 'XXX',
	    -db_version => $dbversion,
	    -port => $port
	    );
    } else {
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org',
	    -user => 'anonymous',
	    -db_version => $dbversion,
	    -port => $port
	    );
    }

    $registry->set_reconnect_when_lost();

    $log  = scalar(localtime())." : Connecting to Ensembl version ".$dbversion."\nSample = ".$sample_data->{result_dir}."\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    my $gene_adaptor = $registry->get_adaptor($organism, 'core', 'gene');
    my $slice_adaptor = $registry->get_adaptor($organism, 'core', 'slice');
    my $variation_adaptor = $registry->get_adaptor($organism, 'variation', 'variation');
    my $variation_feature_adaptor = $registry->get_adaptor($organism, 'variation', 'variationfeature');

    # Input file
    $log = scalar(localtime())." : Open $allele_file  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    open IN, $allele_file;

    while (my $line =  <IN>) {
	chomp $line;
	next if (($line =~ /^[\#|\;]/) || ($line eq ""));
	my ($chromosome, $start, $ref_allele, $allele1, $allele2, $hh, $depth, $score) = split (/\t/, $line);
	my $tiny_slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome, $start, $start);
	my $slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome);

	# Check if alleles are different from ref
	if ($allele1 ne $ref_allele) {
	    my $id = "snp_".$chromosome."_".$start."_".$allele1;
	    my $allele = $ref_allele."/".$allele1;
	    &variation_features($tiny_slice, $slice, $id, $chromosome, $start, $allele, $dbh, $table, $variation_feature_adaptor, $gene_adaptor);
	}

	if ($allele2 ne $ref_allele) {
	    my $id = "snp_".$chromosome."_".$start."_".$allele2;
	    my $allele = $ref_allele."/".$allele2;
	    &variation_features($tiny_slice, $slice, $id, $chromosome, $start, $allele, $dbh, $table, $variation_feature_adaptor, $gene_adaptor);
	}
    }

    # Close filehandles
    close IN;
    $dbh->disconnect();
    print scalar(localtime())." : ensembl SNP annot from $allele_file done\n";
#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
}

##################################################################################
sub variation_features {
    my ($tiny_slice, $slice, $id, $chromosome, $start, $allele, $dbh, $table, $variation_feature_adaptor, $gene_adaptor) = @_;

    # Check if variation is already in our database
    my $req = "select Local_ID from $table where " .
	" chrom = \"$chromosome\" " .
	" and Start = \"$start\" " .
	" and Local_alleles = \"$allele\" ";
    my $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
    $sth->execute()||  print ("\nExecution requete: $req echec\n");

    my $n = $sth->rows;

    if ($n == 0) {
	# Check if variation is already annotated in EnsEMBL
	my $variation_features = $variation_feature_adaptor->fetch_all_by_Slice($tiny_slice);

	foreach my $vf (@{$variation_features}){
	    if ($vf->allele_string eq $allele) {
		&get_ConsequenceTypes($table,$gene_adaptor,$vf,$allele,$id);
	    } else {
		my $new_variation_feature = Bio::EnsEMBL::Variation::VariationFeature->new(
		    -start => $start,
		    -end => $start,
		    -slice => $slice,           # the variation must be attached to a slice
		    -allele_string => $allele,    # the first allele should be the reference allele
		    -strand => 1,
		    -map_weight => 1,
		    -adaptor => $variation_feature_adaptor,           # we must attach a variation feature adaptor
		    -variation_name => "_".$vf->variation_name,
		    );
		&get_ConsequenceTypes($table,$gene_adaptor,$new_variation_feature,$allele,$id,$vf);
	    }
	}

	# Otherwise get consequence types for new variations
	unless (scalar @{$variation_features} > 0) {
	    # create a new VariationFeature object
	    my $new_variation_feature = Bio::EnsEMBL::Variation::VariationFeature->new(
		-start => $start,
		-end => $start,
		-slice => $slice,           # the variation must be attached to a slice
		-allele_string => $allele,    # the first allele should be the reference allele
		-strand => 1,
		-map_weight => 1,
		-adaptor => $variation_feature_adaptor,           # we must attach a variation feature adaptor
		-variation_name => $id,
		);
	    &get_ConsequenceTypes($table,$gene_adaptor,$new_variation_feature,$allele,$id);
	}
    }
}

##################################################################################
# Get the consequence types
sub get_ConsequenceTypes {

    my ($table,$gene_adaptor,$variation_feature,$allele,$id,$known_vf) = @_;
    my ($annot_allele, $seq_region_name, $seq_region_start, $seq_region_end, $consequence, $transcript, $transcript_id, $swiss_id, $aa_change, $gene_id, $gene_name, $biotype, $description, $refseq_mrna, $refseq_peptide);

    my $variation_name = $variation_feature->variation_name;

    if ($known_vf) {
	$annot_allele = $known_vf->allele_string;
    } else {
	$annot_allele = $variation_feature->allele_string;
    }

    $seq_region_name = $variation_feature->seq_region_name;
    $seq_region_start = $variation_feature->seq_region_start;
    $seq_region_end = $variation_feature->seq_region_end;

#   my $consequence_type = $variation_feature->get_consequence_type();
    my $consequence_type = $variation_feature->consequence_type();

    if ($consequence_type =~ /ARRAY/) {
	$consequence = join(",",@$consequence_type);
    } else {
	$consequence = $consequence_type;
    }

#    if ($consequence =~ 'INTERGENIC') {
    if ($consequence =~ 'intergenic_variant') {
	$transcript_id = 'NA';
	$swiss_id = "NA";
	$aa_change = "NA";
	$gene_id = "NA";
	$gene_name = "NA";
	$biotype = "NA";
	$description = "NA";
	$refseq_mrna = "NA";
	$refseq_peptide = "NA";

	# Insert into our SNP database
	my $dbh = DBconnect::connect_db("SNP_annot");
	my $req = "insert into $table set " .
#	    " Local_ID = \"$id\", " .
	    " Annotated_ID = \"$variation_name\", " .
	    " Local_alleles = \"$allele\", " .
	    " Annotated_alleles = \"$annot_allele\", " .
	    " chrom = \"$seq_region_name\", " .
	    " Start = \"$seq_region_start\", " .
	    " End = \"$seq_region_end\", " .
	    " Consequence = \"$consequence\", " .
	    " Transcript = \"$transcript_id\", " .
	    " Swissprot_ID = \"$swiss_id\", " .
	    " AA_change = \"$aa_change\", " .
	    " Gene_ID = \"$gene_id\", " .
	    " Gene_name = \"$gene_name\", " .
	    " Biotype = \"$biotype\", " .
	    " Gene_description = \"$description\", " .
	    " RefSeq_mRNA = \"$refseq_mrna\", " .
	    " RefSeq_peptide = \"$refseq_peptide\" ";

	my $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	$sth->execute()||  print ("\nExecution requete: $req echec\n");

    } else {
	my $transcript_consequences = $variation_feature->get_all_TranscriptVariations();
#	if (scalar(@{$transcript_consequences}) > 0) {
	foreach my $tr_consequence(@{$transcript_consequences}) {
	    $swiss_id = "NA";
	    $aa_change = "NA";
	    $gene_id = "NA";
	    $gene_name = "NA";
	    $biotype = "NA";
	    $description = "NA";
	    $refseq_mrna = "NA";
	    $refseq_peptide = "NA";
	    foreach my $feature_type(@{$tr_consequence->consequence_type}) {
		$consequence = $feature_type;
		$transcript_id = $tr_consequence->transcript->stable_id;
		$transcript = $tr_consequence->transcript();

		if (defined $tr_consequence->pep_allele_string) {
		    $aa_change = $tr_consequence->pep_allele_string;
#		    if (($feature_type eq 'NON_SYNONYMOUS_CODING') || ($feature_type eq 'STOP_GAINED')){
		    if (($feature_type eq 'missense_variant') || ($feature_type eq 'initiator_codon_variant') || ($feature_type eq 'stop_gained')){
			my @aas = split (/\//, $aa_change);
			my $strand = $transcript->strand();
			my $translation = $transcript->translation();

			my @swissprot_entries = @{$translation->get_all_DBEntries('Uniprot/SWISSPROT')};

			foreach my $swissprot_entry (@swissprot_entries) {
			    $swiss_id = $swissprot_entry->primary_id();
			}
#			unless ($swiss_id) {$swiss_id = 'NA';}

			my @pep_coords = $transcript -> genomic2pep($variation_feature->seq_region_start, $variation_feature->seq_region_end, $strand);
			my $pep_pos;
			foreach my $pep_coord (@pep_coords) {
			    $pep_pos = $pep_coord->{'start'};
			}
			$aa_change = $aas[0].$pep_pos.$aas[1];
#		    } else {
#			$swiss_id = 'NA';
#			unless ($aa_change) {$aa_change = 'NA';}
		    }
		}

		my $gene = $gene_adaptor->fetch_by_transcript_id($tr_consequence->transcript->dbID);
		if (defined $gene->stable_id) {
		    $gene_id = $gene->stable_id;
		}
		if (defined $gene->external_name) {
		    $gene_name = $gene->external_name;
		}
		if (defined $gene->biotype) {
		    $biotype = $gene->biotype;
		}
		if (defined $gene->description) {
		    $description = $gene->description;
		}

#		my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
		my @dbentries = @{$transcript->get_all_DBLinks('RefSeq%')};
		foreach my $dbentry(@dbentries) {
		    if ($dbentry->dbname() eq 'RefSeq_mRNA') {
			$refseq_mrna = $dbentry->display_id;
		    } elsif ($dbentry->dbname() eq 'RefSeq_peptide') {
			$refseq_peptide = $dbentry->display_id;
		    }
		}

		# Insert into our SNP database
		my $dbh = DBconnect::connect_db("SNP_annot");
		my $req = "insert into $table set " .
#		    " Local_ID = \"$id\", " .
		    " Annotated_ID = \"$variation_name\", " .
		    " Local_alleles = \"$allele\", " .
		    " Annotated_alleles = \"$annot_allele\", " .
		    " chrom = \"$seq_region_name\", " .
		    " Start = \"$seq_region_start\", " .
		    " End = \"$seq_region_end\", " .
		    " Consequence = \"$consequence\", " .
		    " Transcript = \"$transcript_id\", " .
		    " Swissprot_ID = \"$swiss_id\", " .
		    " AA_change = \"$aa_change\", " .
		    " Gene_ID = \"$gene_id\", " .
		    " Gene_name = \"$gene_name\", " .
		    " Biotype = \"$biotype\", " .
		    " Gene_description = \"$description\", " .
		    " RefSeq_mRNA = \"$refseq_mrna\", " .
		    " RefSeq_peptide = \"$refseq_peptide\" ";

		my $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
		$sth->execute()||  print ("\nExecution requete: $req echec\n");
	    }
	}
    }
}

##################################################################################
### Convert Grch38 to hg19 and annotate with hg19 positions
sub launch_hg19_pos_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);

    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_hg19_pos';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

### Ajout mehdi 15/11/18 pour ajouter colonne pos_hg19 dans le script "start_analysis"
### Ajout mehdi 05/12/18 pour ajouter colonne pos_hg19 dans le script "start_analysis"

    if ($fam == 0){
  $fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
  $fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
  AnnotSNP::hg19_pos_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
  my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
  my @filesIn = @$filesIn_ref;
  my @filesOut = @$filesOut_ref;
  if ($combinaison == 0){
    $fileIn = $filesIn[0];
    $fileOut = $filesOut[0];
    AnnotSNP::hg19_pos_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
  } elsif ($combinaison == 1) {
      for (my $i = 0; $i < scalar @filesIn; $i++) {
    AnnotSNP::hg19_pos_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
    }
}
  }
    print "hg19 annot snp from $from DONE\n";

}

##################################################################################
### Annotate with hg19 position
sub hg19_pos_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $chrom, $pos, $hg19_pos, @tab);

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    my $dbh = DBconnect::connect_db("SNP_annot");
    my $table = "mut_convert_pos";
    my ($req,$sth,$req2,$sth2);

    $log ="open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

#    $log = "check for position in hg19";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\thg19_pos\n";
	    next;
	}
	@tab = split(/\t/, $line);

	# Check for position in hg19
	$chrom = $tab[0];
	$pos = $tab[1];


	### check locally
	$req = "select hg19_pos from $table where chrom = \"$chrom\" and hg38_pos = \"$pos\" ";
	$sth = $dbh->prepare($req) ||  print ("\nPreparation requete $req: echec\n");
	$sth->execute() ||  print ("\nExecution requete: $req echec\n");

	## if exists, print
	if ($sth->rows > 0) {
	    for ( my @data = $sth->fetchrow_array ) {
		print FHO "$line\t$data[0]\n";
	    }
	} else {
	    ## else, request and insert
	    my $http = HTTP::Tiny->new();

	    my $server = 'http://rest.ensembl.org';
	    my $ext = "/map/human/GRCh38/${chrom}:${pos}..${pos}:1/GRCh37?";
	    my $response = $http->get($server.$ext, {
		headers => { 'Content-type' => 'application/json' }
				      });

	    $log .= "$chrom: $pos not found in hg19 !\n" unless $response->{success};

	    ###
	    my $status = $response->{status};
	    my $remain = $response->{'headers'}->{'x-ratelimit-remaining'};
	    my $sleeptime = $response->{'headers'}->{'x-ratelimit-reset'};
	    if ($remain < 10){
		print "# remain $remain requests... sleep for $sleeptime seconds.\n";
		sleep($sleeptime);
		$response = $http->get($server.$ext, {
		    headers => { 'Content-type' => 'application/json' }
				       });
	    }
	    if ($status == 429){
		print "status $status sleep for  $sleeptime seconds.\n";
		sleep($sleeptime);
		$response = $http->get($server.$ext, {
		    headers => { 'Content-type' => 'application/json' }
				       });
	    }

	    if(length $response->{content}) {
		my $hash = decode_json($response->{content});
		$hg19_pos = $hash->{'mappings'}[0]{'mapped'}{'start'};
		unless (defined $hg19_pos) {
		    $hg19_pos = "NA";
		}
	    } else {
		$hg19_pos = "NA";
	    }

	    print FHO "$line\t$hg19_pos\n";

	    $req2 = "insert into $table set chrom = \"$chrom\", hg38_pos = \"$pos\", hg19_pos = \"$hg19_pos\" ";
	    $sth2 = $dbh->prepare($req2)||  print ("\nPreparation requete $req2: echec\n");
	    $sth2->execute()||  print ("\nExecution requete: $req2 echec\n");

	}
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);

    }
    #$sth->finish();

    close(FHI);
    close(FHO);
    $log = "$n_out mutation(s) not found in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    #$dbh->disconnect();

}

##################################################################################
### Annotate with control covers when available (otherwise 0)
sub launch_cover_other_snp {
    my $run_name = shift;
    my $analysis_data = shift;
    my $combinaison = shift;
    my $recess = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my ($log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $run_name;
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $from_short = "F_coding";
    $from_short =~ s/F_//;

    my $cover_other_folder = $result_folder."/SNP_filter/F_cover_other";
    system("mkdir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/F_cover_other") unless (-d "$result_folder/SNP_filter/F_cover_other");
    system("mkdir $result_folder/SNP_filter/F_cover_other_combi") unless (-d "$result_folder/SNP_filter/F_cover_other_combi");

    # On parcourt tous les samples to get file names
    my @dna_cases;
    my @dna_controls;
    my @parent_dna_controls;
    my @other_dna_controls;
    my $family;
    my $sickness;
    foreach my $k (keys %$analysis_data) {
	if ($$analysis_data{$k}{'atteint'} == 1) {
	    push @dna_cases, $$analysis_data{$k}{'result_dir'};
	}
	if ($$analysis_data{$k}{'atteint'} == 0) {
		push @dna_controls, $$analysis_data{$k}{'result_dir'};
	    if (($$analysis_data{$k}{'parent'} eq 'M') || ($$analysis_data{$k}{'parent'} eq 'F')) {
		push @parent_dna_controls, $$analysis_data{$k}{'result_dir'};
	    } else {
		push @other_dna_controls, $$analysis_data{$k}{'result_dir'};
	    }
	}
	$family = $$analysis_data{$k}{'pop'};
	$sickness = $$analysis_data{$k}{'maladie'};
    }
    my $cases = join('_',  sort @dna_cases);
    my ($controls);

    my $fileOut;
    my $fileIn;

    my $suffix;
    if ($recess == 1) {
	$suffix = "_recess";
    } else {
	$suffix = "";
    }

    #Cas normal : cas versus controles : ecrit dans "cover_other"
    if (@parent_dna_controls && @other_dna_controls) {
	$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @other_dna_controls);
    } elsif (@parent_dna_controls && !(@other_dna_controls)) {
	$controls = join('_', sort @parent_dna_controls);
    } elsif (!(@parent_dna_controls) && @other_dna_controls) {
	$controls = join('_', sort @other_dna_controls);
    }

    if(@dna_controls) {

#   if(defined(@dna_controls)) {
    # if(@dna_controls) {
    # 	if ($recess == 1) {
    # 	    my @non_parent_dna_controls = @dna_controls;
    # 	    splice(@non_parent_dna_controls, 0, 2);
    # 	    my @parent_dna_controls = ($dna_controls[0], $dna_controls[1]);
    # 	    if (@non_parent_dna_controls) {
    # 		$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @non_parent_dna_controls);
    # 	    } else {
    # 		$controls = join('_', sort @parent_dna_controls);
    # 	    }
    # 	} else {
    # 	    $controls = join('_',  sort @dna_controls);
    # 	}
	$fileOut = $cover_other_folder."/F_cover_other_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	$fileIn = $result_folder."/SNP_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	cover_other_snp($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 0, "none", \@dna_controls);
    } else {
	$fileOut = $cover_other_folder."/F_cover_other_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	$fileIn = $result_folder."/SNP_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	cover_other_snp($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 0, "none", "none");
    }

    #Autre cas qui teste toutes les autres combinaisons
    if ($combinaison == 1) {
    	my @bitMask = ();
    	#On cherche tous les subsets possibles pour la liste d'adn
    	while (FilterSNP::generate_mask(\@bitMask, \@dna_cases)) {
    	    my $d_cases = FilterSNP::get_one_subset(\@bitMask, \@dna_cases);

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
    		my @cases_as_controls = @d_controls;
    		#On ajoute nos "vrais" controles
    		# a cette liste de controles
    		if(@dna_controls) {
    		    push @d_controls, @dna_controls;
    		}
    		my $ca = join('_', sort @$d_cases);
    		my $co = join('_', sort @d_controls);
		$fileOut = $cover_other_folder."_combi/F_cover_other_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		$fileIn = $result_folder."/SNP_filter/${from}_combi/${from}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		cover_other_snp($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 1, \@cases_as_controls, \@dna_controls);
    	    }
    	}
    }
}

##################################################################################
sub cover_other_snp {
    my ($logFolder, $logFile, $result_folder, $fileIn, $fileOut, $filter, $ref_cases, $ref_controls) = @_;
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my ($log, $n_in, $n_out, $chrom, $pos, @tab);
    my ($controls_cover, $cover, $flag, $cases_cover);

    my @dna_cases;
    if ($ref_cases eq "none") {
	$dna_cases[0] = "none";
    } else {
	@dna_cases = @{$ref_cases};
    }
    my @dna_controls;
    if ($ref_controls eq "none") {
	$dna_controls[0] = "none";
    } else {
	@dna_controls = @{$ref_controls};
    }

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	#print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log ="open $fileIn\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;
    my @threads;
    my %results;

    $log = "";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    next;
	}

	if ($n_in % 10 == 0) {
#	    print "n : $n_in \n";
	    for my $thr (@threads) {
		my $res = $thr->join();
#		print "@$res\n";
		$results{$$res[0], $$res[1]} = $$res[2];
	    }
	    undef(@threads);
	}

	@tab = split(/\t/, $line);
	$chrom = $tab[0];
	$pos = $tab[1];

	# for controls
	$controls_cover = '';
	if ($dna_controls[0] ne "none") {
	    foreach my $control (@dna_controls) {
		my $file = "$result_folder/SEQ_cover/all/$control/chr$chrom.cov";
		my @arg = ($file, $pos, $n_in, $control);
		my $thread = threads->create('grep_cover', @arg);
		push @threads, $thread;
	    }
	}

	# for cases used as controls (combination cases)
	if ($dna_cases[0] ne "none") {
	    foreach my $case (@dna_cases) {
		my $file = "$result_folder/SEQ_cover/all/$case/chr$chrom.cov";
		my @arg = ($file, $pos, $n_in, $case);
		my $thread = threads->create('grep_cover', @arg);
		push @threads, $thread;
	    }
	}
    }
    close(FHI);

    for my $thr (@threads) {
#	  print "wait $thr \n";
	  my $res = $thr->join();
#	  print "@$res\n";
	  $results{$$res[0], $$res[1]} = $$res[2];
    }

    $n_in = 0;
#   print "open again $fileIn dna=$dna\n";
    open (FHI2, $fileIn) or die("cannot open $fileIn");

    while (my $line = <FHI2>) {
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    print FHO "$line\tcover_other\n";
	    next;
	}
	$n_in++;

	# for controls
	$controls_cover = '';
	unless ($dna_controls[0] eq "none") {
	    foreach my $control (@dna_controls) {
		$cover = $results{$n_in, $control};
		$controls_cover .= "$control:$cover;";
	    }
	}


	$controls_cover =~ s/;$//;

	# for cases used as controls (combination cases)
	if ($dna_cases[0] ne "none") {
	    $cases_cover = '';
	    $flag = 0;
	    foreach my $case (@dna_cases) {
		$cover = $results{$n_in, $case};
		$cases_cover .= "$case:$cover;";
		if($filter eq '1') {
		    if ($cover >= 8) {
			# Pas de flag si on veut garder tous les SNPs
			$flag = 1;
		    }
		}
	    }
	    $cases_cover =~ s/;$//;

	    unless ($flag == 1) {
		if ($dna_controls[0] eq "none") {
		    print FHO "$line\t$cases_cover\n";
		} else {
		    print FHO "$line\t$controls_cover;$cases_cover\n";
		}
		$n_out++;
	    }
	} else {
	    print FHO "$line\t$controls_cover\n";
	    $n_out++;
	}
    }

#   LogOut::printLog($log, $current_sub, $analysis_data->{pop}, $analysis_data->{dna});

    close(FHI2);
    close(FHO);
    if ($filter == 1) {
	$log = "$n_out mutation(s) with insufficient cover in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    } else {
	$log = "$n_out mutation(s) written in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    }
    #print $log;
   LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
sub grep_cover {
    my $file = shift @_;
    my $pos = shift @_;
    my $line = shift @_;
    my $dna = shift @_;

#   print "launch grep $pos $file \n";
    my $cover_line;
    if (-f $file) {
	$cover_line = `grep -w $pos $file`;
    } else {
	my $gz_file = $file.".gz";
	$cover_line = `zgrep -w $pos $gz_file`;
    }
    my @cover_tab = split(/ /, $cover_line);
    my $cover = $cover_tab[1];
    if ($cover) {
	chomp $cover;
    } else {
	$cover = 0;
    }

#   print "finished $dna, $line = $cover\n" if ($cover == 0);

    my @tab = ($line, $dna, $cover);

    return \@tab;

}

##################################################################################
### Prepare for dbSNP_annot_snp
sub launch_dbSNP_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $dbversion = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);

    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_dbSNP';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotSNP::dbSNP_annot_snp($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotSNP::dbSNP_annot_snp($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotSNP::dbSNP_annot_snp($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "dbSNP annot snp from $from DONE\n";
}

##################################################################################
### Annotate with dbSNP ID when available (otherwise 0)
sub dbSNP_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $dbversion = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $start, $end, $req, $sth, @row, @tab);

    my $dbh = DBconnect::connect_db("filtre_seq");

    my $not_in_dbsnp = 0;

    $log = "dbSNP : $dbversion\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log ="open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "check for mutation in dbSNP\n";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\t$dbversion\n";
	    next;
	}
	@tab = split(/\t/, $line);

	# Check for mutation in dbSNP
	if ($tab[18] ne "NA") {
#	    $start = $tab[1] - 1;
	    $start = $tab[18] - 1;  #col with hg19 pos
#	    $end = $tab[1];
	    $end = $tab[18];  #col with hg19 pos
	} else {
	    print FHO "$line\tNA\n";
	    next;
	}

	my @alleles;

	if ($fam == 1) {
	    my @ind_list = split(/\] /, $tab[16]);
	    foreach my $indiv(@ind_list) {
		$indiv =~ s/\]//;
		$indiv =~ s/.*\[//;
		$indiv =~ s/,.*//;
		my @ind_alleles = split(/\//, $indiv);
		foreach my $ia(@ind_alleles) {
		    if ($ia ne $tab[2]) {
			push(@alleles,$ia);
		    }
		}
	    }
	} else {
	    my @two_alleles = ($tab[3], $tab[4]);
	    foreach my $ia(@two_alleles) {
		if ($ia ne $tab[2]) {
		    push(@alleles,$ia);
		}
	    }
	}

	# Remove duplicate alleles
	my %seen=();
	my @uniq_alleles = grep { ! $seen{$_} ++ } @alleles;

	my $in_dbsnp;

	$req = "select observed, strand, refNCBI, name, alleleFreqs from $dbversion where " .
	      " chrom = \"chr$tab[0]\" " .
	      " and chromEnd = \"$end\" " .
	      " and class = \'single\' ".
	      " and refNCBI != \'-\' ";

	$sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	$sth->execute()||  print ("\nExecution requete: $req echec\n");

	while ( my @data = $sth->fetchrow_array ) {
	    my $observed = $data[0];
	    my $strand = $data[1];
	    my $ref_ncbi = $data[2];
	    my $name = $data[3];
	    my $freq = $data[4];
	    chop($freq);
	    unless ($freq) {
		$freq = 'NA';
	    }

	    my $ref_eq_obs = 0;

	    my @observed = split('/', $observed);
	    my @reversed_observed;

	    if ($strand eq '-') {
		# For each observed allele
		for my $obs (@observed) {
		    my $seq_obj = Bio::Seq->new(-seq => $obs, -alphabet => 'dna' );
		    my $reversed_obj = $seq_obj->revcom;
		    my $reversed_obs = $reversed_obj->seq;
		    push @reversed_observed, $reversed_obs;
		    # Check that ref_NCBI corresponds to at least one observed. If it is the case, one allele must correspond to at least one observed
		    if($ref_ncbi eq $reversed_obs){
			$ref_eq_obs = 1;
		    }
		}
	    } else {
		for my $obs (@observed) {
		    if ($ref_ncbi eq $obs){
			$ref_eq_obs = 1;
		    }
		}
	    }

	    # If ref equals none of the observed alleles: error and test the next
	    if ($ref_eq_obs == 0) {
		next;
	    }

	    for (my $i = 0; $i < scalar(@uniq_alleles); $i ++) {
		my $allele = $uniq_alleles[$i];
		if ($strand eq '-') {
		    # For each observed allele
		    for my $obs (@reversed_observed) {
			if ($obs eq $allele) {
			    $in_dbsnp = $in_dbsnp . "$allele:$name($freq);";
			    # Allele removed from table
			    splice (@uniq_alleles, $i);
			    last;
			}
		    }
		} else {
		    for my $obs (@observed) {
			if ($obs eq $allele) {
			    $in_dbsnp = $in_dbsnp . "$allele:$name($freq);";
			    # Allele removed from table
			    splice (@uniq_alleles, $i);
			    last;
			}
		    }
		}
	    }
	}

	# Then set all the ones not found to 0
	foreach my $allele (@uniq_alleles) {
	    $not_in_dbsnp = 1;
	    $in_dbsnp = $in_dbsnp . "$allele:0;";
	    $n_out++;
	}

	chop($in_dbsnp);

	if ($filter == 1) {
	    if ($not_in_dbsnp == 0) {
		next;
	    } else {
		$not_in_dbsnp = 0;
	    }
	}

	print FHO "$line\t$in_dbsnp\n";

    }

    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

    close(FHI);
    close(FHO);
    $dbh->disconnect();
    $log = "$n_out mutation(s) not found in db $dbversion in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
### Prepare for CG_54genomes_annot_snp
sub launch_CG_54genomes_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_54genomes';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotSNP::CG_54genomes_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotSNP::CG_54genomes_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotSNP::CG_54genomes_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "54genomes annot snp from $from DONE\n";
}

##################################################################################
### Annotate with observed frequency when present in Complete Genomics 54 genomes (otherwise 0)
sub CG_54genomes_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $start, $end, $req, $sth, @row, @tab);

    my $dbh = DBconnect::connect_db("SNP_annot");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
#	print $log;
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
#	print $log;
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "check for mutation in 54genomes\n";
    while (my $line = <FHI>) {
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\t54_genomes\n";
	    next;
	}
	@tab = split(/\t/, $line);

	# Check for mutation in 54 genomes
	if ($tab[18] ne "NA") {
#	    $start = $tab[1];
	    $start = $tab[18];  #col with hg19 pos
	} else {
	    print FHO "$line\tNA\n";
	    next;
	}


	my @alleles;

	if ($fam == 1) {
	    my @ind_list = split(/\] /, $tab[16]);
	    foreach my $indiv(@ind_list) {
		$indiv =~ s/\]//;
		$indiv =~ s/.*\[//;
		$indiv =~ s/,.*//;
		my @ind_alleles = split(/\//, $indiv);
		foreach my $ia(@ind_alleles) {
		    if ($ia ne $tab[2]) {
			push(@alleles,$ia);
		    }
		}
	    }
	} else {
	    my @two_alleles = ($tab[3], $tab[4]);
	    foreach my $ia(@two_alleles) {
		if ($ia ne $tab[2]) {
		    push(@alleles,$ia);
		}
	    }
	}

	# Remove duplicate alleles
	my %seen=();
	my @uniq_alleles = grep { ! $seen{$_} ++ } @alleles;

	my $freq_alleles;

	foreach my $allele(@uniq_alleles) {
	    $req = "select frequency from 54genomes where " .
		" chrom = \"chr$tab[0]\" " .
		" and pos = \"$start\" " .
		" and ref = \"$tab[2]\" ".
		" and allele = \"$allele\" ";

	    $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	    $sth->execute()||  print ("\nExecution requete: $req echec\n");

	    if ($sth->rows > 0) {
		# There is at least one line so the mutation is found in at least one of the 54 genomes
		while (my @func = $sth->fetchrow_array) {
		    $freq_alleles = $freq_alleles ."$allele:$func[0];";
		}
		$sth->finish();
	    } else {
		# Unknown mutation, set to 0
		$freq_alleles = $freq_alleles ."$allele:0;";
		$n_out++;
	    }
	}
	chop($freq_alleles);
	print FHO "$line\t$freq_alleles\n";

    }

    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

    close(FHI);
    close(FHO);
    $log  = "$n_out mutation(s) not found in 54 genomes in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
#    print $log;
    #$dbh->disconnect();

}

##################################################################################
### Prepare for expr_beta_annot_snp
sub launch_expr_beta_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $data_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_expr';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotSNP::expr_beta_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotSNP::expr_beta_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotSNP::expr_beta_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "expr_beta annot snp from $from DONE\n";
}

##################################################################################
### Annotate with 1 or 0 depending on expression in beta cell (lab data)
sub expr_beta_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $data_folder = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $expr, @tab);

    my $expr_file = "$data_folder/Genes_Expr_Beta_union.txt";

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

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
	    print FHO "$line\tExpr\n";
	    next;
	}
	@tab = split(/\t/, $line);

	# Check if gene is expressed, based on gene symbol
	if ($fam == 0) {
	    $expr = `grep -w $tab[10] $expr_file | wc -l`;
	} else {
	    $expr = `grep -w $tab[5] $expr_file | wc -l`;
	}

	if ($expr > 0) {
	    # There is at least one line so the gene is expressed
	    print FHO "$line\t1\n";
	    $n_out++;
	} else {
	    # No data for this gene, set to 0
	    print FHO "$line\t0\n";
	}
    }

#    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) in genes expressed in beta cell in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
### Prepare for good_annot_snp
sub launch_good_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $data_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_good';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotSNP::good_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotSNP::good_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotSNP::good_annot_snp($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "GOOD annot snp from $from DONE\n";
}

##################################################################################
### Annotate with GOOD diabetes or obesity ranking of genes where snps are found
sub good_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $data_folder = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $good_file, %good_data, $good_line, @good_tab, $good_score, @tab);

    my $disease;
    if ($fam == 0) {
	if ($sample_data->{maladie} eq 'diab'){
	    $disease = 'diab';
	    $good_file = "$data_folder/GOOD_no_plus.txt";
	} elsif ($sample_data->{maladie} eq 'obes'){
	    $disease = 'obes';
	    $good_file = "$data_folder/OB_GOOD_ranks.txt";
	}
    } elsif ($fam == 1) {
	foreach my $k (keys %$sample_data) {
	    if ($$sample_data{$k}->{maladie} eq 'diab') {
		$disease = 'diab';
		$good_file = "$data_folder/GOOD_no_plus.txt";
	    } elsif ($$sample_data{$k}->{maladie} eq 'obes'){
		$disease = 'obes';
		$good_file = "$data_folder/OB_GOOD_ranks.txt";
	    }
	}
    }

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    if (($sample_data->{maladie} eq 'diab') || ($sample_data->{maladie} eq 'obes')){
	open (FHG, $good_file) or die ("Cannot open $good_file");
	while ($good_line = <FHG>) {
	    $good_line =~ s/\012//;
	    $good_line =~ s/\015//;
	    if ($good_line =~ /^\#/) {
		next;
	    }
	    @good_tab = split(/\t/, $good_line);
#	    if ($sample_data->{maladie} eq 'diab'){
	    if ($disease eq 'diab'){
		$good_data{$good_tab[0]} = $good_tab[1]."\t".$good_tab[2];
#	    } elsif ($sample_data->{maladie} eq 'obes'){
	    } elsif ($disease eq 'obes'){
	    $good_data{$good_tab[0]} = $good_tab[1];
	    }
	}
	close (FHG);
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
	    print FHO "$line\tGOOD\n";
	    next;
	}

	if ($sample_data->{maladie} eq 'unkn'){
	    print FHO "$line\tNA\n";
	    next;
	} elsif ($sample_data->{maladie} eq 'lipo'){
	    print FHO "$line\tNA\n";
	    next;
	} elsif ($sample_data->{maladie} eq 'neph'){
	    print FHO "$line\tNA\n";
	    next;
	} elsif ($sample_data->{maladie} eq 'anor'){
	    print FHO "$line\tNA\n";
	    next;
	}

	@tab = split(/\t/, $line);

	# Check whether gene is ranked by gene symbol in GOOD
	if ($fam == 0) {
	    $good_line = $good_data{$tab[10]};
	} else {
	    $good_line = $good_data{$tab[5]};
	}

	if ($good_line) {
	    # There is at least one line so the gene is ranked in GOOD
#	    if ($sample_data->{maladie} eq 'diab'){
	    if ($disease eq 'diab'){
		@good_tab = split(/\t/, $good_line);
		unless ($good_tab[1] =~ 'VALEUR') {
		    $good_score = $good_tab[1];
		} else {
		    $good_score = $good_tab[0];
		}
#	    } elsif ($sample_data->{maladie} eq 'obes'){
	    } elsif ($disease eq 'obes'){
		$good_score = $good_line;
	    }
	    $good_score =~ s/\012//;
	    $good_score =~ s/\015//;
	    print FHO "$line\t$good_score\n";
	    $n_out++;
	} else {
	    # Check if gene is ranked by its Ensembl transcript ID
	    if ($fam == 0) {
		$good_line = $good_data{$tab[8]};
	    } else {
		$good_line = $good_data{$tab[3]};
	    }

	    if ($good_line) {
		# There is at least one line so the gene is ranked in GOOD
#		if ($sample_data->{maladie} eq 'diab'){
		if ($disease eq 'diab'){
		    @good_tab = split(/\t/, $good_line);
		    unless ($good_tab[1] =~ 'VALEUR') {
			$good_score = $good_tab[1];
		    } else {
			$good_score = $good_tab[0];
		    }
#		} elsif ($sample_data->{maladie} eq 'obes'){
		} elsif ($disease eq 'obes'){
		    $good_score = $good_line;
		}
		$good_score =~ s/\012//;
		$good_score =~ s/\015//;
		print FHO "$line\t$good_score\n";
		$n_out++;
	    } else {
		# Gene is not ranked, set to NA
		print FHO "$line\tNA\n";
	    }
	}
    }

#    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) in genes ranked in GOOD in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
}

##################################################################################
### Prepare for dbnsfp_annot_snp
sub launch_dbnsfp_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    if ($fam == 0) {
	if ($sample_data->{maladie} eq 'unkn'){
	    $from = "F_54genomes";
	}
    }

    my $to = 'F_dbNSFP';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	#AnnotSNP::dbnsfp_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	AnnotSNP::dbnsfp_annot_snp_mongo($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    #AnnotSNP::dbnsfp_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	    AnnotSNP::dbnsfp_annot_snp_mongo_fam($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		#AnnotSNP::dbnsfp_annot_snp($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
		AnnotSNP::dbnsfp_annot_snp_mongo_fam($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "dbNSFP annot snp from $from DONE\n";
}

##################################################################################
### Annotate with results of SIFT, PolyPhen2, LRT, MutationTaster, MutationAssessor and FATHMM functional predictions,
### PhyloP conservation score, 1000 genomes and ESP MAFs (from dbNSFP)
sub dbnsfp_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $pos, $line, $req, $sth, @tab);

    my $dbh = DBconnect::connect_db("SNP_annot");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

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
	    ###print FHO "$line\tPPh2_Humdiv\tPPh2_Humvar\tSIFT\tMutationTaster\tLRT\tMutationAssessor\tFATHMM\tPhyloP\tFreq_1KG\tFreq_African_1KG\tFreq_European_1KG\tFreq_American_1KG\tFreq_Asian_1KG\tFreq_African_American_ESP6500\tFreq_European_American_ESP6500\tAncestral_allele\tPathway\tFunction\tDisease\tMIM_id\tMIM_disease\tTrait_association\tExpression_egenetics\tExpression_GNF_Atlas\tInteractions_IntAct\tInteractions_BioGRID\n";
	    print FHO "$line\tPPh2_Humdiv\tPPh2_Humvar\tSIFT\tMutationTaster\tLRT\tMutationAssessor\tFATHMM\tPhyloP\tFreq_1KG\tFreq_African_1KG\tFreq_European_1KG\tFreq_American_1KG\tFreq_East_Asian_1KG\tFrequence_South_Asian_1KG\tFreq_African_American_ESP6500\tFreq_European_American_ESP6500\tAncestral_allele\tPathway\tFunction\tDisease\tMIM_id\tMIM_disease\tTrait_association\tExpression_egenetics\tExpression_GNF_Atlas\tInteractions_IntAct\tInteractions_BioGRID\n";
	    next;
	}

	@tab = split(/\t/, $line);

	my $alleles = '';
	if ($fam == 0) {
	    $alleles = $tab[3].$tab[4];
	} else {
 #modif mehdi 03/10/18 :

	    my @individus = split(/\] /, $tab[16]);
 #
	    foreach my $individu(@individus) {
		$individu =~ s/.*\///; # mutated allele is always the second for heterozygotes
		$individu =~ s/,.*//;
		if ($individu ne $alleles) {
		    $alleles .= $individu;


		}
	    }
	}

	$req = "select alt, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred , SIFT_pred, 1000Gp3_AF, 1000Gp3_AFR_AF, 1000Gp3_EUR_AF, 1000Gp3_AMR_AF, 1000Gp3_EAS_AF, 1000Gp3_SAS_AF, Ancestral_allele, genename, MutationAssessor_pred, FATHMM_pred, PhyloP7way_vertebrate, ESP6500_AA_AF, ESP6500_EA_AF from dbNSFP_v3_0b2a_variants where chr = \"$tab[0]\" and pos = \"$tab[1]\"";
	###$req = "select alt, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred , SIFT_score, 1000Gp1_AF, 1000Gp1_AFR_AF, 1000Gp1_EUR_AF, 1000Gp1_AMR_AF, 1000Gp1_ASN_AF, Ancestral_allele, genename, MutationAssessor_pred, FATHMM_score, PhyloP, ESP6500_AA_AF, ESP6500_EA_AF from dbNSFP_v2_0_variants where chr = \"$tab[0]\" and pos = \"$tab[1]\"";

	$sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	$sth->execute()||  print ("\nExecution requete: $req echec\n");

	if ($sth->rows > 0) {
	    # There is at least one line, the mutation is found in the dbNSFP_variants table
	    my ($Alt, $PPh2_hdiv, $PPh2_hvar, $LRT, $MT, $SIFT, $G1000_AF, $G1000_AFR, $G1000_EU, $G1000_AM, $G1000_EAS, $G1000_SAS, $anc_all, $genename, $MutAss, $FATHMM, $PhyloP, $ESP_AA, $ESP_EA);
	    while (my @funcpreds = $sth->fetchrow_array) {
		$Alt = $funcpreds[0];
		$genename = $funcpreds[12];
		if ($alleles =~ /$Alt/) {
		    if ($funcpreds[1] ne '.') { # missing data : "."
			$PPh2_hdiv .= $funcpreds[1];
		    }
		    if ($funcpreds[2] ne '.') {
			$PPh2_hvar .= $funcpreds[2];
		    }
		    if ($funcpreds[3] ne '.') {
			$LRT .= $funcpreds[3];
		    }
		    if ($funcpreds[4] ne '.') {
			$MT .= $funcpreds[4];
		    }
		    #if ($funcpreds[5] ne '.') {
			#if ($funcpreds[5] < 0.05) {
			    #$SIFT = 'D';
			#} else {
			    #$SIFT = 'T';
			#}
		    #}
		    if ($funcpreds[5] ne '.') {
			$SIFT .= $funcpreds[5];
		    }
		    if ($funcpreds[6] ne '.') {
		       $G1000_AF .= $funcpreds[6];
		    }
		    if ($funcpreds[7] ne '.') {
		        $G1000_AFR .= $funcpreds[7];
		    }
		    if ($funcpreds[8] ne '.') {
		        $G1000_EU .= $funcpreds[8];
		    }
		    if ($funcpreds[9] ne '.') {
		        $G1000_AM .= $funcpreds[9];
		    }
		    if ($funcpreds[10] ne '.') {
		        $G1000_EAS .= $funcpreds[10];
		    }
		    if ($funcpreds[11] ne '.') {
		        $G1000_SAS .= $funcpreds[11];
		    }
		    if ($funcpreds[12] ne '.') {
		        $anc_all .= $funcpreds[12];
		    }
		    if ($funcpreds[14] ne '.') {
			unless ((defined($MutAss)) && ($funcpreds[14] =~ $MutAss)) {
			    $MutAss .= $funcpreds[14];
			}
		    }
		    if ($funcpreds[15] ne '.') {
			unless ((defined($FATHMM)) && ($funcpreds[15] =~ $FATHMM)) {
			    #if ($funcpreds[15] < -1.5) {
				#$FATHMM .= 'D';
			    #} else {
				#$FATHMM .= 'T';
			    #}
			    $FATHMM .= $funcpreds[15];
			}
		    }
		    if ($funcpreds[16] ne '.') {
			unless ((defined($PhyloP)) && ($funcpreds[16] =~ $PhyloP)) {
			    $PhyloP .= $funcpreds[16];
			}
		    }
		    if ($funcpreds[17] ne '.') {
		        $ESP_AA .= $funcpreds[17];
		    }
		    if ($funcpreds[18] ne '.') {
		        $ESP_EA .= $funcpreds[18];
		    }
		}
	    }

	    unless (defined($PPh2_hdiv)) {
		$PPh2_hdiv = 'NA';
	    }
	    unless (defined($PPh2_hvar)) {
		$PPh2_hvar = 'NA';
	    }
	    unless (defined($LRT)) {
		$LRT = 'NA';
	    }
	    unless (defined($MT)) {
		$MT = 'NA';
	    }
	    unless (defined($SIFT)) {
		$SIFT = 'NA';
	    }
	    unless (defined($G1000_AF)) {
		$G1000_AF = 'NA';
	    }
	    unless (defined($G1000_AFR)) {
		$G1000_AFR = 'NA';
	    }
	    unless (defined($G1000_EU)) {
		$G1000_EU = 'NA';
	    }
	    unless (defined($G1000_AM)) {
		$G1000_AM = 'NA';
	    }
	    unless (defined($G1000_EAS)) {
		$G1000_EAS = 'NA';
	    }
	    unless (defined($G1000_SAS)) {
		$G1000_SAS = 'NA';
	    }
	    unless (defined($anc_all)) {
		$anc_all = 'NA';
	    }
	    unless (defined($MutAss)) {
		$MutAss = 'NA';
	    }
	    unless (defined($FATHMM)) {
		$FATHMM = 'NA';
	    }
	    unless (defined($PhyloP)) {
		$PhyloP = 'NA';
	    }
	    unless (defined($ESP_AA)) {
		$ESP_AA = 'NA';
	    }
	    unless (defined($ESP_EA)) {
		$ESP_EA = 'NA';
	    }

	    $sth->finish();

	    ###$req = "select Pathway, Function_description, Disease_description, MIM_phenotype_id, MIM_disease, Trait_association_GWAS, Expression_epigenetics, Expression_GNF_Atlas, Interactions_IntAct, Interactions_BioGRID from dbNSFP_v2_0_genes where Gene_name = \"$genename\"";
	    $req = "select Pathway_Uniprot, Function_description, Disease_description, MIM_phenotype_id, MIM_disease, Trait_association_GWAS, Expression_epigenetics, Expression_GNF_Atlas, Interactions_IntAct, Interactions_BioGRID from dbNSFP_v3_0b2a_genes where Gene_name = \"$genename\"";

	    $sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	    $sth->execute()||  print ("\nExecution requete: $req echec\n");

	    if ($sth->rows > 0) {
		# There is at least one line, the mutation is found in the dbNSFP_genes table
		my ($pathway, $function, $disease, $MIM_id, $MIM_disease, $trait_assoc, $expr_egenetics, $expr_GNF, $intact, $biogrid);
		while (my @gene_info = $sth->fetchrow_array) {
		    if ($gene_info[0] ne '.') {
			$pathway .= $gene_info[0];
		    }
		    if ($gene_info[1] ne '.') {
			$function .= $gene_info[1];
		    }
		    if ($gene_info[2] ne '.') {
			$disease .= $gene_info[2];
		    }
		    if ($gene_info[3] ne '.') {
			$MIM_id .= $gene_info[3];
		    }
		    if ($gene_info[4] ne '.') {
			$MIM_disease .= $gene_info[4];
		    }
		    if ($gene_info[5] ne '.') {
			$trait_assoc .= $gene_info[5];
		    }
		    if ($gene_info[6] ne '.') {
			$expr_egenetics .= $gene_info[6];
		    }
		    if ($gene_info[7] ne '.') {
			$expr_GNF .= $gene_info[7];
		    }
		    if ($gene_info[8] ne '.') {
			$intact .= $gene_info[8];
		    }
		    if ($gene_info[9] ne '.') {
			$biogrid .= $gene_info[9];
		    }
		}

		unless (defined($pathway)) {
		    $pathway = 'NA';
		}
		unless (defined($function)) {
		    $function = 'NA';
		}
		unless (defined($disease)) {
		    $disease = 'NA';
		}
		unless (defined($MIM_id)) {
		    $MIM_id = 'NA';
		}
		unless (defined($MIM_disease)) {
		    $MIM_disease = 'NA';
		}
		unless (defined($trait_assoc)) {
		    $trait_assoc = 'NA';
		}
		unless (defined($expr_egenetics)) {
		    $expr_egenetics = 'NA';
		}
		unless (defined($expr_GNF)) {
		    $expr_GNF = 'NA';
		}
		unless (defined($intact)) {
		    $intact = 'NA';
		}
		unless (defined($biogrid)) {
		    $biogrid = 'NA';
		}

		$sth->finish();

		print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$anc_all\t$pathway\t$function\t$disease\t$MIM_id\t$MIM_disease\t$trait_assoc\t$expr_egenetics\t$expr_GNF\t$intact\t$biogrid\n";
		$n_out++;
	    } else {
		print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$anc_all\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
		$n_out++;
	    }
	} else {
	    # Unknown mutation => NA
	    print FHO "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
    }

#    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) found in dbNSFP in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

sub dbnsfp_annot_snp_mongo {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $line, $req, $sth);
    my ($l, $header, @tab, $alleles, @individus, $individu, $search_annot, $search_gene, $count_annot, $count_gene, $genename, $alt_allele, $search, $annot, $next, $genes, $next2, $mim, $chrom, $pos);
    my ($PPh2_hdiv, $PPh2_hvar, $LRT, $MT, $SIFT, $G1000_AF, $G1000_AFR, $G1000_EU, $G1000_AM, $G1000_EAS, $G1000_SAS, $anc_all, $MutAss, $FATHMM, $PhyloP, $ESP_AA, $ESP_EA, $ExAC_AF, $ExAC_AFR_AF, $ExAC_AMR_AF, $ExAC_NFE_AF, $ExAC_EAS_AF, $ExAC_SAS_AF);
    my ($pathway, $function, $disease, $MIM_id, $MIM_disease, $trait_assoc, $expr_egenetics, $expr_GNF, $intact, $biogrid);

    my $client = MongoDB->connect('mongodb-calcul.egid.local:27017');
    #my $client = MongoDB->connect('10.220.223.138:27017');
    my $db   = $client->get_database( 'dbNSFPv3_4a' );
    my $var = $db->get_collection('dbNSFP_variant');
    my $gen = $db->get_collection('dbNSFP_gene');
    my $snv = $db->get_collection('dbscSNV');

    my @term = ("chr","pos_1-based","alt");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	chomp $line;
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^#/) {
	    $n_in--;
	    print FHO $line."\tPPh2_Humdiv\tPPh2_Humvar\tSIFT\tMutationTaster\tLRT\tMutationAssessor\tFATHMM\tPhyloP\tFreq_1KG\tFreq_African_1KG\tFreq_European_1KG\tFreq_American_1KG\tFreq_East_Asian_1KG\tFrequence_South_Asian_1KG\tFreq_African_American_ESP6500\tFreq_European_American_ESP6500\tExAC_AF\tExAC_AFR_AF\tExAC_NFE_AF\tExAC_AMR_AF\tExAC_EAS_AF\tExAC_SAS_AF\tAncestral_allele\tPathway\tFunction\tDisease\tMIM_id\tMIM_disease\tTrait_association\tExpression_egenetics\tExpression_GNF_Atlas\tInteractions_IntAct\tInteractions_BioGRID\n";
	    next;
	}

	@tab = split(/\t/, $line);
	$chrom = $tab[0];
	$pos = $tab[1];


	my $alleles = '';
	if ($fam == 0) {
	    $alleles = $tab[3].$tab[4];
	} else {
 print "individus : $tab[16] !!!";
	    my @individus = split(/\] /, $tab[16]);
	    foreach my $individu(@individus) {
		$individu =~ s/.*\///; # mutated allele is always the second for heterozygotes
		$individu =~ s/,.*//;
		if ($individu ne $alleles) {
		    $alleles .= $individu;
  print "alleles : $alleles !!!";
		}
	    }
	}

	$search_annot = {
	    $term[0]=>$tab[0],
		$term[1]=>$tab[1],
		$term[2]=>$tab[4]
	};

	$count_annot = $var->count($search_annot);
	if ($count_annot == 0){
	    print FHO "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	    $n_out++;
	} else {
	    $annot = $var->find($search_annot)->fields({"chr"=>1,"hg19_chr"=>1,"pos_1-based"=>1,"hg19_pos_1-based"=>1,"_id"=>0,"alt"=>1,"ref"=>1,"genename"=>1,"Polyphen2_HDIV_pred"=>1,"Polyphen2_HVAR_pred"=>1,"LRT_pred"=>1,"MutationTaster_pred"=>1,"SIFT_pred"=>1,"1000Gp3_AF"=>1,"1000Gp3_AFR_AF"=>1,"1000Gp3_EUR_AF"=>1,"1000Gp3_AMR_AF"=>1,"1000Gp3_EAS_AF"=>1,"1000Gp3_SAS_AF"=>1,"Ancestral_allele"=>1,"MutationAssessor_pred"=>1,"FATHMM_pred"=>1,"phyloP100way_vertebrate"=>1,"ESP6500_AA_AF"=>1,"ESP6500_EA_AF"=>1,"ExAC_AF"=>1,"ExAC_AFR_AF"=>1,"ExAC_NFE_AF"=>1,"ExAC_AMR_AF"=>1,"ExAC_EAS_AF"=>1,"ExAC_SAS_AF"=>1}); ### add columns
	    while ($next = $annot->next) {
		$alt_allele = $next->{"alt"};
		$genename = $next->{"genename"};
		$PPh2_hdiv = $next->{"Polyphen2_HDIV_pred"};
		$PPh2_hvar = $next->{"Polyphen2_HVAR_pred"};
		$LRT = $next->{"LRT_pred"};
		$MT = $next->{"MutationTaster_pred"};
		$SIFT = $next->{"SIFT_pred"};
		$G1000_AF = $next->{"1000Gp3_AF"};
		$G1000_AFR = $next->{"1000Gp3_AFR_AF"};
		$G1000_EU = $next->{"1000Gp3_EUR_AF"};
		$G1000_AM = $next->{"1000Gp3_AMR_AF"};
		$G1000_EAS = $next->{"1000Gp3_EAS_AF"};
		$G1000_SAS = $next->{"1000Gp3_SAS_AF"};
		$anc_all = $next->{"Ancestral_allele"};
		$MutAss = $next->{"MutationAssessor_pred"};
		$FATHMM = $next->{"FATHMM_pred"};
		$PhyloP = $next->{"phyloP100way_vertebrate"};
		$ESP_AA = $next->{"ESP6500_AA_AF"};
		$ESP_EA = $next->{"ESP6500_EA_AF"};
		$ExAC_AF = $next->{"ExAC_AF"};
		$ExAC_AFR_AF = $next->{"ExAC_AFR_AF"};
		$ExAC_NFE_AF = $next->{"ExAC_NFE_AF"};
		$ExAC_AMR_AF = $next->{"ExAC_AMR_AF"};
		$ExAC_EAS_AF = $next->{"ExAC_EAS_AF"};
		$ExAC_SAS_AF = $next->{"ExAC_SAS_AF"};

		$search_gene = {
		    "Gene_name"=>$genename,
		};
		if ($alleles =~ /$alt_allele/){
		    $count_gene = $gen->count($search_gene);
		    if($count_gene == 0){
			print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$ExAC_AF\t$ExAC_AFR_AF\t$ExAC_NFE_AF\t$ExAC_AMR_AF\t$ExAC_EAS_AF\t$ExAC_SAS_AF\t$anc_all\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
			$n_out++;
		    } else {
			$genes = $gen->find($search_gene)->fields({"_id"=>0,"Pathway_Uniprot"=>1,"Function_description"=>1,"Disease_description"=>1,"MIM_phenotype_id"=>1,"MIM_disease"=>1,"Trait_association_GWAS"=>1,"Expression_egenetics"=>1,"Expression_GNF-Atlas"=>1,"Interactions_IntAct"=>1,"Interactions_BioGRID"=>1});
			while ($next2 = $genes->next) {
			    $pathway = $next2->{"Pathway_Uniprot"};
			    $function = $next2->{"Function_description"};
			    $disease = $next2->{"Disease_description"};
			    $MIM_id = $next2->{"MIM_phenotype_id"};
			    $MIM_disease = $next2->{"MIM_disease"};
			    $trait_assoc = $next2->{"Trait_association_GWAS"};
			    $expr_egenetics = $next2->{"Expression_egenetics"};
			    $expr_GNF = $next2->{"Expression_GNF-Atlas"};
			    $intact = $next2->{"Interactions_IntAct"};
			    $biogrid = $next2->{"Interactions_BioGRID"};
			    print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$ExAC_AF\t$ExAC_AFR_AF\t$ExAC_NFE_AF\t$ExAC_AMR_AF\t$ExAC_EAS_AF\t$ExAC_SAS_AF\t$anc_all\t$pathway\t$function\t$disease\t$MIM_id\t$MIM_disease\t$trait_assoc\t$expr_egenetics\t$expr_GNF\t$intact\t$biogrid\n";
			    $n_out++;
			}
		    }
		} else {
		    #print $alt_allele."\t".$alleles."\t*** ".$l."\n";
		}
	    }
	}
    }

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) found in dbNSFP in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    #    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}


#####################################################################################################################
sub dbnsfp_annot_snp_mongo_fam {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $line, $req, $sth);
    my ($l, $header, @tab, $alleles, @individus, $individu, $search_annot, $search_gene, $count_annot, $count_gene, $genename, $alt_allele, $search, $annot, $next, $genes, $next2, $mim, $chrom, $pos);
    my ($PPh2_hdiv, $PPh2_hvar, $LRT, $MT, $SIFT, $G1000_AF, $G1000_AFR, $G1000_EU, $G1000_AM, $G1000_EAS, $G1000_SAS, $anc_all, $MutAss, $FATHMM, $PhyloP, $ESP_AA, $ESP_EA, $ExAC_AF, $ExAC_AFR_AF, $ExAC_AMR_AF, $ExAC_NFE_AF, $ExAC_EAS_AF, $ExAC_SAS_AF);
    my ($pathway, $function, $disease, $MIM_id, $MIM_disease, $trait_assoc, $expr_egenetics, $expr_GNF, $intact, $biogrid);

    my $client = MongoDB->connect('XXX:XXX');
    my $db   = $client->get_database( 'XXX' );
    my $var = $db->get_collection('dbNSFP_variant');
    my $gen = $db->get_collection('dbNSFP_gene');
    my $snv = $db->get_collection('dbscSNV');

    my @term = ("chr","pos_1-based","alt");

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	chomp $line;
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^#/) {
	    $n_in--;
	    print FHO $line."\tPPh2_Humdiv\tPPh2_Humvar\tSIFT\tMutationTaster\tLRT\tMutationAssessor\tFATHMM\tPhyloP\tFreq_1KG\tFreq_African_1KG\tFreq_European_1KG\tFreq_American_1KG\tFreq_East_Asian_1KG\tFrequence_South_Asian_1KG\tFreq_African_American_ESP6500\tFreq_European_American_ESP6500\tExAC_AF\tExAC_AFR_AF\tExAC_NFE_AF\tExAC_AMR_AF\tExAC_EAS_AF\tExAC_SAS_AF\tAncestral_allele\tPathway\tFunction\tDisease\tMIM_id\tMIM_disease\tTrait_association\tExpression_egenetics\tExpression_GNF_Atlas\tInteractions_IntAct\tInteractions_BioGRID\n";
	    next;
	}

	@tab = split(/\t/, $line);
	$chrom = $tab[0];
	$pos = $tab[1];


	my $alleles = '';
	if ($fam == 0) {
	    $alleles = $tab[3].$tab[4];
	} else {
 #modif mehdi 03/10/18 :
 print "individus : $tab[16] !!!";
	    my @individus = split(/\] /, $tab[16]);
 #
	    foreach my $individu(@individus) {
		$individu =~ s/.*\///; # mutated allele is always the second for heterozygotes
		$individu =~ s/,.*//;
		if ($individu ne $alleles) {
		    $alleles .= $individu;
  print "alleles : $alleles !!!";
		}
	    }
	}

	$search_annot = {
	    $term[0]=>$tab[0],
		$term[1]=>$tab[1],
		$term[2]=>$alleles
	};

	$count_annot = $var->count($search_annot);
	if ($count_annot == 0){
	    print FHO "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	    $n_out++;
	} else {
	    $annot = $var->find($search_annot)->fields({"chr"=>1,"hg19_chr"=>1,"pos_1-based"=>1,"hg19_pos_1-based"=>1,"_id"=>0,"alt"=>1,"ref"=>1,"genename"=>1,"Polyphen2_HDIV_pred"=>1,"Polyphen2_HVAR_pred"=>1,"LRT_pred"=>1,"MutationTaster_pred"=>1,"SIFT_pred"=>1,"1000Gp3_AF"=>1,"1000Gp3_AFR_AF"=>1,"1000Gp3_EUR_AF"=>1,"1000Gp3_AMR_AF"=>1,"1000Gp3_EAS_AF"=>1,"1000Gp3_SAS_AF"=>1,"Ancestral_allele"=>1,"MutationAssessor_pred"=>1,"FATHMM_pred"=>1,"phyloP100way_vertebrate"=>1,"ESP6500_AA_AF"=>1,"ESP6500_EA_AF"=>1,"ExAC_AF"=>1,"ExAC_AFR_AF"=>1,"ExAC_NFE_AF"=>1,"ExAC_AMR_AF"=>1,"ExAC_EAS_AF"=>1,"ExAC_SAS_AF"=>1}); ### add columns
	    while ($next = $annot->next) {
		$alt_allele = $next->{"alt"};
		$genename = $next->{"genename"};
		$PPh2_hdiv = $next->{"Polyphen2_HDIV_pred"};
		$PPh2_hvar = $next->{"Polyphen2_HVAR_pred"};
		$LRT = $next->{"LRT_pred"};
		$MT = $next->{"MutationTaster_pred"};
		$SIFT = $next->{"SIFT_pred"};
		$G1000_AF = $next->{"1000Gp3_AF"};
		$G1000_AFR = $next->{"1000Gp3_AFR_AF"};
		$G1000_EU = $next->{"1000Gp3_EUR_AF"};
		$G1000_AM = $next->{"1000Gp3_AMR_AF"};
		$G1000_EAS = $next->{"1000Gp3_EAS_AF"};
		$G1000_SAS = $next->{"1000Gp3_SAS_AF"};
		$anc_all = $next->{"Ancestral_allele"};
		$MutAss = $next->{"MutationAssessor_pred"};
		$FATHMM = $next->{"FATHMM_pred"};
		$PhyloP = $next->{"phyloP100way_vertebrate"};
		$ESP_AA = $next->{"ESP6500_AA_AF"};
		$ESP_EA = $next->{"ESP6500_EA_AF"};
		$ExAC_AF = $next->{"ExAC_AF"};
		$ExAC_AFR_AF = $next->{"ExAC_AFR_AF"};
		$ExAC_NFE_AF = $next->{"ExAC_NFE_AF"};
		$ExAC_AMR_AF = $next->{"ExAC_AMR_AF"};
		$ExAC_EAS_AF = $next->{"ExAC_EAS_AF"};
		$ExAC_SAS_AF = $next->{"ExAC_SAS_AF"};

		$search_gene = {
		    "Gene_name"=>$genename,
		};
		if ($alleles =~ /$alt_allele/){
		    $count_gene = $gen->count($search_gene);
		    if($count_gene == 0){
			print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$ExAC_AF\t$ExAC_AFR_AF\t$ExAC_NFE_AF\t$ExAC_AMR_AF\t$ExAC_EAS_AF\t$ExAC_SAS_AF\t$anc_all\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
			$n_out++;
		    } else {
			$genes = $gen->find($search_gene)->fields({"_id"=>0,"Pathway_Uniprot"=>1,"Function_description"=>1,"Disease_description"=>1,"MIM_phenotype_id"=>1,"MIM_disease"=>1,"Trait_association_GWAS"=>1,"Expression_egenetics"=>1,"Expression_GNF-Atlas"=>1,"Interactions_IntAct"=>1,"Interactions_BioGRID"=>1});
			while ($next2 = $genes->next) {
			    $pathway = $next2->{"Pathway_Uniprot"};
			    $function = $next2->{"Function_description"};
			    $disease = $next2->{"Disease_description"};
			    $MIM_id = $next2->{"MIM_phenotype_id"};
			    $MIM_disease = $next2->{"MIM_disease"};
			    $trait_assoc = $next2->{"Trait_association_GWAS"};
			    $expr_egenetics = $next2->{"Expression_egenetics"};
			    $expr_GNF = $next2->{"Expression_GNF-Atlas"};
			    $intact = $next2->{"Interactions_IntAct"};
			    $biogrid = $next2->{"Interactions_BioGRID"};
			    print FHO "$line\t$PPh2_hdiv\t$PPh2_hvar\t$SIFT\t$MT\t$LRT\t$MutAss\t$FATHMM\t$PhyloP\t$G1000_AF\t$G1000_AFR\t$G1000_EU\t$G1000_AM\t$G1000_EAS\t$G1000_SAS\t$ESP_AA\t$ESP_EA\t$ExAC_AF\t$ExAC_AFR_AF\t$ExAC_NFE_AF\t$ExAC_AMR_AF\t$ExAC_EAS_AF\t$ExAC_SAS_AF\t$anc_all\t$pathway\t$function\t$disease\t$MIM_id\t$MIM_disease\t$trait_assoc\t$expr_egenetics\t$expr_GNF\t$intact\t$biogrid\n";
			    $n_out++;
			}
		    }
		} else {
		    #print $alt_allele."\t".$alleles."\t*** ".$l."\n";
		}
	    }
	}
    }

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) found in dbNSFP in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    #    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}


##################################################################################
### Prepare for indiv_annot_snp
sub launch_indiv_annot_snp {
    my $run_name = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;
    my $count = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_mut_indiv';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");
    system("mkdir $result_folder/SNP_filter/${to}_combi") unless (-d "$result_folder/SNP_filter/${to}_combi");

print "indiv_annot snp from $from LAUNCHED\n";

    if ($fam == 0){
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotSNP::indiv_annot_snp_redis($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
      print "indiv_annot snp for FAM LAUNCHED\n";
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotSNP::coseg_file_snp($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotSNP::indiv_annot_snp_redis($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotSNP::indiv_annot_snp_redis($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "indiv_annot snp from $from DONE\n";
}

##################################################################################
### Annotate with the list of samples where the mutation is found
sub indiv_annot_snp {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;
    my $count = shift;

print "indiv_annot snp normal lancement\n";

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, $start, $end, $req, $sth, @row, @tab, $req2, $sth2, $count_atteint, $count_control, %tmp);

    my $dbh = DBconnect::connect_db("SNP_annot");
    my $dbh2 = DBconnect::connect_db("RUNS");

    # Load indiv and their status
    $req2 = "select name_sample,patient_control from SAMPLE";
    $sth2 = $dbh2->prepare($req2)||  print ("\nPreparation requete $req2: echec\n");
    $sth2->execute()||  print ("\nExecution requete: $req2 echec\n");

    if ($sth2->rows > 0) {
	while (my @status = $sth2->fetchrow_array) {
	    $tmp{$status[0]} = $status[1];
	}
    }
    $sth2->finish();

    $dbh2->disconnect();

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	$count_atteint = 0;
	$count_control = 0;
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    if ($count == 1) {
		print FHO "$line\tN_atteint\tN_control\n";
	    } else {
		print FHO "$line\tIndiv\n";
	    }
	    next;
	}
	@tab = split(/\t/, $line);

	# Check for mutation in mut_indiv table
	$req = "select Individu,ARA1,ARA2,Score from mut_indiv where " .
	    " Chrom = \"$tab[0]\" " .
	    " and Pos = \"$tab[1]\" ";
	$sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	$sth->execute()||  print ("\nExecution requete: $req echec\n");

	if ($sth->rows > 0) {
	    # There is at least one line so the mutation is found in at least one of our samples
	    my ($individu, $ara1, $ara2, $score);
	    my $list = '';
	    while (my @individus = $sth->fetchrow_array) {
		$individu = $individus[0];
#		next if (($individu =~ /hg19/) && ($assembly eq 'hg18'));
#		next if (($individu !~ /hg19/) && ($assembly eq 'hg19'));
		next if ($individu !~ /hg19/);
		next if ($individu =~ /Sample/); # sauf exception
		next if ($individu =~ /_S/);
		next if ($individu =~ /test/);
		next if ($individu =~ /uniq/);
		next if ($individu =~ /Next/);
		$ara1 = substr($individus[1],2,1);
		$ara2 = substr($individus[2],2,1);
		$score = $individus[3];
		$list .= "$individu($ara1/$ara2,$score),";
		if ($count == 1) {
		    # Check status of indiv
		    my $indiv = $individu;
		    $indiv =~s/_hg19//;
#		    $req2 = "select patient_control from SAMPLE where name_sample = \"$indiv\"";
#		    $sth2 = $dbh2->prepare($req2)||  print ("\nPreparation requete $req2: echec\n");
#		    $sth2->execute()||  print ("\nExecution requete: $req2 echec\n");

#		    if ($sth2->rows > 0) {
#			my $count_atteint_sample = 0;
#			my $count_control_sample = 0;
#			while (my @status = $sth2->fetchrow_array) {
		    if (defined $tmp{$indiv}) {
			    if ($tmp{$indiv} == 1) { # est un controle si patient_control = 1
				$count_control++;
			    } else {
				$count_atteint++;
			    }
		    }
#			}
#			if ($count_atteint_sample ne 0) {
#			    $count_atteint++;
#			}
#			if ($count_control_sample ne 0) {
#			    $count_control++;
#			}
#			$sth2->finish();
#		    }
		}
	    }

	    $sth->finish();
	    $list =~ s/,$//;
#	    print "BANG ! TOT=$n_rows_tot NE $count_atteint + $count_control" if (($count_atteint + $count_control) ne $n_rows_tot);  # table SAMPLE existe depuis moins longtemps que mut_indiv !
	    if ($count == 1) {
		print FHO "$line\t$count_atteint\t$count_control\n";
	    } else {
		print FHO "$line\t$list\n";
	    }
	} else {
	    # Unknown mutation, set to 0
	    if ($count == 1) {
		print FHO "$line\t0\t0\n";
	    } else {
		print FHO "$line\t0\n";
	    }
	    $n_out++;
	}
    }

#    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
    $dbh->disconnect();
    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) not found in our samples in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}


##################################################################################
### Annotate with the number of samples (ill and control) where the mutation is found
sub indiv_annot_snp_redis {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;

print "indiv_annot snp redis lancement\n";

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, @tab, $count_atteint, $count_control, $count_atteint_hg38, $count_control_hg38, $chr, $pos, $list_name);

    my $redis = Redis->new(
	server => 'XXX:XXX',
	);

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
	LogOut::printLog($log, $current_sub, $logFolder, $logFile);
	return;
    }

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;

    $log = "";
    while (my $line = <FHI>) {
	$count_atteint = 0;
	$count_control = 0;
	$count_atteint_hg38 = 0;
	$count_control_hg38 = 0;
	$n_in++;
	$log .= "$n_in $fileIn\n" if (($n_in % 10000) == 0);
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line =~ /^\#/) {
	    $n_in--;
	    print FHO "$line\tN_atteint_hg19\tN_control_hg19\tN_atteint_hg38\tN_control_hg38\n";
	    next;
	}

	@tab = split(/\t/, $line);

        $chr = $tab[0];
        $pos = $tab[18];

	$list_name = "$chr"."-"."$pos"."-cas";
	$count_atteint = $redis->scard( $list_name );

	$list_name = "$chr"."-"."$pos"."-control";
        $count_control = $redis->scard( $list_name );

	print FHO "$line\t$count_atteint\t$count_control\t";

        $pos = $tab[1];

	$list_name = "$chr"."-"."$pos"."-cas-hg38";
	$count_atteint_hg38 = $redis->scard( $list_name );

	$list_name = "$chr"."-"."$pos"."-control-hg38";
        $count_control_hg38 = $redis->scard( $list_name );

	print FHO "$count_atteint_hg38\t$count_control_hg38\n";

	if (($count_atteint == 0) && ($count_control == 0) && ($count_atteint_hg38 == 0) && ($count_control_hg38 == 0)) {
	    $n_out++;
	}
    }

    close(FHI);
    close(FHO);

    $log = "$n_out mutation(s) not found in our samples in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}


##################################################################################
### Generate combinations of samples and associated file names
sub coseg_file_snp {
    my $analysis_data = shift;
    my $combinaison = shift;
    my $result_folder = shift;
    my $from = shift;
    my $to = shift;
    my $recess = shift;

    my $from_short = "F_coding";
    $from_short =~ s/F_//;

    my (@filesIn, @filesOut, @case_bag, @control_bag);

    # On parcourt tous les samples to get file names
    my @dna_cases;
    my @dna_controls;
    my @parent_dna_controls;
    my @other_dna_controls;
    my $family;
    my $sickness;
    foreach my $k (keys %$analysis_data) {
	if ($$analysis_data{$k}{'atteint'} == 1) {
	    push @dna_cases, $$analysis_data{$k}{'result_dir'};
	}
	if ($$analysis_data{$k}{'atteint'} == 0) {
	    push @dna_controls, $$analysis_data{$k}{'result_dir'};
	    if (($$analysis_data{$k}{'parent'} eq 'M') || ($$analysis_data{$k}{'parent'} eq 'F')) {
		push @parent_dna_controls, $$analysis_data{$k}{'result_dir'};
	    } else {
		push @other_dna_controls, $$analysis_data{$k}{'result_dir'};
	    }
	}
	$family = $$analysis_data{$k}{'pop'};
	$sickness = $$analysis_data{$k}{'maladie'};
    }
    my $cases = join('_',  sort @dna_cases);
    my ($controls);

    my $fileOut;
    my $fileIn;

    my $suffix;
    if ($recess == 1) {
	$suffix = "_recess";
    } else {
	$suffix = "";
    }

    #Cas normal : cas versus controles : ecrit dans "cover_other"
    if (@parent_dna_controls && @other_dna_controls) {
	$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @other_dna_controls);
    } elsif (@parent_dna_controls && !(@other_dna_controls)) {
	$controls = join('_', sort @parent_dna_controls);
    } elsif (!(@parent_dna_controls) && @other_dna_controls) {
	$controls = join('_', sort @other_dna_controls);
#	if ($recess == 1) {
#	    die "What the fuck : recess analysis without parents !!!";
#	}
#    } elsif (!(@parent_dna_controls) && !(@other_dna_controls)) {
#	if ($recess == 1) {
#	    die "What the fuck : recess analysis without controls !!!";
#	}
    }

    if($controls) {
S
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	$fileIn = $result_folder."/SNP_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
    } else {
	$fileOut = $result_folder."/SNP_filter/${to}/${to}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	$fileIn = $result_folder."/SNP_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
    }
    push @filesIn, $fileIn;
    push @filesOut, $fileOut;
    push @case_bag, 'none';
    push @control_bag, \@dna_controls;

    # Autre cas qui teste toutes les autres combinaisons
    if ($combinaison == 1) {
    	my @bitMask = ();
    	#On cherche tous les subsets possibles pour la liste d'adn
    	while (FilterSNP::generate_mask(\@bitMask, \@dna_cases)) {
    	    my $d_cases = FilterSNP::get_one_subset(\@bitMask, \@dna_cases);

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
    		my @cases_as_controls = @d_controls;
    		#On ajoute nos "vrais" controles
    		# a cette liste de controles
    		if(@dna_controls) {
    		    push @d_controls, @dna_controls;
    		}
    		my $ca = join('_', sort @$d_cases);
    		my $co = join('_', sort @d_controls);
		$fileOut = $result_folder."/SNP_filter/${to}_combi/${to}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		$fileIn = $result_folder."/SNP_filter/${from}_combi/${from}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		push @filesIn, $fileIn;
		push @filesOut, $fileOut;
		push @case_bag, \@cases_as_controls;
		push @control_bag, \@dna_controls;
    	    }
    	}
    }
    return (\@filesIn, \@filesOut, \@case_bag, \@control_bag);
}

##################################################################################
###
sub launch_base_count_snp {
    my $run_name = shift;
    my $run_disk = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($chr, $vcf_file, $log, $fileIn, $fileOut, $log_folder, $log_file, $position, @line, @line_vcf, @count, @position, $position_vcf, $count_format, $ref_alt, %allele_count, %f_coding, $head);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_basecount';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");

    $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/Sample_".$sample_data->{dna}."/".$sample_data->{dna}."_snps.vcf";
    $fileOut = $result_folder."/SNP_filter/".$to."/F_basecount_".$sample_data->{result_dir}.".txt";
    $fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
	return;
    }
    if (! -f $vcf_file) {
	$log = scalar(localtime())."  !!!missing file $vcf_file\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
#	print $log;
	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
	return;
    }

    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    # recuperer comptage bases dans les fichiers vcf
    $log = "open $vcf_file  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
    open (VCF, $vcf_file) or die("cannot open $vcf_file");

    while (my $l=<VCF>){
	chomp $l;
	next if $l=~/^#/;
	@line_vcf = split("\t",$l);
	$chr = $line_vcf[0];
	$chr =~ s/chr//;
	$position_vcf = $chr.":".$line_vcf[1];
	next if $line_vcf[9]=~/^.\/.$/;
	@count = split(":",$line_vcf[9]);
	$count_format= $count[1];
	$count_format=~s/,/\//;
	$ref_alt = $line_vcf[3]."/".$line_vcf[4];
	$allele_count{$position_vcf} = $ref_alt.":".$count_format;

    }
    close VCF;

    # recuperer lignes dans F_mut_indiv
    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    while (my $l=<FHI>){
	chomp $l;
	if ($l=~/^#/){
	    $head = $l;
	    print FHO "$head\tgenotype_count\n";
	} else {
	    @line = split("\t", $l);
	    $position = $line[0].":".$line[1];
	    if ($allele_count{$position}){
		print FHO "$l\t$allele_count{$position}\n";
	    } else {
		print FHO "$l\tNA\n";
	    }
	}
    }
    close FHI;
    close FHO;
    system("gzip $fileOut");
}


##################################################################################
###
sub launch_disease_inher_snp {
    my $run_name = shift;
    my $run_disk = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $fileIn, $fileOut, $log_folder, $log_file, @tab, $gene, $n_in, $n_out);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

    my $to = 'F_disease_inher';

    system("mdkir $result_folder/SNP_filter") unless (-d "$result_folder/SNP_filter");
    system("mkdir $result_folder/SNP_filter/${to}") unless (-d "$result_folder/SNP_filter/${to}");

    $fileOut = $result_folder."/SNP_filter/".$to."/F_disease_inher_".$sample_data->{result_dir}.".txt";
    $fileIn = $result_folder."/SNP_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

    my $fileIn_gz = $fileIn . ".gz";
    if ((-e $fileIn_gz) && (! -e $fileIn)) {
	system("gunzip $fileIn_gz");
    }

    if (! -f $fileIn) {
	$log = scalar(localtime())."  !!!missing file $fileIn\n";
	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
	return;
    }

    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists, replacing it\n";
	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
#	return;
    }

    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    # Get diseases and inheritance from file
    my $gene_file = "/media/Data/genes_raindance_3.tab";
    open (FHG, $gene_file) or die ("cannot open $gene_file");
    my (%gene_hash, @gene_tab);
    while (my $l = <FHG>) {
	$l =~ s/\012//;
	$l =~ s/\015//;
	@gene_tab = split(/\t/, $l);
	$gene_hash{$gene_tab[0]} = [$gene_tab[1],$gene_tab[2]];
    }

    # Get lines from previous step and add new annot
    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
    $log = "create $fileOut [dna]=$sample_data->{dna}\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
    open (FHI, $fileIn) or die("cannot open $fileIn");

    $n_in = 0;
    $n_out = 0;

    while (my $line=<FHI>){
	$n_in++;
	$line =~ s/\012//;
	$line =~ s/\015//;
	if ($line=~/^#/){
	    $n_in--;
	    print FHO "$line\tDisease\tInheritance\n";
	    next;
	}
	@tab = split("\t", $line);
	$gene = $tab[10];
	if (exists($gene_hash{$gene})) {
	    print FHO "$line\t$gene_hash{$gene}[0]\t$gene_hash{$gene}[1]\n";
	} else {
	    # Unlisted gene, set to NA
	    print FHO "$line\tNA\tNA\n";
	    $n_out++;
	}
    }
    close FHI;
    close FHO;

    unless (-e $fileIn_gz) {
	system("gzip $fileIn");
    }
    system("gzip $fileOut");

    $log = "$n_out mutation(s) not found in raindance 3 gene panel in $fileOut, $n_in read from $fileIn\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
}


1;
