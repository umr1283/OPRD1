package AnnotINDEL;

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
### From a tab delimited input file containing chromosome, start, end, CIGAR, allele, etc.
### this module retrieves from Ensembl Local ID, Annotated ID, Local alleles,
### Annotated alleles, Chromosome, Start, End, Consequence, Transcript, AA change,
### Gene ID, and Gene name and writes it in a db table
sub ensembl_annot_indel {
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

    my $F_score_folder = $result_folder."/INDEL_filter/F_score";
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
	my ($chromosome, $start, $CIGAR, $ref_allele, $allele1, $allele2, $hh, $depth, $score) = split (/\t/, $line);

	my $allele = $ref_allele.'/'.$allele2;
	my $id = 'indel_'.$chromosome.'_'.$start.'_'.$CIGAR;
	my $size = length ($allele2);
	my $type = substr $CIGAR,-1;
	my $end;
	if (($CIGAR eq 'BP_RIGHT') || ($CIGAR eq 'BP_LEFT')) {
	    $end = $start + $size;
	} else {
	    if ($type eq 'D') {
		$end = $start + $size -1;
	    } else {
		$end = $start + $size;
	    }
	}

	my $tiny_slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome, $start, $end);
	my $slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome);

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
    close IN;
    $dbh->disconnect();
    print scalar(localtime())." : ensembl INDEL annot from $allele_file done\n";
#   LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
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

    if ($seq_region_start > $seq_region_end) {
	my $tmp_start = $seq_region_start;
	$seq_region_start = $seq_region_end;
	$seq_region_end = $tmp_start;
    }

    my $consequence_type = $variation_feature->consequence_type();

    if ($consequence_type =~ /ARRAY/) {
	$consequence = join(",",@$consequence_type);
    } else {
	$consequence = $consequence_type;
    }

#   if ($consequence =~ 'INTERGENIC') {
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

	# Insert into our mutations database
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
		$refseq_mrna = "NA";
		$refseq_peptide = "NA";
                my @dbentries = @{$transcript->get_all_DBLinks('RefSeq%')};
                foreach my $dbentry(@dbentries) {
                    if ($dbentry->dbname() eq 'RefSeq_mRNA') {
                        $refseq_mrna = $dbentry->display_id;
                    } elsif ($dbentry->dbname() eq 'RefSeq_peptide') {
                        $refseq_peptide = $dbentry->display_id;
                    }
                }

		# Insert into our mutations database
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
sub launch_hg19_pos_annot_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    ### Ajout mehdi 15/01/18 pour ajouter colonne pos_hg19 dans le script "start_analysis"


    if ($fam == 0){
    $fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
    $fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

    AnnotINDEL::hg19_pos_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
  my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
  my @filesIn = @$filesIn_ref;
  my @filesOut = @$filesOut_ref;
  if ($combinaison == 0){
    $fileIn = $filesIn[0];
    $fileOut = $filesOut[0];
    AnnotINDEL::hg19_pos_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
  } elsif ($combinaison == 1) {
      for (my $i = 0; $i < scalar @filesIn; $i++) {
    AnnotINDEL::hg19_pos_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
    }
  }
  }

    AnnotINDEL::hg19_pos_annot_indel($log_folder, $log_file, $sample_data, $fileIn, $fileOut);


    print "hg19 annot indel from $from DONE\n";

}

##################################################################################
### Annotate with hg19 position
sub hg19_pos_annot_indel {
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
	    my $status = $response->{'status'};
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
	    ###

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
    #$dbh->disconnect();
    close(FHI);
    close(FHO);
    $log = "$n_out mutation(s) not found in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
### Annotate with control covers when available (otherwise 0)
sub launch_cover_other_indel {
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

    my $cover_other_folder = $result_folder."/INDEL_filter/F_cover_other";
    system("mkdir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/F_cover_other") unless (-d "$result_folder/INDEL_filter/F_cover_other");
    system("mkdir $result_folder/INDEL_filter/F_cover_other_combi") unless (-d "$result_folder/INDEL_filter/F_cover_other_combi");

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


    if (@dna_controls) {
	# if ($recess == 1) {
	#     my @non_parent_dna_controls = @dna_controls;
	#     splice(@non_parent_dna_controls, 0, 2);
	#     my @parent_dna_controls = ($dna_controls[0], $dna_controls[1]);
	#     if (@non_parent_dna_controls) {
	# 	$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @non_parent_dna_controls);
	#     } else {
	# 	$controls = join('_', sort @parent_dna_controls);
	#     }
	# } else {
	#     $controls = join('_',  sort @dna_controls);
	# }
	$fileOut = $cover_other_folder."/F_cover_other_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	$fileIn = $result_folder."/INDEL_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	cover_other_indel($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 0, "none", \@dna_controls);
    } else {
	$fileOut = $cover_other_folder."/F_cover_other_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	$fileIn = $result_folder."/INDEL_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	cover_other_indel($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 0, "none", "none");
    }
#   cover_other_snp($result_folder, $fileIn, $fileOut, 0, \@dna_cases, \@dna_controls);


    #Autre cas qui teste toutes les autres combinaisons
    if($combinaison == 1) {
    	my @bitMask = ();
    	#On cherche tous les subsets possibles pour la liste d'adn
    	while (FilterINDEL::generate_mask(\@bitMask, \@dna_cases)) {
    	    my $d_cases = FilterINDEL::get_one_subset(\@bitMask, \@dna_cases);

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
    		if (@dna_controls) {
    		    push @d_controls, @dna_controls;
    		}
    		my $ca = join('_', sort @$d_cases);
    		my $co = join('_', sort @d_controls);
		$fileOut = $cover_other_folder."_combi/F_cover_other_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		$fileIn = $result_folder."/INDEL_filter/${from}_combi/${from}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		cover_other_indel($log_folder, $log_file, $result_folder, $fileIn, $fileOut, 1, \@cases_as_controls, \@dna_controls);
    	    }
    	}
    }
    print "cover_other indel from $from DONE\n";
}

##################################################################################
sub cover_other_indel {
    my ($logFolder, $logFile, $result_folder, $fileIn, $fileOut, $filter, $ref_cases, $ref_controls) = @_;
    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)
    my ($log, $n_in, $n_out, $chrom, $pos, $size, @tab, $cover_line, @cover_tab);
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
	LogOut::printLog($log, $current_sub,  $logFolder, $logFile);
	return;
    }
    if (-f $fileOut) {
	$log = scalar(localtime())."  $fileOut already exists\n";
	#print $log;
	LogOut::printLog($log, $current_sub,  $logFolder, $logFile);
	return;
    }

    $log ="open $fileIn\n";
    LogOut::printLog($log, $current_sub,  $logFolder, $logFile);
    open (FHI, $fileIn) or die("cannot open $fileIn");
    $log = "create $fileOut\n";
    LogOut::printLog($log, $current_sub,  $logFolder, $logFile);
    open (FHO, ">$fileOut") or die("cannot create $fileOut");

    $n_in = 0;
    $n_out = 0;
    my (@control_threads, @case_threads);
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
	    for my $thr (@control_threads) {
		my $res = $thr->join();
#		print "control : @$res\n";
		$results{$$res[0], $$res[1], 'control'} = $$res[2];
	    }
	    undef(@control_threads);

	    for my $thr (@case_threads) {
		my $res = $thr->join();
#		print "case : @$res\n";
		$results{$$res[0], $$res[1], 'case'} = $$res[2]."\t".$$res[3];
	    }
	    undef(@case_threads);
	}

	@tab = split(/\t/, $line);
	$chrom = $tab[0];
	$pos = $tab[1];
	$size = length($tab[3]);

	# for controls
	$controls_cover = '';
	if ($dna_controls[0] ne "none") {
	    foreach my $control (@dna_controls) {
		my $file = "$result_folder/SEQ_cover/all/$control/chr$chrom.cov";
		my @arg = ($file, $pos, $n_in, $control);
		my $thread = threads->create('grep_cover_control', @arg);
		push @control_threads, $thread;
	    }
	}

	# for cases used as controls (combination cases)
	if ($dna_cases[0] ne "none") {
	    foreach my $case (@dna_cases) {
		my $file = "$result_folder/SEQ_cover/all/$case/chr$chrom.cov";
		my @arg = ($file, $pos, $size, $filter, $n_in, $case);
		my $thread = threads->create('grep_cover_case', @arg);
		push @case_threads, $thread;
	    }
	}
    }
    close(FHI);

    for my $thr (@control_threads) {
#	  print "wait $thr \n";
	  my $res = $thr->join();
#	  print "control : @$res\n";
	  $results{$$res[0], $$res[1], 'control'} = $$res[2];
    }

    undef(@control_threads);

    for my $thr (@case_threads) {
	  my $res = $thr->join();
#	  print "case : @$res\n";
	  $results{$$res[0], $$res[1], 'case'} = $$res[2]."\t".$$res[3];
    }
    undef(@case_threads);

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
		$cover = $results{$n_in, $control, 'control'};
		$controls_cover .= "$control:$cover;";
	    }
	}
	$controls_cover =~ s/;$//;

	# for cases used as controls (combination cases)
	if ($dna_cases[0] ne "none") {
	    $cases_cover = '';
	    $flag = 0;
	    foreach my $case (@dna_cases) {
		my $res =  $results{$n_in, $case, 'case'};
		my @res = split("\t", $res);
		$flag = $res[1];
		$cases_cover .= "$case:$res[0]";
		if ($#dna_cases > 0) {
		    $cases_cover .= ";";
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
    $log = "$n_out mutation(s) with insufficient cover in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    #print $log;
    LogOut::printLog($log, $current_sub,  $logFolder, $logFile);

}

##################################################################################
sub grep_cover_control {
  my $file = shift @_;
  my $pos = shift @_;
  my $line = shift @_;
  my $dna = shift @_;

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

#  print "finished $dna, $line = $cover\n";

  my @tab = ($line, $dna, $cover);

  return \@tab;

}

##################################################################################
sub grep_cover_case {
  my $file = shift @_;
  my $pos = shift @_;
  my $size = shift @_;
  my $filter = shift @_;
  my $line = shift @_;
  my $dna = shift @_;

  my ($cover_res);
  my $flag = 0;

  my $cover_line = `grep -w $pos $file`;
  my @cover_tab = split(/ /, $cover_line);
  my $cover = $cover_tab[1];
  if ($cover) {
      chomp $cover;
      $cover_res = $cover;
      if ($cover >= 8) {
	  if($filter == 1) {
	      $flag = 1;
	  }
	  if ($size > 1) {
	      for (my $i = 1; $i < $size; $i++) {
		  my $posi = $pos + $i;
		  $cover_line = `grep -w $posi $file`;
		  @cover_tab = split(/ /, $cover_line);
		  $cover = $cover_tab[1];
		  if ($cover) {
			chomp $cover;
		       if ($cover < 8) {
				$flag = 0;
			}
		      $cover_res .= "-$cover";
		  } else {
		      $flag = 0;
		      $cover_res .= "-0";
		  }
	      }
	  }
      }
  } else {
      $cover_res = 0;
  }

#  print "finished $dna, $line = $cover_res, $flag\n";

  my @tab = ($line, $dna, $cover_res, $flag);

  return \@tab;
}

##################################################################################
### Prepare for dbSNP_annot_indel
sub launch_dbSNP_annot_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotINDEL::dbSNP_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotINDEL::dbSNP_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotINDEL::dbSNP_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "dbSNP annot indel from $from DONE\n";
}

##################################################################################
### Annotate with dbSNP ID when available (otherwise 0)
sub dbSNP_annot_indel {
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
#   print $log;
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

    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
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
	if ($tab[19] ne "NA") {
#	    $start = $tab[1] - 1;
	    $start = $tab[19] - 1;  #col with hg19 pos
	    if (($tab[2] eq 'BP_RIGHT') || ($tab[2] eq 'BP_LEFT')) {
		$end = $start + length ($tab[3]);
	    } elsif ($tab[2] =~ /M/) {
		$end = $start;
	    } elsif (($tab[2] =~ /D/) && ($tab[2] !~ /I/)) {
		$end = $start + length ($tab[3]) - 1;
	    } elsif (($tab[2] !~ /D/) && ($tab[2] =~ /I/)) {
		$end = $start;
	    } elsif (($tab[2] =~ /D/) && ($tab[2] =~ /I/)) {
		$end = $start + length ($tab[3]) - 1;
	    }
	} else {
	    print FHO "$line\tNA\n";
	    next;
	}

	my @alleles;

	if ($fam == 1) {
	    my @ind_list = split(/\] /, $tab[17]);
	    foreach my $indiv(@ind_list) {
		$indiv =~ s/\]//;
		$indiv =~ s/.*\[//;
		$indiv =~ s/,.*//;
		my @ind_alleles = split(/\//, $indiv);
		foreach my $ia(@ind_alleles) {
		    if ($ia ne $tab[3]) {
			push(@alleles,$ia);
		    }
		}
	    }
	} else {
	    my @two_alleles = ($tab[4], $tab[5]);
	    foreach my $ia(@two_alleles) {
		if ($ia ne $tab[3]) {
		    push(@alleles,$ia);
		}
	    }
	}

	# Remove duplicate alleles
	my %seen=();
	my @uniq_alleles = grep { ! $seen{$_} ++ } @alleles;

	my $in_dbsnp;

	$req = "select observed, strand, refNCBI, name, class, alleleFreqs from $dbversion where " .
	      " chrom = \"chr$tab[0]\" " .
	      " and chromEnd = \"$end\" ";

	$sth = $dbh->prepare($req)||  print ("\nPreparation requete $req: echec\n");
	$sth->execute()||  print ("\nExecution requete: $req echec\n");
	while ( my @data = $sth->fetchrow_array ) {
	    my $observed = $data[0];
	    my $strand = $data[1];
	    my $ref_ncbi = $data[2];
	    my $name = $data[3];
	    my $class = $data[4];
	    my $freq = $data[5];
	    chop($freq);
	    unless ($freq) {
		$freq = 'NA';
	    }

	    my $ref_eq_obs = 0;

	    my @observed;

	    if ($class eq 'microsatellite') {
		@observed = ($observed);
	    } else {
		@observed = split('/', $observed);
	    }
	    my @reversed_observed;

	    if ($strand eq '-') {
		#print "- strand \n";
		# For each observed allele
		for my $obs (@observed) {
		    my $reversed_obs;
		    if ($obs =~ /LARGE/) {
			$reversed_obs = $obs;
		    } elsif ($obs =~ /INSERTION/) {
			$reversed_obs = $obs;
		    } elsif ($obs =~ /DELETION/) {
			$reversed_obs = $obs;
		    } elsif ($obs =~ /Long/) {
			$reversed_obs = $obs;
		    } elsif ($obs =~ /\//) {
			$reversed_obs = $obs;
		    } else {
		    my $seq_obj = Bio::Seq->new(-seq => $obs, -alphabet => 'dna' );
		    my $reversed_obj = $seq_obj->revcom;
		    $reversed_obs = $reversed_obj->seq;
		    }
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
### Prepare for CG_54genomes_annot_indel
sub launch_CG_54genomes_annot_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotINDEL::CG_54genomes_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotINDEL::CG_54genomes_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotINDEL::CG_54genomes_annot_indel($log_folder, $log_file, $sample_data, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "54genomes annot indel from $from DONE\n";
}

##################################################################################
### Annotate with observed frequency when present in Complete Genomics 54 genomes (otherwise 0)
sub CG_54genomes_annot_indel {
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
	if ($tab[19] ne "NA") {
#	    $start = $tab[1];
	    $start = $tab[19];  #col with hg19 pos
	} else {
	    print FHO "$line\tNA\n";
	    next;
	}

	my @alleles;

	if ($fam == 1) {
	    my @ind_list = split(/\] /, $tab[17]);
	    foreach my $indiv(@ind_list) {
		$indiv =~ s/\]//;
		$indiv =~ s/.*\[//;
		$indiv =~ s/,.*//;
		my @ind_alleles = split(/\//, $indiv);
		foreach my $ia(@ind_alleles) {
		    if ($ia ne $tab[3]) {
			push(@alleles,$ia);
		    }
		}
	    }
	} else {
	    my @two_alleles = ($tab[4], $tab[5]);
	    foreach my $ia(@two_alleles) {
		if ($ia ne $tab[3]) {
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
    ##$dbh->disconnect();
    $log = "$n_out mutation(s) not found in 54 genomes in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
#    print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
### Prepare for expr_beta_annot_indel
sub launch_expr_beta_annot_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotINDEL::expr_beta_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotINDEL::expr_beta_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotINDEL::expr_beta_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "expr_beta annot indel from $from DONE\n";
}

##################################################################################
### Annotate with 1 or 0 depending on expression in beta cell (lab data)
sub expr_beta_annot_indel {
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

    $log ="open $fileIn  [dna]=$sample_data->{dna}\n";
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
	    $expr = `grep -w $tab[11] $expr_file | wc -l`;
	} else {
	    $expr = `grep -w $tab[6] $expr_file | wc -l`;
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
### Prepare for good_annot_indel
sub launch_good_annot_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotINDEL::good_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotINDEL::good_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $fileIn, $fileOut);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotINDEL::good_annot_indel($log_folder, $log_file, $sample_data, $data_folder, $fam, $filter, $filesIn[$i], $filesOut[$i]);
	    }
	}
    }
    print "GOOD annot indel from $from DONE\n";
}

##################################################################################
### Annotate with GOOD diabetes or obesity ranking of genes where indels are found
sub good_annot_indel {
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

    $log ="open $fileIn  [dna]=$sample_data->{dna}\n";
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
	    if ($disease eq 'diab'){
		$good_data{$good_tab[0]} = $good_tab[1]."\t".$good_tab[2];
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
	    $good_line = $good_data{$tab[11]};
	} else {
	    $good_line = $good_data{$tab[6]};
	}

	if ($good_line) {
	    # There is at least one line so the gene is ranked in GOOD
	    if ($disease eq 'diab'){
		@good_tab = split(/\t/, $good_line);
		unless ($good_tab[1] =~ 'VALEUR') {
		    $good_score = $good_tab[1];
		} else {
		    $good_score = $good_tab[0];
		}
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
		$good_line = $good_data{$tab[9]};
	    } else {
		$good_line = $good_data{$tab[4]};
	    }

	    if ($good_line) {
		# There is at least one line so the gene is ranked in GOOD
		if ($disease eq 'diab'){
		    @good_tab = split(/\t/, $good_line);
		    unless ($good_tab[1] =~ 'VALEUR') {
			$good_score = $good_tab[1];
		    } else {
			$good_score = $good_tab[0];
		    }
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
#   print $log;
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);

}

##################################################################################
### Generate combinations of samples and associated file names
sub coseg_file_indel {
    my $analysis_data = shift;
    my $combinaison = shift;
    my $result_folder = shift;
    my $from = shift;
    my $to = shift;
    my $recess = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

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


    if (@dna_controls) {
#	if ($recess == 1) {
#	    my @non_parent_dna_controls = @dna_controls;
#	    splice(@non_parent_dna_controls, 0, 2);
#	    my @parent_dna_controls = ($dna_controls[0], $dna_controls[1]);
#	    if (@non_parent_dna_controls) {
#		$controls = join('_', sort @parent_dna_controls)."_".join('_',  sort @non_parent_dna_controls);
#	    } else {
#		$controls = join('_', sort @parent_dna_controls);
#	    }
#	} else {
#	    $controls = join('_',  sort @dna_controls);
#	}
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
	$fileIn = $result_folder."/INDEL_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}_vs_${controls}${suffix}.txt";
    } else {
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
	$fileIn = $result_folder."/INDEL_filter/${from}/${from}_from_${from_short}_in_family_${family}_adn_${cases}${suffix}.txt";
    }
    push @filesIn, $fileIn;
    push @filesOut, $fileOut;
    push @case_bag, 'none';
    push @control_bag, \@dna_controls;

    #Autre cas qui teste toutes les autres combinaisons
    if($combinaison == 1) {
    	my @bitMask = ();
    	#On cherche tous les subsets possibles pour la liste d'adn
    	while (FilterINDEL::generate_mask(\@bitMask, \@dna_cases)) {
    	    my $d_cases = FilterINDEL::get_one_subset(\@bitMask, \@dna_cases);

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
    		if (@dna_controls) {
    		    push @d_controls, @dna_controls;
    		}
    		my $ca = join('_', sort @$d_cases);
    		my $co = join('_', sort @d_controls);
		$fileOut = $result_folder."/INDEL_filter/${to}_combi/${to}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		$fileIn = $result_folder."/INDEL_filter/${from}_combi/${from}_combi_from_${from_short}_in_family_${family}_adn_${ca}_vs_${co}.txt";
		push @filesIn, $fileIn;
		push @filesOut, $fileOut;
		push @case_bag, \@cases_as_controls;
		push @control_bag, \@dna_controls;
    	    }
    	}
    }
    return (\@filesIn, \@filesOut, \@case_bag, \@control_bag);
}

sub launch_base_count_indel {
    my $run_name = shift;
    my $run_disk = shift;
    my $sample_data = shift;
    my $result_folder = shift;
    my $from = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($chr, $pos, $vcf_file, $log, $fileIn, $fileOut, $log_folder, $log_file, $position, @line, @line_vcf, @count, @position, $position_vcf, $count_format, $ref_alt, %allele_count, %f_coding, $head);
    if ($run_name=~/^run/){
	$log_folder = $run_name;
	$log_file = $sample_data->{dna};
    } else {
	$log_folder = $run_name;
	$log_file = $run_name;
    }

#    if ($sample_data->{maladie} eq 'unkn'){
#	$from = "F_maf";
#    }

    my $to = 'F_basecount';

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");

    $vcf_file = $run_disk."/".$run_name."/variation/vcf_gatk/Sample_".$sample_data->{dna}."/".$sample_data->{dna}."_indels.vcf";
    $fileOut = $result_folder."/INDEL_filter/".$to."/F_basecount_".$sample_data->{result_dir}.".txt";
    $fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

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
	$pos = $line_vcf[1] + 1;
	$position_vcf = $chr.":".$pos;
	next if $line_vcf[9]=~/^.\/.$/;
	@count = split(":",$line_vcf[9]);
	$count_format= $count[1];
	$count_format=~s/,/\//;
	$ref_alt = $line_vcf[3]."/".$line_vcf[4];
	$allele_count{$position_vcf} = $ref_alt.":".$count_format;

    }
    close VCF;

    # recuperer lignes dans F_coding
    $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
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
}

# sub launch_inCNV_annot_indel {
#     my $run_name = shift;
#     my $sample_data = shift;
#     my $result_folder = shift;
#     my $from = shift;
#     my $cnv_chr = shift;

#     my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

#     my @chr = (1..22);
#     my ($id_cnv, $chr, $pos, %cnv_chr, %lines, $start, $stop, $interval, %out, $head, $head2, $log_folder, $log_file, $fileIn, $fileOut, $log);
#     if ($run_name=~/^run/){
# 	$log_folder = $run_name;
# 	$log_file = $sample_data->{dna};
#     } else {
# 	$log_folder = $run_name;
# 	$log_file = $run_name;
#     }

#     my $to = 'F_incnv';

#     system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
#     system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");

#     $fileOut = $result_folder."/INDEL_filter/".$to."/F_incnv_".$sample_data->{result_dir}.".txt";
#     $fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

#     if (! -f $fileIn) {
# 	$log = scalar(localtime())."  !!!missing file $fileIn\n";
# #	print $log;
# 	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
# 	return;
#     }

#     if (-f $fileOut) {
# 	$log = scalar(localtime())."  $fileOut already exists\n";
# #	print $log;
# 	LogOut::printLog($log, $current_sub, $log_folder, $log_file);
# 	return;
#     }

#     $log = "open $fileIn  [dna]=$sample_data->{dna}\n";
#     LogOut::printLog($log, $current_sub, $log_folder, $log_file);

#     open (IN,"$fileIn") || die "Cannot open $fileIn : $!\n";
#     while (my $l=<IN>){
# 	chomp $l;
# 	if ($l =~ /^#/){
# 	    $head = $l;
# 	} elsif ($l =~ /(^[\d]+)\s+([\d]+)/){
# 	    $pos = $2;
# 	    $chr = $1;
# 	    $lines{$chr}{$pos} = $l;
# 	}
#     }
#     close IN;

#     $head2 = $head."\t1KG_CNV_id";

#     for my $chr (@chr){
# 	for my $mut (keys $lines{$chr}){
# 	    for my $interval (keys $$cnv_chr{$chr}){
# 		($start,$stop)= split("-",$interval);
# 		if ($mut ~~ [$start..$stop]){
# 		    my $x = $start."-".$stop;
# 		    $out{"$chr:$mut"} = $lines{$chr}{$mut}."\t".$$cnv_chr{$chr}{$x};
# 		    last;
# 		} else {
# 		    $out{"$chr:$mut"} = $lines{$chr}{$mut}."\tNA";
# 		}
# 	    }
# 	}
#     }

#     $log = "create $fileOut [dna]=$sample_data->{dna}\n";
#     LogOut::printLog($log, $current_sub, $log_folder, $log_file);
#     open (OUT,">$fileOut") || die "Cannot open $fileOut : $!\n";
#     print OUT "$head2\n";
#     foreach my $i (sort keys %out){
# 	print OUT "$out{$i}\n";
#     }
#     close OUT;
# }

##################################################################################
### Prepare for maf_annot_indel
sub launch_maf_annot_indel {
    my $run_name = shift;
    my $sample_data = shift;
    my $dbversion = shift;
    my $result_folder = shift;
    my $from = shift;
    my $fam = shift;
    my $filter = shift;
    my $combinaison = shift;
    my $recess = shift;
    my $local_server = shift;

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
	    $from = "F_expr";
	}
    }

    my $to = 'F_maf';

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");
    system("mkdir $result_folder/INDEL_filter/${to}_combi") unless (-d "$result_folder/INDEL_filter/${to}_combi");

    if ($fam == 0){
	$fileOut = $result_folder."/INDEL_filter/${to}/${to}_".$sample_data->{result_dir}.".txt";
	$fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";
	AnnotINDEL::maf_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut, $local_server);
    } elsif ($fam == 1) {
	my ($filesIn_ref, $filesOut_ref, $case_bag_ref, $control_bag_ref) = AnnotINDEL::coseg_file_indel($sample_data, $combinaison, $result_folder, $from, $to, $recess);
	my @filesIn = @$filesIn_ref;
	my @filesOut = @$filesOut_ref;
	if ($combinaison == 0){
	    $fileIn = $filesIn[0];
	    $fileOut = $filesOut[0];
	    AnnotINDEL::maf_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $fileIn, $fileOut, $local_server);
	} elsif ($combinaison == 1) {
	    for (my $i = 0; $i < scalar @filesIn; $i++) {
		AnnotINDEL::maf_annot_indel($log_folder, $log_file, $sample_data, $dbversion, $fam, $filter, $filesIn[$i], $filesOut[$i], $local_server);
	    }
	}
    }
    print "MAF annot indel from $from DONE\n";
}

##################################################################################
### Annotate with MAF queried from Ensembl
sub maf_annot_indel {
    my $logFolder = shift;
    my $logFile = shift;
    my $sample_data = shift;
    my $dbversion = shift;
    my $fam = shift;
    my $filter = shift;
    my $fileIn = shift;
    my $fileOut = shift;
    my $local_server = shift;

    my $current_sub = (caller(0))[3]; # renvoie le nom de la subroutine actuelle (pour les fichiers de log)

    my ($log, $n_in, $n_out, @tab);

    my $port;
    if (($dbversion < 48) || ($local_server == 1)) {
	$port = 'XXX';
    } else {
	$port = 'XXX';
    }

    my $organism = 'homo_sapiens';

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

    # Get registry
    my $registry = 'Bio::EnsEMBL::Registry';

    if ($local_server == 1) {
	$registry->load_registry_from_db(
	    -host => 'XXX',
	    -user => 'XXX',
	    -pass => "XXX",
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

    my $slice_adaptor = $registry->get_adaptor($organism, 'core', 'slice');
    my $variation_feature_adaptor = $registry->get_adaptor($organism, 'variation', 'variationfeature');

    $n_in = 0;
    $n_out = 0;

    while (my $line =  <FHI>) {
	$n_in++;
	chomp $line;
	if ($line =~ /^\#/){
	    $n_in--;
	    print FHO "$line\t1KG_ALL\t1KG_AFR\t1KG_AMR\t1KG_ASN\t1KG_EUR\tESP6500_African_American\tESP6500_European_American\n";
	    next;
	}
#	my ($chromosome, $position, $CIGAR, $ref_allele, $allele1, $allele2, @rest) = split (/\t/, $line);
	@tab = split(/\t/, $line);
	my ($chromosome, $position, $allele2, @alleles);
	$chromosome = $tab[0];
	$position = $tab[1];

	if ($fam == 1){
	    my @ind_list = split(/\] /, $tab[17]);
	    foreach my $indiv(@ind_list) {
		$indiv =~ s/\]//;
		$indiv =~ s/.*\[//;
		$indiv =~ s/,.*//;
		my @ind_alleles = split(/\//, $indiv);
		foreach my $ia(@ind_alleles) {
		    if ($ia ne $tab[3]) {
			push(@alleles,$ia);
		    }
		}
	    }
	} else {
	    my @two_alleles = ($tab[4], $tab[5]);
	    foreach my $ia(@two_alleles) {
		if ($ia ne $tab[3]) {
		    push(@alleles,$ia);
		}
	    }
	}

	# Remove duplicate alleles
	my %seen=();
	my @uniq_alleles = grep { ! $seen{$_} ++ } @alleles;
	$allele2 = $uniq_alleles[0];
	if ($#uniq_alleles > 0) {
	    $log  = scalar(localtime())." : More than one non ref allele at $chromosome - $position, only $allele2 treated \nSample = ".$sample_data->{result_dir}."\n";
	    LogOut::printLog($log, $current_sub, $sample_data->{run}, $sample_data->{dna});
	}

	if ($allele2 =~ /-/){
	    $allele2 = "-";
	}

	my $slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome, $position, $position);

	# Check if maf is annotated in EnsEMBL
	my $variation_features = $variation_feature_adaptor->fetch_all_by_Slice($slice);
	my $KG_ALL_freq = 2; # Initialize with impossible value
	my $KG_AFR_freq = 2; # Initialize with impossible value
	my $KG_AMR_freq = 2; # Initialize with impossible value
	my $KG_ASN_freq = 2; # Initialize with impossible value
	my $KG_EUR_freq = 2; # Initialize with impossible value
	my $ESP_AA_freq = 2; # Initialize with impossible value
	my $ESP_EA_freq = 2; # Initialize with impossible value

	foreach my $vf (@{$variation_features}){
	    my @db_alleles = @{$vf->variation->get_all_Alleles()};
	    my $id = $vf->variation->name;
	    foreach my $allele (@db_alleles) {
		if ($allele->population && ($allele->allele eq $allele2)) {
		    if ($allele->population->name =~ /1000GENOMES.*ALL/) {
			my $KG_ALL_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($KG_ALL_freq_tmp < $KG_ALL_freq){
			    $KG_ALL_freq = $KG_ALL_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /1000GENOMES.*AFR/) {
			my $KG_AFR_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($KG_AFR_freq_tmp < $KG_AFR_freq){
			    $KG_AFR_freq = $KG_AFR_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /1000GENOMES.*AMR/) {
			my $KG_AMR_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($KG_AMR_freq_tmp < $KG_AMR_freq){
			    $KG_AMR_freq = $KG_AMR_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /1000GENOMES.*ASN/) {
			my $KG_ASN_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($KG_ASN_freq_tmp < $KG_ASN_freq){
			    $KG_ASN_freq = $KG_ASN_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /1000GENOMES.*EUR/) {
			my $KG_EUR_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($KG_EUR_freq_tmp < $KG_EUR_freq){
			    $KG_EUR_freq = $KG_EUR_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /ESP6500:African_American/) {
			my $ESP_AA_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($ESP_AA_freq_tmp < $ESP_AA_freq){
			    $ESP_AA_freq = $ESP_AA_freq_tmp;
			}
		    }
		    if ($allele->population->name =~ /ESP6500:European_American/) {
			my $ESP_EA_freq_tmp = (defined($allele->frequency) ? $allele->frequency : "0");
			if ($ESP_EA_freq_tmp < $ESP_EA_freq){
			    $ESP_EA_freq = $ESP_EA_freq_tmp;
			}
		    }
		}
	    }
	}
	my $n_out_all = 0;
	if ($KG_ALL_freq == 2) { # No other value found, freq = 0
	    $KG_ALL_freq = 0;
	    $n_out_all++;
	}
	if ($KG_AFR_freq == 2) { # No other value found, freq = 0
	    $KG_AFR_freq = 0;
	    $n_out_all++;
	}
	if ($KG_AMR_freq == 2) { # No other value found, freq = 0
	    $KG_AMR_freq = 0;
	    $n_out_all++;
	}
	if ($KG_ASN_freq == 2) { # No other value found, freq = 0
	    $KG_ASN_freq = 0;
	    $n_out_all++;
	}
	if ($KG_EUR_freq == 2) { # No other value found, freq = 0
	    $KG_EUR_freq = 0;
	    $n_out_all++;
	}
	if ($ESP_AA_freq == 2) { # No other value found, freq = 0
	    $ESP_AA_freq = 0;
	    $n_out_all++;
	}
	if ($ESP_EA_freq == 2) { # No other value found, freq = 0
	    $ESP_EA_freq = 0;
	    $n_out_all++;
	}

	if ($n_out_all == 0){
	    $n_out++;
	}

	print FHO "$line\t$KG_ALL_freq\t$KG_AFR_freq\t$KG_AMR_freq\t$KG_ASN_freq\t$KG_EUR_freq\t$ESP_AA_freq\t$ESP_EA_freq\n";
    }
    close FHI;
    close FHO;
    $log = "$n_out mutation(s) not found in Ensembl $dbversion in $fileOut, $n_in read from $fileIn\n".scalar(localtime())."\n";
    LogOut::printLog($log, $current_sub, $logFolder, $logFile);
}

##################################################################################
###
sub launch_disease_inher_indel {
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

    system("mdkir $result_folder/INDEL_filter") unless (-d "$result_folder/INDEL_filter");
    system("mkdir $result_folder/INDEL_filter/${to}") unless (-d "$result_folder/INDEL_filter/${to}");

    $fileOut = $result_folder."/INDEL_filter/".$to."/F_disease_inher_".$sample_data->{result_dir}.".txt";
    $fileIn = $result_folder."/INDEL_filter/".$from."/".$from."_".$sample_data->{result_dir}.".txt";

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
	$gene = $tab[11];
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

    $log = "$n_out mutation(s) not found in raindance 3 gene panel in $fileOut, $n_in read from $fileIn\n";
    LogOut::printLog($log, $current_sub, $log_folder, $log_file);
}

1;
