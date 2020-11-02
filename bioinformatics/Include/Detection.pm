package Detection;

use Sys::Hostname;
use Data::Dumper;

sub read_sample {
    my ($disque,$name) = @_;
    open(FH, "$disque/$name/samplesheet.csv");
    while ($tmp = <FH>) {
	@col=split(/\,/, $tmp);
	$echant = $col[2];
	next if $col[2] eq "SampleID";
	push(@echant, $echant);
    }
    close(FH);

    foreach $objet (@echant){
	push(@adn, $objet) unless $deja_vu{$objet}++;
    }
    return @adn;
}

sub gatk_mutation {
    my $disque = shift;
    my $run_name = shift;
    my $adn = shift ;
    my $config = shift;
    my $run_data = shift;
    my $ref_genome = shift;

    my ($sample_data, $output_dir, $capture, $bedfile, $recal_file, $snps_vcf, $indels_vcf, $HC_vcf, $bed_option);

    my $gatk_launch;
    my $host = hostname;


    $gatk_launch = "gatk3";

    my $log_folder = "$disque/$run_name/log";
    my $gatk_alignement_dir = "$disque/$run_name/alignement/gatk";
    my $variation_dir = "$disque/$run_name/variation";
    my $gatk_variation_dir = "$variation_dir/vcf_gatk";
    unless (-d "$variation_dir") { system("mkdir $variation_dir"); }
    unless (-d "$gatk_variation_dir") { system("mkdir $gatk_variation_dir"); }

    foreach my $dna (@{ $adn }){
	$dna_name = $dna; ###### a enlever en fonction erreur noms samples (BDD ou bio)
	$sample_data = $$run_data->{$dna_name};

	$output_dir = "$gatk_variation_dir/Sample_$dna";
	unless (-d "$output_dir") { system("mkdir $output_dir"); }

	$capture = $sample_data->{capture};
	$bedfile = $$config{$capture};
	if ($bedfile eq "/media/Data/bed/capture/Twist/Exome_V_1_3_0_RefSeq_Gencode_Probe_Targets_hg38.bed"){
		$bedfile = "/media/Data/bed/capture/Twist/Exome+50bp/Exome_V_1_3_0_RefSeq_Gencode_Probe_Targets_hg38_+50bp.bed";
	}
	if ($bedfile eq "WGS"){
	    $bed_option = " ";
	} else {
	    $bed_option = "-L $bedfile";
	}
	$recal_file = "$gatk_alignement_dir/Sample_$dna/$dna\_RG.dedup.sorted.realigned.recal.bam";

	if (!(-f $recal_file)){
	    print "recalibrated BAM file $recal_file does not exist!!!\n";
	    next;
	}

	$snps_vcf = "$output_dir/$dna\_snps.vcf";
	$indels_vcf = "$output_dir/$dna\_indels.vcf";
	$HC_vcf = "$output_dir/HC.raw.vcf";

	### GATK : HaplotypeCaller (plus performant mais plus lent)
	 system ("java -Xms24g -Xmx48g -jar $$config{$gatk_launch} -T HaplotypeCaller -allowPotentiallyMisencodedQuals -I  $recal_file -R $ref_genome $bed_option --dbsnp $$config{hg38_dbsnp_file} -stand_call_conf 10.0 -l INFO -log $log_folder/$dna\_HC.log --disable_auto_index_creation_and_locking_when_reading_rods -o $HC_vcf");

	### GATK : SelectVariants : select SNPs
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T SelectVariants -R $ref_genome --variant $HC_vcf -selectType SNP -log $log_folder/$dna\_selectSNP.log --disable_auto_index_creation_and_locking_when_reading_rods -o $snps_vcf -select \"DP > 8.0\"");

	### GATK : SelectVariants : select indels
	system ("java -Xms12g -Xmx24g -jar $$config{$gatk_launch} -T SelectVariants -R $ref_genome --variant $HC_vcf -selectType INDEL -log $log_folder/$dna\_SelectIndels.log --disable_auto_index_creation_and_locking_when_reading_rods -o $indels_vcf -select \"DP > 8.0\"");

    }
}

sub vcf_bwa{
    my ($disque,$name,$samtools,$bcftools,$dna,$genome_bwa,$bed) = @_;
    my $alignement_dir = "$disque/$name/alignement";
    my $bwa_dir = $alignement_dir."/bwa";
    my $variation_dir = "$disque/$name/variation";
    my $variation_bwa_dir = "$variation_dir/bwa";
    unless ( -d "$variation_dir" ) { system("mkdir $variation_dir"); }
    unless ( -d "$variation_bwa_dir" ) { system("mkdir $variation_bwa_dir"); }
    print $dna."\n";
    my $sample_ali_dir = "$bwa_dir/Sample_$dna";
    my $sample_var_dir = "$variation_bwa_dir/Sample_$dna";
    unless ( -d "$sample_var_dir" ) { system("mkdir $sample_var_dir"); }
    system("$samtools mpileup -l $bed -u -g -B -m 3 -C 50 -d 1000000 -L 1000000 -F 0.0002 -Q 0 -f $genome_bwa $sample_ali_dir/$dna.sorted.bam | $bcftools view -v -c -g -p 0.20 - >  $sample_var_dir/$dna.raw.vcf");
    system("$samtools depth $sample_ali_dir/$dna.sorted.bam > $sample_var_dir/$dna.cov");
}


1;
