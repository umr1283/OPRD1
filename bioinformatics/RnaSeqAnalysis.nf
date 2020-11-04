#!/usr/bin/env nextflow
// ###################################################################################################################################
// script to analyse paired and single files
// Adapted from Rnaseq.pm of pipeline anges to analyze RNASeq data
// Autors: Amanzougarene Souhila
// Version 2

// General architecture of a nextflow run
// 
// run_xxx
// |-- Input
// |   |-- ...
// |   |-- Raw
// |   `-- Reference
// |-- Logs
// |   |-- bcl2fastq
// |   |-- RunParameters
// |   |-- MultiQC
// |   `-- Steps
// |-- Scripts
// |-- Temporary
// |-- Tools
// `-- Output
//     |-- ...
//     |-- BAM
//     |-- FastQ
//     |-- RSEM
//     |-- QC
//     |   |-- ...
//     |   `-- BAM
//     `-- Reporting 
//               |-- Stats
// 
// ###################################################################################################################################
///////////////////////////////////////////////////////////////
//                                        _                  //
//   _ __   __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___   //
//  | '_ \ / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|  //
//  | |_) | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \  //
//  | .__/ \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/  //
//  |_|                                                      //
//////////////////////////////////////////////////////////////

if (params.help) { 
    log.info "----------------------------------------------------"
    log.info "| 				      USAGE                      |"
    log.info "----------------------------------------------------"
    log.info ''
    log.info 'nextflow run RnaSeqAnalysis.nf'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info "~~~~~~~~~~~~~~~~~~~~"
    log.info ""
    log.info ""
    log.info '--params.useBasesMask STRING 	The parameter for bcl2fastq2 to be used if the flowcell contains different types of sequenced data (e.g. Paired-End with dual-index and single-index)'
    log.info ""
    log.info '--params.userName		STRING		The name of the user who launched the run analysis'
    log.info '--params.run_no		STRING		This is the run number. (e.g. run_xxx)'
    log.info '--params.output		STRING		The default output folder location: '
    log.info '--params.genome		STRING 		The reference genome to be used'
    log.info ""
    log.info 'Optional arguments:'
    log.info "~~~~~~~~~~~~~~~~~~~"
    log.info '-with-report			STRING		The HTML report that will be created ; a single document which includes many useful metrics about a workflow execution. The report is organised in the three main sections: Summary, Resources and Tasks'
    log.info '-with-trace			STRING		The tracing file that will be created : it  contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used'
    log.info '-with-timeline		STRING		Nextflow can render an HTML timeline for all processes executed in your pipeline'
    log.info '-with-dag				STRING		It creates a file named dag.dot containing a textual representation of the pipeline execution graph in the DOT format'
    log.info ''
    log.info ''
    log.info 'Example command: nohup nextflow run RnaSeqAnalysis.nf --run_no run_111 --sequencer Miniseq --folder 181105_NS500777 --sampleSheet SampleSheet.csv --genome mm10 -w work --userName Souhila --output Run/ -bg > RnaSeqAnalysis_run_111_230419.log'
    exit 0
}else {
  params.help = null
  log.info ''
  log.info "For more information about the Pipeline, type nextflow run RnaSeqAnalysis.nf --help"
  log.info ''
}


//#############################################
// Pipeline parameters params.nameOfVariable
//#############################################

// Name of the user who launched the run analysis 
params.userName = null

// This is the run number. (e.g. run_111)
params.run_no = null

// The default output folder location.
params.output = "/Run/"

// The reference genome, ncRNA to be used
params.genome = null

switch(params.genome) {
	case "hg38":  //human
		genome = "GRCh38.fa"
		genomeDA = "/hg38"
		genomeStar = "/genome_STAR-2.7.3/GenCode/"
		genomebamToCram = "GRCh38.fa"
		genomeBamTranscritToCram = "hg38.transcripts.fa"
		ncrna = "Homo_sapiens.GRCh38.98.ncrna"
		host = "jan2020.archive.ensembl.org"
		dataset = "hsapiens_gene_ensembl"
		knownSites1 = "1000G_omni2.5.hg38.vcf"
		knownSites2 = "Mills_and_1000G_gold_standard.indels.hg38.vcf"
		dbSNP = "dbsnp_146.hg38.vcf"
		
		break
	case "hg19":
		genomeDA = "/Homo_sapiens"
		ncrna = "/Homo_sapiens.GRCh37.75.ncrna"
		host = "feb2014.archive.ensembl.org"
		dataset = "hsapiens_gene_ensembl"
		break
	case "mm10": //mouse mm10
		genome = "genome_mm10_GenCode/mm10.fa"
		genomeDA = "V1.3.1/mm10"
		genomeStar = "genome_STAR-2.7.3/"
		genomebamToCram = "genome_mm10_GenCode/mm10.fa"
		genomeBamTranscritToCram = "V1.3.1/mm10.transcripts.fa"
		ncrna = "Mus_musculus.GRCm38.98.ncrna"
		host = "jan2020.archive.ensembl.org"
		dataset = "mmusculus_gene_ensembl"
		dbSNP ="dbSNP/mus_musculus_snp_indels.vcf"
		break
}

//transcriptome

params.transcriptome = null

switch(params.transcriptome) {
	case "hg19":
		transcriptome = "hg19/RNAseq/transcriptome"
		break
	case "mm10": //mouse mm10
		transcriptome = "mm10/RNAseq/transcriptome"
		break
}

params.ncrna = null

switch(params.ncrna) {
	case "hg38":  //human
		ncrna = "Homo_sapiens.GRCh38.89.ncrna"
		break
	case "hg19":
		ncrna = "Homo_sapiens.GRCh37.75.ncrna"
		break
	case "mm10": //mouse mm10
		ncrna = "Mus_musculus.GRCm38.ncrna"
		break
	case "sscrofa": //pig Sscrofa11.1
		ncrna = "Sus_scrofa.Sscrofa11.1.ncrna"
		break
}


//#############################################
// Check if command line arguments are null
//#############################################

if (params.userName == null) 
{println " \n **************** ERROR : params.userName is empty !!! ****************** \n"
		System.exit(0) //Arréter l'analyse
		}

println "***************** INFOS : The " + params.run_no + " has been launched by : " + params.userName + "*******************"
//#############################################
// Absolute Path 
//#############################################

// Path of directory to copy files TO
to = params.output + params.run_no
// ------------------RUN FOLDER------------------------

// Path to Input folder 
toInput = to + "/Input"

// Path to Scripts folder
toScripts = to + "/Scripts"

// Path to Logs folder
toLogs = to + "/Logs"

// Path to Output folder
toOutput = to + "/Output"

// Path to Temporary folder
toTMP = to + "/Temporary"

// Path to Tools folder
toTools = to + "/Tools"

// -----------------INPUT FOLDER-------------------------

// Path to Raw folder in Input folder
toRaw = toInput + "/Raw"

// Path to Reference folder in Input folder
toRef = toInput + "/Reference"


// -------------------OUTPUT FOLDER-----------------------

//  Path to FastQ folder
toFastQ = toOutput + "/FastQ"

//  Path to FastQ folder after they are grouped
toFastq = toFastQ + "/Fastq"

// Path to FastQ_Cut folder in Output/FastQ/ folder, it contains the fastqc files after cutadapt 
toFastQCut = toFastQ + "/FastQ_Cut"

// Path to FastQ_NCRNA folder in Output/FastQ/ folder, it contains the fastqc files after cutadapt
toFastQncrna = toFastQ + "/FastQ_NCRNA"

// Path to RSEM folder in Output folder, it contains the files .results of RSEM 
toRSEM = toOutput + "/RSEM"

// Path to BAM folder in Output folder after they are being sorted
toBAM = toOutput + "/BAM"

// Path to QC folder in Output folder, il contains files html of fastqc, multiQC and file htm of demultiplexage
toQC = toOutput + "/QC"

// Demultiplexing reports (html file)
toQCDx = toQC + "/1.Demultiplexe"

// FastQC reports (html file)
toFastQC = toQC + "/2.FastQC"

// MultiQC reports (html file)
toMultiFQC = toQC + "/3.MultiQC"

// Path to Reporting folder in Output folder, it contains the stats files
toReporting = toOutput + "/Reporting"

// Path to fastqc_dup folder in Reporting folder, it contains the stats files of fastqc; % of duplicates reads
toReportingFastqc = toReporting + "/Fastqc_dup"

toReportingBowtie = toReporting + "/Bowtie_ncRNA"

toReportingRSEM = toReporting + "/Rsem_mapp"

toReportingDemultpx = toReporting + "/Demultiplexage"

// ------------------LOGS FOLDER------------------------

// Path to run parameters
toRunParameters = toLogs + "/RunParameters"

// Path to bcl2fastq2 log in Logs folder : log of demultiplexage
toLogsBCL = toLogs + "/Bcl2fastq2"

// Path to analysis steps in Logs folder
toLogsSteps = toLogs + "/Steps"

// Path to MultiQC in Logs folder : log of MultiQC
toLogsMultiQC = toLogs + "/MultiQC"

toLogsCutadapt = toLogs + "/Cutadapt"

toLogsFastQC = toLogs + "/FastQC"

toLogsRSEM = toLogs + "/Rsem"

toLogsStar = toLogs + "/starAlign"

toLogsMarkDuplicates = toLogs + "/MarkDuplicates"

toLogsSplitNCigar = toLogs + "/SplitNCigar"

toLogsBQSR = toLogs + "/BQSR"

toLogsHP = toLogs + "/HC"

toLogsVEP = toLogs + "/Vep"

// Path to logs generated by the process reorderBAM
toLogsReorderBam = toLogs + "/ReorderBAM"

// ------------------TEMPORARY FOLDER------------------------

// Path to FastQ folder before they are grouped
toPreFSQ = toTMP + "/1.FastQ"

// Path to BAM folder after they are being aligned
toPreBAM = toTMP + "/2.PreBAM"

// ------------------SNP FOLDER------------------------

toVCF = toOutput + "/VCF"

toVCFOriginal = toVCF + "/Original"

toVCFOriginalsnp = toVCFOriginal + "/SNP"

toVCFOriginalIndel = toVCFOriginal + "/INDEL"

toVCFFilter = toVCF + "/Filtered"

toVCFFilterSNP = toVCFFilter + "/SNP"

toVCFFilterIndel = toVCFFilter + "/INDEL"

toVCFAnnot = toVCF + "/Annotated"

toVCFAnnotSNP = toVCFAnnot + "/SNP"

toVCFAnnotINDEL = toVCFAnnot + "/INDEL"

toBAMstar = toOutput + "/BAM_align"

toDedupBAM = toTMP + "/3.DedupBAM"

toReorderBAM = toTMP + "/4.ReorderBAM"

toRawAnnotatdVCF = toTMP + "/5.AnnotatedRawVCF"

toRawAnnotatdSNPvcf = toRawAnnotatdVCF + "/SNP"

toRawAnnotatdINDELvcf = toRawAnnotatdVCF + "/INDEL"


// This checks if the workflow has to resume the previous run or not
// in case of an unexpected error by checking if the 'resume' flag
// is given in the command line or not
// Then we use this to define the 'to' variable below
isResume = "${workflow.resume}"

if ( isResume == false && new File(toInput).exists() == true ) {
	println("""
   WARNING
   WARNING
   WARNING

   The folder $to exists and the parameter '-resume' was not specified in the command line !
   This means that the Run will be re-analyzed. Are you sure ?

   WARNING
   WARNING
   WARNING
""")
	for (i in 1..10) {
		sleep(1000)
		println("$i/180")
	} 
} else if ( isResume == true && new File(toInput).exists() == true){
	println("\n\n  Resuming ${params.run_no} ...\n\n")
} else {
	println("\n\n This looks like a new run ... \n\n")
}


///////////////////////////////////////////////////
//   __                  _   _                   //
//  / _|_   _ _ __   ___| |_(_) ___  _ __  ___   //
// | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|  //
// |  _| |_| | | | | (__| |_| | (_) | | | \__ \  //
// |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/  //
///////////////////////////////////////////////////

switch (params.level) {
	case "1" :
	level = 1
	listPaths = [to, toInput, toScripts, toLogs, toOutput, toTools, toRaw, toFastQ, toFastq, toFastQC, toQC, toLogsBCL, toLogsSteps, toPreFSQ, toQCDx, toReporting, toReportingDemultpx, toReportingFastqc, toMultiFQC, toLogsFastQC, toLogsMultiQC ]
	break
	case "2" :
	level = 2
	listPaths = [to, toInput, toScripts, toLogs, toOutput, toTMP, toTools, toRaw, toRef, toFastQ,toFastq, toFastQC, toBAM, toQC, toLogsBCL, toLogsSteps, toPreFSQ, toFastQCut, toFastQncrna, toPreBAM, toQCDx, toRSEM, toLogsFastQC, toLogsMultiQC, toLogsCutadapt, toLogsRSEM, toLogsStar, toLogsMarkDuplicates, toLogsSplitNCigar, toLogsBQSR, toLogsHP, toReporting, toReportingDemultpx, toReportingFastqc, toReportingBowtie, toReportingRSEM, toMultiFQC]
	break
	case "3" :
	level = 3
	listPaths = [to, toInput, toScripts, toLogs, toOutput, toTMP, toTools, toRaw, toRef, toFastQ,toFastq, toFastQC, toBAM, toQC, toLogsBCL, toLogsSteps, toPreFSQ, toFastQCut, toFastQncrna, toPreBAM, toQCDx, toRSEM, toLogsFastQC, toLogsMultiQC, toLogsCutadapt, toLogsRSEM, toLogsStar, toLogsMarkDuplicates, toLogsSplitNCigar, toLogsBQSR, toLogsHP, toLogsVEP, toReporting, toReportingDemultpx, toReportingFastqc, toReportingBowtie, toReportingRSEM, toMultiFQC, toBAMstar, toDedupBAM, toLogsReorderBam, toReorderBAM,toVCF,,toVCFOriginal,toVCFOriginalsnp,toVCFOriginalIndel,toVCFFilter,toVCFFilterSNP,toVCFFilterIndel,toVCFAnnot,toVCFAnnotSNP,toVCFAnnotINDEL, toRawAnnotatdVCF, toRawAnnotatdSNPvcf, toRawAnnotatdINDELvcf]
	
}


myFileChannel = Channel.fromPath( 'Path To .fastq file' )
myFileChannelRSEM = Channel.fromPath( 'Path To .fastq file' )
myFileChannelbowtie = Channel.fromPath( 'Path To .fastq file' )

switch(params.genome) {
	case "hg38":
		samplesList = ["GSM4618781", "GSM4618783", "GSM4618785", "GSM4618787", "GSM4747383", "GSM4747385", "GSM4747387", "GSM4747389", "GSM4747391", "GSM4747393", "GSM4747395", "GSM4705022", "GSM4705023", "GSM4705024", "GSM4705025", "GSM4705026", "GSM4705027"]
		break
	case "hg19": //hg19
		samplesList = ["GSM4618781", "GSM4618783", "GSM4618785", "GSM4618787", "GSM4747383", "GSM4747385", "GSM4747387", "GSM4747389", "GSM4747391", "GSM4747393", "GSM4747395", "GSM4705022", "GSM4705023", "GSM4705024", "GSM4705025", "GSM4705026", "GSM4705027"]

		break
	case "mm10": //mouse mm10
		samplesList = ["GSM4618386", "GSM4618390", "GSM4618391", "GSM4618392", "GSM4618393", "GSM4503630", "GSM4503631", "GSM4503632", "GSM4503633", "GSM4503634", "GSM4503635", "GSM4503636", "GSM4503637", "GSM4503638"]
}

sample_size = samplesList.size()
samplesListBash = samplesList.join(' ')

println "SamplesList = " + samplesList
println "sample_size = " + sample_size
println "samplesListBash : " + samplesListBash


read1 = "75"
read2 = "75"
//////////////////////////////////////////////////////////
//  _              _                        _   _       //
// | |_ ___   ___ | |___        _ __   __ _| |_| |__    //
// | __/ _ \ / _ \| / __|      | '_ \ / _` | __| '_ \   //
// | || (_) | (_) | \__ \      | |_) | (_| | |_| | | |  //
//  \__\___/ \___/|_|___/      | .__/ \__,_|\__|_| |_|  //
//                             |_|                      //
//////////////////////////////////////////////////////////


samtools = "/opt/samtools/bin/samtools"
rsem = "/usr/local/bin/"
star = "/usr/local/bin/"
rsem_addgenename = "/home/rnaseq_add_genenameOld.R" 
bowtie2 = "/usr/local/bin/bowtie2"
multiqcconf = "multiqc_config.yaml"
picard220 = "/opt/picard-2.20.1/picard.jar"
starAlign = "/usr/local/bin/star_2-7-3a/STAR"
gatk = "/opt/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"
polyXjar = "/jvarkit/dist/vcfpolyx.jar"

// Check the analysis level 
if (params.level == 1)
{

	println ("**************** INFOS : Running Level 1 analysis : (Demultiplexer + Stats) *******************")

}

else if (params.level == 2)
{

	println ("**************** INFOS : Running Level 2 analysis : (Demultiplexer + Counts + Stats) *******************")
}

else if (params.level == 3)
{

	println ("**************** INFOS : Running Level 3 analysis : (Demultiplexer + Counts + SNP + Stats) *******************")
}


///////////////////////////////////////////////////
//  _ __  _ __ ___   ___ ___  ___ ___  ___  ___  //
// | '_ \| '__/ _ \ / __/ _ \/ __/ __|/ _ \/ __| //
// | |_) | | | (_) | (_|  __/\__ \__ \  __/\__ \ //
// | .__/|_|  \___/ \___\___||___/___/\___||___/ //
// |_|                                           //
///////////////////////////////////////////////////

// // A general process:
// process < name > {
//    [ directives ]
//    input:
//     < process inputs >
//    output:
//     < process outputs >
//    when:
//     < condition >
//    [script|shell|exec]:
//    < user script to be executed >
// }

log.info ""
log.info ""
log.info "~~~~~~~ Nextflow - ANGES Pipeline ~~~~~~~"
log.info "* User name:       ${params.userName}"
log.info "* Run ID:          ${params.run_no}"
log.info "* Read type :      Paired_End Sequencing"
log.info "* Sequencer dir:   ${params.sequencer}"
log.info "* Output dir:      ${toOutput}"
log.info "* scriptFile:      ${workflow.scriptFile}"
log.info "* Launch dir:      ${workflow.launchDir}"
log.info "* Work dir:        ${workflow.workDir}"
log.info " "
log.info "* Cmd line:        $workflow.commandLine"
log.info ""
log.info ""
log.info ""

// Create the Arborescence
process createDirs {
	cache "lenient"

	publishDir "${toLogsSteps}/", pattern: "*.ok"

	input:
		val name from listPaths.join(' ')
	
	output:
		file 'createDirs.ok' into (directories_created1, directories_created2, directories_created3)
	
	script:
		"""
		mkdir -p ${name}

		touch createDirs.ok

		"""
}

process getInfo {
	cache "lenient"

	tag "systemInfo"
	publishDir "${toLogsSteps}/", pattern: "*.ok"
	publishDir "${toLogs}/", pattern: "*.txt"
	
	output:
		file('getInfo.ok') into system_information
	
	script:
		"""
		uname -a > systemInfo.txt
		lscpu >> systemInfo.txt

		touch getInfo.ok
		"""
}

process fastqc {
	echo true
	cache "lenient"
	tag "${sampleID}"
	publishDir "${toTools}/", pattern: "*_version.txt"
	publishDir "${toLogsSteps}/", pattern: "*fastQC.ok"
	publishDir "${toFastQC}/", pattern: "*fastqc.zip"
	publishDir "${toFastQC}/", pattern: "*fastqc.html"
	publishDir "${toLogsFastQC}/", pattern: "*fastQC.log"
	

	input:
	file ok from directories_created1
	val sampleID from samplesList
	file fastq from myFileChannel.collect()

	output:
		file("${sampleID}_fastQC.ok") into fastQCOK
		file("*fastqc.html") into (fastQC_html, fastQC_html)
		file("*fastqc.zip") into (fastQC_zip, fastQC_zip)
		file("*.log") into (fastQC_log, fastQC_log)
		file ("fastqc_version.txt") into (fastQC_version, fastQC_version)
		file("*unzip_r*.out") into (unzip_r1_log)

	script:
	if (params.level == 1 && read1 == read2)

		"""
		fastqc *fastq.bz2 --outdir=. -t ${task.cpus} &> fastQC.log
		fastqc --version &> fastqc_version.txt
		
		unzip ${sampleID}_*R1_fastqc.zip > ${sampleID}_unzip_r1.out
		unzip ${sampleID}_*R2_fastqc.zip > ${sampleID}_unzip_r2.out
		
		mv -f ${sampleID}_*R1_fastqc ${toReportingFastqc}
		mv -f ${sampleID}_*R2_fastqc ${toReportingFastqc}
		touch ${sampleID}_fastQC.ok

		"""

	else if (params.level == 1 && read1 != read2)

		"""
		fastqc *fastq.bz2 --outdir=. -t ${task.cpus} &> fastQC.log
		fastqc --version &> fastqc_version.txt
		
		unzip ${sampleID}_*R1_fastqc.zip > ${sampleID}_unzip_r1.out
		
		mv -f ${sampleID}_*R1_fastqc ${toReportingFastqc}
		touch ${sampleID}_fastQC.ok

		"""

	else if (params.level == 2 || 3 && read1 == read2)

		"""
		fastqc ${sampleID}*.fastq --outdir=. -t ${task.cpus} &> fastQC.log
		fastqc --version &> fastqc_version.txt
		
		unzip ${sampleID}_R1_fastqc.zip > ${sampleID}_unzip_r1.out
		unzip ${sampleID}_R2_fastqc.zip > ${sampleID}_unzip_r2.out
		
		mv -f ${sampleID}_R1_fastqc ${toReportingFastqc}
		mv -f ${sampleID}_R2_fastqc ${toReportingFastqc}
		touch ${sampleID}_fastQC.ok

		"""
	else if (params.level == 2 || 3 && read1 != read2)

		"""
		fastqc *Cut_Pair.fastq.bz2 --outdir=. -t ${task.cpus} &> fastQC.log
		fastqc --version &> fastqc_version.txt
		
		unzip ${sampleID}_R1_Cut_Pair_fastqc.zip > ${sampleID}_unzip_r1.out
		
		mv -f ${sampleID}_R1_Cut_Pair_fastqc ${toReportingFastqc}
		touch ${sampleID}_fastQC.ok

		"""
}

///Process Align
process rsem {
	cache "lenient"
	tag "${sampleID}"
	publishDir "${toLogsSteps}/", pattern: "*rsem.ok"
	publishDir "${toTools}/", pattern: "*version.txt"
	publishDir "${toRef}/", pattern: "GenomeRef.txt"
	publishDir "$toPreBAM/", pattern: "*.genome.cram*"
	publishDir "$toBAM/", pattern: "*sorted.cram*"
	publishDir "${toRSEM}/", pattern: "*.results"
	publishDir "${toLogsRSEM}/", pattern: "rsem.command*"
	publishDir "${toLogsRSEM}/", pattern: "*biomart_genename.log"

	when :
	(params.level == 2 || params.level == 3)
		
	input:
		val sampleID from samplesList
		file fileCutRsemR1 from myFileChannelRSEM.collect()
		file ok from directories_created2
	
	output:
		file("*.results*") into resultsRSEMaddGenome
		file("*.cram*") into resultsRSEMbam
		file("*.cram.crai") into cram_sort_index
		file("GenomeRef.txt") into genomeRefVersion
		file("*rsem.ok") into (rsemMultiQCOK, rsemStatsCOK, rsemMultiSexQCOK)
		file("rsem.command*") into rsemCommandLogs
		file("*.log") into biomartLog
		file("*_version.txt") into toolsVersionRsem

	script:
	if (read1 == read2)

		"""
		$rsem/rsem-calculate-expression -p ${task.cpus} --star --star-path $star --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end ${sampleID}_R1.fastq ${sampleID}_R2.fastq $genomeDA ${sampleID}

		$samtools view -T $genomebamToCram ${sampleID}.genome.sorted.bam -C -S -@ 5 -o ${sampleID}.genome.sorted.cram
		$samtools view -T $genomeBamTranscritToCram ${sampleID}.transcript.sorted.bam -C -S -@ 5 -o ${sampleID}.transcript.sorted.cram

		$samtools index ${sampleID}.genome.sorted.cram
		$samtools index ${sampleID}.transcript.sorted.cram

		rm -rf ${sampleID}.genome.bam ${sampleID}.genome.sorted.bam ${sampleID}.genome.sorted.bam.bai ${sampleID}.transcript.bam ${sampleID}.transcript.sorted.bam ${sampleID}.transcript.sorted.bam.bai

		yes | cp -ru ${sampleID}.temp/${sampleID}Log.final.out ${toReportingRSEM}

		awk '{split (\$1,a,"."); print a[1] "\\t" \$0}' ${sampleID}.genes.results > ${sampleID}.genes2.results && mv ${sampleID}.genes2.results ${sampleID}.genes.results
		awk '{split (\$1,a,"."); print a[1] "\\t" \$0}' ${sampleID}.isoforms.results > ${sampleID}.isoforms2.results && mv ${sampleID}.isoforms2.results ${sampleID}.isoforms.results
		
		Rscript --vanilla $rsem_addgenename ${sampleID} $host $dataset &> rsem_biomart_genename.log	
		echo $genome > GenomeRef.txt

		$rsem/rsem-calculate-expression --version &> rsem_version.txt
		$star/STAR --version &> starRSEM_version.txt
		samtools --version > samtools_version.txt
		Rscript --version > Rscript_version.txt

		echo $host > Ensembl_Biomart_version.txt

		mv .command.log rsem.command.log
		mv .command.err rsem.command.err
		touch ${sampleID}_rsem.ok

		"""
	else

		"""
		fnameR1=\$(echo ${fileFastqRsemR1} | cut -d'_' -f1)

		$rsem/rsem-calculate-expression -p ${task.cpus} --star --star-path $star --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files ${sampleID}_R1.fastq $genomeDA ${sampleID}

		$samtools view -T $genomebamToCram ${sampleID}.genome.sorted.bam -C -S -@ 5 -o ${sampleID}.genome.sorted.cram
		$samtools view -T $genomeBamTranscritToCram ${sampleID}.transcript.sorted.bam -C -S -@ 5 -o ${sampleID}.transcript.sorted.cram

		$samtools index ${sampleID}.genome.sorted.cram
		$samtools index ${sampleID}.transcript.sorted.cram

		rm -rf ${sampleID}.genome.bam ${sampleID}.genome.sorted.bam ${sampleID}.genome.sorted.bam.bai ${sampleID}.transcript.bam ${sampleID}.transcript.sorted.bam ${sampleID}.transcript.sorted.bam.bai

		yes | cp -ru ${sampleID}.temp/${sampleID}Log.final.out ${toReportingRSEM}

		awk '{split (\$1,a,"."); print a[1] "\\t" \$0}' ${sampleID}.genes.results > ${sampleID}.genes2.results && mv ${sampleID}.genes2.results ${sampleID}.genes.results
		awk '{split (\$1,a,"."); print a[1] "\\t" \$0}' ${sampleID}.isoforms.results > ${sampleID}.isoforms2.results && mv ${sampleID}.isoforms2.results ${sampleID}.isoforms.results
		#Rscript --vanilla $rsem_addgenename ${sampleID} $host $dataset &> rsem_biomart_genename.log

		echo $genome > GenomeRef.txt

		$rsem/rsem-calculate-expression --version &> rsem_version.txt
		$star/STAR --version &> starRSEM_version.txt
		samtools --version > samtools_version.txt
		Rscript --version > Rscript_version.txt

		echo $host > Ensembl_Biomart_version.txt

		mv .command.log rsem.command.log
		mv .command.err rsem.command.err
		touch ${sampleID}_rsem.ok

		"""
		}

// launch mapping VS ncRNA
process bowtie2 {
	cache "lenient"
	tag "${sampleID}"
	publishDir "${toTools}/", pattern: "*version.txt"
	publishDir "${toLogsSteps}/", pattern: "bowtie2.ok"
	publishDir "${toFastQncrna}/", pattern: "*.minusNCRNA.fastq.bz2"

	when :
	(params.level == 2 || params.level == 3)
	
	input:
		file ok from directories_created3
		val sampleID from samplesList
		file fileCutBowtieR1 from myFileChannelbowtie.collect()

	output:
		file("*_R1.minusNCRNA.fastq.bz2") into fastq_minusNCRNA_r1
		file("*_R2.minusNCRNA.fastq.bz2") optional true into fastq_minusNCRNA_r2
		file('bowtie2.ok') into (bowtieOK, bowtieMQOK, bowtieStatOK)
		file("bowtie2_version.txt") into bowtie2Version

	script:
	if (read1 == read2)
		"""

		$bowtie2 -x $ncrna -U ${sampleID}_R1.fastq -S ${sampleID}_R1_NCRNA.sam \
		--un ${sampleID}_R1.minusNCRNA.fastq > ${sampleID}_R1_mapNCRNA_summary.txt 2>&1
		
		$bowtie2 -x $ncrna -U ${sampleID}_R2.fastq -S ${sampleID}_R2_NCRNA.sam \
		--un ${sampleID}_R2.minusNCRNA.fastq > ${sampleID}_R2_mapNCRNA_summary.txt 2>&1
		
		lbzip2 -f -n ${task.cpus} --best ${sampleID}*.fastq

		rm -rf ${sampleID}_R1_NCRNA.sam ${sampleID}_R2_NCRNA.sam
		mv -f *_mapNCRNA_summary.txt ${toReportingBowtie}
		$bowtie2 --version &> bowtie2_version.txt
		touch bowtie2.ok

		"""
	else
		"""
		fnameR1=\$(echo ${fileCutBowtieR1} | cut -d'_' -f1)

		$bowtie2 -x $ncrna -U ${sampleID}_R1.fastq -S ${sampleID}_R1_NCRNA.sam \
		--un ${sampleID}_R1.minusNCRNA.fastq > ${sampleID}_R1_mapNCRNA_summary.txt 2>&1
		
		lbzip2 -f -n ${task.cpus} --best ${sampleID}*.fastq

		rm -rf ${sampleID}_R1_NCRNA.sam
		mv -f *_mapNCRNA_summary.txt ${toReportingBowtie}
		$bowtie2 --version &> bowtie2_version.txt
		touch bowtie2.ok

		"""
}


process multiQC {
	tag "${params.run_no}"
	cache "lenient"
	tag " MultiQC "
	publishDir "${toLogsSteps}/", pattern: "*.ok"
	publishDir "${toMultiFQC}/", pattern: "*.html"
	publishDir "${toTools}/", pattern: "*.txt"
	publishDir "${toLogsMultiQC}/", pattern: "*.log"

	when : 
	(params.level == 1 || params.level == 2 || params.level == 3)
	
	input:
		file fastQCok from fastQCOK.collect().ifEmpty([])
		file("rsem.ok") from rsemMultiQCOK.collect().ifEmpty([])
		file("bowtie2.ok") from bowtieOK.collect().ifEmpty([])
		file("VEP_filter.ok") from filterVepAnnotSnpOK.collect().ifEmpty([])
		file("VEP_INDELs_filter.ok") from filterVepAnnotIndelOK.collect().ifEmpty([])

	output:
		file("*.html") into multiQCComplete
		file("multiqc.ok") into multiQCOK
		file("*.log") into multiQCLog

	script:
		"""
		multiqc ${to} -n ${params.run_no} -c ${multiqcconf}  > multiQC.log 2>&1

		multiqc --version > multiqc.txt
		touch multiqc.ok
		"""
}

process stats {
	tag "${params.run_no}"
	cache "lenient"
	echo = true
	publishDir "${toLogsSteps}/", pattern: "*.ok", mode : "copy"

	input:
		file ("multiqc.ok") from multiQCOK.ifEmpty([])
		file ('bowtie2.ok') from bowtieStatOK.collect().ifEmpty([])
		file ok from rsemStatsCOK.collect().ifEmpty([])

	output:	
		file('stats.ok') optional true into statsOK_L1
		file('stats.ok') optional true into (statsOKSexVerif, statsOKReplaceLink)

	script:

	if (params.level == 1 && read1 == read2)

		"""
		echo -e sample\ttotal_reads\t%R1_dup\t%R2_dup > ${toReporting}/${params.run_no}_stats.csv
		echo $samplesList
		for sample in $samplesListBash
		do

		if [ -f "${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt" ] && [ -f "${toReportingFastqc}/\$sample"_R2_fastqc"/fastqc_data.txt" ]; then
		echo "All files of sample \$sample for the stats exist and not empty"

		total_reads=\$(grep -w "Total Sequences" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{print \$3*2}')

		r1_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')
		r2_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R2_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')

		echo -e \$sample\t\$total_reads\t\$r1_dup\t\$r2_dup>> ${toReporting}/${params.run_no}_stats.csv
		else
		echo One or more files for the stats are not exist or empty
		exit
		fi
		done

		touch statsL1.ok

		"""
	else if (params.level == 1 && read1 != read2)

		"""
		echo -e sample"\t"total_reads"\t"%R1_dup > ${toReporting}/${params.run_no}_stats.csv
		echo $samplesList
		for sample in $samplesListBash
		do

		if [ -f "${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt" ]; then
		echo "All files of sample \$sample for the stats exist and not empty"

		total_reads=\$(grep -w "Total Sequences" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{print \$3}')

		r1_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')

		echo -e \$sample "\t"\$total_reads"\t"\$r1_dup >> ${toReporting}/${params.run_no}_stats.csv
		else
		echo One or more files for the stats are not exist or empty
		exit
		fi
		done

		touch stats.ok

		"""

 	else if (params.level == 2 || 3 && read1 == read2)

		"""
		echo -e sample"\\t"total_reads"\\t"%R1_dup"\\t"%R2_dup"\\t"total_map_reads"\\t"%R1_map_ncRNA"\\t"%R2_map_ncRNA > ${toReporting}/${params.run_no}_stat.tsv
		echo $samplesList
		for sample in $samplesListBash
		do

		if [ -f "${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt" ] && [ -f "${toReportingFastqc}/\$sample"_R2_fastqc"/fastqc_data.txt" ] && [ -f ${toReportingRSEM}/\$sample"Log.final.out" ] && [ -f ${toReportingBowtie}/\$sample"_R1_mapNCRNA_summary.txt" ] && [ -f ${toReportingBowtie}/\$sample"_R2_mapNCRNA_summary.txt" ]; then
		echo "All files of sample \$sample for the stats exist and not empty"

		total_reads=\$(grep -w "Total Sequences" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{print \$3*2}')

		r1_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R1_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')
		r2_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R2_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')

		count1_map_RNA=\$(grep -w "Uniquely mapped reads %" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$2}')
		count2_map_RNA=\$(grep -w "% of reads mapped to multiple loci" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$8}')
		count3_map_RNA=\$(grep -w "% of reads mapped to too many loci |" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$9}')
		total_map_reads=\$(echo "(\$count1_map_RNA+\$count2_map_RNA+\$count3_map_RNA)" | bc)

		r1_map_ncRNA=\$(grep -w "overall alignment rate" ${toReportingBowtie}/\$sample"_R1_mapNCRNA_summary.txt" | awk '{split(\$0,a, "%"); print a[1]}')
		r2_map_ncRNA=\$(grep -w "overall alignment rate" ${toReportingBowtie}/\$sample"_R2_mapNCRNA_summary.txt" | awk '{split(\$0,a, "%"); print a[1]}')

		echo -e \$sample "\\t"\$total_reads"\\t"\$r1_dup "\\t" \$r2_dup "\\t" \$total_map_reads"\\t" \$r1_map_ncRNA "\\t" \$r2_map_ncRNA >> ${toReporting}/${params.run_no}_stat.tsv
		else
		echo " One or more files for the stats are not exist or empty"
		exit
		fi
		done

		touch stats.ok

		"""
	else if (params.level == 2 || 3 && read1 != read2)
		
		"""
		echo -e sample"\\t"total_reads"\\t"%R1_dup"\\t"total_map_reads"\\t"%R1_map_ncRNA > ${toReporting}/${params.run_no}_stat.tsv
		echo $samplesList
		for sample in $samplesListBash
		do

		if [ -f "${toReportingFastqc}/\$sample"_R1_Cut_fastqc"/fastqc_data.txt" ] && [ -f ${toReportingRSEM}/\$sample"Log.final.out" ] && [ -f ${toReportingBowtie}/\$sample"_R1_mapNCRNA_summary.txt" ]; then
		echo "All files of sample \$sample for the stats exist and not empty"

		total_reads=\$(grep -w "Total Sequences" ${toReportingFastqc}/\$sample"_R1_Cut_fastqc"/fastqc_data.txt | awk '{print \$3}')

		r1_dup=\$(grep -w "#Total Deduplicated Percentage" ${toReportingFastqc}/\$sample"_R1_Cut_fastqc"/fastqc_data.txt | awk '{printf "%.2f\\n", 100-\$4}')

		count1_map_RNA=\$(grep -w "Uniquely mapped reads %" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$2}')
		count2_map_RNA=\$(grep -w "% of reads mapped to multiple loci" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$8}')
		count3_map_RNA=\$(grep -w "% of reads mapped to too many loci |" ${toReportingRSEM}/\$sample"Log.final.out" | awk '{split(\$0,a, "%"); print a[2]}' | awk '{print \$9}')
		total_map_reads=\$(echo "(\$count1_map_RNA+\$count2_map_RNA+\$count3_map_RNA)" | bc)

		r1_map_ncRNA=\$(grep -w "overall alignment rate" ${toReportingBowtie}/\$sample"_R1_mapNCRNA_summary.txt" | awk '{split(\$0,a, "%"); print a[1]}')

		echo -e \$sample "\\t"\$total_reads"\\t"\$r1_dup "\\t" \$total_map_reads"\\t" \$r1_map_ncRNA >> ${toReporting}/${params.run_no}_stat.tsv
		else
		echo " One or more files for the stats are not exist or empty"
		exit
		fi
		done

		touch stats.ok

		"""
}

process replaceSymLinks {
	cache "lenient"
	tag "${params.run_no}"
	tag "Replace all Symbolic Links with files"
	
	when : 
	(params.level == 1 || params.level == 2 || params.level == 3)

	input:
	 	file("statsL1.ok") from statsOK_L1.ifEmpty([])
	 	file("stats.ok") from statsOKReplaceLink.ifEmpty([])

	output:
		file('removeLinks.ok') into replaceLinksOKstat

	script:
	"""
	for link in \$(find ${to} -type l)
	do
	  loc="\$(dirname "\$link")"
	  dir="\$(readlink "\$link")"
	  mv -f "\$dir" "\$loc"
	done

	
	touch removeLinks.ok
	cp -rp removeLinks.ok ${toLogsSteps}/
	"""
}

workflow.onComplete {

    sendMail {
    to 'souhila.amanzougarene@cnrs.fr'
    subject "The job of ${params.run_no} is finished -> Success : ${workflow.success}"
    """
    Pipeline execution summary
    ****************************************************************
    Run ID 			:  ${params.run_no}
    Sample size 	:  ${sample_size}
    User name 		:  ${params.userName}
    Sequencer dir 	:  ${params.sequencer}
    Output dir 		:  ${toOutput}
    scriptFile 		:  ${workflow.scriptFile}
    Launch dir 		:  ${workflow.launchDir}
    Work dir 		:  ${workflow.workDir}
    Completed at 	:  ${workflow.complete}
    Duration 		:  ${workflow.duration}
    Success 		:  ${workflow.success}
    workDir 		:  ${workflow.workDir}
    exit status 	:  ${workflow.exitStatus}
    error message 	:  ${workflow.errorMessage}

    Cmd line 		:  $workflow.commandLine

    """
    }
}