#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

	===================================================================
	 L A R G E   I N D E L S   P I P E L I N E  - I R Y C I S    v 0.1
	===================================================================

	Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow_LargeIndels.nf [OPTIONS]

    Options:
      --indir [DIR]                           The input directory, all fastq files or csv files in this directory will be processed. (default: "data")
      --outdir [DIR]                          The output directory where the results will be saved (default: "my-results")
      --genome <GRCh37 | GRCh38 | [FILE]>     Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
      --adapter_file [FILE]                   Adapter file to trimm reads by. (Trimmomatic adapters provided in $baseDir/adapter)
      --region_intervals [BED FILE]           Complete path to specific genomic region in .list format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)
      --paired <true | false>                 Execute pipleine in single-end or paired-end mode. If "--paired true" then all fastq files in $params.indir will be processed as samples from the same experiment.
                                              If "--paired false" a csv with the single-end files path and their IDs will be used to identify the fasq files. Options: true, false (default: true)
      --min_alt_fraction [NUM]                Freebayes specific option, minimumn threshold at which allele frequency is considered real. (default: 0.2)

      --reads [GLOB]                          Glob pattern to identify paired-end reads or the single-end csv file. All data must be compressed. (default: "$baseDir/data/*R{1,2}*.fastq.gz (if paired) "$baseDir/data/*.csv" (if single-end))
                                              CSV format: SampleID, [path/to/read].fastq
      --remove_duplicates <true | false>      Remove marked as duplicated reads. Options: true, false (default: false)

      """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if(params.paired){

  //Reads must be read twice, both as a trupple and as an array
   Channel
      .fromFilePairs(params.reads, size: 2,  checkIfExists: true, flat:true)
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .into{ [ch_pre_samples, ch_FlowCell_lane] }


   //Make ch_pre_samples a tuple to process SampleId in processes more easily and to include it in ReadGroup (RG) info.
   ch_pre_samples.map{it -> new Tuple(it[0].split("_")[0],it[1,2])}.set{ ch_samples_with_id }

   //Identify fastq data for Read group definition "ID". RG:ID = {flowcell}.{lane}.{uniqueId} [must be unique despite documentation stating otherwise]
   ch_FlowCell_lane.splitFastq(record: true, pe: true, limit: 1).map{it -> new Tuple(it[0].split("_")[0], it[1].readHeader.split(":")[2,3,9].join("."))}.set{ch_RG_ID}

}else{

   Channel
      .fromPath(params.reads, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> [row.sampleId, [row.read]] }
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into{ [ch_samples, ch_sampleName, ch_FlowCell_lane] }

   ch_FlowCell_lane.splitFastq(record: true, pe: true, limit: 1).map{it -> new Tuple(it[0].split("_")[0], it[1].readHeader.split(":")[2,3,9].join("."))}.set{ch_RG_ID}

 }
//Fuse Ids and samples to manage SampleId as a tuple together with the samples they identify.
ch_RG_ID.concat(ch_samples_with_id).groupTuple().map{ it -> [[it[0],it[1][0]],it[1][1]] }.set{ ch_samples }

ch_dbSNP = file(params.dbSNP)

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
================================================================
genome               : $params.genome
reads                : $params.reads
adapters             : $params.adapters

region_intervals     : $params.region_intervals
dbSNP                : ${params.dbSNP}

paired               : $params.paired
aligner              : minimap2
variant_caller       : transindel
remove_duplicates    : $params.remove_duplicates

read_directory       : $params.indir
results              : $params.outdir
===============================================================
"""

//------------------------------------------------------------Genome Reference Handling-------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

if(params.genome != "GRCh37" && params.genome != "GRCh38"){


  ch_reference = file(params.genome, checkIfExists: true)


  process Indexing_custom_genome {
    tag "Indexes supplied reference FASTA file (Samtools)"
    label 'med_mem'

    publishDir "data/custom_reference", mode: 'copy'


    input:
    file reference_file from ch_reference

    output:
    file("${reference_file[0]}")
    file '*.fai'
    
    """
    samtools faidx ${reference_file[0]}
    """
    
  }

   process Building_genome_dictionary {
    tag "Creating Dictionary file (GATK)"
    label 'med_mem'

    publishDir "data/custom_reference", mode: 'copy'

    input:
    file reference_file from ch_reference

    output:
    file '*.dict'

    script:
    
    """
    gatk CreateSequenceDictionary -R ${reference_file[0]}
    """
    
  }
}

  // Extract file name:
  custom_reference = file("$params.genome")
  prefixRef = custom_reference.name.take(custom_reference.name.lastIndexOf('.'))

def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/../data/custom_reference/${prefixRef}.fasta": file("${params.indexRef}")

//------------------------------------------------------------Trimming-----------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

def adapter_trimm = params.adapters != 'NO_FILE' ? "ILLUMINACLIP:${params.adapters}:2:30:10" : ''

process FASTQ_Trimming {

  tag "Quality checks and trimms reads using a trimmomatic"
  label 'med_mem'

  publishDir "$params.outdir/trimmed_data"

   input:
   set sampleId, file(samples) from ch_samples

   output:
   set sampleId, file('*.fastq.gz') into ch_alignment

   script:

   //WARNING: Trimmomatic does not mind if the adapter file is not found
   //it will continue without processing it and without a warning.
   if(params.paired){
      // def adapter_trimm does not work with template command.
      //template 'trimmomatic/trimmomatic_PE_adapter_test'

      """
      trimmomatic PE -threads ${params.threads} ${samples[0]} ${samples[1]} ${sampleId[0]}_R1.fastq.gz bad_1 ${sampleId[0]}_R2.fastq.gz bad_2 ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """

   }else{

      """
      trimmomatic SE -threads ${params.threads} ${samples[0]} ${sampleId[0]}_trimmed.fastq.gz ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """
   }
}

//------------------------------------------------------------Alignment-----------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

process Alignment {

tag "Aligns reads to a reference genome using minimap2"
label 'med_mem'
 
  publishDir "$params.outdir/alignment/minimap"

 input:
 set sampleId, file(fastq_file) from ch_alignment

 output:
 set sampleId, file('*.sam') into ch_sam_to_bam

 script:  

 if(params.paired){

    """
    minimap2 -t 4 -Y -ax sr ${reference} ${fastq_file[0]} ${fastq_file[1]} > ${sampleId[0]}.TraIn.sam
    """
  }else{

    """
    minimap2 -t 4 -Y -ax ${reference} ${fastq_file[0]} > ${sampleId[0]}.TraIn.sam
    """
  }
}

process SAM_to_BAM{

  tag "Converts SAM file to BAM file using samtools view"
  label 'med_mem'

  publishDir "$params.outdir/alignment/transindel"

  input:
  set sampleId, file(sam_file) from ch_sam_to_bam

  output:
  set sampleId, file('*sort.bam*'), file('*.bam.bai') into ch_bam_sorting

  script:

  """
    samtools view --threads ${params.threads} -b -o ${sampleId[0]}.TraIn.bam ${sam_file[0]}
    samtools sort --threads ${params.threads} -o ${sampleId[0]}.TraIn.sort.bam ${sampleId[0]}.TraIn.bam

    samtools index ${sampleId[0]}.TraIn.sort.bam    
  """
} 

def min_alt_fraction_var = params.min_alt_fraction == '' ? 0.2:"${params.min_alt_fraction}"

process transIndel_build_DNA {
  tag "transindel_build"
  publishDir "$params.outdir/alignment/transindel"

  input:
    set sampleId, file(bam_files) from ch_bam_sorting
  output:
    set sampleId, file('*.OutStep1.bam'), file('*.OutStep1.bam.bai') into ch_transindel_call

  script:
  """
    transIndel_build_DNA -i ${bam_files[0]} -o ${sampleId[0]}.OutStep1.bam

  """
}

process transIndel_call {
  tag "transindel_call"
  publishDir "$params.outdir/raw_variant_calling_files/transindel"

  input:
    set sampleId, file(bam_files) from ch_transindel_call
  output:
    set sampleId, file('*indel.vcf') into ch_transindel

  script:
  """
    transIndel_call -f ${min_alt_fraction_var} -i ${bam_files[0]} -o ${sampleId[0]}.OutStep2
  """
}


def correction = params.position_correction == 'default' ? 118:params.position_correction

process pos_corrector {
  tag "Corrects the position of a vcf file to match that of a custom reference. Default: 118"
  publishDir "$params.outdir/Large_indels_variant_calling_files", mode: 'copy'

  input:
    file(vcf_file) from ch_transindel.collect()
  output:
    file("*fix.vcf") into ch_pos_corrector

  script:
  """
  $baseDir/fix_pos_isoforma_corta.R ${vcf_file} ${correction}
  """
}

def annot = params.exon_annotation == '' ? '': file(params.exon_annotation)
if(annot != ''){

process exon_mapper {
  tag "Maps to the position of a variant to an exon"
  publishDir "$params.outdir/Large_indels_variant_calling_files/exons", mode: 'copy'

  input:
    file(vcf_file) from ch_pos_corrector.collect()
  output:
    file('*.exons.vcf')

  script:
  """
  $baseDir/custom_ref_exon_mapper.R ${vcf_file} ${annot}
  """
  }
}