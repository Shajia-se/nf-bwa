#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def bwa_output = params.bwa_output ?: "bwa_output"


process bwa_index {
  tag "${params.organism}.${params.release}"
  stageInMode  'symlink'
  stageOutMode 'move'

  output:
    val(true)

  script:
  """
  set -eux
  index_folder=${params.genomes}/${params.organism}/${params.release}/toplevel_bwa
  ref_fa=${params.reference_fasta}
  index_fa="\$index_folder/index.fa"
  mkdir -p "\$index_folder"
  cd "\$index_folder"

  if [[ -f index.fa.bwt ]]; then
    echo "BWA index already exists in \$index_folder, skipping indexing."
    exit 0
  fi

  [[ -s "\$ref_fa" ]] || { echo "Reference FASTA not found: \$ref_fa" >&2; exit 1; }
  if [[ "\$ref_fa" != "\$index_fa" ]]; then
    ln -sf "\$ref_fa" index.fa
  fi

  bwa index -a bwtsw -p index.fa index.fa
  """
}



process mapping {
  tag "${pair_id}"
  stageInMode  'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${bwa_output}", mode: 'copy'

  input:
    tuple val(pair_id), path(reads)
    val(ready_signal)

  output:
    path "${pair_id}.sam"
    path "${pair_id}.bam"
    path "${pair_id}.bam.stat"
    path "${pair_id}.sorted.bam"
    path "${pair_id}.sorted.bam.bai"

  script:
    def ref = "${params.genomes}/${params.organism}/${params.release}/toplevel_bwa/index.fa"

    if( reads instanceof Path ) {
      """
      set -eux

      bwa mem -t ${task.cpus} -M "${ref}" ${reads} > ${pair_id}.sam

      samtools view -bS ${pair_id}.sam > ${pair_id}.bam
      samtools flagstat ${pair_id}.bam > ${pair_id}.bam.stat
      samtools sort -@ ${task.cpus} -o ${pair_id}.sorted.bam ${pair_id}.bam
      samtools index ${pair_id}.sorted.bam
      """
    } else {
      """
      set -eux

      bwa mem -t ${task.cpus} -M "${ref}" ${reads[0]} ${reads[1]} > ${pair_id}.sam

      samtools view -bS ${pair_id}.sam > ${pair_id}.bam
      samtools flagstat ${pair_id}.bam > ${pair_id}.bam.stat
      samtools sort -@ ${task.cpus} -o ${pair_id}.sorted.bam ${pair_id}.bam
      samtools index ${pair_id}.sorted.bam
      """
    }
}

workflow {

  def outdir = "${params.project_folder}/${bwa_output}"
  def pattern = params.bwa_pattern ?: "*_R{1,2}.fastp.trimmed.fastq.gz"

  def read_pairs = Channel
    .fromFilePairs("${params.bwa_raw_data}/${pattern}", checkIfExists: true)
    .ifEmpty { exit 1, "ERROR: No FASTQ pairs found for pattern: ${params.bwa_raw_data}/${pattern}" }
    .filter { pair_id, reads ->
      ! file("${outdir}/${pair_id}.sorted.bam.bai").exists()
    }

  def index_ready = bwa_index()
  mapping(read_pairs, index_ready)
}
