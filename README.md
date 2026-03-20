# nf-bwa

`nf-bwa` is a Nextflow DSL2 module for aligning trimmed FASTQ pairs with BWA-MEM and generating sorted/indexed BAM files.

## What This Module Does

1. Ensures reference index exists in a versioned reference directory such as `.../bwa/GRCm39_vM27/primary_bwa/` (`bwa index`).
2. Reads paired FASTQ files (priority: `samples_master` > `bwa_raw_data + bwa_pattern`).
3. Aligns each pair with `bwa mem -M`.
4. Converts, summarizes, sorts, and indexes alignments:
   - `samtools view`
   - `samtools flagstat`
   - `samtools sort`
   - `samtools index`
5. Skips samples that already have `${sample}.sorted.bam.bai`.

## Input

Preferred mode (`samples_master`):
- `params.samples_master`: CSV with at least `sample_id` (optional `enabled`)
- Module resolves input FASTQ as:
  - `${bwa_raw_data}/${sample_id}_R1.fastp.trimmed.fastq.gz`
  - `${bwa_raw_data}/${sample_id}_R2.fastp.trimmed.fastq.gz`

Fallback mode (pattern-based):
- Directory: `params.bwa_raw_data`
- Pattern: `params.bwa_pattern` (default: `*_R{1,2}.fastp.trimmed.fastq.gz`)
- Expected reads: paired-end trimmed FASTQ

Reference-related inputs:
- `params.reference_root`
- `params.genome_release`
- `params.gencode_version`
- `params.reference_fasta`
- `params.index_dir`

## Output

Under `${project_folder}/${bwa_output}`:
- `${sample}.sam`
- `${sample}.bam`
- `${sample}.bam.stat`
- `${sample}.sorted.bam`
- `${sample}.sorted.bam.bai`

Downstream modules typically use:
- `${sample}.sorted.bam`
- `${sample}.sorted.bam.bai`

## Key Parameters

- `bwa_raw_data`: input FASTQ folder
- `samples_master`: preferred sample list input
- `bwa_pattern`: paired FASTQ matching pattern
- `bwa_output`: output folder name
- `reference_root`, `genome_release`, `gencode_version`: versioned reference layout
- `reference_fasta`: FASTA used to prepare/link `index.fa`
- `index_dir`: where BWA index files are created and reused
- `cpus`, `memory`, `time`: process resources

## Run

```bash
nextflow run main.nf -profile local
```

```bash
nextflow run main.nf -profile hpc
```

Custom example:
Recommended run (`samples_master`):
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bwa_raw_data /path/to/nf-fastp/fastp_output
```

Fallback run (pattern):

```bash
nextflow run main.nf -profile hpc \
  --bwa_raw_data /path/to/fastp_output \
  --bwa_pattern "*_R{1,2}.fastp.trimmed.fastq.gz" \
  --reference_fasta /ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/bwa/GRCm39_vM27/GRCm39.primary_assembly.genome.fa \
  --index_dir /ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/bwa/GRCm39_vM27/primary_bwa
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```

## Notes

- First run can be slow due to index creation.
- If `index.fa.bwt` already exists, index step is skipped.
- Recommended layout is versioned, for example:
  - FASTA: `/ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/bwa/GRCm39_vM27/GRCm39.primary_assembly.genome.fa`
  - index: `/ictstr01/groups/idc/projects/uhlenhaut/jiang/reference/bwa/GRCm39_vM27/primary_bwa/`
- Keep treatment and control samples processed with the same upstream settings before MACS3.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
