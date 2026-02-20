# nf-bwa

`nf-bwa` is a Nextflow DSL2 module for aligning trimmed FASTQ pairs with BWA-MEM and generating sorted/indexed BAM files.

## What This Module Does

1. Ensures reference index exists in `.../toplevel_bwa/` (`bwa index`).
2. Reads paired FASTQ files from `bwa_raw_data` using `bwa_pattern`.
3. Aligns each pair with `bwa mem -M`.
4. Converts, summarizes, sorts, and indexes alignments:
   - `samtools view`
   - `samtools flagstat`
   - `samtools sort`
   - `samtools index`
5. Skips samples that already have `${sample}.sorted.bam.bai`.

## Input

- Directory: `params.bwa_raw_data`
- Pattern: `params.bwa_pattern` (default: `*_R{1,2}.fastp.trimmed.fastq.gz`)
- Expected reads: paired-end trimmed FASTQ

Reference-related inputs:
- `params.genomes`
- `params.organism`
- `params.release`
- `params.reference_fasta`

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
- `bwa_pattern`: paired FASTQ matching pattern
- `bwa_output`: output folder name
- `genomes`, `organism`, `release`: reference directory layout
- `reference_fasta`: FASTA used to prepare/link `index.fa`
- `cpus`, `memory`, `time`: process resources

## Run

```bash
nextflow run main.nf -profile local
```

```bash
nextflow run main.nf -profile hpc
```

Custom example:

```bash
nextflow run main.nf -profile hpc \
  --bwa_raw_data /path/to/fastp_output \
  --bwa_pattern "*_R{1,2}.fastp.trimmed.fastq.gz" \
  --reference_fasta /path/to/Mus_musculus.GRCm39.dna.toplevel.fa
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```

## Notes

- First run can be slow due to index creation.
- If `index.fa.bwt` already exists, index step is skipped.
- Keep treatment and control samples processed with the same upstream settings before MACS3.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
