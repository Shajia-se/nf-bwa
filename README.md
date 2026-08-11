# nf-bwa

`nf-bwa` aligns trimmed paired-end FASTQ files to a reference genome with BWA-MEM.

This module is normally the third step in the ChIP-seq pipeline:

```text
nf-fastqc -> nf-fastp -> nf-bwa -> nf-picard
```

`nf-fastp` produces trimmed FASTQ files. `nf-bwa` maps those reads and produces sorted/indexed BAM files for downstream BAM QC and filtering.

## What It Does

1. Checks that a reference FASTA was provided.
2. Creates or reuses a BWA index.
3. Finds trimmed paired-end FASTQ files from the `nf-fastp` output folder.
4. Runs `bwa mem -M`.
5. Converts alignments to BAM, sorts, indexes, and writes a mapping summary.
6. Saves only final BAM/BAI/QC outputs, not large SAM intermediate files.

## Before You Run

You need:

- Nextflow installed
- Trimmed paired-end FASTQ files from `nf-fastp`
- A reference FASTA file
- Write permission to the BWA index folder, unless the index already exists
- For local runs: Docker available
- For HPC runs: Slurm and Singularity available
- For HPC runs: the BWA/Samtools Singularity image path in `configs/slurm.config` must exist

Default HPC notification email:

```text
molendo.hpc@gmail.com
```

## Input Options

The input FASTQ files should usually come from `nf-fastp`.

Expected filenames:

```text
sample_R1.fastp.trimmed.fastq.gz
sample_R2.fastp.trimmed.fastq.gz
```

### Option 1: Samples Master CSV Recommended

Use this when running the full ChIP-seq pipeline.

The CSV should contain at least:

```csv
sample_id,enabled
WT_rep1,true
WT_rep2,true
```

Notes:

- `sample_id` must match the fastp output prefix.
- For sample `WT_rep1`, this module expects:
  - `WT_rep1_R1.fastp.trimmed.fastq.gz`
  - `WT_rep1_R2.fastp.trimmed.fastq.gz`
- Rows with `enabled=false` are skipped.
- Empty `enabled` values are treated as enabled.

Run with:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bwa_raw_data /path/to/fastp_output \
  --reference_fasta /path/to/reference.fa \
  --project_folder /path/to/output_project \
  --bwa_output bwa_output
```

### Option 2: FASTQ Folder + Pair Pattern

Use this for quick module testing.

```bash
nextflow run main.nf -profile hpc \
  --bwa_raw_data /path/to/fastp_output \
  --bwa_pattern "*_R{1,2}.fastp.trimmed.fastq.gz" \
  --reference_fasta /path/to/reference.fa \
  --project_folder /path/to/output_project \
  --bwa_output bwa_output
```

If `--index_dir` is not provided, the module creates or reuses an index folder next to the reference FASTA:

```text
/path/to/primary_bwa/
```

To use a specific index directory:

```bash
--index_dir /path/to/reference_bwa_index
```

## Output

Results are written to:

```text
${project_folder}/${bwa_output}/
```

Example:

```text
/path/to/output_project/bwa_output/
  WT_rep1.bam.stat
  WT_rep1.sorted.bam
  WT_rep1.sorted.bam.bai
```

Output files:

- `*.sorted.bam`: coordinate-sorted alignment file
- `*.sorted.bam.bai`: BAM index
- `*.bam.stat`: `samtools flagstat` mapping summary

Large intermediate SAM files are not published.

## Recommended HPC Run

From inside the `nf-bwa` folder:

```bash
cd /path/to/nf-bwa

nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bwa_raw_data /path/to/fastp_output \
  --reference_fasta /path/to/reference.fa \
  --project_folder /path/to/output_project \
  --bwa_output bwa_output
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc -resume \
  --samples_master /path/to/samples_master.csv \
  --bwa_raw_data /path/to/fastp_output \
  --reference_fasta /path/to/reference.fa \
  --project_folder /path/to/output_project \
  --bwa_output bwa_output
```

Override the HPC notification email:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --bwa_raw_data /path/to/fastp_output \
  --reference_fasta /path/to/reference.fa \
  --project_folder /path/to/output_project \
  --mail_user molendo.hpc@gmail.com
```

## Local Test Run

Use local mode only for small test data:

```bash
nextflow run main.nf -profile local \
  --bwa_raw_data /path/to/test_fastp_output \
  --bwa_pattern "*_R{1,2}.fastp.trimmed.fastq.gz" \
  --reference_fasta /path/to/reference.fa \
  --project_folder ./test_output
```

## Key Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `samples_master` | CSV used to select enabled sample IDs | `null` |
| `bwa_raw_data` | Input folder containing fastp-trimmed FASTQ pairs | `null` |
| `bwa_pattern` | Paired FASTQ filename pattern | `*_R{1,2}.fastp.trimmed.fastq.gz` |
| `reference_fasta` | Reference genome FASTA | `null` |
| `index_dir` | BWA index directory | next to reference FASTA |
| `project_folder` | Base output folder | current folder |
| `bwa_output` | BWA output subfolder | `bwa_output` |
| `cpus` | CPUs per mapping task | `8` |
| `memory` | Memory per task | `16GB` |
| `time` | Runtime limit per task | `6h` |
| `mail_user` | HPC notification email | `molendo.hpc@gmail.com` |

## Existing Results Are Skipped

For each sample, the module checks whether all three expected output files already exist:

```text
sample.sorted.bam
sample.sorted.bam.bai
sample.bam.stat
```

If all three exist, that sample is skipped. If any file is missing, BWA runs again for that sample.

## How To Check Results

Check the mapping summary:

```bash
less ${project_folder}/${bwa_output}/sample.bam.stat
```

Important values to inspect:

- total reads
- mapped reads and mapping percentage
- properly paired reads
- singleton reads

The sorted BAM and BAI files are the inputs for the next module, usually `nf-picard`.

## Troubleshooting

If the run fails:

1. Check the main Nextflow log:

```bash
less .nextflow.log
```

2. Check the failed task error file:

```bash
less work/<hash>/.command.err
```

3. Common problems:

- `--bwa_raw_data` does not point to the `nf-fastp` output folder.
- `sample_id` in `samples_master` does not match the fastp output prefix.
- `--reference_fasta` is missing or points to the wrong genome build.
- The BWA index directory is not writable during first-time index creation.
- `configs/slurm.config` points to a missing Singularity image.
- The HPC bind path in `extra_mounts` does not include the FASTQ, reference, index, or output location.
- Docker is not running for local mode.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
