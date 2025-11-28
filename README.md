# nf-bwa

A simple, portable BWA-MEM alignment pipeline using Nextflow.
Designed for aligning trimmed FASTQ files (e.g., from nf-fastp)

---

## 📁 Demo data

* Loads trimmed FASTQs: `*1.fastp.trimmed.fastq.gz`, `**1.fastp.trimmed.fastq.gz`


---

## 🧬 Reference download (Ensembl - primary_assembly or toplevel)

**Mouse GRCm39 & GRCm38**

```bash
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
```

**Human GRCh38**

```bash
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

Unzip:

```bash
gunzip *.fa.gz
```

---

## 🚀 Run locally (Docker)

```bash
nextflow run main.nf -profile local
```

**Requirements**

* Nextflow
* Docker
* A combined BWA+samtools image, e.g.:

```bash
docker pull shajiase/bwa-sam:1.16.1
```

* `configs/local.config`

---

## 🚀 Run on HPC (Slurm + Singularity)

```bash
nextflow run main.nf -profile hpc
```

**Requirements**

* Slurm
* Singularity available
* Singularity image, e.g.:

```bash
singularity pull bwa-sam.sif docker://shajiase/bwa-sam:1.16.1
```

* `configs/slurm.config`

---

## 📤 Output

Results are written to:

```
${project_folder}/${bwa_output}/
```

Each input pair produces:

* `SAMPLE.sam`
* `SAMPLE.bam`
* `SAMPLE.bam.stat`
* `SAMPLE.sorted.bam`
* `SAMPLE.sorted.bam.bai`

The files needed for downstream analysis are:

* `SAMPLE.sorted.bam`
* `SAMPLE.sorted.bam.bai`

---

## 📂 Project structure

```
main.nf
nextflow.config
configs/
├── local.config
└── slurm.config
test_data/
```

---

## ✔️ What this pipeline does

* Builds BWA index if not present

* Runs:

  ```
  bwa mem -M
  samtools view
  samtools flagstat
  samtools sort
  samtools index
  ```

* Skips samples that already have `*.sorted.bam.bai`
