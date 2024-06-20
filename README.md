# Genome_assembly
This repository contains bash scripts used for the genome assembly of different strains of *Microctonus aethiopoides*. The scripts facilitate various stages of the genome assembly process, including quality control, read trimming, assembly, and post-assembly analysis.

# NANOPORE BASECALLING & DEMULTIPLEXING
Eight *M. aethiopoides* strains were sequenced using both ONT and Illumina sequencing platforms. ONT was used for genome assembly, while Illumina was used for polishing. For ONT, two libraries were made, resulting in two datasets (R0119, R0120). Basecalling is done individually on each library. Conversly, two samples (MO_04 and MO_05) did not produce good quality data, so a replacement library (R0149) with two new samples was sequenced. In the following steps, MO_04 and MO_05 will not be represented in the codes.


1. Files received after sequencing are POD5 files (similar to fast5 files).
   - Fast5 = OLD
   - POD5 = NEW

## STEP 1: Basecalling
Tool used is Dorado with simplex "basecaller" super accuracy model (sup).

### Sample Sheet Format
Dorado requires a sample sheet, which has to be in the following format. Note that the "alias" column information will be used to produce output filenames for each sample in a library.

#### Sample Sheet 1: R0119

| kit            | experiment_id | flow_cell_id | barcode  | alias |
| -------------- | ------------- | ------------ | -------- | ----- |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode01 | MO_03 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode02 | MO_04 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode03 | MO_05 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode04 | MO_07 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode05 | MO_08 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode06 | MO_13 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode07 | MO_16 |
| SQK-NBD114-24  | R0119         | PAQ57576     | barcode08 | MO_19 |

#### Sample Sheet 2: R0120

| kit            | experiment_id | flow_cell_id | barcode  | alias |
| -------------- | ------------- | ------------ | -------- | ----- |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode01 | MO_03 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode02 | MO_04 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode03 | MO_05 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode04 | MO_07 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode05 | MO_08 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode06 | MO_13 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode07 | MO_16 |
| SQK-NBD114-24  | R0120         | PAQ57576     | barcode08 | MO_19 |

#### Sample Sheet 3: R0149

| kit            | experiment_id | flow_cell_id | barcode  | alias |
| -------------- | ------------- | ------------ | -------- | ----- |
| SQK-NBD114-24  | R0149         | PAW14881     | barcode09 | MO_06 |
| SQK-NBD114-24  | R0149         | PAW14881     | barcode19 | MO_40 |

Note: The information in the 'alias' column will be used to produce output filenames for each sample in a library.

### Dorado Basecaller Syntax

```bash
### Library R0119
dorado basecaller sup /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/R0119/LX0030/20240115_1457_1D_PAQ57576_8b278f23/pod5_pass --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0119.csv > /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/BAM_basecall/R0119_sup_calls.bam

### Library R0120
dorado basecaller sup /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/R0120/LX0030/20240117_1251_1D_PAQ57576_b43c06c8/pod5_pass --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0120.csv > /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/BAM_basecall/R0120_sup_calls.bam

### Library R0149
dorado basecaller sup /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/01_raw/03_R0149/LX0038/20240402_1733_1F_PAW14881_72e2c281/pod5_pass  --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0149.csv > /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/01_BAM_basecall/R0149_sup_calls.bam
```
Note: In my case, two unmapped BAM files are produced after basecalling, one for each of the two libraries. Since running in slurm gpu node will be activated **"#SBATCH --gpus-per-node=A100:1"**

## STEP 2: Demultiplexing

Tool used is Dorado "demux". This process separates individual samples into unmapped BAM files within a library, using the names provided in the "alias" column of the previous sample sheet.

```bash
#### Library : R0119
dorado demux --output-dir demux_sample --no-classify /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/BAM_basecall/R0119_sup_calls.bam

#### Library : R0120
dorado demux --output-dir demux_sample --no-classify /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/BAM_basecall/R0120_sup_calls.bam

#### Library : R0149
dorado demux --output-dir /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/02_demux_sample/R0149 --no-classify /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/01_BAM_basecall/R0149_sup_calls.bam
```
### STEP 3: BAM to Fastq conversion

Tool used here is SAMtools "fastq". Unmapped individual BAM files are transformed into FASTQ files during this step.

```bash
#### Library : R0119
for i in 03 04 05 07 08 13 16 19; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
### Library : R0120
for i in 03 04 05 07 08 13 16 19; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
### Library : R0149
for i in 06 40; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
```
