# Exploring *Microctonus aethiopoides*: Genome Assembly and Annotation of Eight Ecotypes
This repository contains bash scripts used for the genome assembly and annotation of different strains of *Microctonus aethiopoides*. The scripts facilitate various stages of the genome assembly process, including basecalling, demultiplexing, quality control, read trimming, assembly, and post-assembly analysis.

# Genome Assembly
Eight *M. aethiopoides* strains were sequenced using both ONT and Illumina sequencing platforms. ONT was used for genome assembly, while Illumina was used for polishing. For ONT, two libraries were made, resulting in two datasets (R0119, R0120). Basecalling is done individually on each library. Conversly, two samples (MO_04 and MO_05) did not produce good quality data, so a replacement library (R0149) with two new samples was sequenced. In the following steps, MO_04 and MO_05 will not be represented in the codes.


Files received after sequencing are POD5 files (similar to fast5 and blow5 files). It's like:
   - Fast5, BLOW5 = OLD versions
   - POD5 = NEW version

## STEP 1: ONT reads Basecalling
Tool used is Dorado with simplex "basecaller" super accuracy model (sup).

### Sample Sheet Format
Dorado requires a sample sheet, which has to be in the following format.
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
Note: In my case, two unmapped BAM files are produced after basecalling, one for each of the two libraries. Since running in slurm gpu node is activated using following in config section **"#SBATCH --gpus-per-node=A100:1"**

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
## STEP 3: BAM to Fastq conversion

Tool used here is SAMtools "fastq". Unmapped individual BAM files are transformed into FASTQ files during this step.

```bash
#### Library : R0119
for i in 03 07 08 13 16 19; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
### Library : R0120
for i in 03 07 08 13 16 19; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
### Library : R0149
for i in 06 40; do samtools fastq MO_${i}.bam > MO_${i}.fastq ; done;
```
## Step 4: Concatenate R0119 and R0120
During this step, the samples from libraries R0119 and R0120 are combined:

```bash
cat ../04_demux_sample/R0119/01_fastq_files/MO_03.fastq ../04_demux_sample/R0120/01_fastq_files/MO_03.fastq > MO_03_cat.fastq
cat ../04_demux_sample/R0119/01_fastq_files/MO_07.fastq ../04_demux_sample/R0120/01_fastq_files/MO_07.fastq > MO_07_cat.fastq
cat ../04_demux_sample/R0119/01_fastq_files/MO_08.fastq ../04_demux_sample/R0120/01_fastq_files/MO_08.fastq > MO_08_cat.fastq
cat ../04_demux_sample/R0119/01_fastq_files/MO_13.fastq ../04_demux_sample/R0120/01_fastq_files/MO_13.fastq > MO_13_cat.fastq
cat ../04_demux_sample/R0119/01_fastq_files/MO_16.fastq ../04_demux_sample/R0120/01_fastq_files/MO_16.fastq > MO_16_cat.fastq
cat ../04_demux_sample/R0119/01_fastq_files/MO_19.fastq ../04_demux_sample/R0120/01_fastq_files/MO_19.fastq > MO_19_cat.fastq
```
This concatenation process merges the corresponding samples from libraries R0119 and R0120 into single FASTQ files (e.g., MO_03_cat.fastq, MO_07_cat.fastq, etc.).

## STEP 5: QC of fastq files

Tool used is NanoPlot & BBMap.

#### NanoPlot Syntax:

```bash
for i in 03 06 07 08 13 16 19 40;
do 
NanoPlot --verbose -t 8 --fastq MO_${i}_cat.fastq -o 01_QC/MO_${i} ;
done;
```
#### BBMap Syntax:

```bash
for i in 03 06 07 08 13 16 19 40;
do 
readlength.sh in=${i}_cat_fil.fastq out=${i}_histogram.txt
done;
```
This step involves using NanoPlot to perform quality checks on the FASTQ files. The BBMap tool is also utilized to check the read length statistics for specific files.

## STEP 6: Filtering reads using chopper

This step involves using chopper to filter reads based on quality (minimum quality score of 10) and read length (minimum length of 500 bases). The --headcrop 10 parameter is used to remove the first 10 bases from each read. The filtered reads are then saved to new FASTQ files (e.g., MO_${i}_cat_fil.fastq) for each sample.

```bash
for i in 03 06 07 08 13 16 19 40;
do 
chopper --threads $SLURM_CPUS_PER_TASK -q 10 -l 500 --headcrop 10 < ../MO_${i}_cat.fastq > MO_${i}_cat_fil.fastq ;
done;
```

## STEP 7: Repeat Quality Check (NanoPlot, BBMap, FastQC)
Repeat previous QC steps on filtered fastq files.

## STEP 8: Genome Assembly using FLYE
In this step, FLYE is used for genome assembly using the cleaned reads (MO_${i}cat_fil.fastq). The assembly results are saved to the directory (./MO${i}_flye). The --threads 16 option is used to specify the number of threads , and -g 129m specifies an estimated genome size of 129 megabases (`varies based on species`).

```bash
for i in 03 06 07 08 13 16 19 40; do
    flye --nano-hq /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq --out-dir ./MO_${i}_flye --threads 16 -g 129m;
done
```
## STEP 9: Polishing using MEDAKA
MEDAKA is used to polish the assembled genome using the filtered reads (MO_${i}cat_fil.fastq). The polished assembly results are saved to the directory ./MO${i}_flye_medaka.

```bash
for i in 03 04 05 07 08 13 16 19; do
    medaka_consensus -t $SLURM_CPUS_PER_TASK -i /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq -d MO_${i}_flye/assembly.fasta -o MO_${i}_flye_medaka ;
done
```
## STEP 10: Remove haplotigs and contig overlaps using purge_dups

For this step, the output from MEDAKA is used. The following code demonstrates the process for one sample, which is similarly applied to other samples.

#### 1. Convert the BAM file produced in MEDAKA to PAF file

**BAM to SAM:**

```bash
samtools view -@ 8 -h calls_to_draft.bam > calls_to_draft.sam
```
**SAM to PAF using paftools.js:**
```bash
paftools.js sam2paf calls_to_draft.sam > MO_40.paf
```
#### 2. Load purge_dups in NeSI and execute the following commands
`Generate statistics for PAF file:`
```bash
pbcstat MO_40.paf
```
`Calculate cutoffs:`
```bash
calcuts PB.stat > cutoffs 2> calcuts.log
```
`Split the consensus fasta file:`
```bash
split_fa ../consensus.fasta > con_split.fa
```
`Generate self-mapping PAF file:`
```bash
minimap2 -xasm5 -DP con_split.fa con_split.fa | gzip -c - > con_split.self.paf.gz
```
`Purge duplicates:`
```bash
purge_dups -2 -T cutoffs -c PB.base.cov con_split.self.paf.gz > dups.bed 2> purge_dups.log
```
`Extract sequences:`
```bash
get_seqs -e dups.bed ../consensus.fasta
```
Two files are produced: purged.fa and hap.fa. The former is used for further analysis.

Refer to this link for more information: [purge_dups GitHub](https://github.com/dfguan/purge_dups) & [purge_dups academic paper](https://academic.oup.com/bioinformatics/article-pdf/36/9/2896/48986490/btaa025.pdf).

## STEP 11: QC using CompleASM and QUAST

### Compleasm
In this step, BUSCO assessment is performed on the polished assemblies using CompleASM. Here I am not downloading the specific lineage database. Instead, I have mentioned the lineage database using the “-l” flag, which automatically downloads the required files.

```bash
for i in 03 06 07 08 13 16 19 40; do
    compleasm.py run -a MO_${i}_purged/purged.fasta -o MO_${i}_compleasm -l hymenoptera_odb10 -t 12 ;
done
```
### QUAST
QUAST (Quality Assessment Tool for Genome Assemblies) is a software tool used for evaluating the quality of genome assemblies. It provides various metrics such as N50, number of contigs, genome size, and misassemblies that is used to assess and compare the accuracy and completeness of assembled genomes. QUAST also generates informative visualizations to facilitate the interpretation of assembly results.

**Individual quast**
Here I used the closet reference genome(optional) for the comparative evaluation of individual genome
```bash
for i in 03 06 07 08 13 16 19 40; do
quast.py --fragmented -t 16 -o ./MO_${i}_purged/01_quast -r ../../../06_reference/ncbi_dataset/data/GCA_030347275.1/GCA_030347275.1_UoO_Maeth_IR_genomic.fna MO_${i}_purged/purged.fasta;
done
```
**Comparative quast**
Compare stats of all the assembly produced
```bash
quast.py -t 16 -o ./01_quast_compare -l 'MO_03,MO_06,MO_07,MO_08,MO_13,MO_16,MO_19,MO_40'  MO_03_purged/purged.fasta MO_06_purged/purged.fasta MO_07_purged/purged.fasta MO_08_purged/purged.fasta MO_13_purged/purged.fasta MO_16_purged/purged.fasta MO_19_purged/purged.fasta MO_40_purged/purged.fasta
```

# Processing Illumina Reads
## STEP 12: Concatebating Fastq files
The samples in Novoseq are run in two lanes. Therefore, once the data is received, samples from replicate lanes are combined using the `cat` command.

**Script for Concatenating Illumina FASTQ Reads from Two Lanes**
```bash
# Loop through all FASTQ.gz files recursively in the current directory
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
do 
    echo "Merging R1"
    cat "$i"_L00*_R1_001.fastq.gz > /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data/"$i"_R1.fastq.gz
    echo "Merging R2"
    cat "$i"_L00*_R2_001.fastq.gz > /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data/"$i"_R2.fastq.gz
done
```
## STEP 13: Quality check using Fastqc & Multiqc




