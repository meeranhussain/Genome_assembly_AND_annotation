# Exploring *Microctonus aethiopoides*: Genome Assembly and annotation of Eight Ecotypes (ONPROGRESS)
This repository contains bash scripts used for the genome assembly of different strains of *Microctonus aethiopoides*. The scripts facilitate various stages of the genome assembly process, including basecalling, demultiplexing, quality control, read trimming, assembly, and post-assembly analysis.

# Genome Assembly
Eight *M. aethiopoides* strains were sequenced using both ONT and Illumina sequencing platforms. ONT was used for genome assembly, while Illumina was used for polishing. For ONT, two libraries were made, resulting in two datasets (R0119, R0120). Basecalling is done individually on each library. Conversly, two samples (MO_04 and MO_05) did not produce good quality data, so a replacement library (R0149) with two new samples was sequenced. In the following steps, MO_04 and MO_05 will not be represented in the codes.


Files received after sequencing are POD5 files (similar to fast5 and blow5 files). It's like:
   - Fast5, BLOW5 = OLD versions
   - POD5 = NEW version

## 📑 Table of Contents

- [STEP 1: ONT reads Basecalling](#step-1-ont-reads-basecalling)
- [STEP 2: Demultiplexing](#step-2-demultiplexing)
- [STEP 3: BAM to Fastq conversion](#step-3-bam-to-fastq-conversion)
- [STEP 4: Concatenate R0119 and R0120](#step-4-concatenate-r0119-and-r0120)
- [STEP 5: QC of fastq files](#step-5-qc-of-fastq-files)
- [STEP 6: Filtering reads using chopper](#step-6-filtering-reads-using-chopper)
- [STEP 7: Repeat Quality Check (NanoPlot, BBMap, FastQC)](#step-7-repeat-quality-check)
- [STEP 8: Genome Assembly using FLYE](#step-8-genome-assembly-using-flye)
- [STEP 9: Polishing using MEDAKA](#step-9-polishing-using-medaka)
- [STEP 10: Remove haplotigs and contig overlaps using purge_dups](#step-10-remove-haplotigs-and-contig-overlaps-using-purge_dups)
- [STEP 11: QC using CompleASM and QUAST](#step-11-qc-using-compleasm-and-quast)
- [STEP 12: Concatenating Fastq files](#step-12-concatenating-fastq-files)
- [STEP 13: Quality check using Fastqc & Multiqc](#step-13-quality-check-using-fastqc-and-multiqc)
- [STEP 14: Illumina data Filtration using TrimGalore](#step-14-illumina-data-filtration-using-trimgalore)
- [STEP 15: Nextpolish Genome Polishing using Illumina Reads](#step-15-nextpolish-genome-polishing-using-illumina-reads)
- [STEP 16: Perform QC on polished assembly using Compleasm and Quast](#step-16-perform-qc-on-polished-assembly-using-compleasm-and-quast)
- [STEP 17: Kmer analysis using Illumina reads & Genome quality completeness estimation using Merqury](#step-17-kmer-analysis-using-illumina-reads-and-genome-quality-completeness-estimation-using-merqury)
- [STEP 18: KAT analysis (similar to merqury, uses Illumina reads)](#step-18-kat-analysis)
- [STEP 19: Scaffolding Using RagTag](#step-19-scaffolding-using-ragtag)
- [STEP 20: Quality Control (QC) Steps - Merqury, KAT Analysis, Qualimap, Quast, Compleasm](#step-20-quality-control-steps)
- [STEP 21: BlobToolKit for QC and Detecting Contamination in Assembly](#step-21-blobtoolkit-for-qc-and-detecting-contamination-in-assembly)
- [STEP 22: Filtration of Scaffolded Genome Assembly Using Custom Script](#step-22-filtration-of-scaffolded-genome-assembly-using-custom-script)
- [STEP 23: Repeat Annotation](#step-23-repeat-annotation)
- [STEP 24: Gene Annotation](#step-24-gene-annotation)
- [STEP 25: Annotation Quality Control (QC)](#step-25-annotation-quality-control-qc)
- [STEP 26: Functional Annotation](#step-26-functional-annotation)


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
dorado basecaller sup /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/R0119/LX0030/20240115_1457_1D_PAQ57576_8b278f23/pod5_pass --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0119.csv > /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/BAM_basecall/R0119_sup_calls.bam

### Library R0120
dorado basecaller sup /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/R0120/LX0030/20240117_1251_1D_PAQ57576_b43c06c8/pod5_pass --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0120.csv > /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/BAM_basecall/R0120_sup_calls.bam

### Library R0149
dorado basecaller sup /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/03_R0149/LX0038/20240402_1733_1F_PAW14881_72e2c281/pod5_pass  --recursive --device 'cuda:all' --kit-name SQK-NBD114-24 --sample-sheet sample_sheet_R0149.csv > /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/01_BAM_basecall/R0149_sup_calls.bam
```
Note: In my case, two unmapped BAM files are produced after basecalling, one for each of the two libraries. Since running in slurm gpu node is activated using following in config section **"#SBATCH --gpus-per-node=A100:1"**

## STEP 2: Demultiplexing

Tool used is Dorado "demux". This process separates individual samples into unmapped BAM files within a library, using the names provided in the "alias" column of the previous sample sheet.

```bash
#### Library : R0119
dorado demux --output-dir demux_sample --no-classify /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/BAM_basecall/R0119_sup_calls.bam

#### Library : R0120
dorado demux --output-dir demux_sample --no-classify /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/BAM_basecall/R0120_sup_calls.bam

#### Library : R0149
dorado demux --output-dir /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/02_demux_sample/R0149 --no-classify /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/01_BAM_basecall/R0149_sup_calls.bam
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

## STEP 7: Repeat Quality Check
Repeat previous QC steps  (NanoPlot, BBMap, FastQC) on filtered fastq files.

## STEP 8: Genome Assembly using FLYE
In this step, FLYE is used for genome assembly using the cleaned reads (MO_${i}cat_fil.fastq). The assembly results are saved to the directory (./MO${i}_flye). The --threads 16 option is used to specify the number of threads , and -g 129m specifies an estimated genome size of 129 megabases (`varies based on species`).

```bash
for i in 03 06 07 08 13 16 19 40; do
    flye --nano-hq /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq --out-dir ./MO_${i}_flye --threads 16 -g 129m;
done
```
## STEP 9: Polishing using MEDAKA
MEDAKA is used to polish the assembled genome using the filtered reads (MO_${i}cat_fil.fastq). The polished assembly results are saved to the directory ./MO${i}_flye_medaka.

```bash
for i in 03 04 05 07 08 13 16 19; do
    medaka_consensus -t $SLURM_CPUS_PER_TASK -i /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq -d MO_${i}_flye/assembly.fasta -o MO_${i}_flye_medaka ;
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

**Individual Quast**
Here I used the closet reference genome(optional) for the comparative evaluation of individual genome
```bash
for i in 03 06 07 08 13 16 19 40; do
quast.py --fragmented -t 16 -o ./MO_${i}_purged/01_quast -r ../../../06_reference/ncbi_dataset/data/GCA_030347275.1/GCA_030347275.1_UoO_Maeth_IR_genomic.fna MO_${i}_purged/purged.fasta;
done
```
**Comparative Quast**
Compare stats of all the assembly produced
```bash
quast.py -t 16 -o ./01_quast_compare -l 'MO_03,MO_06,MO_07,MO_08,MO_13,MO_16,MO_19,MO_40'  MO_03_purged/purged.fasta MO_06_purged/purged.fasta MO_07_purged/purged.fasta MO_08_purged/purged.fasta MO_13_purged/purged.fasta MO_16_purged/purged.fasta MO_19_purged/purged.fasta MO_40_purged/purged.fasta
```

# Processing Illumina Reads
## STEP 12: Concatenating Fastq files
The samples in Novoseq are run in two lanes. Therefore, once the data is received, samples from replicate lanes are combined using the `cat` command.

**Script for Concatenating Illumina FASTQ Reads from Two Lanes**
```bash
# Loop through all FASTQ.gz files recursively in the current directory
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)
do 
    echo "Merging R1"
    cat "$i"_L00*_R1_001.fastq.gz > /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data/"$i"_R1.fastq.gz
    echo "Merging R2"
    cat "$i"_L00*_R2_001.fastq.gz > /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data/"$i"_R2.fastq.gz
done
```
## STEP 13: Quality check using Fastqc and Multiqc

```bash
**Fastqc**
fastqc -t 8 -o 00_QC/fastqc/ /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data/*.fastq.gz

**MultiQC**
multiqc /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/02_cat_data
```

## STEP 14: Illumina data Filtration using TrimGalore
I performed quality trimming of bases with a Phred score of 20 and below.

```bash
for i in MI-19 MI-03 MI-07 MI-08 MI-13 MI-16 MI-40  MI-06;
do
    trim_galore -q 20 --paired --fastqc --cores 20 ${i}_R1.fastq.gz ${i}_R2.fastq.gz -o ../03_fil_data;
done
```
## STEP 15: Nextpolish Genome Polishing using Illumina Reads

NextPolish is a tool used for polishing genomes using Illumina reads. In the script, NextPolish is applied to the purged genome fasta file obtained from a previous step. The process involves the following steps:

### Process
-   **Alignment:** The purged genome fasta file is indexed and aligned to filtered Illumina paired-end reads using `BWA-MEM`.
-   **Alignment Processing:** The aligned reads are processed using `SAMtools` to remove duplicate reads, sort the alignments, and mark duplicates.
-   **Polishing Round 1:** `NextPolish` is used with specific parameters (-t 1) to polish the genome based on the alignments from the first round. This step generates a temporary polished genome file (genome.polishtemp.fa).
-   **Polishing Round 2:** The temporary polished genome file from Round 1 is indexed and aligned again to the Illumina reads. `NextPolish` is then applied again (-t 2) to perform a second round of polishing, resulting in the final polished genome file **(genome.nextpolish.fa)**.

```bash

for sam in 03 07 08 13 16 19 06 40;
do
    ##########
    # PARAMS #
    INDIR=/nesi/nobackup/acc_name/PX024_Parasitoid_wasp/01_raw/04_LIC_Data/03_fil_data/
    read1=${INDIR}MI_${sam}_R1.fq.gz
    read2=${INDIR}MI_${sam}_R2.fq.gz
    OUTDIR=/nesi/nobackup/acc_name/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/05_nextpolish/Maethio_${sam}
    REFDIR=/nesi/nobackup/acc_name/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/04_purg_dups/02_af_medaka/MO_${sam}/
    REF=${REFDIR}MO_${sam}_purged.fasta
    round=2
    threads=20
    NEXTPOLISH=/nesi/project/acc_name/softwares/NextPolish/lib/nextpolish1.py
    #####################################################################
    mkdir -p $OUTDIR
    cd $OUTDIR
    for ((i=1; i<=${round};i++)); do
        # Step 1:
        echo index the genome file and do alignment $i;
        bwa index ${REF};
        bwa mem -t ${threads} ${REF} ${read1} ${read2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3 - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam;
        echo index bam and genome files $i;
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${REF};
        echo polish genome file $i;
        python ${NEXTPOLISH} -g ${REF} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
        REF=genome.polishtemp.fa;
        echo round $i complete;
        
        # Step 2:
        echo index genome file and do alignment $i;
        bwa index ${REF};
        bwa mem -t ${threads} ${REF} ${read1} ${read2} | samtools view --threads 3 -F 0x4 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 2g --threads 5 - | samtools markdup --threads 5 -r - sgs.sort.bam;
        echo index bam and genome files $i;
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${REF};
        echo polish genome file $i;
        python ${NEXTPOLISH} -g ${REF} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
        REF=genome.nextpolish.fa;
        echo round $i complete;
    done
    echo genome polished for Maethio_${sam}
done
# Finally polished genome file: genome.nextpolish.fa
```
## STEP 16: Perform QC on polished assembly using Compleasm and Quast
Repeat Step 11

## Step 17: Kmer analysis using Illumina reads and Genome quality completeness estimation using Merqury

1. **Create kmer database of high-quality Illumina reads using meryl**
   - Steps to follow: [Merqury Meryl DB Preparation](https://github.com/marbl/merqury/wiki/1.-Prepare-meryl-dbs)
  
2. **Create histogram plot in Genomescope**
   Generate histogram from meryl data
   Example:
   - `meryl histogram MO_03_k18.meryl > MO_03_k18.hist`

   **Note:** The `.hist` file contains two columns separated by `\t` (TAB) which needs to be replaced by a space.
   - `tr '\t' ' ' <MO_03_k18.hist > MO_03_k18_s.hist`
   - Run GenomeScope by loading the `.hist` file at [GenomeScope](http://genomescope.org/)
   - Adjust kmer length accordingly

4. **Run merqury**
```slurm
#!/bin/bash -e
#SBATCH --account=acc_name
#SBATCH --job-name=merqury
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=26G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output merqury_%j.out    # save the output into a file
#SBATCH --error merqury_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfare
module purge
## load module 
module load Merqury/1.3-Miniconda3
module load R
# Need argparse, ggplot2, scales
Rscript -e 'install.packages(c("argparse", "ggplot2", "scales"),repos = "http://cran.us.r-project.org")'

module load SAMtools/1.19-GCC-12.3.0
module load BEDTools/2.30.0-GCC-11.3.0
module load Java/20.0.2
module load IGV/2.16.1

for i in 03 06 07 08 13 16 19 40;
do
mkdir Maethio_${i}/00_QC/02_merqury
cd Maethio_${i}/00_QC/02_merqury
#### Link files
#meryl
ln -s /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/05_ncgenome/02_fil_Illumina/00_QC/02_meryl/MI_${i}_k18.meryl
#assembly
ln -s /nesi/nobackup/acc_name/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/05_nextpolish/Maethio_${i}/genome.nextpolish.fa
####meryl count
/nesi/project/acc_name/softwares/merqury-master/merqury.sh MI_${i}_k18.meryl genome.nextpolish.fa Maethio_${i}_nxtplsh_mqry
cd ../../../
done
```

## Step 18: KAT analysis 
Similar to merqury, uses Illumina reads
```slurm
#!/bin/bash -e
#SBATCH --account=acc_name
#SBATCH --job-name=kat_analysis
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=2
#SBATCH --mem=26G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output kat_analysis_%j.out    # save the output into a file
#SBATCH --error kat_analysis_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfere
module purge
ml KAT/2.4.2-gimkl-2018b-Python-3.7.3

########## Loop ##############
for i in 03 06 07 08 13 16 19 40;
do
    mkdir -p Maethio_${i}/00_QC/04_kat_results
    cd Maethio_${i}/00_QC/04_kat_results
    #### Link files
    #assembly
    #ln -s /nesi/nobackup//PX024_Parasitoid_wasp/05_ncgenome/01_assembly/05_nextpolish/Maethio_${i}/genome.nextpolish.fa
    ####kat command
    kat comp -t 16 -o Mei_${i}_kat “/nesi/nobackup/acc_name/PX024_Parasitoid_wasp/05_ncgenome/02_fil_Illumina/MI_06_R?.fq.gz” genome.nextpolish.fa
    cd ../../../
done
```

## Step 19: Scaffolding Using RagTag

In this step, the polished genome is scaffolded using a reference genome with the **RagTag** tool.

---

### About RagTag

**RagTag** is a genome scaffolding and assembly tool designed to use one genome assembly (usually a high-quality reference) to improve the structure of another (usually a draft genome). The main input and output for RagTag scaffolding are as follows:

- **Input**:
  - **Reference genome**: The high-quality genome assembly used to scaffold the draft genome.
  - **Query genome**: The draft genome assembly to be scaffolded.

- **Output**:
  - **Scaffolded genome**: A new version of the draft genome with an improved structure based on the reference genome.

---

### SLURM Job Script for RagTag Scaffolding

The following SLURM job script runs RagTag to scaffold multiple draft genomes using a reference genome:

```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=ragtag_scfld
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=2
#SBATCH --mem=26G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output ragtag_scfld_%j.out    # save the output into a file
#SBATCH --error ragtag_scfld_%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfere
module purge

# Load necessary modules
ml Python/3.11.6-foss-2023a
ml minimap2/2.28-GCC-12.3.0
ml unimap/0.1-GCC-11.3.0
pip install RagTag

# Run RagTag scaffolding
for i in 03 06 07 08 13 16 19 40;
do
    mkdir -p 01_scaffold
    mkdir -p 01_scaffold/Maethio_${i}
    ragtag.py scaffold /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/06_reference/01_Maethio_IR/GCA_030347275.1_UoO_Maeth_IR_genomic.fna \
    /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/05_nextpolish/Maethio_${i}/genome.nextpolish.fa \
    -o 01_scaffold/Maethio_${i}/
done
```

## Step 20: Quality Control Steps 
Merqury, KAT Analysis, Qualimap, Quast, Compleasm
Perform QC on scaffolded genome using the tools mentioned below. Since BAM files are not generated in the previous scaffolding step, additional steps are required to create them for **Qualimap** analysis.

---

### QC Tools Overview (Steps on using these tools are described already, check above)

- **Merqury**: Evaluates the quality of genome assemblies using k-mers from short reads to estimate assembly completeness and accuracy.
- **KAT (K-mer Analysis Toolkit)**: Provides various tools for k-mer analysis of sequencing data, useful for genome assembly QC.
- **Qualimap**: Performs quality control on BAM files, providing insights into coverage, GC content, and other sequencing metrics.
- **Quast**: Assesses the accuracy and completeness of genome assemblies.
- **Compleasm**: Used for completeness assessment in assemblies, helping verify if all genomic regions are covered.

---

### Step-by-Step Instructions for BAM File Generation and Qualimap Analysis

To run **Qualimap** on scaffolded genome assemblies, align short reads to the scaffolded genome, generate BAM files, and perform **Qualimap bamqc** on the aligned files.

```bash
for i in 07 08 13 16 19 40;
do
    cd Maethio_${i};

    # Index the scaffolded genome for alignment
    bwa index ragtag.scaffold.fasta;

    # Align reads to the scaffolded genome and process BAM file
    bwa mem -t 24 ragtag.scaffold.fasta \
    /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/02_fil_Illumina/MI_${i}_R1.fq.gz \
    /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/02_fil_Illumina/MI_${i}_R2.fq.gz | \
    samtools view --threads 16 -F 0x4 -b - | \
    samtools fixmate -m --threads 16 - - | \
    samtools sort -m 2g --threads 16 - | \
    samtools markdup --threads 16 -r - 00_QC/04_qualimap/Maethio_${i}_sort_sfld.bam;

    # Index the BAM file
    samtools index -@ 24 00_QC/04_qualimap/Maethio_${i}_sort_sfld.bam -o 00_QC/04_qualimap/Maethio_${i}_sort_sfld.bam.bai;

    # Perform Qualimap bamqc on the BAM file
    /nesi/project/uow03744/softwares/qualimap_v2.3/qualimap bamqc \
    -bam 00_QC/04_qualimap/Maethio_${i}_sort_sfld.bam \
    -outdir 00_QC/04_qualimap \
    -nt 16 --java-mem-size=8G;

    cd ../;
done
```

## Step 21: BlobToolKit for QC and Detecting Contamination in Assembly

BlobToolKit is used here for quality control and to detect potential contamination in the genome assembly. Several files need to be generated for BlobToolKit to function effectively:
  - **Blast output** (to identify taxonomic affiliations)
  - **Aligned Nanopore reads** (BAM file with `.csi` index)
  - **Assembly** (FASTA format)
  - **BUSCO report** (CSV format)

---

## Steps for Generating Required Files

### a) BLAST the Assembly Against the NT Database

This step uses BLAST to align the scaffolded assembly to the nucleotide (nt) database, identifying likely taxonomic sources.

#### SLURM Job Script for BLAST

```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=blastn
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=46
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output Blastn_%j.out
#SBATCH --error Blastn_%j.err

# Purge any pre-loaded modules
module purge

# Load BLAST modules
module load BLASTDB/2024-07 BLAST/2.13.0-GCC-11.3.0

# Run BLAST for each scaffolded genome
for i in 06 07 08 13 16 19 40; do 
    blastn \
    -query ./Maethio_${i}/ragtag.scaffold.fasta \
    -task megablast \
    -db nt \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -culling_limit 10 \
    -num_threads 42 \
    -evalue 1e-3 \
    -out ./Maethio_${i}/Maethio_${i}_megablast.out;
done
```
#### a2) Alternative BLAST with GNU Parallel

This method splits the assembly into smaller chunks and runs BLAST on each file in parallel, speeding up the process.

```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=blastn
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output Blastn_%j.out
#SBATCH --error Blastn_%j.err

# Purge all loaded modules
module purge

# Load BLAST and Parallel modules
module load BLASTDB/2024-07 BLAST/2.13.0-GCC-11.3.0 
module load Parallel/20220922 Python

# Install biopython (required for splitting FASTA files)
pip install biopython

# Split the FASTA file by contigs
python split_fasta_by_contigs.py ragtag.scaffold.fasta query_chunk 10

# Run BLAST in parallel on each chunk
ls query_chunk* | parallel "blastn -query {} -task megablast -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 10 -num_threads 4 -evalue 1e-3 -out {.}.out"

# Combine all individual BLAST outputs into a single file
cat query_chunk_*.out > Maethio_08_megablast.out

```

### b) Minimap2 to Align Long Reads to Scaffold

This step uses Minimap2 to align Nanopore long reads to the scaffolded assembly. The genome results from previous QC steps are used for filtration.

#### SLURM Job Script for Minimap2

```bash
#!/bin/bash -e
#SBATCH --job-name=minimap
#SBATCH --account=uow03744
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output minimap_%j.out
#SBATCH --error minimap_%j.err

# Load required modules
ml SAMtools
ml minimap2/2.28-GCC-12.3.0

# Align long reads for each sample
for i in 03 06 07 08 13 16 19 40; do
    cd Maethio_${i}
    
    # Index the scaffolded genome
    minimap2 -t 24 -x map-ont ragtag.scaffold.fasta -d ragtag.scaffold.mmi
    
    # Map long reads to the scaffold and create sorted BAM files
    minimap2 -t 24 -ax map-ont ragtag.scaffold.mmi /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/02_Basecaller_sup/03_fastqfiles/02_filtered/MO_${i}_cat_fil.fastq | \
    samtools view --threads 24 -F 0x4 -b - | \
    samtools fixmate -m --threads 24 - - | \
    samtools sort -m 2g --threads 24 - | \
    samtools markdup --threads 24 -r - 00_QC/06_nanoreads_map/Maethio_${i}_sort.bam

    # Index the sorted BAM file
    samtools index -bc -@ 24 00_QC/06_nanoreads_map/Maethio_${i}_sort.bam -o 00_QC/06_nanoreads_map/Maethio_${i}_sort.bam.csi
    
    cd ../
done
```
### c) BUSCO Analysis for Completeness

In this step, BUSCO is run separately to assess genome completeness. This analysis checks for the presence of conserved single-copy orthologs, indicating the completeness of the genome assembly. Modifying the BUSCO output format may make it compatible with BlobToolKit.

#### SLURM Job Script for BUSCO

```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=BUSCO
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output BUSCO_%j.out
#SBATCH --error BUSCO_%j.err

# Purge any pre-loaded modules
module purge

# Load BUSCO module
ml BUSCO/5.6.1-gimkl-2022a

# Run BUSCO for each scaffolded genome
for i in 03 06 07 08 13 16 19 40; do
    busco -i Maethio_${i}/ragtag.scaffold.fasta -m genome -l hymenoptera_odb10 -c 24 -o Maethio_${i}/Maethio_${i}_busco
done
```

### d) BlobToolKit Analysis

BlobToolKit is used to create a dataset for visualizing genome assembly metrics such as coverage, GC content, and possible contamination. Refer to the [BlobToolKit NeSI guide](https://nesi.github.io/blobtools-jupyter-vdt/) for specific usage on the NeSI platform.

#### Steps for BlobToolKit

1. **Create Blob Dataset**

   Use the following command to create a BlobToolKit dataset with BUSCO, FASTA, coverage, and BLAST hit data.

   ```bash
   blobtools create --busco Maethio_03_busco/run_hymenoptera_odb10/full_table.tsv --fasta ragtag.scaffold.fasta --cov Maethio_03_sort.bam --hits Maethio_03_megablast.out --replace --taxdump ../00_taxon Ma_USA
   ```
2. **View Blob Dataset**

   To view the BlobToolKit dataset, use the host command to start a local server:

   ```bash
   blobtools host --port 8081 --api-port 9001 --hostname localhost Ma_USA/
   ```

   - **Download the Taxonomy Table: After starting the server, download the taxonomy table as "Ma_USA.csv"**
   - **Explore Plots: Review plots for coverage, GC content, and taxonomic distribution; download any required visualizations.**

# Step 22: Filtration of Scaffolded Genome Assembly Using Custom Script 

The `Filter_assembly.sh` script (Provided in this repository) is utilized to remove contaminants and filter out contigs that are less than 1000 bp in length from the scaffolded genome assembly. **If you need to remove contigs of greater length, modify script as needed**

### Usage

```bash
./Filter_assembly.sh -i <fasta_file> -b <busco_file> -a <asm_stats_file> -f <blob_file> -o <output_name>
```
### Options

- `-i <fasta_file>`: Path to the FASTA file.
- `-b <busco_file>`: Path to the BUSCO file (use the `full_table.tsv` from BUSCO).
- `-a <asm_stats_file>`: Path to the assembly stats file (e.g., `genome_results.txt` from Qualimap).
- `-f <blob_file>`: Path to the blob file (table from BlobToolKit).
- `-o <output_name>`: Desired name for the output file.
- `-h`: Display this help message.

> **Note:** It is recommended to run the script from the directory where you want to store the filtered FASTA file. Future improvements will aim to enhance the efficiency of this script.

### Perform Quality Check of filtered assembly using previously mentioned QC tools
After filtration, conduct quality control using tools like Qualimap and Compleasm to ensure that the genome assembly meets the desired quality standards.


# ANNOTATION
## STEP 23: Repeat Annotation

## Rename Contig IDs

**Purpose**:  
EDTA requires the input FASTA file to have contig IDs that are **no longer than 13 characters**. To ensure compatibility, rename the contigs in a sequential format such as `contig_01`, `contig_02`, ..., `contig_n`.

**Usage**:
```bash
./contig_name.sh input_fasta output_fasta
```
This script replaces original contig names in your input FASTA with a sequential naming scheme that adheres to EDTA's requirements.

## Transposable Element (TE) Annotation using EDTA
This step performs repeat identification and annotation using EDTA (Extensive de-novo TE Annotator), followed by repeat masking using RepeatMasker.

 **Note:** In this workflow, EDTA was run without the --anno 1 flag. If you want EDTA to output a masked FASTA file directly, run it with --anno 1. In this case, a custom repeat library TE.lib.fa generated by EDTA was used for masking with RepeatMasker instead.

### Repeat Annotation with EDTA
#### SLURM Job Script:
```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=EDTA_03
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=EDTA_03_%j.out
#SBATCH --error=EDTA_03_%j.err

# Purge conflicting modules
module purge

# Load EDTA module
ml EDTA/2.1.0

# Run EDTA
EDTA.pl --genome /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/07_rm_cont/Maethio_03/01_filtered_scfld/Maethio_03_scfld_fil_mod.fasta --threads 32 --sensitive 1
```
This step builds a TE library by identifying and classifying repeats in the genome. The output includes a custom .TElib.fa file that can be used for masking.

### Repeat Masking with RepeatMasker
RepeatMasker was used separately for masking using the .TElib.fa file generated by EDTA.

#### SLURM Job Script:
```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=Repeatmasker
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=36
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=Repeatmasker_%j.out
#SBATCH --error=Repeatmasker_%j.err

# Purge conflicting modules
module purge

# Load RepeatMasker
ml RepeatMasker/4.1.0-gimkl-2020a

# Loop through selected strains and run RepeatMasker
for i in 03 06 07 08 13 16 19 40; do
    cd Maethio_${i}
    RepeatMasker -pa 32 -xsmall -lib Maethio_${i}_scfld_fil_mod.fasta.mod.EDTA.TElib.fa Maethio_${i}_scfld_fil_mod.fasta
    cd ../
done
```
##### Key Parameters Explained
- `-pa 32`: Runs 32 threads in parallel

- `-xsmall`: Converts masked regions to lowercase in the output FASTA

- `-lib`: Specifies a custom library for masking, generated by EDTA

## STEP 24: Gene Annotation

**Purpose**:  
To annotate gene models on the repeat-masked genome using homology-based evidence from protein databases, leveraging **BRAKER3**  pipelines.

**Preparation**:
- A soft link was created to the masked genome file `Maethio_03_scfld_fil_mod.fasta.masked` and renamed as `Maethio_03_masked.fasta` within the respective ecotype annotation folder.
- Protein evidence used for annotation was a concatenated database comprising:
  - `Uniprot-SwissProt` proteins
  - `OrthoDB` Arthropoda proteins

```bash
# Combine protein databases
cat Arthropoda.fasta Uniprot_sprot.fasta > proteins.fasta
```

##  BRAKER3-based Annotation

#### SLURM Job Script:
```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=braker_03
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=braker_%j.out
#SBATCH --error=braker_%j.err

# Clean environment
module purge

# Load Singularity module
ml Singularity/3.11.3

# Link the masked genome and protein database
ln -s /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/08_EDTA/Maethio_03/Maethio_03_scfld_fil_mod.fasta.masked ./Maethio_03_masked.fasta
ln -s ../proteins.fasta

# Run BRAKER3
singularity exec /nesi/project/uow03744/softwares/braker3.sif braker.pl \
  --threads=12 \
  --genome=Maethio_03_masked.fasta \
  --prot_seq=proteins.fasta \
  --species=Ma_usa_braker1 \
  --gff3 \
  --AUGUSTUS_ab_initio \
  --crf
```
##### Key Parameters Explained

- `--genome`: The repeat-masked genome file used as input.
- `--prot_seq`: The combined protein evidence in FASTA format (e.g., SwissProt + OrthoDB).
- `--species`: A species identifier used for AUGUSTUS model training.
- `--gff3`: Ensures that the output is in standard GFF3 format.
- `--AUGUSTUS_ab_initio`: Enables ab initio gene prediction using AUGUSTUS.
- `--crf`: Activates Conditional Random Fields (CRF) for improved gene model prediction accuracy.

## STEP 25: Annotation Quality Control (QC)

###  BUSCO Assessment Using Protein Mode (`-p`)

To evaluate the completeness of gene annotations generated by BRAKER3, **BUSCO** was run in **protein mode** using the `compleasm.py` wrapper. This assesses how well the predicted proteins match conserved orthologs from the Hymenoptera lineage.

**Lineage dataset used**: `hymenoptera_odb10`  
**Mode**: Protein  
**Threads**: 12

**Command:**
```bash
# Run BUSCO on BRAKER3-predicted protein sequences for each strain
for i in 03 06 07 08 13 16 19 40 FRN MOR; do
  compleasm.py protein \
    -p Maethio_${i}/braker/braker.aa \
    -o Maethio_${i}/00_QC/01_compleasm_proteincheck \
    -l hymenoptera_odb10 \
    -t 12
done
```
This step outputs BUSCO scores indicating the proportion of: Complete (C), Fragmented (F), Missing (M) genes in the annotated protein set, providing an estimate of annotation quality.


### GeneValidator Assessment

**Purpose**:  
GeneValidator is used to assess the quality and biological plausibility of predicted protein-coding genes. It evaluates gene models based on metrics such as alignment quality, completeness, length distribution, and homology to reference proteins.

**Reference Database Used**:  
- `uniprot_sprot.fasta` (Swiss-Prot curated protein database)

> **Tip**: Before running GeneValidator, ensure your reference protein file is formatted as a BLAST database using `makeblastdb`.

---

#### SLURM Job Script for GeneValidator

```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=genevalidator
#SBATCH --time=110:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=genevalidator_%j.out
#SBATCH --error=genevalidator_%j.err

# Clean environment
module purge

# Load BLAST module (required for homology search)
ml BLAST/2.16.0-GCC-12.3.0

# Optional: Create BLAST database (do once before running the loop)
# makeblastdb -in ./uniprot_sprot.fasta -dbtype prot -parse_seqids

# Run GeneValidator for each ecotype
for i in 03 06 07 08 13 16 19 40 FRN MOR; do
  /nesi/project/uow03744/softwares/genevalidator/genevalidator/bin/genevalidator \
    -d ./uniprot_sprot.fasta \
    --num_threads 24 \
    -m 8 \
    -o Maethio_${i}/00_QC/02_genevalidator_uni_swissprot \
    Maethio_${i}/braker/braker.aa
done
```

##### Key Parameters Explained

- `-d ./uniprot_sprot.fasta`: Specifies the reference protein database for homology comparison (Swiss-Prot).
- `--num_threads 24`: Runs 24 threads for faster performance.
- `-m 8`: Sets the maximum number of genes to display per gene cluster (optional for reducing HTML report size).
- `-o`: Output directory where GeneValidator reports (HTML, TSV) will be saved.
- *Final argument*: The input FASTA file containing predicted protein sequences from BRAKER (`braker.aa`).

## STEP 26: Functional Annotation

**Tool Used**: [EggNOG-mapper v2.1.12](http://eggnog-mapper.embl.de/)  
**Purpose**: Assign functional annotations (GO terms, KEGG pathways, COG categories, PFAM domains, etc.) to predicted proteins from structural annotation.

### Database Setup and Annotation Pipeline

#### SLURM Job Script:
```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=eggnog
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=eggnog_%j.out
#SBATCH --error=eggnog_%j.err

# Clean environment
module purge

# Load EggNOG-mapper module
module load eggnog-mapper/2.1.12-gimkl-2022a

# Loop through ecotypes for functional annotation
for i in 03 06 07 08 13 16 19 40 FRN MOR; do
  cd Maethio_${i}

  # Optional: Create softlink to protein FASTA file
  # ln -s /path/to/braker.aa Maethio_${i}_protein.fa

  # Copy GFF3 file for functional annotation decoration
  cp /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/09_annotation/01_str_anno/01_BRAKER/Maethio_${i}/braker/braker.gff3 ./Maethio_${i}_braker.gff3

  # Run EggNOG-mapper
  emapper.py \
    -i Maethio_${i}_protein.fa \
    -o Maethio_${i}_eggnog \
    --data_dir /nesi/nobackup/uow03744/PX024_Parasitoid_wasp/05_ncgenome/01_assembly/09_annotation/03_eggnog/00_eggnog_db \
    --decorate_gff Maethio_${i}_braker.gff3 \
    --cpu 22 \
    --resume

  cd ../
done
```
##### Key Parameters Explained

- `-i Maethio_${i}_protein.fa`: Input FASTA file of predicted protein sequences.
- `-o Maethio_${i}_eggnog`: Output prefix for EggNOG-mapper results.
- `--data_dir`: Path to the locally downloaded EggNOG-mapper database (recommended for speed and offline use).
- `--decorate_gff`: Annotates the GFF3 file with functional terms (e.g., GO terms, KEGG pathways, COGs).
- `--cpu 22`: Number of threads used for parallel processing.
- `--resume`: Resume the job if previously interrupted or partially completed.

## Step 27: Submit InterProScan

### Remove Asterisk (`*`) Stop Codon Symbols before Interproscan 

BRAKER3 outputs protein sequences with an asterisk (`*`) at the end of each protein sequence to indicate a stop codon. These symbols must be removed before running InterProScan to prevent errors or misinterpretation of sequences.

#### Command to Remove `*` from Sequences:
```bash
# Remove "*" from all sequences in the input FASTA file
sed -e 's/*//g' Maethio_03_braker.aa_modified.faa > Maethio_03_final_fixed.faa
```

### Run InterProScan 
The following SLURM batch script runs InterProScan with selected applications and includes gene ontology and pathway annotations. Ensure the cleaned protein file (*_final_fixed.faa) is used.
```bash
#!/bin/bash -e
#SBATCH --account=uow03744
#SBATCH --job-name=interproscan_06
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meeranhussain1996@gmail.com
#SBATCH --output=eggnog_%j.out    # Save standard output to file
#SBATCH --error=eggnog_%j.err     # Save error output to file

# Clean environment
module purge

# Load required modules
module load InterProScan/5.66-98.0-gimkl-2022a Perl-5.34.1 Python-3.11.3 PCRE2/10.42-GCCcore-12.3.0 GCC/11.3.0

# Run InterProScan
interproscan.sh \
  -i Maethio_06_final_fixed.faa \
  -t p \
  -dp \
  -appl Pfam,SMART,PrositeProfiles,Gene3D,SUPERFAMILY,PANTHER \
  --goterms \
  --pathways \
  --iprlookup \
  -cpu 22 \
  -T temp/
```
##### Key Parameters Explained

 - `-i`          : Input FASTA file (cleaned protein file)                            
 - `-t p`        : Input type is protein                                              
 - `-appl`       : Specifies which InterPro applications/databases to run             
 - `--goterms`   : Include Gene Ontology term mappings                                
 - `--pathways`  : Include pathway annotations (e.g., KEGG)                           
 - `--iprlookup` : Enable InterPro cross-reference lookup                             
 - `-cpu`        : Number of threads to use (adjust based on HPC allocation)          
 - `-T temp/`    : Temporary working directory for InterProScan output and processing 
