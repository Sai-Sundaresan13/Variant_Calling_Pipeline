# Running HISAT2 Alignment for Trimmed Files
High-quality reads are aligned to the reference genome using HISAT2. This step produces SAM files representing mapped read positions. The reference genome used for this project was Human Genome GRCh 38.

## 1. Preparing the Environment

Open the terminal and navigate to the directory containing all the trimmed FASTQ files:

```
cd /path/to/trimmed_output
```

## 2. Test if HISAT2 is Working

Run the following command to check if HISAT2 is functioning correctly on a single file:

```
hisat2 -x /home/sundar/project_PU/trimmed_output_301/GRCh38_index \
-U SRR31632778_trimmed.fastq.gz \
-S output.sam \
--threads 8
```

## 3. Running HISAT2 on Multiple Single-End Files

To align multiple single-end FASTQ files using a loop:

```
for file in *.fastq.gz; do
    base=$(basename $file .fastq.gz)
    hisat2 -x /home/sundar/project_PU/trimmed_output_301/GRCh38_index \
    -U $file \
    -S ${base}.sam
done
```

## 4. Running HISAT2 on Paired-End Files

Use the following command to align paired-end reads:

```
hisat2 -x grch38_index \
-1 ERR10851963_1_paired.fq.gz \
-2 ERR10851963_2_paired.fq.gz \
--threads 4 \
-S sample.bam
```
