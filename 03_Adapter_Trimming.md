# Running Trimmomatic for FastQC Results
Adapters and low-quality bases are removed using Trimmomatic. This step improves read quality and enhances the accuracy of alignment.

## 1. Running Trimmomatic on Single-End Files

Basic usage for trimming a single-end FASTQ file:

```
java -jar /opt/trimmomatic/trimmomatic-0.39.jar SE -phred33 filename.fastq.gz output_trimmed.fastq.gz \
ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```


Batch-processing all `.fastq.gz` files in a directory (used in this project):

```
mkdir trimmed_output

for file in *.fastq.gz; do
    base=$(basename "$file" .fastq.gz)
    trimmomatic SE -phred33 "$file" "trimmed_output/${base}_trimmed.fastq.gz" \
    ILLUMINACLIP:/path/to/adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
done
```

## 2. Running Trimmomatic on Paired-End Files

```
mkdir -p trimmed_output  # Create output directory

for file in *.fastq.gz; do
    base=$(basename "$file" .fastq.gz)  # Extract sample name

    trimmomatic PE -threads 4 -phred33 \
        "${base}.fastq.gz" "${base}.fastq.gz" \
        "trimmed_output/${base}_paired.fq.gz" "trimmed_output/${base}_unpaired.fq.gz" \
        "trimmed_output/${base}_paired.fq.gz" "trimmed_output/${base}_unpaired.fq.gz" \
        ILLUMINACLIP:/home/sundar/trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:20 MINLEN:36
done
```


