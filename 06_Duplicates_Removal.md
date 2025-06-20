# Removing Duplicates from BAM Files After Alignment

Duplicate reads are marked and removed using Picard to reduce PCR bias. This ensures more reliable variant calling in downstream steps.

## 1. Tool Used

Duplicate reads in BAM files are removed using **Picard**'s `MarkDuplicates` tool.

## 2. Command to Remove Duplicates and Index BAM Files

Use the following loop to process all `.bam` files in the directory:

```
for file in *.bam; do
    picard MarkDuplicates \
        I="$file" \
        O="${file%.bam}_dedup.bam" \
        M="${file%.bam}_metrics.txt" \
        REMOVE_DUPLICATES=true

    samtools index "${file%.bam}_dedup.bam"
done
```

> This loop removes duplicates from each BAM file, generates a deduplicated output, writes duplication metrics to a text file, and creates an index for the deduplicated BAM file.
