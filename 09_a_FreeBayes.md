# Running FreeBayes on BAM Files

FreeBayes is a haplotype-based variant detector, used to call genetic variants (SNPs and indels) from aligned BAM files.

## 1. Run FreeBayes on a Single BAM File

Use the following command to call variants from a single sample:

```
freebayes -f reference.fasta SRR31632778_trimmed_dedup_RG_recalibrated.bam > raw_variants.vcf
```

## 2. Run FreeBayes on Multiple BAM Files

To perform variant calling on multiple BAM files in a loop:

```
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    freebayes -f reference.fasta "$bam" > "${sample}.vcf"
done
```

> Ensure that the reference genome file (`reference.fasta`) matches the one used during alignment and BQSR.
