# Running Base Quality Score Recalibration (BQSR)

Base Quality Score Recalibration (BQSR) is performed using GATK. It corrects systematic errors made by the sequencer during base calling.

## 1. Generating Recalibration Tables

Use the following loop to generate recalibration tables for all BAM files with read groups (`*_RG.bam`):

```
for bam in *_RG.bam; do
    gatk BaseRecalibrator \
        -R nw/hg38.fa \
        -I "$bam" \
        --known-sites nw/Homo_sapiens_assembly38.dbsnp138.vcf \
        -O "${bam%.bam}_recal.table"
done
```

## 2. Applying BQSR

### Applying BQSR to a Single File

```
gatk ApplyBQSR \
    -I SRR31632778_trimmed_dedup_RG.bam \
    -R reference.fasta \
    --bqsr-recal-file SRR31632778_trimmed_dedup_RG_recal_data.table \
    -O SRR31632778_recalibrated.bam
```

### Applying BQSR to All Files

To apply BQSR to multiple BAM files in a batch:

```
for bam in *.bam; do
    gatk ApplyBQSR \
        -I "$bam" \
        -R reference.fasta \
        --bqsr-recal-file "${bam%.bam}_recal_data.table" \
        -O "${bam%.bam}_recalibrated.bam"
done
```

> Ensure that the reference genome and known-sites VCF file paths match your environment.
