# Running GATK HaplotypeCaller for Variant Calling

GATK's HaplotypeCaller is used to perform variant calling on aligned and preprocessed BAM files. It provides robust SNP and indel detection with built-in error modeling.

## Command to Run HaplotypeCaller on Multiple BAM Files

Use the following loop to call variants from each BAM file and generate a compressed VCF (`.vcf.gz`) file:

```
for bam in *.bam; do
    sample_name=$(basename "$bam" .bam)
    gatk HaplotypeCaller \
        -R reference.fasta \
        -I "$bam" \
        -O "${sample_name}.vcf.gz" \
        --native-pair-hmm-threads 4
done
```

> Ensure that `reference.fasta` is the same reference genome used during alignment and preprocessing.
