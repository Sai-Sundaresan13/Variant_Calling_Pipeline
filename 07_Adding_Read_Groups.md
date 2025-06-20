# Adding Read Groups to BAM Files Prior to Base Quality Recalibration

Read groups must be added to BAM files before running tools like GATK BaseRecalibrator. This is done using **Picard's `AddOrReplaceReadGroups`** tool.

## Command

Use the following loop to add read groups to all BAM files:

```
for bam in *.bam; do  
    picard AddOrReplaceReadGroups \
        I="$bam" \
        O="${bam%.bam}_RG.bam" \
        RGID="${bam%.bam}" \
        RGLB="lib1" \
        RGPL="ILLUMINA" \
        RGPU="unit1" \
        RGSM="${bam%.bam}"
done
```

> This will generate a new BAM file (`*_RG.bam`) with read group information required for downstream tools like GATK.
