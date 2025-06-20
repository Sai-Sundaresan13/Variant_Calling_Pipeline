# Annotating chr17 VCF Files Against dbSNP to Add Gene Information

SnpSift is used to annotate VCF files with gene and variant information from a reference dbSNP VCF file. Gene-related information is added to each variant for biological interpretation.


## 1. Command for Annotating a Single File

```
java -jar /home/sundar/snpEff/SnpSift.jar annotate \
    filtered_chr17_dbsnp_146.hg38.vcf.gz \
    filtered_chr17_SRR31632778_trimmed_dedup_RG_recalibrated_filtered.vcf.gz \
    > annotated_dbSNP.vcf
```

## 2. Command for Annotating Multiple Files

Use this loop to annotate all `chr17`-filtered VCF files in the directory:

```
for vcf in *.vcf.gz; do
    output="$(basename "$vcf" .vcf.gz)_annotated.vcf"
    java -jar /home/sundar/snpEff/SnpSift.jar annotate \
        filtered_chr17_dbsnp_146.hg38.vcf.gz "$vcf" > "$output"
done
```

> Ensure the dbSNP file (`filtered_chr17_dbsnp_146.hg38.vcf.gz`) is indexed and properly formatted for annotation.
