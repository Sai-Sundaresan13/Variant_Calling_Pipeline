# Separating SNPs from Indel Variants Using BCFtools

SNPs were separated and retrieved from the variants list. Two separate files were generated — one containing only **SNPs** and the other containing only **Indels** — for downstream analysis.

## Command

Use the following loop to process each `.vcf` file:

```
for file in *.vcf; do
    bcftools view -v snps "$file" -o "${file%.vcf}_snps.vcf"
    bcftools view -v indels "$file" -o "${file%.vcf}_indels.vcf"
done
```

This command will generate two output files per sample:
- `${sample}_snps.vcf` containing only SNPs
- `${sample}_indels.vcf` containing only Indels

> Make sure the input VCFs are properly formatted and indexed if needed.
