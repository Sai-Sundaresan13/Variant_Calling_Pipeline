# 1. Filtering Only Chromosome 17 Variants from All Samples

Filtration was done to extract only the variants located on chromosome 17 (`chr17`) from each VCF file using `bcftools`.

## Command

Use the following loop to filter all `.vcf.gz` files and retain only `chr17` variants:

```
for file in *.vcf.gz; do
    bcftools view -r chr17 "$file" -o "filtered_chr17_${file}" -O z
done
```

> Make sure your reference and VCF files use `chr17` naming (not just `17`) to match correctly.


# 2. Filtering Variants for MAPT Gene Using BCFtools

Filtration was done to extract only those variants from `chr17` VCF files that are associated with the **MAPT** gene, based on the `GENEINFO` field in the INFO column.

## Command to Filter Multiple Files

Use the following loop to filter for variants containing `"MAPT"` in the `INFO/GENEINFO` field:

```
for file in filtered_chr17_*.vcf; do
    bcftools view -i 'INFO/GENEINFO~"MAPT"' "$file" > "${file%.vcf}_MAPT.vcf"
done
```

> Make sure that the input VCF files contain the `GENEINFO` tag and that gene annotations have been added (e.g., using SnpSift).


# 3. Filtering Variants Based on Chromosomal Position Using BCFtools

Filtration was done to extract variants that lie within a specific region of a chromosome (e.g., chr17:40,000,000–45,000,000) from compressed and indexed VCF files.

---

## For a Single VCF File

```
bcftools view -r chr17:40000000-45000000 -Oz -o output_filtered.vcf.gz input.vcf.gz
```

**Explanation:**
- `-r chr17:40000000-45000000` → Filters the region from 40M to 45M on chromosome 17
- `-Oz` → Outputs a bgzipped VCF
- `-o` → Specifies the output filename

---

## For Multiple VCF Files

```
for vcf in *.vcf.gz; do
    output="${vcf%.vcf.gz}_chr17_40M_45M.vcf.gz"
    bcftools view -r chr17:40000000-45000000 -Oz -o "$output" "$vcf"
    tabix -p vcf "$output"
done
```

This loop:
- Filters each `.vcf.gz` file for chr17:40,000,000–45,000,000
- Writes the filtered variants to a new compressed VCF
- Indexes the new file using `tabix`
