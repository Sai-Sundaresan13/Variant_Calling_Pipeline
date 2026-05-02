# ========== USER SETTINGS ==========
ADAPTER_PATH="path/to/TruSeq3-PE-2.fa"                        # Path to Trimmomatic adapter file
REF="path/to/hg38.fa"                                         # Path to reference FASTA
KNOWN_SITES="path/to/Homo_sapiens_assembly38.dbsnp138.vcf"    # Path to known sites VCF for BQSR
SNPSIFT_JAR="path/to/snpEff/SnpSift.jar"                      # Path to SnpSift jar
DBSNP="path/to/dbsnp_146.hg38.vcf.gz"                         # Path to dbSNP VCF
THREADS=8                                                     # Number of threads
# ====================================


# 1. FastQC - Quality Control
mkdir -p 1_fastqc_output
fastqc * -o 1_fastqc_output/


# 2. MultiQC - Aggregate QC Report
mkdir -p 2_multiqc_output
multiqc 1_fastqc_output/ --interactive -o 2_multiqc_output/


# 3. Trimmomatic - Adapter Trimming (Paired-End)
mkdir -p 3_trimmed_output
for file in *_R1.fastq.gz; do
    base=$(basename "$file" _R1.fastq.gz)
    trimmomatic PE -threads "${THREADS}" -phred33 \
        "${base}_R1.fastq.gz" "${base}_R2.fastq.gz" \
        "3_trimmed_output/${base}_R1_paired.fq.gz" "3_trimmed_output/${base}_R1_unpaired.fq.gz" \
        "3_trimmed_output/${base}_R2_paired.fq.gz" "3_trimmed_output/${base}_R2_unpaired.fq.gz" \
        ILLUMINACLIP:"${ADAPTER_PATH}":2:30:10 \
        SLIDINGWINDOW:4:20 MINLEN:36
done


# 4. HISAT2 - Alignment
cd 3_trimmed_output
# download the reference file and index it
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
hisat2-build hg38.fa hg38_index

# run the command for alignment. sort and index the .bam files
mkdir -p 4_alignment_output
for file in *_R1_paired.fq.gz; do
    base=$(basename "$file" _R1_paired.fq.gz)
    hisat2 -x hg38_index \
        -1 "${base}_R1_paired.fq.gz" \
        -2 "${base}_R2_paired.fq.gz" \
        --threads "${THREADS}" | \
    samtools view -bS | \
    samtools sort -o "4_alignment_output/${base}.bam" && \
    samtools index "4_alignment_output/${base}.bam"
done


# 5. Picard - Duplicate Removal
cd 4_alignment_output
mkdir -p 5_picard_output
for file in *.bam; do
    base=$(basename "$file" .bam)
    picard MarkDuplicates \
        I="$file" \
        O="5_picard_output/${base}_dedup.bam" \
        M="5_picard_output/${base}_metrics.txt" \
        REMOVE_DUPLICATES=true
    samtools index "5_picard_output/${base}_dedup.bam"
done


# 6. Picard - Add Read Groups
mkdir -p 6_readgroups_output
for bam in 5_picard_output/*_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    picard AddOrReplaceReadGroups \
        I="$bam" \
        O="6_readgroups_output/${base}_RG.bam" \
        RGID="${base}" \
        RGLB="lib1" \
        RGPL="ILLUMINA" \
        RGPU="unit1" \
        RGSM="${base}"
done


# 7. GATK - Base Quality Score Recalibration (BQSR)
mkdir -p 7_bqsr_output
for bam in 6_readgroups_output/*_RG.bam; do
    base=$(basename "$bam" _RG.bam)
    gatk BaseRecalibrator \
        -R "${REF}" \
        -I "$bam" \
        --known-sites "${KNOWN_SITES}" \
        -O "7_bqsr_output/${base}_recal.table"
done

for bam in 6_readgroups_output/*_RG.bam; do
    base=$(basename "$bam" _RG.bam)
    gatk ApplyBQSR \
        -I "$bam" \
        -R "${REF}" \
        --bqsr-recal-file "7_bqsr_output/${base}_recal.table" \
        -O "7_bqsr_output/${base}_recalibrated.bam"
done


# 8. GATK - Variant Calling with HaplotypeCaller
mkdir -p 8_haplotypecaller_output
for bam in 7_bqsr_output/*_recalibrated.bam; do
    base=$(basename "$bam" _recalibrated.bam)
    gatk HaplotypeCaller \
        -R "${REF}" \
        -I "$bam" \
        -O "8_haplotypecaller_output/${base}.vcf.gz" \
        --native-pair-hmm-threads "${THREADS}"
done


# 9. GATK - Variant Filtration
mkdir -p 9_variantfiltration_output
for vcf in 8_haplotypecaller_output/*.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    gatk VariantFiltration \
        -R "${REF}" \
        -V "$vcf" \
        -O "9_variantfiltration_output/${base}_filtered.vcf.gz" \
        --filter-name "QD_filter" --filter-expression "QD < 2.0" \
        --filter-name "FS_filter" --filter-expression "FS > 60.0" \
        --filter-name "MQ_filter" --filter-expression "MQ < 40.0"
done


# 10. SnpSift - Variant Annotation
mkdir -p 10_annotation_output
for vcf in 9_variantfiltration_output/*_filtered.vcf.gz; do
    base=$(basename "$vcf" _filtered.vcf.gz)
    java -jar "${SNPSIFT_JAR}" annotate \
        "${DBSNP}" "$vcf" \
        > "10_annotation_output/${base}_annotated.vcf"
done


# 11. BCFtools - Separate SNPs and Indels
mkdir -p 11_snps_output
mkdir -p 11_indels_output
for file in 10_annotation_output/*.vcf; do
    base=$(basename "$file" .vcf)
    bcftools view -v snps "$file" -o "11_snps_output/${base}_snps.vcf" # filter SNPs
    bcftools view -v indels "$file" -o "11_indels_output/${base}_indels.vcf" # filter Indels
done
