ADAPTER_PATH="/home/sundar/trimmomatic_adapters/TruSeq3-PE-2.fa"  # Change this to your adapter file path for trimmomatic runs
THREADS=8 
# run this script after entering the directory containing your files

# 1. Fastqc analysis for quality control
mkdir 1_fastqc_output # run only once to create the folder
fastqc * -o 1_fastqc_output/ # saves the output files into the folder

# 2. Multiqc
mkdir 2_multiqc_output
multiqc * --interactive -o 2_multiqc_output/

# 3. Adapter trimming using trimmomatic (paired end data) - change parameters as per usage
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

# 4. Alignment using HISAT2
cd 3_trimmed_output
# download the reference fasta file and index it
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
hisat2-build GRCh38.fasta GRCh38_index # indexing the reference

mkdir 4_alignment_output
# run the alignment save the .bam files in an output folder
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

# 5. Duplicate removal using PICARD
cd  4_alignment_output
mkdir 5_picard_output

for file in *.bam; do
    base=$(basename "$file" .bam)
    picard MarkDuplicates \
        I="$file" \
        O="5_picard_output/${base}_dedup.bam" \
        M="5_picard_output/${base}_metrics.txt" \
        REMOVE_DUPLICATES=true

    samtools index "5_picard_output/${base}_dedup.bam"
done

# 6. Add read groups
# Step 1: Add Read Groups
mkdir 6_readgroups_output

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

REF="path/to/hg38.fa"                          # <-- Change to your reference path
KNOWN_SITES="path/to/Homo_sapiens_assembly38.dbsnp138.vcf"

# Step 2: Generate Recalibration Tables
mkdir -p 7_bqsr_output

for bam in 6_readgroups_output/*_RG.bam; do
    base=$(basename "$bam" _RG.bam)
    gatk BaseRecalibrator \
        -R "${REF}" \
        -I "$bam" \
        --known-sites "${KNOWN_SITES}" \
        -O "7_bqsr_output/${base}_recal.table"
done


# Step 3: Apply BQSR
for bam in 6_readgroups_output/*_RG.bam; do
    base=$(basename "$bam" _RG.bam)
    gatk ApplyBQSR \
        -I "$bam" \
        -R "${REF}" \
        --bqsr-recal-file "7_bqsr_output/${base}_recal.table" \
        -O "7_bqsr_output/${base}_recalibrated.bam"
done

#################################################
REF="path/to/hg38.fa"  # <-- Change to your reference path
THREADS=8               # <-- Change to your number of threads
# ====================================

mkdir 8_haplotypecaller_output

for bam in 7_bqsr_output/*_recalibrated.bam; do
    base=$(basename "$bam" _recalibrated.bam)
    gatk HaplotypeCaller \
        -R "${REF}" \
        -I "$bam" \
        -O "8_haplotypecaller_output/${base}.vcf.gz" \
        --native-pair-hmm-threads "${THREADS}"
done
############################################
REF="path/to/hg38.fa"  # <-- Change to your reference path
# ====================================

mkdir  9_variantfiltration_output

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

#################################################
SNPSIFT_JAR="path/to/snpEff/SnpSift.jar"      # <-- Change to your SnpSift jar path
DBSNP="path/to/dbsnp_146.hg38.vcf.gz"         # <-- Change to your dbSNP VCF path
# ====================================

mkdir 10_annotation_output

for vcf in 9_variantfiltration_output/*_filtered.vcf.gz; do
    base=$(basename "$vcf" _filtered.vcf.gz)
    java -jar "${SNPSIFT_JAR}" annotate \
        "${DBSNP}" "$vcf" \
        > "10_annotation_output/${base}_annotated.vcf"
done

##########################################
mkdir -p 11_snps_output
mkdir -p 11_indels_output

for file in 10_annotation_output/*.vcf; do
    base=$(basename "$file" .vcf)
    bcftools view -v snps "$file" -o "11_snps_output/${base}_snps.vcf"
    bcftools view -v indels "$file" -o "11_indels_output/${base}_indels.vcf"
done
