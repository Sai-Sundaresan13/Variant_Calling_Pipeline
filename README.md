# Variant_Calling_Pipeline
This repository consists of the bash commands corresponding to Variant Calling, used in my master's thesis titled "A Study on the Impact of Variants on the Structure of Tau protein filaments in Frontotemporal Dementia". 

Here is a basic summary of all the steps used in this Variant Calling Pipeline: 

## 1. Data Collection:
Raw sequencing data is collected from public databases or sequencing facilities. This step ensures all input FASTQ files are organized and ready for downstream analysis.

## 2. Quality Control:
FASTQ files are evaluated using FastQC to assess read quality, GC content, and adapter contamination. Results are summarized with MultiQC for a comprehensive overview.

## 3. Adapter Trimming:
Adapters and low-quality bases are removed using Trimmomatic. This step improves read quality and enhances the accuracy of alignment.

## 4. Alignment:
High-quality reads are aligned to the reference genome using HISAT2. This step produces SAM files representing mapped read positions.

## 5. SAM to BAM:
SAM files are converted to sorted and indexed BAM files using SAMtools. BAM files are more efficient for storage and further processing.

## 6. Duplicates Removal:
Duplicate reads are marked and removed using Picard to reduce PCR bias. This ensures more reliable variant calling in downstream steps.

## 7. Adding Read Groups:
Read groups are added using Picard AddOrReplaceReadGroups. These are required for tools like GATK to differentiate samples and lanes.

## 8. Base Quality Score Recalibration:
Base Quality Score Recalibration (BQSR) is performed using GATK. It corrects systematic errors made by the sequencer during base calling.

## 9. Variant Calling:
   a. FreeBayes: Variants are called using FreeBayes, a haplotype-based variant detector. Output is a VCF file containing SNPs and indels.
   
   b. GATK: An alternative variant calling method using GATK HaplotypeCaller. It provides robust SNP and indel detection with built-in error modeling.

## 10. Quality Filtration:
Variants are filtered based on quality metrics like QD, FS, and MQ using GATK VariantFiltration. This step improves confidence in variant calls.

## 11. Annotation:
VCF files are annotated against databases like dbSNP using SnpSift. Gene-related information is added to each variant for biological interpretation.

## 12. Filtration:
Variants are filtered based on chromosome, gene (e.g., MAPT), and other criteria using BCFtools. This refines the dataset for focused analysis.

## 13. SNP Retrieval:
SNPs are separated from indels and further filtered by position or region. This results in a final high-confidence SNP dataset for interpretation.


## Directory Structure

```
Pipeline for variant calling/
│
├── 01_Data_Collection
├── 02_Quality_Control
├── 03_Adapter_Trimming
├── 04_Alignment
├── 05_SAM_to_BAM
├── 06_Duplicates_Removal
├── 07_Adding_Read_Groups
├── 08_Base_Quality_Score_Recalibration
├── 09_a_Freebayes_VC
├── 09_b_GATK_VC
├── 10_Quality_Filtration
├── 11_Annotation
├── 12_Filtration
├── 13_SNP_Retrieval

```
