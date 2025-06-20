# Running FastQC for All Files and MultiQC
Quality control is performed to assess the read quality, GC content, Per base sequence content, adapter content, etc. To perform quality control using FastQC and summarize the results using MultiQC, follow these steps:


Step 1: Navigate to the directory containing the sequencing files
```
cd /path/to/your/sequencing/files
```
Step 2: Run FastQC on all files in the directory
```
fastqc *
```

Step 3: Run MultiQC to aggregate FastQC results interactively
```
multiqc * --interactive
```

