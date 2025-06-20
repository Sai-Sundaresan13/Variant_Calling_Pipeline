# Converting SAM Files to BAM, Sorting, and Indexing
SAM files are converted to sorted and indexed BAM files using SAMtools. BAM files are more efficient for storage and further processing.

## 1. Navigate to the Directory

Open the terminal and navigate to the directory containing all your `.sam` files:

```
cd /path/to/sam_files
```


## 2. One-Liner Command: Convert, Sort, and Index in One Go

To perform all three steps (convert to BAM, sort, and index) in a single command:

```
for file in *.sam; do
    samtools view -bS "$file" | samtools sort -o "${file%.sam}.bam"
    samtools index "${file%.sam}.bam"
done
```

## 3. Running Commands Separately

If you prefer to do each step individually, use the following:

### a. Convert SAM to BAM:

```
for file in *.sam; do
    samtools view -bS "$file" -o "${file%.sam}.bam"
done
```

### b. Sort BAM Files:

```
for file in *.bam; do
    samtools sort "$file" -o "${file%.bam}_sorted.bam"
done
```

### c. Index Sorted BAM Files:

```
for file in *_sorted.bam; do
    samtools index "$file"
done
```
