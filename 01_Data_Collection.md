# Downloading Files from ENA

## Description
This step involves downloading `.fastq` or `.sra` files from the European Nucleotide Archive (ENA) using `wget` links.

## Requirements
- A file (e.g., `download_links.sh`) containing all `wget` links.
- Executable permission on the shell script.

## Steps

1. **Move to the directory** containing the download script:
   ```
   cd /path/to/your/directory
   ```

2. **Make the script executable**:
   ```
   chmod +x download_links.sh
   ```

3. **Run the script to start downloading**:
   ```
   ./download_links.sh
   ```

## Notes
- Make sure your script uses valid `wget` links for ENA `.fastq.gz` or `.sra` files.
- You can monitor progress using tools like `htop` or `screen` for long downloads.

