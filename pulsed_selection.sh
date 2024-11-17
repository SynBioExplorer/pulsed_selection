#!/bin/bash

# Define paths and filenames
WORKDIR="/Volumes/Felix_SSD/PhD/Sequencing_Pulse_selection/N2413980_30-1043295145_Reseq_2024-08-01/raw_data/final"
ANALYSIS_DIR="$WORKDIR/analysis"
REFERENCE="$WORKDIR/BY4742_pRS416.fasta"
SNPEFF_DIR="/Users/MQ46079823/snpEff"
export JAVA_HOME="/opt/homebrew/opt/openjdk@21"
export PATH="/opt/homebrew/opt/openjdk@21/bin:$PATH"

# Check software dependencies
dependencies=("java" "fastp" "bwa" "samtools" "bcftools" "tabix")
for cmd in "${dependencies[@]}"; do
    command -v $cmd >/dev/null 2>&1 || { echo "Error: $cmd is not installed." >&2; exit 1; }
done

# Check Java version
java -version || { echo "Error: Java is not configured properly." >&2; exit 1; }

# Create main analysis directory and subdirectories
mkdir -p $ANALYSIS_DIR/{trimmed,alignments,variants,annotations,logs}
cd $WORKDIR

# Check and download SnpEff database if necessary
if [ ! -d "$SNPEFF_DIR/data/R64-1-1.99" ]; then
    echo "Downloading R64-1-1.99 database for snpEff..."
    java -jar $SNPEFF_DIR/snpEff.jar download R64-1-1.99 || { echo "Error: Failed to download SnpEff database." >&2; exit 1; }
fi

# Index the reference genome (only needed once)
if [ ! -f "$REFERENCE.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index $REFERENCE || { echo "Error: Failed to index reference genome with BWA." >&2; exit 1; }
    samtools faidx $REFERENCE || { echo "Error: Failed to index reference genome with samtools." >&2; exit 1; }
fi

# Process each paired-end sample
for file1 in *L2_1.fq.gz; do
    file2="${file1%_1.fq.gz}_2.fq.gz"
    sample="${file1%_L2_1.fq.gz}"  # Extract sample name

    echo "Processing sample: $sample" | tee "$ANALYSIS_DIR/logs/${sample}_log.txt"

    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "Error: Missing paired files for sample $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
        continue
    fi

    # Quality Control & Trimming
    fastp -i "$file1" -I "$file2" \
          -o "$ANALYSIS_DIR/trimmed/${sample}_trimmed_1.fq.gz" \
          -O "$ANALYSIS_DIR/trimmed/${sample}_trimmed_2.fq.gz" \
          -h "$ANALYSIS_DIR/trimmed/${sample}_fastp_report.html" \
          -j "$ANALYSIS_DIR/trimmed/${sample}_fastp_report.json" || {
              echo "Error: fastp failed for $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
              continue
          }

    # Alignment with BWA
    bwa mem $REFERENCE "$ANALYSIS_DIR/trimmed/${sample}_trimmed_1.fq.gz" "$ANALYSIS_DIR/trimmed/${sample}_trimmed_2.fq.gz" > "$ANALYSIS_DIR/alignments/${sample}_aligned.sam" || {
        echo "Error: Alignment failed for $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
        continue
    }

    # Convert SAM to BAM, Sort and Index
    samtools view -Sb "$ANALYSIS_DIR/alignments/${sample}_aligned.sam" > "$ANALYSIS_DIR/alignments/${sample}_aligned.bam"
    samtools sort "$ANALYSIS_DIR/alignments/${sample}_aligned.bam" -o "$ANALYSIS_DIR/alignments/${sample}_sorted.bam"
    samtools index "$ANALYSIS_DIR/alignments/${sample}_sorted.bam"

    # Remove intermediate SAM file to save space
    rm "$ANALYSIS_DIR/alignments/${sample}_aligned.sam"

    # Variant Calling with bcftools
    bcftools mpileup -f $REFERENCE "$ANALYSIS_DIR/alignments/${sample}_sorted.bam" | \
    bcftools call -mv -Oz -o "$ANALYSIS_DIR/variants/${sample}_variants.vcf.gz" || {
        echo "Error: Variant calling failed for $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
        continue
    }

    # Index the VCF file
    tabix -p vcf "$ANALYSIS_DIR/variants/${sample}_variants.vcf.gz"

    # Variant Filtering with bcftools
    bcftools filter -e 'INFO/DP < 10 || QUAL < 20' "$ANALYSIS_DIR/variants/${sample}_variants.vcf.gz" -Oz -o "$ANALYSIS_DIR/variants/${sample}_filtered_variants.vcf.gz" || {
        echo "Error: Variant filtering failed for $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
        continue
    }

    # Annotation with snpEff
    java -jar $SNPEFF_DIR/snpEff.jar R64-1-1.99 "$ANALYSIS_DIR/variants/${sample}_filtered_variants.vcf.gz" > "$ANALYSIS_DIR/annotations/${sample}_annotated_variants.vcf" || {
        echo "Error: Annotation failed for $sample." | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
        continue
    }

    echo "Completed processing for sample: $sample" | tee -a "$ANALYSIS_DIR/logs/${sample}_log.txt"
done

echo "All samples processed."
