
# RNA-seq Analysis Pipeline (Appendix)
# ----------------------------------------------------------
# NOTE: RNA-seq workflow was not completed and is not 
# part of the main thesis. It is included here as 
# supplementary code to support the appendix B part in main thesis.

# ----------------------------------------------------------
# Step 1: Build STAR genome index
# Original script: star genome.sh
# Prepares genome for alignment with STAR
# ----------------------------------------------------------
#SBATCH --job-name=star_index_build
#SBATCH --output=star_index_build_%j.out
#SBATCH --error=star_index_build_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00

module load STAR/2.7.0f

genome_fa="/Data3/samiksha/Haem_PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa"
gff3="/Data3/samiksha/Haem_PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3"
outdir="/Data3/samiksha/Haem_PRJEB506/star_index"

mkdir -p "$outdir"

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$outdir" \
     --genomeFastaFiles "$genome_fa" \
     --sjdbGTFfile "$gff3" \
     --sjdbGTFfeatureExon exon \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbOverhang 99


# ----------------------------------------------------------
# Step 2: Align RNA-seq reads with STAR
# Original script: rna align.sh
# Maps the RNA-seq reads to genome to generate sorted BAM files.
# ----------------------------------------------------------
#SBATCH --job-name=rnaseq_align_star
#SBATCH --output=star_align_%A_%a.out
#SBATCH --error=star_align_%A_%a.err
#SBATCH --array=0-17 
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00

module load STAR/2.7.0f
module load samtools/1.12

# Sample directories
samples=($(ls -d /Data3/samiksha/Trimmed/Sample_*))
sample_dir="${samples[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$sample_dir")

r1=$(find "$sample_dir" -name '*_R1*.fastq.gz' | head -n 1)
r2=$(find "$sample_dir" -name '*_R2*.fastq.gz' | head -n 1)

# STAR genome index
genome_dir="/Data3/samiksha/Haem_PRJEB506/star_index"

# Output directory
outdir="/Data3/samiksha/RNAseq_alignments_star"
mkdir -p "$outdir"

# Run STAR
echo "Running STAR alignment for sample $sample_name"

STAR --runThreadN 8 \
     --genomeDir "$genome_dir" \
     --readFilesIn "$r1" "$r2" \
     --readFilesCommand zcat \
     --outFileNamePrefix "$outdir/${sample_name}." \
     --outSAMtype BAM SortedByCoordinate

# Index the sorted BAM
samtools index "$outdir/${sample_name}.Aligned.sortedByCoord.out.bam"


# ----------------------------------------------------------
# Step 3: Quantify gene counts with HTSeq
# Original script: rna_seq_counts.sh
# Quantifies read counts per gene from aligned BAM files and outputs a counts file for each sample.
# ----------------------------------------------------------
#SBATCH --job-name=htseq_counts
#SBATCH --output=htseq_counts_%j.out
#SBATCH --error=htseq_counts_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4

module load htseq-count/0.11.0 samtools/1.17

# Paths
gff3=/Data3/samiksha/RNAseq_alignments_star/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3
bamdir=/Data3/samiksha/RNAseq_alignments_star
outdir=/Data3/samiksha/RNAseq_counts_htseq

mkdir -p $outdir

# Loop over BAMs
for bam in $bamdir/*.Aligned.sortedByCoord.out.bam
do
    sample=$(basename $bam .Aligned.sortedByCoord.out.bam)

    echo "Processing sample: $sample"
    htseq-count -f bam -r pos -s no -t exon -i Parent $bam $gff3 > $outdir/${sample}.counts.txt
done
