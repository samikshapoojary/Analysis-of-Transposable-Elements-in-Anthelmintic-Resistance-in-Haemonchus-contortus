#!/bin/bash

############################################################
# Master Bash Script: TE Analysis in Haemonchus contortus
# All Bash scripts used in the analysis
############################################################


# -----------------------------------------------------------
# Step 1: FASTQ header / format fixing
# Original script: fix_headers.sh
# Fixes inconsistent FASTQ headers by replacing spaces and 
# ensuring proper read numbering. Runs in parallel batches.
# -----------------------------------------------------------

#SBATCH --job-name=fix_fastq_headers
#SBATCH --output=fix_fastq_headers_%j.log
#SBATCH --error=fix_fastq_headers_%j.err
#SBATCH --cpus-per-task=36

mkdir -p fixed_fastq

fix_headers() {
  fq1="$1"
  base="${fq1%_R1.fastq.gz}"
  fq2="${base}_R2.fastq.gz"

  out1="fixed_fastq/${base}_R1.fixed.fastq.gz"
  out2="fixed_fastq/${base}_R2.fixed.fastq.gz"

  echo "Fixing headers in $fq1 and $fq2"

  zcat "$fq1" | awk 'NR%4==1 {
    split($0, a, " ");
    if (length(a) > 1) {
      readnum = substr(a[2], 1, 1);
      print a[1] "/" readnum;
    } else {
      print $0;
    }
    next;
  } { print }' | gzip > "$out1" &

  zcat "$fq2" | awk 'NR%4==1 {
    split($0, a, " ");
    if (length(a) > 1) {
      readnum = substr(a[2], 1, 1);
      print a[1] "/" readnum;
    } else {
      print $0;
    }
    next;
  } { print }' | gzip > "$out2" &
}

max_jobs=18
current_jobs=0
for fq1 in *_R1.fastq.gz; do
  fix_headers "$fq1"
  ((current_jobs++))
  if (( current_jobs >= max_jobs )); then
    wait
    current_jobs=0
  fi
done
wait


# -----------------------------------------------------------
# Step 2: EarlGrey TE library generation
# Original script: earlgrey.sh
# Runs EarlGrey to annotate transposable elements in genome.
# -----------------------------------------------------------

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --job-name=earlGrey

module load earlGrey/6.0.3_dfam39

genome_file="../Haem_PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa"
annot_file="../Haem_PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3"

output_dir="../Haem_PRJEB506/earlgrey_output2"
mkdir -p "$output_dir"

rm -f genome.fa annotations.gff3
ln -s "$genome_file" genome.fa
ln -s "$annot_file" annotations.gff3

earlGrey -g "$genome_file" -s Hcontortus -t 24 -e yes -c yes -o "$output_dir"


# -----------------------------------------------------------
# Step 3: Filter and prepare TE-only annotations
# Original script: filter_only_te.sh
# Removes non-TE features to retain only TE annotations.
# -----------------------------------------------------------

#SBATCH --job-name=filter_TEs
#SBATCH --output=filter_TEs.out
#SBATCH --error=filter_TEs.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4

awk '$3 != "Satellite" && $3 != "Simple_repeat" && $3 != "Low_complexity" && $3 != "Unknown"' \
  Hcontortus.filteredRepeats.gff > Hcontortus.TEonly.gff


# -----------------------------------------------------------
# Step 4: Merge genome + TE FASTA
# Original script: new_only_te_merged.sh
# Creates merged reference combining masked genome and TE sequences.
# -----------------------------------------------------------

#!/bin/bash
#SBATCH --job-name=merged_te_process
#SBATCH --output=merged_only_te_%j.log
#SBATCH --error=merged_only_te_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4

module load bedtools/2.29.2
module load samtools/1.9

genome="/Data3/samiksha/Haem_PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa"
gff="/Data3/samiksha/Hcontortus_EarlGrey/Hcontortus_summaryFiles/Hcontortus.filteredRepeats.gff"
outdir="/Data3/samiksha/Hcontortus_EarlGrey/merged"
mkdir -p "${outdir}"

filtered_gff="${outdir}/Hcontortus.TEonly.gff"
filtered_bed="${outdir}/Hc_TE_annotation.bed"
filtered_bed_1kb="${outdir}/Hc_TE_annotation.1kb.bed"
te_fasta="${outdir}/hc_only_te_seqs.1kb.fasta"
masked_fasta="${outdir}/hc.only_te.masked_genome.fasta"
merged_fasta="${outdir}/hc.only_te.temergedref.fasta"

awk '$3 != "Satellite" && $3 != "Simple_repeat" && $3 != "Low_complexity" && $3 != "Unknown"' "${gff}" > "${filtered_gff}"
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $3, $6, $7}' "${filtered_gff}" > "${filtered_bed}"
awk '{ if ($3 - $2 >= 1000) print $0 }' "${filtered_bed}" > "${filtered_bed_1kb}"

if [[ "${genome}.fai" -ot "${genome}" ]]; then
    samtools faidx "${genome}"
fi

bedtools getfasta -s -name -fi "${genome}" -bed "${filtered_bed_1kb}" -fo "${te_fasta}"
bedtools maskfasta -fi "${genome}" -bed "${filtered_bed_1kb}" -fo "${masked_fasta}"
cat "${masked_fasta}" "${te_fasta}" > "${merged_fasta}"


# -----------------------------------------------------------
# Step 5: Create TE hierarchy file
# Original script: hier_only_te.sh
# Builds hierarchy file (id/family/order) from TE FASTA headers.
# -----------------------------------------------------------

#!/bin/bash
#SBATCH --job-name=make_te_hierarchy
#SBATCH --output=make_te_hierarchy_%j.out
#SBATCH --error=make_te_hierarchy_%j.err
#SBATCH --cpus-per-task=4

cd /Data3/samiksha/Hcontortus_EarlGrey/merged || exit 1

merged_fasta="hc.only_te.temergedref.fasta"
output="Hcontortus_only_TEhierarchy.txt"

(
    echo -e "id\tfamily\torder"
    grep '^>' "$merged_fasta" | sed 's/^>//' | grep '::' | \
    awk -F'::' '{
        id = $0
        split($1, parts, "/")
        if (length(parts)==2) {
            order = parts[1]
            family = parts[2]
        } else {
            order = parts[1]
            family = parts[1]
        }
        print id "\t" family "\t" order
    }'
) > "$output"


# -----------------------------------------------------------
# Step 6: PoPoolationTE2 individual sample analysis
# Original scripts: new bowtie.sh + bowtie popte.sh
# Maps reads with Bowtie2, processes BAMs, runs PoPoolationTE2.
# -----------------------------------------------------------

#!/bin/bash
module load samtools/1.9
module load bowtie2/2.3.4.3

samtools=$(which samtools)
bowtie2=$(which bowtie2)
popte2="/Data3/samiksha/tools/popte2-v1.10.03.jar"
refg="/Data3/samiksha/Hcontortus_EarlGrey/merged/hc.only_te.temergedref.fasta"
hier="/Data3/samiksha/Hcontortus_EarlGrey/merged/Hcontortus_only_TEhierarchy.txt"
bowtie2_index="/Data3/samiksha/Hcontortus_EarlGrey/merged/hc_onlyte_bt2_index"

if [ ! -f "${bowtie2_index}.1.bt2" ]; then
    bowtie2-build "$refg" "$bowtie2_index"
fi

file1=$1
file2=$2
of=$3
mincount=$4

mkdir -p "$of"

$bowtie2 -x "$bowtie2_index" -1 "$file1" -2 "$file2" -p $SLURM_CPUS_PER_TASK | \
$samtools view -Sb - | \
$samtools sort -@ $SLURM_CPUS_PER_TASK -o "$of/map-pe.sort.bam"

$samtools index "$of/map-pe.sort.bam"

java -jar $popte2 ppileup --bam "$of/map-pe.sort.bam" --map-qual 15 --hier "$hier" --output "$of/pp.gz"
java -jar $popte2 identifySignatures --ppileup "$of/pp.gz" --mode separate --min-count "$mincount" --output "$of/te.signatures"
java -jar $popte2 updatestrand --bam "$of/map-pe.sort.bam" --signature "$of/te.signatures" --output "$of/testrand.signatures" --hier "$hier" --map-qual 15 --max-disagreement 0.4
java -jar $popte2 frequency --ppileup "$of/pp.gz" --signature "$of/testrand.signatures" --output "$of/te.freqsignatures"
java -jar $popte2 pairupsignatures --signature "$of/te.freqsignatures" --ref-genome "$refg" --hier "$hier" --min-distance -200 --max-distance 300 --output "$of/tes.finalresult.txt"

#SBATCH --job-name=bowtie_popte2
#SBATCH --array=1-9
#SBATCH --cpus-per-task=32
#SBATCH --output=onlyte_bowtie_%A_%a.out
#SBATCH --error=onlyte_bowtie_%A_%a.err

module load samtools/1.9

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples_bowtie2.txt)
set -- $line
fastq1="$1"
fastq2="$2"
outdir="$3"
mincount="$4"
./new_bowtie2.sh "$fastq1" "$fastq2" "$outdir" "$mincount"


# -----------------------------------------------------------
# Step 7: PoPoolationTE2 joint analysis
# Original scripts: bamlist.sh + final joint.sh
# Prepares BAMs, builds bamlist, runs pooled joint analysis.
# -----------------------------------------------------------

#!/bin/bash
#SBATCH --job-name=rename_bams
#SBATCH --output=rename_bams.out
#SBATCH --error=rename_bams.err
#SBATCH --cpus-per-task=8

final="/Data3/samiksha/mapped_new/final/bowtie"
final_out="/Data3/samiksha/mapped_new/final_joint_bams"

mkdir -p  "$final_out"

for dir in "$final"/*; do
    if [[ -d "$dir" ]]; then
        sample=$(basename "$dir")
        bam="$dir/map-pe.sort.bam"
        bai="$dir/map-pe.sort.bam.bai"
        if [[ -f "$bam" ]]; then cp "$bam" "$final_out/${sample}.bam"; fi
        if [[ -f "$bai" ]]; then cp "$bai" "$final_out/${sample}.bam.bai"; fi
    fi
done

#!/bin/bash
#SBATCH --job-name=final_joint_popte2
#SBATCH --output=final_joint_popte2.out
#SBATCH --error=final_joint_popte2.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32

module load samtools/1.9
module load bowtie2/2.3.4.3

popte2="/Data3/samiksha/tools/popte2-v1.10.03.jar"
refg="/Data3/samiksha/Hcontortus_EarlGrey/merged/hc.only_te.temergedref.fasta"
hier="/Data3/samiksha/Hcontortus_EarlGrey/merged/cleaned_te_hierarchy.txt"
bamdir="/Data3/samiksha/mapped_new/final_joint_bams"
outdir="/Data3/samiksha/joint_popte2_final_output"

mkdir -p "$outdir"

cd $bamdir
bam_args=$(ls *.bam | awk '{print "--bam "$1}' | tr '\n' ' ')

java -jar $popte2 ppileup $bam_args --map-qual 15 --hier $hier --output "$outdir/joint.pp.gz"
java -jar $popte2 identifySignatures --ppileup "$outdir/joint.pp.gz" --mode joint --min-count 2 --signature-window fix100 --min-valley fix100 --output "$outdir/joint.te.signatures"
java -jar $popte2 updatestrand $bam_args --signature "$outdir/joint.te.signatures" --output "$outdir/joint.testrand.signatures" --hier $hier --map-qual 15 --max-disagreement 0.4
java -jar $popte2 frequency --ppileup "$outdir/joint.pp.gz" --signature "$outdir/joint.testrand.signatures" --output "$outdir/joint.te.freqsignatures"
java -jar $popte2 pairupsignatures --signature "$outdir/joint.te.freqsignatures" --ref-genome $refg --hier $hier --min-distance -200 --max-distance 300 --output "$outdir/joint.tes.finalresult.txt"
