## QC

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P9/raw_reads
for ea in *.fastq.gz; do fastqc $ea --outdir fastqc; done;
multiqc fastqc

## Merge fastqs across lanes

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/raw_reads
ids=$(ls -1 *R1*.gz | awk -F '_' '{print $4"_"$5"_"$6}' | sort | uniq)
for ea in $ids; do
echo GH_CHIP_8881_$ea\_L001_R1_001.fastq.gz GH_CHIP_8881_$ea\_L002_R1_001.fastq.gz \> GH_CHIP_8881_$ea\_R1.fastq.gz
cat GH_CHIP_8881_$ea\_L001_R1_001.fastq.gz GH_CHIP_8881_$ea\_L002_R1_001.fastq.gz > merged/GH_CHIP_8881_$ea\_R1.fastq.gz
echo GH_CHIP_8881_$ea\_L001_R2_001.fastq.gz GH_CHIP_8881_$ea\_L002_R2_001.fastq.gz \> GH_CHIP_8881_$ea\_R1.fastq.gz
cat GH_CHIP_8881_$ea\_L001_R2_001.fastq.gz GH_CHIP_8881_$ea\_L002_R2_001.fastq.gz > merged/GH_CHIP_8881_$ea\_R2.fastq.gz
done;

## STAR & FEATURECOUNTS

## Already made index:
## REMAKE INDEX WITH OVERHANG = 49

## Align (v. 2.7.8a)

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/raw_reads/merged
ref="/home/shared/hg_align_db/GRCm39_gencode_primary/star_index"
outdir="/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/STAR_aligned"
for R1 in GH_CHIP*R1*fastq.gz; do 
  R2=$(echo $R1 | sed "s/R1/R2/")
  STAR \
    --runThreadN 15 \
    --genomeDir $ref  \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn $R1 $R2 \
    --outFilterType BySJout \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMax 1000000 \
    --outFilterMultimapNmax 20 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outSAMmultNmax 1 \
    --sjdbOverhang 75 \
    --outFileNamePrefix ${outdir}/$(echo $R1 | sed "s/R1_001.fastq.gz//")
done;
multiqc .

## Align summary:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/STAR_aligned
annot="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
bam_list=$(ls *bam)
/home/shared/programs/qualimap_v2.2.1/./qualimap rnaseq \
  -a proportional \
  -outdir qualimap \
  -bam $bam_list \
  -pe -p strand-specific-forward \
  --java-mem-size=20G \
  -gtf $annot

## Read quant:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/P5/STAR_aligned
annot="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
featurecounts="/opt/subread-2.0.3-Linux-x86_64/bin/featureCounts"
for ea in *bam; do 
  $featurecounts \
    -T 15 -a $annot \
    -o ../featureCounts/$(echo $ea | sed "s/_Aligned.sortedByCoord.out.bam//") \
    -p -s 2 -M -O --fraction \
    --countReadPairs \
    --extraAttributes gene_name \
    $ea
done;
multiqc ../featureCounts
