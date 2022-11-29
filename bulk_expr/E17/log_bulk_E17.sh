## QC:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads
for ea in *.fastq.gz
do 
  sem -j 8 fastqc $ea --outdir fastqc
done;
cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/fastqc
multiqc .

## Trim:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads
ref="/home/shared/programs/bbmap/resources/adapters.fa"
outdir="/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed"
for R1 in XZ*R1*fastq.gz
  do
  R2=$(echo $R1 | sed "s/R1/R2/")
  /home/shared/programs/bbmap/bbduk.sh -Xmx1g \
  in1=$R1 \
  in2=$R2 \
  out1=${outdir}/$(echo $R1 | sed "s/.fastq/_trimmed.fastq/") \
  out2=${outdir}/$(echo $R2 | sed "s/.fastq/_trimmed.fastq/") \
  ref=$ref t=10 ktrim=r minlength=60 entropy=0.01 maq=10 tpe tbo
done;

## QC:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed
for ea in *fastq.gz
do
  sem -j 10 fastqc $ea --outdir fastqc
done;
cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed/fastqc
multiqc .

## STAR & FEATURECOUNTS 

## Make index:

cd /home/shared/hg_align_db/GRCm39_gencode_primary
STAR --runMode genomeGenerate \
--runThreadN 15 \
--genomeDir star_index \
--genomeFastaFiles GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile gencode.vM28.annotation.gtf \
--sjdbOverhang 75

## Align (v 2.7.8a)

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed
ref="/home/shared/hg_align_db/GRCm39_gencode_primary/star_index"
outdir="/mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned"
for R1 in XZ*R1*fastq.gz
  do R2=$(echo $R1 | sed "s/R1/R2/")
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
    --outFileNamePrefix ${outdir}/$(echo $R1 | sed "s/R1_001_trimmed.fastq.gz/trimmed_/")
done;

## Align summary:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
annot="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
for ea in *bam
do
  /home/shared/programs/qualimap_v2.2.1/./qualimap rnaseq \
  -a proportional \
  -outdir qualimap \
  -bam $ea \
  -pe \
  -p strand-specific-forward \
  --java-mem-size=20G \
  -gtf $annot
done;
cd /home/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
multiqc .

## Read quant:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/STAR_aligned
annot="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
featurecounts="/opt/subread-2.0.3-Linux-x86_64/bin/featureCounts"
for ea in *bam
  do $featurecounts \
    -T 15 \
    -a $annot \
    -o ../featureCounts/$(echo $ea | sed "s/_Aligned.sortedByCoord.out.bam//") \
    -p \
    -s 2 \
    -M -O --fraction \
    --countReadPairs \
    --extraAttributes gene_name \
    $ea
done;
cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/featureCounts
multiqc .

## SALMON

## Make decoy for indexing:

## Cat genome to end of reference transcriptome:

cd /home/shared/hg_align_db/GRCm39_genome_plus_transcripts
cat ../GRCm39_gencode_transcripts/gencode.vM28.transcripts.fa ../GRCm39_gencode_primary/GRCm39.primary_assembly.genome.fa > GRCm39.vM28.transcript_primary_assembly.genome.fa

## Extract sequence headers from genome:

mkdir salmon_index
cd salmon_index
grep "^>" ../../GRCm39_gencode_primary/GRCm39.primary_assembly.genome.fa | cut -d " " -f 1 > GRCm39.primary_assembly_decoys.txt
sed -i.bak -e 's/>//g' GRCm39.primary_assembly_decoys.txt

## Index:

conda activate salmon
salmon index -t ../GRCm39.vM28.transcript_primary_assembly.genome.fa -i GRCm39.vM28.transcript_primary_assembly.genome_Salmon_index --decoys GRCm39.primary_assembly_decoys.txt -k 29

## Align and quantify:

cd /mnt/bdata/rebecca/collabs/ghashghaei/bulk_expr/raw_reads/trimmed
ref="/home/shared/hg_align_db/GRCm39_genome_plus_transcripts/salmon_index/GRCm39.vM28.transcript_primary_assembly.genome_Salmon_index"
annot="/home/shared/hg_align_db/GRCm39_gencode_primary/gencode.vM28.annotation.gtf"
conda activate salmon
for R1 in X*R1*
  do R2=$(echo $R1 | sed "s/R1/R2/")
  salmon quant -i $ref \
  -l ISR \
  -g $annot \
  -1 $R1 -2 $R2 \
  -p 15 \
  --validateMappings \
  --gcBias \
  --seqBias \
  -o ../../Salmon_aligned/decoy/$(echo $R1 | sed "s/R1_001_trimmed.fastq.gz/trimmed/")_transcripts_quant
done;

