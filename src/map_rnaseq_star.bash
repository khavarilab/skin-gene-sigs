#!/bin/bash

mkdir $2;
STAR --genomeDir $1 \
     --outFileNamePrefix $2/$3 \
     --readFilesIn $5 $6 \
     --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --sjdbScore 1 \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outWigStrand Stranded \
     --twopassMode Basic \
     --twopass1readsN -1 \
     --limitBAMsortRAM 250000000000 \
     --runThreadN $4;
