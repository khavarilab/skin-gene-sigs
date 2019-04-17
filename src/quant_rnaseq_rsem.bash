#!/bin/bash

#module add rsem

# Given a inDir and outDir
inDir=$1
prefix=$2
RSEMrefDir=$3
dataType=$4
nThreadsRSEM=$5
outDir=$6

# # input: gzipped fastq file read1 [read2 for paired-end] 
# #        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
# read1=$1 #gzipped fastq file for read1
# read2=$2 #gzipped fastq file for read1, use "" if single-end
# STARgenomeDir=$3 
# RSEMrefDir=$4
# dataType=$5 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
# nThreadsSTAR=$6 # number of threads for STAR
# nThreadsRSEM=$7 # number of threads for RSEM

RSEM=rsem-calculate-expression

######### RSEM

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
# NOTE: IF USING THIS IN THE FUTURE, USE SAMBAMBA SORT AND USE TMP SSD SPACE
#cd $inDir

#echo "Preparing for RSEM..."
mkdir $outDir

#trBAMsortRAM=12G

##cp ${inDir}/${prefix}_Aligned.toTranscriptome.out.bam ${outDir}/Tr.bam 

#rsync -avzL --progress ${inDir}/${prefix}_Aligned.toTranscriptome.out.bam ${outDir}/
#mv ${outDir}/${prefix}_Aligned.toTranscriptome.out.bam ${outDir}/Tr.bam 

#case "$dataType" in
#str_SE|unstr_SE)
      # single-end data
#      cat <( samtools view -H ${outDir}/Tr.bam ) <( samtools view -@ $nThreadsRSEM ${outDir}/Tr.bam | sort --parallel=$nThreadsRSEM -S $trBAMsortRAM -T ${outDir} ) | samtools view -@ $nThreadsRSEM -bS - > ${outDir}/${prefix}_Aligned.toTranscriptome.out.bam
#      ;;
#str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
#      cat <( samtools view -H ${outDir}/Tr.bam ) <( samtools view -@ $nThreadsRSEM ${outDir}/Tr.bam | awk '{printf $0 " "; getline; print}' | sort --parallel=$nThreadsRSEM -S $trBAMsortRAM -T ${outDir} | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > ${outDir}/${prefix}_Aligned.toTranscriptome.out.bam
#      ;;
#esac

#'rm' ${outDir}/Tr.bam
#echo "Done with RSEM prep!"


# RSEM parameters: common
#RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"
RSEMparCommon="--bam --estimate-rspd --no-bam-output --seed 12345" # Don't calculate confidence intervals! Huge computational cost

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac

###### RSEM command
echo "Running RSEM quantification..."
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType ${inDir}/${prefix}_Aligned.toTranscriptome.out.bam $RSEMrefDir ${outDir}/Quant >& ${outDir}/${prefix}_Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType ${inDir}/${prefix}_Aligned.toTranscriptome.out.bam $RSEMrefDir ${outDir}/Quant >& ${outDir}/${prefix}_Log.rsem
echo "Done running RSEM!"
#'rm' ${outDir}/${prefix}_Aligned.toTranscriptome.out.bam

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model ${outDir}/Quant ${outDir}/${prefix}_Quant.pdf
rsem-plot-model ${outDir}/Quant ${outDir}/${prefix}_Quant.pdf
