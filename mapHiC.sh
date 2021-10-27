cpu=$1
case_to_process=$2
case_replicate=$3
case=$4
pth=$5
name=$6
coreseq_r1=${pth}/${name}_R1_001.fastq.gz
coreseq_r2=${pth}/${name}_R2_001.fastq.gz

r1=${name}_R1_001.fastq.gz
r2=${name}_R2_001.fastq.gz

echo $1 $2 $3 $4 $5 $6 $7 $8


if [[ "$case_to_process" == "$case" ]]
then 
    codeDir='/confidential/FamilyR13/DATA/10x/SCRIPT/local_denovo_asm/scale/'
    HICREPAIE=${codeDir}HiC_repair.py
    ref_h37=/confidential/FamilyR13_data/DATA/10x/case_17-08/local_denovo_asm/dataLinks_WG/Homo_sapiens_assembly19.fa
    BWA=/project/moeinzadeh_tools/TOOLS/miniConda/bin/bwa
    PYTHON=/usr/bin/python
    SAMTOOLS=/usr/local/package/bin/samtools
    PICARD=/confidential/FamilyR13/CODE/Anaconda2/bin/picard

    echo $case_replicate

    if true 
    then
        cp $coreseq_r1 .
        cp $coreseq_r2 .


        $BWA mem -t $cpu -B 8 -M ${ref_h37} $r1 |  $SAMTOOLS view -bS -o - - | $SAMTOOLS sort  -@ $cpu -n -o ${case_replicate}'_'1.bam - # 40min
        $BWA mem -t $cpu -B 8 -M ${ref_h37} $r2 |  $SAMTOOLS view -bS -o - - | $SAMTOOLS sort  -@ $cpu -n -o ${case_replicate}'_'2.bam - 
        $PYTHON $HICREPAIE -b1 ${case_replicate}'_'1.bam -b2 ${case_replicate}'_'2.bam -o ${case_replicate}'_'hic_p.bam -m 20 #2h
        $SAMTOOLS fixmate -@ ${cpu} ${case_replicate}'_'hic_p.bam - | $SAMTOOLS sort  -@ $cpu  -o ${case_replicate}'_'hic_sorted.bam - # 20min
        # $PICARD MarkDuplicates READ_NAME_REGEX= null INPUT= ${case_replicate}'_'hic_sorted.bam OUTPUT= ${case_replicate}_hic_m_dup.bam METRICS_FILE= ${case_replicate}_hic.metrics ASSUME_SORTED= true
        /confidential/FamilyR13/CODE/Anaconda2/bin/java \
        -XX:-UseGCOverheadLimit -Xms512m -Xmx4g -jar \
        /confidential/FamilyR13/CODE/Anaconda2/share/picard-2.20.8-0/picard.jar  \
        MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.1 READ_NAME_REGEX= null \
        INPUT= ${case_replicate}'_'hic_sorted.bam OUTPUT= ${case_replicate}_hic_m_dup.bam METRICS_FILE= ${case_replicate}_hic.metrics ASSUME_SORTED= true 
        $SAMTOOLS index -@ $cpu ${case_replicate}_hic_m_dup.bam

        # cp ${case_replicate}_hic_m_dup.bam /confidential/iGenVar/Chromothripsis/Hi-C-Chromothripsis/alignment
    fi
fi

