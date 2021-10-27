
cr=$1
case=$2
dataFolder=$3
outFolder=$4
variantQ='30'
cpu=5
# r=':27298614-30852969'
r=''

freebayes=/confidential/FamilyR13/CODE/miniconda3/bin/freebayes

ref_h37=/confidential/FamilyR13_data/DATA/10x/case_17-08/local_denovo_asm/Homo_sapiens_assembly19.fa
HAPCUT2='/confidential/FamilyR13/DATA/HAPS/tools/HapCUT2/build/HAPCUT2'
whatshap='/project/moeinzadeh_tools/conda/python3_7conda/bin/whatshap'
extractHAIRS='/confidential/FamilyR13/DATA/HAPS/tools/HapCUT2/build/extractHAIRS'
ref_h37=/confidential/FamilyR13_data/DATA/10x/case_17-08/local_denovo_asm/dataLinks_WG/Homo_sapiens_assembly19.fa

codeDir='/confidential/FamilyR13/DATA/10x/SCRIPT/local_denovo_asm/scale/'


caseDir=$outFolder
mkdir -p $caseDir
cd $caseDir


pacbio_bam=${dataFolder}"pb."${case}".bam"








Illumina_bam=${dataFolder}"bwa."${case}".bam"
hic_bam=${dataFolder}"hic."${case}".bam"
breakfiles=${dataFolder}"bp."${case}".csv"
# vcf=${dataFolder}"gatk_hc."${case}".vcf.gz"
# vcf=${dataFolder}"gatk_hc."${case}".vcf"

vcf='freebayes.'${case}'_chr'${cr}'.GRCh37.vcf'
vcfgz='freebayes.'${case}'_chr'${cr}'.GRCh37.vcf.gz'
vcf_h37='fb.'${case}'_chr'${cr}'.GRCh37.vcf'
vcf_hg19='fb.'${case}'_chr'${cr}'.hg19.vcf'

echo 'Chromosome:' $cr
# variants calling
echo 'variants calling and filtering quality'
#!!!!!!!!!!!
if true
then
    echo "here"
    $freebayes -f $ref_h37 -r $cr$r -u  $Illumina_bam > $vcf
    bgzip -c $vcf > $vcfgz
    tabix -p vcf $vcfgz
    echo "Filter low quality variants (Q = "$variantQ")"
    bcftools view -i '%QUAL>='$variantQ $vcfgz $cr | grep -v "\./." | grep -v "1/1" > $vcf_h37
fi


echo 'haplotyping'
if false
then


    # echo 'Convert vcf to hg19'
    # cat <(bcftools view -h $vcf_h37) <(bcftools view -H $vcf_h37 | awk '$0="chr"$0') >  $vcf_hg19

    echo '  hic'
    samtools view -@ $cpu -b $hic_bam $cr${r} > 'hic.chr'$cr'.bam'
    samtools index -@ $cpu  'hic.chr'$cr'.bam'
    $extractHAIRS --HiC 1 --bam  'hic.chr'$cr'.bam' --vcf $vcf_h37  > chr$cr.fragHiC
fi
if true
then
    samtools view -@ $cpu -b $pacbio_bam $cr${r} > 'pb.'$cr'.bam'
    samtools index -@ $cpu 'pb.'$cr'.bam'
    $extractHAIRS --pacbio 1 --new_format 1 --bam  'pb.'$cr'.bam' --vcf $vcf_h37  --ref $ref_h37 > chr$cr.fragPB


    echo '  cat fragments'
    cat chr${cr}.fragHiC chr${cr}.fragPB > chr${cr}.hicPB

    echo 'hapcut'
    $HAPCUT2 --fragments chr${cr}.hicPB --vcf $vcf_h37 --output chr${cr}.hap --hic 1 --htrans_data_outfile chr${cr}.model_file

    # echo $whatshap hapcut2vcf -o chr${cr}_phased.vcf  $vcf_h37 chr${cr}.hap
    echo '  hap to vcf'
    $whatshap hapcut2vcf -o chr${cr}_phased.vcf  $vcf_h37 chr${cr}.hap

    echo '  make stats'
    # echo  $whatshap stats --gtf chr${cr}.gtf  --tsv chr${cr}.tsv chr${cr}_phased.vcf
    $whatshap stats --gtf chr${cr}.gtf  --tsv chr${cr}.tsv chr${cr}_phased.vcf
fi

### fix PS to avoid identical PSs in different chromosomes
echo 'fix PS to avoid identical PSs in different chromosomes'
if true
then
    /usr/bin/python ${codeDir}modify_PS.py chr${cr}_phased.vcf ${cr} | bgzip -c > chr${cr}_modPS_phased.vcf.gz
    tabix -p vcf chr${cr}_modPS_phased.vcf.gz
fi

### tag pacbio reads
echo 'tag pacbio reads'
if true
then
    # $whatshap haplotag   -r $ref_h37 --ignore-read-groups  >  
    /project/haplotyping/nico/tools/miniconda2/envs/whatshap/bin/whatshap haplotag  -r $ref_h37 --ignore-read-groups -o 'pb.'$cr'.tag.bam'  chr${cr}_modPS_phased.vcf.gz 'pb.'$cr'.bam'
    samtools index 'pb.'$cr'.tag.bam' 
fi



# if true
# then
#     /project/moeinzadeh_tools/conda/python3_7conda/bin/python3.7 ${codeDir}RNAseq_SNPcountRead.py $caseDir $breakfiles $caseOutDir'phased.withbreak.vcf.gz' $case $cr 'exon' $case 2> ${logs}${cr}'_exon.err' 1> ${logs}${cr}'_exon.log' &

#     /project/moeinzadeh_tools/conda/python3_7conda/bin/python3.7 ${codeDir}RNAseq_SNPcountRead.py $caseOutDir $breakfiles $caseOutDir'phased.withbreak.vcf.gz' $case $cr 'gene' $case 2> ${logs}${cr}'_gene.err' 1> ${logs}${cr}'_gene.log' &

#     echo 'sleeping second: '$sleeptime
#     sleep $sleeptime
# fi


touch ${cr}.hapsFinished

echo 'Done'
