from typing import Mapping
import os


configfile: "config.yaml"




sample = config['sample']
male_samples = config['maleSamples']
WorkingDir = config['WorkingDir']
ModifyPhaseSet = os.path.join(config['CodeDir'],config['ModifyPhaseSet'])
break_vs_haplotags = os.path.join(config['CodeDir'],config['break_vs_haplotags'])
breakpoint_to_phasedHaps = os.path.join(config['CodeDir'],config['breakpoint_to_phasedHaps'])
RNAseq_SNPcountRead = os.path.join(config['RNAseq_SNPcountRead'])

### tools
python = config['MyPython']
python_3_7 = config ['MyPython_3_7']
bwa = config['MyBwa']
samtools = config['MySamtools']
picard = config['MyPicard']
java = config['MyJava']
hapcut2 = config['MyHAPCUT2']
whatshap = config['Mywhatshap']
extractHairs = config['MyExtractHAIRS']
freebayes = config['MyFreebayes']
bcftools = config['MyBcftools']
nanocaller = config ['myNanoCaller']
correct_overlapping_smps = config ['correct_overlapping_smps']


### data
#reference
#l_refs = [os.path.join(WorkingDir,'InputDir','Refs',ref) for ref in config['refs']]
# l_refs = config['refs']


l_pairs = ['R1','R2']
l_gene_exon = ['exon','gene']
l_samples = sample
l_refs =['hg19']
l_smp_variant_caller = ['variantCallingFreebayes']

def get_illumina_reads(wildcards):
    return [config['ill_fastq'][wildcards.sample]['R1'], config['ill_fastq'][wildcards.sample]['R2']]




chromosomes = ['chr'+str(i) for i in range(1,23)] + ['chrX']
map_sample2chromosomes = {}
for s in sample:
    map_sample2chromosomes [s] = chromosomes + ['chrY'] if s in male_samples else chromosomes

# print(map_sample2chromosomes)


sampleHiC,repHiC = glob_wildcards(os.path.join(WorkingDir, 'InputDir','Samples','{sampleHiC}','HiC','{repHiC}.R1.fastq.gz'))

map_sample2rep = {s:[] for s in sampleHiC}
for s,r in zip (sampleHiC,repHiC):
    map_sample2rep [s] += [r]





# Generating 'Wild.bam','Chrm.bam','Bridge.bam','Break.bam' is is deactivated 
# break_vs_haplotags_output = ['PS_w.txt','PS_wild.txt','PS_chromo.txt','PS_Conflict_betweenWild_Chrm.txt','conflict_of_halotags.png','PS_final.txt','Wild.bam','Chrm.bam','Bridge.bam','Break.bam']
break_vs_haplotags_output = ['PS_w.txt','PS_wild.txt','PS_chromo.txt','PS_Conflict_betweenWild_Chrm.txt','conflict_of_halotags.png','PS_final.txt']


wildcard_constraints:
    pair = 'R1|R2',
    chr = 'chr[0-9]+|chrX|chrY'


rule all:
    input:
        expand(os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','PhaseRNAReads','RNAcount.{gene_or_exon}.csv'), ref=l_refs, sample = l_samples, gene_or_exon=l_gene_exon, smpVarCaller = l_smp_variant_caller ),

rule bwa_map_hic:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        read=os.path.join(WorkingDir, 'InputDir','Samples','{sample}','HiC','{rep}.{pair}.fastq.gz')
    output: 
        temp(os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.{pair}.bam'))
    log: "log/bwa_map_hic_{ref}_{sample}_{rep}.{pair}.log"
    threads: 16
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell: "(time {bwa} mem -t {threads} -B 8 {input.ref} {input.read}|  {samtools} view -bS -o - - | {samtools} sort  -@ {threads} -n -o {output} -) > {log} 2>&1"

rule hicRepair:
    input:
        expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}' ,'{{sample}}','HiC','{{rep}}.{pair}.bam') ,pair = l_pairs)
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_p.bam')
    log: "log/hic_repair.{ref}_{sample}_{rep}.log"
    params:
        mappingQ = 20 # mappingQ
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {python} Utils/HiC_repair.py -b1 {input[0]} -b2 {input[1]} -o {output} -m {params.mappingQ}) > {log} 2>&1
        '''

rule hicRepairFixmate:
    input: 
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_p.bam')
    output: 
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_sorted.bam')
    log: "log/hic_repairFixmate.{ref}_{sample}_{rep}.log"
    threads: 16
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '(time {samtools} fixmate -@ {threads} {input} - | {samtools} sort  -@ threads  -o {output} -) > {log} 2>&1'

rule hicMarkDup:
    input: os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_sorted.bam')
    output: 
        bam_mdup=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_m_dup.bam'),
        index_bam_mdup=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic_m_dup.bam.bai'),        
        file_metrics=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','{rep}.hic.metrics'),
    log:
        "log/hicMarkDup.{ref}_{sample}_{rep}.log"
    threads:1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {java} \
        -XX:-UseGCOverheadLimit -Xms512m -Xmx4g -jar \
        {picard}  \
        MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.1 READ_NAME_REGEX= null \
        INPUT={input} OUTPUT={output.bam_mdup} METRICS_FILE={output.file_metrics} ASSUME_SORTED= true &&
        {samtools} index -@ {threads} {output.bam_mdup}) > {log} 2>&1
        '''

rule hic_merge_replicate:
    input:
        bam=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'HiC','{rep}.hic_m_dup.bam'), rep=map_sample2rep[wc.sample]),
        index_bam=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'HiC','{rep}.hic_m_dup.bam.bai'), rep=map_sample2rep[wc.sample])
        
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','hic_m_dup.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','hic_m_dup.bam.bai'),
    threads: 8
    log:
        "log/hic_rerge_replicate.{ref}_{sample}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time samtools merge -@ {threads} --output-fmt SAM - {input.bam} | samtools sort -@ {threads} -o {output.bam} &&
        samtools index -@ {threads} {output.bam} ) > {log} 2>&1
        '''

#### -----------------------------------------------------------------------------------------
#### -----------------------------------------HAPLOTYPING-------------------------------------
#### -----------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------



rule pb_index:
    input:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.bam'),
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.bam.bai'),  
    threads:1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    log:"log/pb_index.{ref}_{sample}.log"
    shell:
        '''
        (time {samtools} index -@ {threads} {input}) > {log} 2>&1
        '''










rule ill_map_fixmate:
    input:
        fq = get_illumina_reads,
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa')
    output:
        bam=temp(os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.fixmate.bam'))
    threads: 40
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    log:"log/ill_map_fixmate.{ref}_{sample}.log"
    shell:
        '''
        (time {bwa} mem {input.ref} -t {threads} -v 3 -M -R \'@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA\' {input.fq[0]} {input.fq[1]} | samtools view -b - | samtools fixmate -m - {output.bam} ) > {log} 2>&1
        '''


rule ill_map_markdup:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        bam_fixmate=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.fixmate.bam')
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.markdup.bam')
    threads: 16
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    log:"log/rule_ill_map_markdup.{ref}_{sample}.log"
    shell:
        '''
        (time samtools sort -@ {threads} {input.bam_fixmate} | samtools markdup -@ {threads} - {output.bam}) > {log} 2>&1
        '''

rule ill_index:
    input:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.markdup.bam'),
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.markdup.bam.bai'),  
    threads:1
    log:"log/ill_index.{ref}_{sample}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {samtools} index -@ {threads} {input}) > {log} 2>&1
        '''

rule freebayes:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        illBam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','ILL','ill.markdup.bam'),
    output:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingFreebayes','variants','freebayes.{chr}.vcf'),
        vcf_gz=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingFreebayes','variants','freebayes.{chr}.vcf.gz'),
        vcf_gz_idx=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingFreebayes','variants','freebayes.{chr}.vcf.gz.tbi'),
        vcf_gz_filtered=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingFreebayes','variants','{chr}_filter.vcf'),
    threads: 1
    log:"log/freebayes.{ref}_{sample}_{chr}_variantCallingFreebayes.log"
    params:
        variantQ=30
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {freebayes} -f {input.ref} -r {wildcards.chr} -u  {input.illBam} > {output.vcf} &&
        bgzip -c {output.vcf} > {output.vcf_gz} &&
        tabix -p vcf {output.vcf_gz} &&
        echo "Filter low quality variants (Q = {params.variantQ})" &&
        {bcftools} view -i '%QUAL>='{params.variantQ} {output.vcf_gz} {wildcards.chr} | grep -v "\./." | grep -v "1/1" > {output.vcf_gz_filtered}) > {log} 2>&1
        '''


rule nanocaller:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        pacBam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.bam'),
    output:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}','variant_calls.final.vcf.gz'),
        vcftbi=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}','variant_calls.final.vcf.gz.tbi'),
        phasedBamBySNPs=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}','variant_calls.snps.phased.bam'),
    threads: 20
    log:"log/nanocaller.{ref}_{sample}_{chr}_variantCallingNanoCaller.log"
    conda:
        "/confidential/FamilyR13/DATA/10x/SCRIPT/local_denovo_asm/scale/condaEnvs/nanoCaller.yml"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        outdir=$(dirname {output[0]}) &&
        /home/Cmoeinza/miniconda3/envs/NanoCaller/bin/python {nanocaller} -bam {input.pacBam} -p clr -o $outdir -chrom {wildcards.chr}  -ref {input.ref} -cpu {threads} -keep_bam > {log} 2>&1 
        '''

rule nanocaller_filter:
    input:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}','variant_calls.final.vcf.gz'),
        vcftbi=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}','variant_calls.final.vcf.gz.tbi'),
    output:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','variantCallingNanoCaller','variants','{chr}_filter.vcf'),
    threads: 1
    log:"log/nanocaller_filter.{ref}_{sample}_{chr}_nanocaller_filter.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    params:
        variantQ=30
    shell:
        '''
        ( time echo "Filter low quality variants (Q = {params.variantQ})" &&
        {bcftools} view -i '%QUAL>='{params.variantQ} {input.vcf} {wildcards.chr} | grep -v "\./." | grep -v "1/1" | grep -v "1|1" | grep -v "0/0" | grep -v "0|0" | {python} {correct_overlapping_smps} > {output.vcf} ) > {log} 2>&1
        '''

        # {bcftools} view -i '%QUAL>='{params.variantQ} {input.vcf} {wildcards.chr} | grep -v "\./." | grep -v "1/1" > {output.vcf}) > {log} 2>&1



        # {bcftools} view -i '%QUAL>='{params.variantQ} {output.vcf_gz} {wildcards.chr} | grep -v "\./." | grep -v "1/1" | grep -v "1|1" | grep -v "0/0" | grep -v "0|0" > {output.vcf_gz_filtered}) > {log} 2>&1

# outdir=$(dirname {output[0]}) && /home/Cmoeinza/miniconda3/envs/NanoCaller/bin/python {nanocaller} -bam {input.pacBam} -p clr -o ${outdir} -chrom {wildcards.chr}  -ref {input.ref} -cpu {threads} ) > {log} 2>&1 )

rule split_HiC_by_chromosome:
    input:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','hic_m_dup.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','hic_m_dup.bam.bai'),
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','HiC.{chr}.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','HiC.{chr}.bam.bai'),
    threads: 1
    log:"log/split_HiC_by_chromosome.{ref}_{sample}_{chr}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time echo '  hic' &&
        {samtools} view -@ {threads} -b {input.bam} {wildcards.chr} > {output.bam} &&
        {samtools} index -@ {threads} {output.bam} ) > {log} 2>&1
        '''
        

rule extractHairs_by_chrm_HiC:
    input:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','HiC.{chr}.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','HiC','HiC.{chr}.bam.bai'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','variants','{chr}_filter.vcf'),
    output:
        fragHiC=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','HiC.{chr}.frag'),
    threads: 1
    log:"log/extractHairs_by_chrm_HiC.{ref}_{sample}_{chr}_{smpVarCaller}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {extractHairs} --HiC 1 --bam {input.bam} --vcf {input.vcf}  > {output.fragHiC}) > {log} 2>&1
        '''


rule split_pb_by_chr:
    input:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.bam.bai'),
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam.bai'),
    threads: 1
    log:"log/split_pb_by_chr.{ref}_{sample}_{chr}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time echo 'pb' && {samtools} view -@ {threads} -b {input.bam} {wildcards.chr} > {output.bam} && {samtools} index -@  {threads} {output.bam} ) > {log} 2>&1
        '''


rule svim:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam.bai'),
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','svim','{chr}','variants.vcf'),
    params:
        working_dir=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','svim','{chr}')
    log:"log/svim_{ref}_{sample}_{chr}.log"
    shell:
        '''
        (time svim alignment --sample {wildcards.sample} --read_names {params.working_dir} {input.bam} {input.ref}) > {log} 2>&1
        '''


rule extractHairs_by_chrm_Pb:
    input:
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam.bai'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','variants','{chr}_filter.vcf'),
    output:
        fragPb=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','pb.{chr}.frag'),
    threads: 1
    log:"log/extractHairs_by_chrm_Pb.{ref}_{sample}_{chr}_{smpVarCaller}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {extractHairs} --pacbio 1 --new_format 1 --bam  {input.bam} --vcf {input.vcf}  --ref {input.ref} >  {output.fragPb} ) > {log} 2>&1 
        '''



rule hapcut2_from_fragments:
    input:
        fragHiC=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','HiC.{chr}.frag'),
        fragPb=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','pb.{chr}.frag'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','variants','{chr}_filter.vcf'),
    output:
        Concatenated_fragment=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','pb.{chr}.fragHicPb'),
        hapfile=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.hap'),
        hapModelFile=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.model_file'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.vcf'),
        hapStatgtf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.gtf'),
        hapStattsv=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.tsv'),
    log:
        "log/hapcut2_from_fragments.{ref}_{sample}_{chr}_{smpVarCaller}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time echo 'hapcut2_from_fragments'
        cat {input.fragHiC} {input.fragPb} > {output.Concatenated_fragment} &&
        {hapcut2} --fragments {output.Concatenated_fragment} --vcf {input.vcf} --output {output.hapfile} --hic 1 --htrans_data_outfile {output.hapModelFile} &&
        {whatshap} hapcut2vcf -o {output.vcf} {input.vcf} {output.hapfile} &&
        {whatshap} stats --gtf {output.hapStatgtf}  --tsv {output.hapStattsv} {output.vcf} ) > {log} 2>&1
        '''


rule modify_ps:   #fix PS to avoid identical PSs in different chromosomes
    input:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.vcf'),
    output:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.modPS_phased.vcf.gz'),
        vcftbi=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.modPS_phased.vcf.gz.tbi'),
    log:"log/modidy_PS.{ref}_{sample}_{chr}_{smpVarCaller}.log"
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {python} {ModifyPhaseSet} {input.vcf} {wildcards.chr} | bgzip -c > {output.vcf} &&
        tabix -p vcf {output.vcf} ) > {log} 2>&1
        '''



rule keep_onlyPrimaryreads:
    input:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam'),        
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.bam.bai'),
    output:
        bamTMP=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.onlyPrimaryreads.bam'),
        bamTMPbai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.onlyPrimaryreads.bam.bai'),
    log:
        "log/keep_onlyPrimaryreads.{ref}_{sample}_{chr}.log"
    threads: 1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {samtools} view -h -F 2048 {input.bam} -o {output.bamTMP}  &&
        {samtools} index -@ {threads} {output.bamTMP}) > {log} 2>&1
        '''


rule haplotag_with_new_phase_set:
    input:
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.modPS_phased.vcf.gz'),
        vcftbi=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping','{chr}.modPS_phased.vcf.gz.tbi'),
        ref=os.path.join(WorkingDir,'InputDir','Refs', '{ref}.fa'),
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.onlyPrimaryreads.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','PB','pb.{chr}.onlyPrimaryreads.bam.bai'),
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','PB_tagged','pb.{chr}.tag.bam'),
        bambai=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','PB_tagged','pb.{chr}.tag.bam.bai'),   

    log:
        "log/haplotag_with_new_phase_set.{ref}_{sample}_{chr}_{smpVarCaller}.log"
    threads: 1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {whatshap} haplotag  -r {input.ref} --ignore-read-groups -o {output.bam}  {input.vcf} {input.bam} &&
        {samtools} index -@ {threads} {output.bam} ) > {log} 2>&1
        '''

rule merge_bam_and_vcf:
    input:
        bams=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'{{smpVarCaller}}','PB_tagged','pb.{chr}.tag.bam'),chr=map_sample2chromosomes[wc.sample]), 
        bamsbai=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'{{smpVarCaller}}','PB_tagged','pb.{chr}.tag.bam.bai'),chr=map_sample2chromosomes[wc.sample]), 
        vcfs=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'{{smpVarCaller}}','haplotyping','{chr}.modPS_phased.vcf.gz'),chr=map_sample2chromosomes[wc.sample]), 
        vcfs_index=lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'{{smpVarCaller}}','haplotyping','{chr}.modPS_phased.vcf.gz.tbi'),chr=map_sample2chromosomes[wc.sample]), 
    output:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','pb.tag.bam'),
        bam_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','pb.tag.bam.bai'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','modPS_phased.vcf.gz'),
        vcf_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','modPS_phased.vcf.gz.tbi')
    log:
        "log/merge_bam_and_vcf.{ref}_{sample}_{smpVarCaller}.log"
    threads: 8
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {samtools} merge -f -h {input.bams[0]} -@ {threads} -O BAM {output.bam} {input.bams} &&
        {samtools} index -@ {threads} {output.bam} &&
        {bcftools} merge --force-samples --threads {threads} -o {output.vcf} -O z {input.vcfs} &&
        tabix -p vcf {output.vcf} ) > {log} 2>&1
        '''


rule break_vs_haplotags:
    input:
        bam=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','pb.tag.bam'),
        bam_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','pb.tag.bam.bai'),
        break_file=os.path.join(WorkingDir,'InputDir','BreakFiles','{ref}', '{sample}.breakpoints'),
    output:
        expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}' ,'{{sample}}','{{smpVarCaller}}','haplotyping_with_breakpoints','{out}'),out=break_vs_haplotags_output),
    log:
        err="log/break_vs_haplotags.{ref}_{sample}_{smpVarCaller}.err",
        out="log/break_vs_haplotags.{ref}_{sample}_{smpVarCaller}.out",
    threads: 8
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time outdir=$(dirname {output[0]}) &&
        echo $outdir &&
        {python} {break_vs_haplotags} $outdir {input.bam} {input.break_file} {threads}) 2> {log.err} 1> {log.out}
        '''

rule breakpoint_to_phasedHaps:
    input:
        ps_file=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','PS_w.txt'),
        vcf=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','modPS_phased.vcf.gz'),
        vcf_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','modPS_phased.vcf.gz.tbi')        
    output:
        vcf_withBreak=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf'),
        vcf_onlyBreak=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.onlybreak.vcf'),
    log:
        err="log/breakpoint_to_phasedHaps.{ref}_{sample}_{smpVarCaller}.err",
        out="log/breakpoint_to_phasedHaps.{ref}_{sample}_{smpVarCaller}.out"
    threads: 1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time {python} {breakpoint_to_phasedHaps} {input.ps_file} {input.vcf} {output.vcf_withBreak} {output.vcf_onlyBreak}) 2> {log.err} 1> {log.out}
        '''
rule zip_and_index_break_vcfs:
    input:
        vcf_withBreak=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf'),
        vcf_onlyBreak=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.onlybreak.vcf'),
    output:
        vcf_withBreak_gz=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf.gz'),
        vcf_withBreak_gz_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf.gz.tbi'),
        vcf_onlyBreak_gz=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.onlybreak.vcf.gz'),
        vcf_onlyBreak_gz_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.onlybreak.vcf.gz.tbi')
    log:
        err="log/zip_and_index_break_vcfs.{ref}_{sample}_{smpVarCaller}.err",
        out="log/zip_and_index_break_vcfs.{ref}_{sample}_{smpVarCaller}.out"
    threads: 1
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        bgzip {input.vcf_withBreak} && tabix -p vcf {output.vcf_withBreak_gz} &&
        bgzip {input.vcf_onlyBreak} && tabix -p vcf {output.vcf_onlyBreak_gz}
        '''

rule RNAseq_SNPcountRead:
    input:
        break_file=os.path.join(WorkingDir,'InputDir','BreakFiles','{ref}', '{sample}.breakpoints'),
        vcf_withBreak=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf.gz'),
        vcf_withBreak_index=os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','haplotyping_with_breakpoints','phased.withbreak.vcf.gz.tbi'),
        params_RNAdata=os.path.join(WorkingDir,'InputDir','RNA_params','{ref}', 'params_{sample}.py'),
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','PhaseRNAReads','RNAcount_{chr}.{gene_or_exon}.csv')
    log:
        err='log/rnaLog/RNAseq_SNPcountRead_{ref}_{sample}_{chr}_{smpVarCaller}_{gene_or_exon}.err',
        out='log/rnaLog/RNAseq_SNPcountRead_{ref}_{sample}_{chr}_{smpVarCaller}_{gene_or_exon}.out'
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'
    shell:
        '''
        (time datadir=$(dirname {input.params_RNAdata}) &&
        outdir=$(dirname {output}) &&
        {python_3_7} {RNAseq_SNPcountRead} $outdir {input.break_file} {input.vcf_withBreak} {wildcards.sample} {wildcards.chr} {wildcards.gene_or_exon} $datadir) 2> {log.err} 1> {log.out}
        '''

rule merge_rnaseq_count:
    input:
        lambda wc: expand(os.path.join(WorkingDir,'OutputDir' ,'{{ref}}', wc.sample,'{{smpVarCaller}}','PhaseRNAReads','RNAcount_{chr}.{{gene_or_exon}}.csv'), chr=map_sample2chromosomes[wc.sample])
    output:
        os.path.join(WorkingDir,'OutputDir' ,'{ref}' ,'{sample}','{smpVarCaller}','PhaseRNAReads','RNAcount.{gene_or_exon}.csv')
    log:
        err='log/rnaLog/RNAseq_SNPcountRead_{ref}_{sample}_{smpVarCaller}_{gene_or_exon}.err',
        out='log/rnaLog/RNAseq_SNPcountRead_{ref}_{sample}_{smpVarCaller}_{gene_or_exon}.out',
    resources:
        tmpdir='/scratch/local2/Cmoeinzadeh'

    shell:
        '''
        (time cat  <(head -1 {input[0]})  <(cat {input} | grep -v Dist2BP) > {output}) 2> {log.err} 1> {log.out}
        ##### todo: complete cleaning and merging
        '''
