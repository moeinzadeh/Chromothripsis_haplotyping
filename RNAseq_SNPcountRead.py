import numpy as np
import vcf
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import math
import pysam
import pyranges as pr
import codecs
import os
import sys
# import params_molgen_ossi
# params=params_molgen_ossi

working_dir = sys.argv[1] + '/'

datadir = sys.argv[7]
case = sys.argv[4]
sys.path.append(datadir)
params = __import__('params_'+sys.argv[4])

# import params_molgen
# params=sys.argv[2]

file_breakpoint = sys.argv[2]
vcf_file = sys.argv[3]
individual_name = params.individual_name

case = sys.argv[4]
interval = sys.argv[5]

bo_phaseReadInBam =params.bo_phaseReadInBam 


bams = [pysam.AlignmentFile(params.path + f) for f in params.RNAseqbams]



if sys.argv[6] == 'exon':
    gene_annotation = params.gene_annotation_exon
    out_df = working_dir+'RNAcount_'+interval+'.exon.csv'
    bo_exon = True

    bams_phased_mother, bams_phased_father = [],[]
    if bo_phaseReadInBam:
        ln = len(params.RNAseqbams)
        bams_phased_mother = [pysam.AlignmentFile(working_dir + f[:-29] + '_' + interval + '_exon_Wild_phased.bam','wb',template=bams[0]) for f,inx in zip(params.RNAseqbams,range(ln))]
        bams_phased_father = [pysam.AlignmentFile(working_dir + f[:-29] + '_' + interval + '_exon_Chrm_phased.bam','wb',template=bams[0]) for f,inx in zip(params.RNAseqbams,range(ln))]
        f_bams_phased_mother = [working_dir + f[:-29] + '_' + interval + '_exon_Wild_phased.bam' for f in params.RNAseqbams]
        f_bams_phased_father = [working_dir + f[:-29] + '_' + interval + '_exon_Chrm_phased.bam' for f in params.RNAseqbams]

else:
    gene_annotation = params.gene_annotation_gene
    out_df = working_dir+'RNAcount_'+interval+'.gene.csv'
    bo_exon = False

    bams_phased_mother, bams_phased_father = [],[]
    if bo_phaseReadInBam:
        ln = len(params.RNAseqbams)
        bams_phased_mother = [pysam.AlignmentFile(working_dir + f[:-29] + '_' + interval + '_gene_Wild_phased.bam','wb',template=bams[0]) for f,inx in zip(params.RNAseqbams,range(ln))]
        bams_phased_father = [pysam.AlignmentFile(working_dir + f[:-29] + '_' + interval + '_gene_Chrm_phased.bam','wb',template=bams[0]) for f,inx in zip(params.RNAseqbams,range(ln))]
        f_bams_phased_mother = [working_dir + f[:-29] + '_' + interval + '_gene_Wild_phased.bam' for f in params.RNAseqbams]
        f_bams_phased_father = [working_dir + f[:-29] + '_' + interval + '_gene_Chrm_phased.bam' for f in params.RNAseqbams]


# interval_len_from_breakpoints = params.interval_len_from_breakpoints

print(sys.argv)
# exit()
print(params.RNAseqbams)



for i in range(2):
    bams_phased = []
# vcf_reader = vcf.Reader(open(params.vcf_file))
vcf_reader = vcf.Reader(filename = vcf_file)


# interval_string = params.interval_string
open_left_region = params.open_left_region
open_right_region = params.open_right_region

gene_annotation_all = params.gene_annotation_all




def count_reads(chrm, view_start, view_end, samfiles, vcf_reader,bo_out2phasedBams,fo_mother_bam,fo_father_bam):






    mother = []
    father = []
    coordinates = []
    l_mother_is_ref = []
    replicate = []
    PSs = []

    # print (chrm, view_start, view_end, samfiles, vcf_reader,bo_out2phasedBams,fo_mother_bam,fo_father_bam)

    for Record in vcf_reader.fetch(chrm, view_start, view_end):  # doctest: +SKIP
        genotype = Record.genotype(individual_name)['GT']
        # not phased or homozygous
        if not ("|" == genotype[1]) or genotype[0]==genotype[2]:
            continue

        # multi nucleotide polymorphism
        one_nc_snp = True
        for var in Record.ALT:
            if len(var) > 1:
                one_nc_snp = False
        if len(Record.REF) > 1 or not one_nc_snp:
            continue


        
        snp = Record.POS

        PS = Record.genotype(individual_name)['PS']

        vars = [Record.REF] + [str(v) for v in Record.ALT]
        
        #print (mp_cnt,vars)
        mothers_allele = vars[int(genotype[0])]
        fathers_allele = vars[int(genotype[2])]

        if int(genotype[0]) == 0:
            mother_is_ref = True
        else:
            mother_is_ref = False

        # print (snp)
        bo_first = True
        cnt = 0
        mp_cnt = {} #>>>>
        for rep, samfile in zip(range(len(samfiles)),samfiles):
            # print (snp)
            mp_cnt = {} #<<<
            for pileupcolumn in samfile.pileup(chrm, snp - 1, snp):
                # print(pileupcolumn.pos)
                cnt += 1
                # print(cnt)

                if bo_first:
                    bo_first = False
                    #print(Record.CHROM, Record.POS, Record.ID, Record.REF, Record.ALT, Record.genotype(individual_name)['GT'])
                if pileupcolumn.pos == snp - 1:
                    for pileupread in pileupcolumn.pileups:
                        # print (pileupread.alignment)

                        if not pileupread.is_del and pileupread.alignment.mapping_quality >= 255:
                            if bo_out2phasedBams and PS == 1:
                                # print(pileupread.alignment)
                                if pileupread.alignment.query_sequence[pileupread.query_position] == mothers_allele:
                                    fo_mother_bam[rep].write(pileupread.alignment)
                                if pileupread.alignment.query_sequence[pileupread.query_position] == fathers_allele:
                                    fo_father_bam[rep].write(pileupread.alignment)


                            mp_cnt[pileupread.alignment.query_sequence[pileupread.query_position]] = mp_cnt.get(
                                pileupread.alignment.query_sequence[pileupread.query_position], 0) + 1
            # exit()
            # mapping quality is lower 
            # print (mp_cnt)
            # exit()

            if len(mp_cnt) == 0:
                mother += [0]
                father += [0]
                coordinates += [Record.POS]
                replicate += [rep]
                PSs += [PS]
                l_mother_is_ref += [mother_is_ref]
                continue
            
            coordinates += [Record.POS]

            replicate += [rep]
            PSs += [PS]
            l_mother_is_ref += [mother_is_ref]
        
            if mothers_allele in mp_cnt:
                mother += [mp_cnt[mothers_allele]]
            else:
                mother += [0]

            if fathers_allele in mp_cnt:
                father += [mp_cnt[fathers_allele]]
            else:
                father += [0]
            #print (mother[-1],father[-1])
            
            # if int(genotype[0]) == 0:
            #     mother_is_ref += [True]
            # else:
            #     mother_is_ref += [False]


    m = np.array(mother, dtype=np.float)
    f = np.array(father, dtype=np.float)
    c = np.array(coordinates, dtype=np.int)
    ref = np.array(l_mother_is_ref, dtype=np.bool)
    rep = np.array(replicate,dtype=np.int)
    pss = np.array(PSs,dtype=np.int64)
    # if len(pss) != 0:
    #     print (pss[0])
    
    # print (m, f, c, ref,rep)
    # print(len(pss),len(rep))
    return m, f, c, ref,rep,pss




gr = pr.read_gtf(gene_annotation)



if params.bo_around_breakpoint:
    chm,st,en = [],[],[]
    interval_len = interval_len_from_breakpoints
    for line in open(file_breakpoint):
        if line[0] == 'c':
            continue
        a = line.split()
        c1,p1,c2,p2 = a[0],int(a[1]),a[3],int(a[4])
        chm += [c1]
        st += [ p1 - interval_len ]
        en += [ p1 + '_' + interval_len ]
        chm += [ c2 ]
        st += [ p2 - interval_len ]
        en += [ p2 + '_' + interval_len ]

    gr_interval_of_breakpoints = pr.from_dict({'Chromosome':chm, 'Start':st ,'End':en}).merge()

    # gr = pr.read_gtf('/confidential/FamilyR13_data/DATA/10x/case_17-08/phaseRNA/proteinCoding_2-5-11-16-18.gtf')
    gr = pr.read_gtf(params.gene_annotation)

    gr = gr[['gene_name']]
    s = False
    # for line in params.l_interval_string:
        # a = line.replace(":","-").split('-')
        # chm,st,en = 'chr'+a[0],int(a[1]),int(a[2])
    for k,v in gr_interval_of_breakpoints:
        for i in range(len(v)):
            s |= ( (gr.Chromosome == 'chr'+str(k)) & (gr.Start >= v.Start[i]) & (gr.End <= v.End[i]))
        # break        
    gr_selected_genes = gr[s]
    # print(s)
    # # print(gr)
    print(gr[s])
    df_selected_genes = gr[s].df


    chm, p = [],[]
    interval_len = interval_len_from_breakpoints
    for line in open(file_breakpoint):
        if line[0] == 'c':
            continue
        a = line.split()
        c1,p1,c2,p2 = a[0],int(a[1]),a[3],int(a[4])
        chm += [c1]
        p += [ p1 ]
        chm += [ c2 ]
        p += [ p2 ]
    def dist_to_bp(x):
        prev_min = 2000000
        for c_,p_ in zip(chm,p):
            if 'chr'+str(c_) == x.Chromosome:
                prev_min = min ( abs (x.Start - p_),abs (x.End - p_) , prev_min )
        # print (c_,p_,x.Chromosome,abs (x.Start - p_),abs (x.End - p_),prev_min,x.gene_name)

        return prev_min

    df_selected_genes['Dist2BP'] = df_selected_genes.apply(dist_to_bp,axis=1)
elif params.bo_permutation:
    chm,st,en = [],[],[]
    interval_len = 1000000
    # for line in open(file_breakpoint):
    #     if line[0] == 'c':
    #         continue
    #     a = line.split()
    #     c1,p1,c2,p2 = a[0],int(a[1]),a[3],int(a[4])
    #     chm += [c1]
    #     st += [ p1 - interval_len ]
    #     en += [ p1 + '_' + interval_len ]
    #     chm += [ c2 ]
    #     st += [ p2 - interval_len ]
    #     en += [ p2 + '_' + interval_len ]


    # gr = pr.read_gtf('/confidential/FamilyR13_data/DATA/10x/case_17-08/phaseRNA/proteinCoding_2-5-11-16-18.gtf')
    gr = pr.read_gtf(params.gene_annotation)

    gr = gr[['gene_name']]
    s = False
    for line in l_interval_string:
        a = line.replace(":","-").split('-')
        chm,st,en = ['chr'+a[0]],[int(a[1])],[int(a[2])]
    print(chm,st,en)
    gr_interval_of_breakpoints = pr.from_dict({'Chromosome':chm, 'Start':st ,'End':en}).merge()
    print (gr_interval_of_breakpoints)

    for k,v in gr_interval_of_breakpoints:
        # print(k)
        for i in range(len(v)):
            s |= ( (gr.Chromosome == k) & (gr.Start >= v.Start[i]) & (gr.End <= v.End[i]))
        # break
    # print(gr[s])   
    gr_selected_genes = gr[s]
    # print(s)
    # # print(gr)
    # print(gr[s])
    df_selected_genes = gr[s].df

    df_selected_genes['Dist2BP'] = None
else:
    # print (gr)
    print (interval)
    gr = gr[gr.Chromosome == 'chr'+interval]
    print(gr)
    df_selected_genes = gr.df
    # df_selected_genes = df_selected_genes[df_selected_genes['gene_name']=='IL17RA']
    df_selected_genes['Dist2BP'] = None
    print (df_selected_genes)

df_main_ = pd.DataFrame({'Wild':[], 'Chrm':[],'Pos':[],'WildIsRef':[], 'Rep': [],'Region':[],'Gene':[],'hasSNP':[],'Dist2BP':[]})


open_left_region, open_right_region = 0, 0



#######################################################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cnt = 0
for index, row in df_selected_genes.iterrows():
    chm, st, en,gene,dist2bp = row['Chromosome'], row['Start'], row['End'], row['gene_name'],row['Dist2BP']
    view_start = st - open_left_region
    view_end = en + open_right_region
    m,f,c,ref,rep,pss = count_reads(chm[3:],st,en,bams,vcf_reader,bo_phaseReadInBam,bams_phased_mother,bams_phased_father)
    # print (m,f,c,ref,rep)
    # break
    no = len(m) if len(m) > 0 else 1
    if len(m) > 0:
        d = {'Wild':m, 'Chrm':f,'Pos':c,'WildIsRef':ref, 'Rep': rep,'PS': pss,'Region':chm+':'+str(st)+'-'+str(en),'Gene':gene,'Chromosome':chm, 'hasSNP':'Yes','Dist2BP':dist2bp}
    else:
        d = {'Wild':[None], 'Chrm':[None],'Pos':[None],'WildIsRef':[None], 'Rep': [None],'PS':[None],'Region':chm+':'+str(st)+'-'+str(en),'Gene':gene,'Chromosome':chm, 'hasSNP':'No','Dist2BP':dist2bp }
    df_main_ = df_main_.append(pd.DataFrame(d))
    # print (d)
    
    # print(df_main.shape)
    cnt += 1
    if cnt % 50  == 0:
        # break
        print (cnt)
        # break
if bo_phaseReadInBam:
    for f in bams_phased_mother + bams_phased_father:
        f.close()
    for file_fo in f_bams_phased_mother +f_bams_phased_father:
        pysam.sort  (file_fo,'-o',file_fo[:-4]+'.sort.bam','-O','BAM')
        pysam.index (file_fo[:-4]+'.sort.bam')
        os.system('rm '+file_fo)
        


# df_main_['Ratio'] = df_main_['Mother'] / df_main_['Father']
print ('done')


df_main_['Case'] = case

if params.bo_around_breakpoint:
    df_main_.to_csv(working_dir+'RNAcount_chr_breakpoint'+str(interval_len_from_breakpoints)+'.csv')
else:
    if params.bo_permutation:
        df_main_.to_csv(working_dir+'RNAcount_chr'+interval+'.permutaion.csv')
    else:
        if bo_exon:
           df_main_.to_csv(out_df)
           open(interval+'.exonFinished','a').close()
        else:
            df_main_.to_csv(out_df)
            open(interval+'.geneFinished','a').close()
            


