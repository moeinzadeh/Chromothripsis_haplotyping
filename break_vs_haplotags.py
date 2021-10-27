import pysam
from collections import Counter
import pandas as pd
import seaborn as sn
import os
from bisect import bisect_left 
import sys
import pdb


interval = 90
mappingQual = 20

# working_dir = '/confidential/FamilyR13_data/DATA/10x/case_17-08/local_denovo_asm/extract_sick_region/main_region/cluterPBreads_BridgeBreak_onlyPB/'



propagated = False
primary = False
new_primary = True
debug = False
to_plot=True



working_dir = sys.argv[1]
file_tagged_bam =  sys.argv[2]
file_breakpoints =  sys.argv[3]


print (working_dir,file_tagged_bam,file_breakpoints)

samfile = pysam.AlignmentFile(file_tagged_bam)
# statfile = working_dir + 'phasesetStat_onBreakpoints.tsv'


def cluster_reads(bo_just_PSTable):


    

    l_wildtype = {}
    l_chromocriptic = {}
    l_breakpoint_flancking = []



    bridge_all_readname = set([])
    break_all_readname = set([])

    
    fo_breakPoints_readTags_eachBreak = open(os.path.join(working_dir, 'Breakpoint_haplotype_PS.txt'),'w')
    fo_breakPoints_readTags_eachBreak.write('\t'.join(['chrm','bp', 'interval_size','PS' , 'bridge_tagged_H1', 'bridge_tagged_H2', 'bridge_untagged', 'break_tagged_H1', 'break_tagged_H2', 'break_untagged'])+'\n')


    for line in open (file_breakpoints):
        if line.startswith('chrom1'): # drop first line
            continue

        a = line.split()

        # if a[11] in  ['no' ,'No']:
        #     print (line,'\nbreakpoint not passed')
        #     continue


        chr1, bp1, chr2, bp2 = a[0], int(a[1]), a[3], int(a[4])

        l_loci = []
        if (chr1, bp1) in l_loci or (chr2, bp2) in l_loci:
            print ('---------------------------------')
            continue
        l_loci += [(chr1, bp1),(chr2, bp2)]

        for chrm,bp in [(chr1, bp1),(chr2, bp2)]:

            left  = set ([ read.qname for read in samfile.fetch(chrm, bp - interval, bp - interval + 1  ) if read.mapping_quality >= mappingQual ])
            right = set ([ read.qname for read in samfile.fetch(chrm, bp + interval, bp + interval + 1  ) if read.mapping_quality >= mappingQual ])

            bridge = left & right
            leftonly = left - right
            rightonly = right - left

            breakk = (leftonly | rightonly) - bridge

            if debug:
                fo_debug.write('---start-------------'+'\n')
                for i in bridge:
                    fo_debug.write(i+'\t'+'bridge'+'\t'+ str(chrm)+' ' +str(bp)+ ' ' + str( [(chr1, bp1),(chr2, bp2)])+'\n')
                fo_debug.write('----------------'+'\n')
                for i in breakk:
                    fo_debug.write(i+'\t'+'break'+'\t'+ str(chrm)+' ' +str(bp) +'\n')
                fo_debug.write('---end-------------'+'\n')
                
                print (breakk & bridge)
                print ('********')


            bridge_tagged, left_tagged, right_tagged = 0, 0, 0
            bridge_untagged, left_untagged, right_untagged = 0, 0, 0
            l_bridge_tagged, l_left_tagged, l_right_tagged = [],[],[]

            
            bridge_tagged_H1, bridge_tagged_H2 = 0, 0
            break_tagged_H1, break_tagged_H2, break_untagged = 0, 0, 0
            PS = 'NaN'

            for read in samfile.fetch(chrm, bp - interval, bp + interval + 1  ):
                l_breakpoint_flancking+=[read.qname]
                if read.mapping_quality < mappingQual:
                    continue

                if read.qname in bridge:
                    if read.has_tag('HP'):
                        # print (read)
                        l_bridge_tagged += [(read.get_tag('HP'),read.get_tag('PS') if read.has_tag('PS') else 'NAN')]
                        PS = read.get_tag('PS') if read.has_tag('PS') else 'NAN'
                        HP = read.get_tag('HP')
                        l_wildtype [(PS,HP)] = l_wildtype.get ((PS,HP), 0) + 1
                        if HP == 1:
                           bridge_tagged_H1 += 1
                        if HP == 2:
                            bridge_tagged_H2 += 1  

                        bridge_tagged += 1
                        
                    else:
                        bridge_untagged += 1
                elif read.qname in breakk:
                    if read.qname in leftonly:
                        if read.has_tag('HP'):
                            l_left_tagged += [(read.get_tag('HP'),read.get_tag('PS') if read.has_tag('PS') else 'NAN')]
                            left_tagged += 1
                            PS = read.get_tag('PS') if read.has_tag('PS') else 'NAN'
                            HP = read.get_tag('HP')
                            l_chromocriptic [(PS,HP)] = l_chromocriptic.get ((PS,HP), 0) + 1

                        else:
                            left_untagged += 1

                    elif read.qname in rightonly:
                        if read.has_tag('HP'):
                            l_right_tagged += [(read.get_tag('HP'),read.get_tag('PS') if read.has_tag('PS') else 'NAN')]
                            right_tagged += 1
                            PS = read.get_tag('PS') if read.has_tag('PS') else 'NAN'
                            HP = read.get_tag('HP')
                            l_chromocriptic [(PS,HP)] = l_chromocriptic.get ((PS,HP), 0) + 1
                    
                        else:
                            right_untagged += 1

                    if read.has_tag('HP') and read.get_tag('HP') == 1:
                        break_tagged_H1 += 1
                    elif read.has_tag('HP') and read.get_tag('HP') == 2:
                        break_tagged_H2 += 1    


                    break_untagged = right_untagged + left_untagged


            st_tmp = '\t'.join([str(i) for i in [chrm,bp, interval,PS, bridge_tagged_H1, bridge_tagged_H2, bridge_untagged, break_tagged_H1, break_tagged_H2, break_untagged]]) + '\n'
            fo_breakPoints_readTags_eachBreak.write(st_tmp)
            


            bridge_all_readname |= bridge
            break_all_readname |= breakk
    PSset_c = set([ i[0] for i in l_chromocriptic])
    PSset_w = set([ i[0] for i in l_wildtype])
    print (len(PSset_c | PSset_w))
    print (PSset_c - PSset_w)
    print (PSset_w - PSset_c)

    # print (l_wildtype)
    print ('conflict in haplotags---l_wildtype')
    file_w = os.path.join(working_dir, 'PS_wild.txt')
    fo = open(file_w,'w')
    fo.write('PS\tWHP1\tWHP2\n')
    for ps in set([ i[0] for i in l_wildtype]):
        if ((ps,1) in l_wildtype) and ((ps,2)in l_wildtype):
            fo.write (str(ps)+'\t'+str(l_wildtype [(ps,1)])+'\t'+str( l_wildtype [(ps,2)] )+'\n')
        elif (ps,1) in l_wildtype:
            fo.write (str(ps)+'\t'+str(l_wildtype [(ps,1)])+'\t0'+'\n')
        elif (ps,2) in l_wildtype:
            fo.write (str(ps)+'\t0\t'+str(l_wildtype [(ps,2)])+'\n')
    fo.close()

    print ('conflict in haplotags---l_chromocriptic')            
    file_c  = os.path.join(working_dir, 'PS_chromo.txt')
    fo = open(file_c,'w')
    fo.write('PS\tCHP1\tCHP2\n')
    for ps in set([ i[0] for i in l_chromocriptic]):
        if ((ps,1) in l_chromocriptic) and ((ps,2)in l_chromocriptic):
            fo.write (str(ps)+'\t'+str(l_chromocriptic [(ps,1)])+'\t'+str( l_chromocriptic [(ps,2)] )+'\n')
        elif (ps,1) in l_chromocriptic:
            fo.write (str(ps)+'\t'+str(l_chromocriptic [(ps,1)])+'\t0'+'\n')
        elif (ps,2) in l_chromocriptic:
            fo.write (str(ps)+'\t0\t'+str(l_chromocriptic [(ps,2)])+'\n')
    fo.close()
    print ('conflict in haplotags---inter check')            
    file_conflict_WILDandCHRM = os.path.join(working_dir , 'PS_Conflict_betweenWild_Chrm.txt')
    fo = open(file_conflict_WILDandCHRM,'w')
    fo.write('PS\tCHP\tWHP\n')
    for ps in set([ i[0] for i in l_chromocriptic]):
        if ((ps,1) in l_chromocriptic) and ((ps,1)in l_wildtype):
            fo.write (str(ps)+'\t'+'1_'+str(l_chromocriptic [(ps,1)])+'\t'+'1_'+str( l_wildtype [(ps,1)] )+'\n')
        if ((ps,2) in l_chromocriptic) and ((ps,2)in l_wildtype):
            fo.write (str(ps)+'\t'+'2_'+str(l_chromocriptic [(ps,2)])+'\t'+'2_'+str( l_wildtype [(ps,2)] )+'\n')
    fo.close()


    df_w = pd.read_csv( file_w, sep='\t')
    df_c = pd.read_csv( file_c, sep='\t')
    os.system('cat '+file_conflict_WILDandCHRM)


    # df_all = pd.concat([df_w, df_c], axis=1)
    # df_all = pd.merge(df_w, df_c, left_index=True, right_index=True, how='outer')
    df_all = pd.merge(df_w, df_c, on='PS', how='outer')
    df_all.fillna(0, inplace=True)



    df_all['senario1'] = df_all['WHP1'] + df_all['CHP2']
    df_all['senario2'] = df_all['WHP2'] + df_all['CHP1']
    def decide (x):
        if x['senario1'] > x['senario2']:
            x['decision_W'] = x['PS'] , 1
            x['decision_C'] = x['PS'] , 2
        else:
            x['decision_W'] = x['PS'] , 2
            x['decision_C'] = x['PS'] , 1
        return x
    df_all = df_all.apply(decide,axis=1)


    if to_plot:
        ax = df_all.plot.scatter(x='senario1',
                            y='senario2',
                            alpha=0.3,figsize=(6,6))   
        ax.set_title('Senario 1: W -> PhaseSet1 & C -> PhaseSet2\nSenario 2: W -> PhaseSet2 & C -> PhaseSet1\n distances to the axes mean higher error')
        mx_lim = 1200
        ax.set_xlim((-10,mx_lim))
        ax.set_ylim((-10,mx_lim))
        figure = ax.get_figure()
        figure.savefig(os.path.join(working_dir , 'conflict_of_halotags.png'))



    print (df_c)
    print (df_all)

    l_W_final = [ (int(i[0]),i[1])for i in df_all['decision_W'].tolist()]
    l_C_final = [ (int(i[0]),i[1])for i in df_all['decision_C'].tolist()]

    print(l_W_final)
    print(l_C_final)
    
    fo = open(os.path.join(working_dir , 'PS_final.txt'),'w')
    for i in l_W_final:
        fo.write(str(i[0])+'\t'+str(i[1])+'\tW\n')
    for i in l_C_final:
        fo.write(str(i[0])+'\t'+str(i[1])+'\tC\n')
    fo.close()

    print (l_W_final)

    fo = open(os.path.join(working_dir , 'PS_w.txt'),'w')
    fo.write('Wild\t')
    for i in l_W_final[:-1]:
        fo.write(str(i[0])+','+str(i[1])+' ')
    fo.write(str(l_W_final[-1][0])+','+str(l_W_final[-1][1]))



    if bo_just_PSTable:
        print ('--------------------------')
        return

    

    file_WildReads = os.path.join(working_dir , 'Wild.bam')
    file_ChrmReads = os.path.join(working_dir , 'Chrm.bam')
    fo_Wild = pysam.AlignmentFile(file_WildReads, "wb",template=samfile)
    fo_Chrm = pysam.AlignmentFile(file_ChrmReads, "wb",template=samfile)

    file_bridge = os.path.join(working_dir , 'Bridge.bam')
    file_break  = os.path.join(working_dir , 'Break.bam')
    fo_bridge = pysam.AlignmentFile(file_bridge, "wb",template=samfile)
    fo_break  = pysam.AlignmentFile(file_break, "wb",template=samfile)



    cnt = 1
    cnt_snp_depleated = 1
    set_l_breakpoint_flancking = set(l_breakpoint_flancking)
    # l_breakpoint_flancking.sort() 



    for read in pysam.AlignmentFile(file_tagged_bam, "rb"):
    # for read in pysam.AlignmentFile(working_dir + "tmp_s.bam", "rb"):
        # for read in pysam.AlignmentFile("/confidential/FamilyR13_data/DATA/10x/case_17-08/local_denovo_asm/extract_sick_region/shattered_regions_haplotaged_sorted.bam", "rb"):
        if read.qname in bridge_all_readname:
            fo_bridge.write(read)
            fo_Wild.write(read)
        elif read.qname in break_all_readname:
            fo_break.write(read)
            fo_Chrm.write(read)
        elif read.has_tag('HP') and read.has_tag('PS'):
            if (read.get_tag('PS'),read.get_tag('HP')) in l_W_final:
                fo_Wild.write(read)
                # print  ('W')
            elif (read.get_tag('PS'),read.get_tag('HP')) in l_C_final:
                fo_Chrm.write(read)
                # print  ('C')
        else:
            # if bisect_left(l_breakpoint_flancking, read.qname): 
            if read.qname not in set_l_breakpoint_flancking:
                fo_Wild.write(read)
                fo_Chrm.write(read)
                cnt_snp_depleated += 1

        cnt += 1
        if cnt % 20000 == 0:
            print ('checked reads:',cnt,'\treads from snp deplited regions',cnt_snp_depleated)

    fo_Wild.close()
    fo_Chrm.close()
    fo_break.close()
    fo_bridge.close()

    for f in [file_WildReads,file_ChrmReads,file_bridge,file_break]:
        print ('sorting ',f)

        pysam.sort(f,'-@','10','-o',os.path.join(working_dir ,'tmp.bam'),'-O','BAM')
        os.rename(os.path.join(working_dir ,'tmp.bam',f))
        pysam.index(f)
    os.system('rm ' + os.path.join( working_dir ,'tmp.bam'))









cluster_reads(True)
