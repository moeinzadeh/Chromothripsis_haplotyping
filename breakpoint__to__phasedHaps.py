import gzip
import sys
file_PS = sys.argv[1] 
infile = sys.argv[2]
outfile = sys.argv[3]
outfile_break = sys.argv[4]

fo = open(outfile,'w')
fo_break =  open(outfile_break,'w')

mp_PSW_hap = {}


lines = open(file_PS).readlines()
if len(lines) >= 1:
    for l in lines[0].split():
        print (l)
        if l[0] == 'W':
            continue
        a = l.split(',')
        mp_PSW_hap[a[0]] = a[1]
    # print (mp_PSW_hap)

for l in gzip.open (infile,'rt'):
    if l[0] == '#':
        fo.write(l)
        fo_break.write(l)
        continue
    a = l.split()
    b = a[9].split(":")

    if b[0][1] == '/':
        fo.write(l)
        continue

    H1 = (b[8],b[0][0])
    H2 = (b[8],b[0][1])

    if b[8] in mp_PSW_hap:
        # print (mp_PSW_hap[b[8]])
        if mp_PSW_hap[b[8]] == '2':
            # print ('swap needed in current VCF')
            tmp_GT = b[0][2]+'|'+b[0][0]
            b[0] = tmp_GT

            tmp_AD = b[2].split(',')[1]+','+b[2].split(',')[0]
            b[2] = tmp_AD
        else:
            # print ('matchin with current VCF')
            pass

        b[8] = '1'
        tmp  = ':'.join(b)
        b = tmp
        a[9] = b
        tmp = '\t'.join(a)
        fo.write(tmp+'\n')
        fo_break.write(tmp+'\n')
    else:
        fo.write(l)
    # break
fo.close()
fo_break.close()