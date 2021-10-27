import sys

chr = sys.argv[2]
cr = chr[3:]
if cr == 'X':
    cr = '23'
cnt_ps=int(cr)*100000
mp = {}
for l in open(sys.argv[1]):
    l = l.rstrip()
    if l[0] == '#':
        print(l)
        continue
        
    variant = l.split()
    if variant[9][1] == '|':

        variant_format = variant[8]
        variant_info = variant[9]

        infoTag = variant_format.split(':')
        info = variant_info.split(':')
        mp_info = dict ( zip (infoTag, info ))
        mp_info_inx = dict (zip (infoTag, range(len(infoTag))))
        ps = mp_info['PS']
        ps_inx = mp_info_inx['PS']

        if ps in mp:
            info [ps_inx] = mp[ps]
        else:
            mp[ps] = str(cnt_ps)
            info[ps_inx] = str(cnt_ps)
            cnt_ps += 1
        variant[9] =':'.join(info)
    print ('\t'.join(variant))