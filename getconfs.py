import ektools as E
import ektools.actions as A
import sys


ix = E.index('crimac-scratch/CRIMAC-FM-testdata/2023/T2023002/ACOUSTIC/EK80/EK80_RAWDATA/D20230803-T230004.raw')

p = E.parse(ix[0][3])
d = p['configuration']

p2 = E.parse(ix[1][3])
d2 = p2['initialparameter']

for k in d.keys():
    tmp = d[k].pop('raw_xml')
    with open(k+'.txt', 'w') as sys.stdout:
       A.showdict(d[k])
       print()
       A.showdict(d2[k])

