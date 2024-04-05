import pandas as pd
import os
import subprocess

# NB: Copy the DataSets.csv from https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-data-organisation/blob/main/DataSets.csv

# Read metadata
df = pd.read_csv('DataSets.csv')
testdata = df[(df['pulsetype'] == 'FM') &
              (df['dataset'].str[0] == 'T')]

testdata.to_csv('testdata.csv')

# Call whoami and capture the output
user = subprocess.run(['whoami'], capture_output=True, text=True).stdout.strip()
source = '/data/crimac/'
dest = '/mnt/c/DATAscratch/to_NR'

for _testdata in testdata['dataset']:
    print(_testdata)
    if not os.path.exists(dest + '/' + _testdata[1:5]+'/'):
        mkdir = 'mkdir '+dest + '/' + _testdata[1:5]+'/'
        print(mkdir)
        os.system(mkdir)
    scp = 'scp -r '+user+'@dedun.hi.no:'+source+_testdata[1:5]+'/'+_testdata + ' ' + dest + '/' + _testdata[1:5]+'/'
    print(scp)
    os.system(scp)

