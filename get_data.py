import os
import csv

# this script copy the data from the backed up crimac stire to the S3 via crimac scratch

# Read metadata
testdata = []
with open("testdata.csv", 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)
    for row in csvreader:
        testdata.append(row[1])

# Call whoami and capture the output
source = 'dedun:/data/crimac/'
dest = 'crimac-scratch/CRIMAC-FM-testdata/'

for _testdata in testdata:
    print(_testdata)
    if not os.path.exists(dest + '/' + _testdata[1:5]+'/'):
        mkdir = 'mkdir -p '+dest + '/' + _testdata[1:5]+'/'
        print(mkdir)
        os.system(mkdir)
    scp = 'scp -r '+source+_testdata[1:5]+'/'+_testdata + ' ' + dest + _testdata[1:5]+'/'
    print(scp)
    os.system(scp)

