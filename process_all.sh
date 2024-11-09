#!/bin/bash
# Script to process the files in testdata.csv one by one adn collect output/error

tail -n +2  $CRIMACSCRATCH/CRIMAC-FM-testdata/testdata.csv | cut -d, -f1 | while read dir; do
    echo $dir
    python raw2pc.py $CRIMACSCRATCH/CRIMAC-FM-testdata/*/$dir/ACOUSTIC/EK80/EK80_RAWDATA/ tmp_out/$dir > $dir.OUT 2> $dir.ERR
done
