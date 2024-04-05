# CRIMAC-FM-testdatapaper

## Test data
The test data should be placed under the crimacs-cratch folder. Download the data and 
set the `$CRIMAC` environement variable to point to the top node of the data directory. 
Each test data set is placed under `$CRIMACSCRATCH/{year}/{testdataset}`.

Use `crimac = os.getenv('CRIMAC')` to access the path variable in python.

The file `DataSets.csv` contain the list of test data sets (subset of the `Data_sets.csv` list).

The `get_data.py` downloads the data from the central IMR data repository. Will only work inside IMRs filre wall.

The `check_data.py` parses key diretories and count files by file extension per standard directory.


## Script for test data processing


