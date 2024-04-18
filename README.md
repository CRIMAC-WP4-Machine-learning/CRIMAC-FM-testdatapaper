# CRIMAC-FM-testdatapaper

## Test data
The test data should be placed under the crimacs-cratch folder. Download the data and 
set the `$CRIMAC` environement variable to point to the top node of the data directory. 
Each test data set is placed under `$CRIMACSCRATCH/{year}/{testdataset}`.

Use `crimac = os.getenv('CRIMAC')` to access the path variable in python.

The file `DataSets.csv` contain the list of test data sets (subset of the `Data_sets.csv` list).

The `get_data.py` moves the data from crimac (backed up) to crimac-scratch/CRIMAC-FM-testdata (with S3 access). Must be run on one of IMR servers. 

The `get_data_S3.py` downloads the data from IMR S3 server.

The `check_data.py` parses key diretories and count files by file extension per standard directory.

## Libraries 
Read LSSS output:
https://github.com/CRIMAC-WP4-Machine-learning/lssstools

Run Korona through Python:
https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-KoronaScript



## Script for test data processing

### Data splitting
Code to split data (if necessary).

### raw2.pc.py - Convert raw data to pulse compressed data
Script to convert raw data to pulse compressed data in netcdf format.

The script reads the raw data for each test data at `ACOUSTIC/EK80/EK80_RAWDATA`, run the KoronaModule and convert to pulsecompressed data. The output is saved to `ACOUSTIC/GRIDDED` as net cdf files. A figure is generated from the netcdf file and placed in the same folder.

TODO: Agree on pc data format. We propose to use the variable names from Andersen et al. Change Korona to write this out. Korona must be updated to accomodate for this.

TODO: Need main frequency from metadata

### raw2tracks.py - Tracking using Korona

The script reads the raw data for each test data set at `ACOUSTIC/EK80/EK80_RAWDATA`, the output is stored under `ACOUSTIC/LSSS/KORONA`.

TODO: Set up the tracker

TODO: Need case specific parateters from metadata


### LSSS Manual tracking

Use LSSS for manual tracking. 

TODO: Read results and merge with the Korona tracking results. Store in Netcdf. Where to store the code for reading?


### Image based tracking

Use Ingrids code to track samples belonging to same target across channels. 

TODO: Write annotation data set for single targets, similar to annotation data.

### pc2tsf - Estimate TS(f)

Based on the pulsecompresse data and the track definitions, estiamte TS(f) per target across all channels.





