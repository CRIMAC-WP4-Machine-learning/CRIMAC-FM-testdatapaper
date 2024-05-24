# CRIMAC-FM-testdatapaper

This repository contains code to dowlonad, process and visluaize the IMR test data sets using a combination of Korona and python libraries. 


## Preparations

### Software

You will need a working installation of Korona and python. It is also recommended to use git for obtaining the latest updates.

The required python libraries are listed in the [requirements.txt](requirements.txt) file.

### Environmental variables

You need to set two environmental variables:

export CRIMACSCRATCH="/mnt/c/DATAscratch/crimac-scratch"
export LSSS=~/lsss/lsss-2.17.0-alpha

where the $CRIMACSCRATCH variable points to the root location of your local data storage and the $LSSS variable points to the location of the installed LSSS installation.


## Download test data

The file `DataSets.csv` contain the list of test data sets.

To obtain the test data you need to run `python3 get_data_S3.py`. The data will be downloaded from the IMR S3 server and placed under `${CRIMACSCRATCH}/CRIMAC-FM-testdata`. Each individual test data set will be placed under `${CRIMACSCRATCH}/{year}/{testdataset}`.

Note that the `get_data.py` moves the data from the backed up crimac storage to the S3 bucket, and must be run on one of IMR servers.

The `check_data.py` parses key diretories and count files by file extension per standard directory.


## Script for test data processing


### raw2pc.py - Convert raw data to pulse compressed data

Script to convert raw data to pulse compressed data in netcdf format.

The script reads the raw data for each test data at `ACOUSTIC/EK80/EK80_RAWDATA`, extract metadata from the raw files using ektools, and run the KoronaModule to convert to pulsecompressed data. The output is saved to `ACOUSTIC/GRIDDED` as net cdf files. A figure is generated from the netcdf file and placed in the same folder. When multiple ping groups are present in the data, a separate folder is generated for each ping group.

### raw2tracks.py - Tracking using Korona

The script reads the raw data for each test data set at `ACOUSTIC/EK80/EK80_RAWDATA`, the output is stored under `ACOUSTIC/LSSS/KORONA`.


### LSSS Manual tracking

TODO: Read results and merge with the Korona tracking results. Store in Netcdf. Where to store the code for reading?


### Image based tracking

Use Ingrids code to track samples belonging to same target across channels. 

TODO: Write annotation data set for single targets, similar to annotation data.

### pc2tsf - Estimate TS(f)

Based on the pulsecompresse data and the track definitions, estiamte TS(f) per target across all channels.


