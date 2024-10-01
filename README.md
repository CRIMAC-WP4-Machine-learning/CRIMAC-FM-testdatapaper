# CRIMAC-FM-testdatapaper

This repository contains code to dowlonad, process and visluaize the IMR test data sets using a combination of Korona and python libraries. 


## Preparations

### Software

It is recommended to use git for obtaining the latest updates for the code.

You will need a working installation of Korona and python. The code has been tested on Python 3.8.

The required python libraries for python 3.8 are listed in the [requirements.txt](requirements.txt) file. To install the packages you can run `pip install -r requirements.txt` from the terminal.

You need to install LSSS/Korona, follow the instructions here: [link](https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-KoronaScript/blob/master/README.md#install-lssskorona)

You also need the OpenGL libraries for running the tracking like so: `sudo apt install libgl1-mesa-glx`

It is recommended to use python environments. You can generate an environment by `python3 -m venv CRIMAC_testdata`, where `CRIMAC_testdata` is the name of the environment. Note that the environment will be installed in the current directory. It is recommented to create a separate directory to organize different environments, e.g. by `mkdir pyvenv` and enter the directory by `cd pyenv` before creating it. 

After generating the environment you need to activate it by typing `source ./CRIMAC_testdata/bin/activate`. 

### Environmental variables

You need to set two environmental variables:

export CRIMACSCRATCH="/mnt/c/DATAscratch/crimac-scratch"
export LSSS=~/lsss/lsss-2.17.0-alpha

where the $CRIMACSCRATCH variable points to the root location of your local data storage and the $LSSS variable points to the location of the installed LSSS installation.


## Download test data

The file `DataSets.csv` contain the list of test data sets.

To obtain the test data you need to run `python3 get_data_S3.py`. The data will be downloaded from the IMR S3 server and placed under `${CRIMACSCRATCH}/CRIMAC-FM-testdata`. Each individual test data set will be placed under `/{year}/{testdataset}`.

Note that the `get_data.py` moves the data from the backed up crimac storage to the S3 bucket, and must be run on one of IMR servers.

The `check_data.py` parses key diretories and count files by file extension per standard directory.


## Script for test data processing


### raw2pc.py - Convert raw data to pulse compressed data

Script to convert raw data to pulse compressed data in netcdf format.

The script reads the raw data for each test data at `ACOUSTIC/EK80/EK80_RAWDATA`, extract metadata from the raw files using ektools, and run the KoronaModule to convert to pulsecompressed data. The output is saved to `ACOUSTIC/GRIDDED/pc_{i}` as net cdf files corresponding to the raw files, where {i} is the ping group number. When multiple ping groups are present in the data, {i} denotes the ping group, otherwise i=1. A figure is generated from the netcdf file and placed in the same folder. .

### raw2tracks.py - Tracking using Korona

The script reads the raw data for each test data set at `ACOUSTIC/EK80/EK80_RAWDATA`, the output is stored under `ACOUSTIC/LSSS/KORONA`.


### pc2annotations.py - Tracking using image based methods

Use Ingrids code to track samples belonging to same target across channels. 

TODO: The script reads the pulse compressed data for each test data set at `ACOUSTIC/GRIDDED/pc_{i}`, run the tracking code and stores the output under `ACOUSTIC/GRIDDED/pc_{i}/tracks.csv`.

The tracking code can be found here: https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-fm-sed

Frå Ingrid: 

Vedlagt er et eksempel på hvordan jeg ser for meg at output fra tracking-algoritmen kan skrives, basert på tidligere diskusjoner. 
En kort beskrivelse av innholdet finnes i README her: https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-fm-sed : Single Echo Detection algorithm on multi channel frequency modulated echosounder data (github.com).

Jeg har også pushet eksempelkode for hvordan csv-fila kan plottes sammen med data: https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-fm-sed/scripts/visualization/visualize_track_csv.py at main. Koden skal produsere det vedlagte bildet. 

Si gjerne ifra hvis det er noe dere ønsker å endre på i output-formatet. 
Fila inneholder ikke gode tracks – det er mye rom for forbedringer i algoritmen som produserte dette. 


### pc2tsf - Estimate TS(f)

Based on the pulsecompressed data and the track definitions, estimate TS(f) per target across all channels.

The tracking code can be found here: ?

TODO: The script reads the pulse compressed data and track definitions for each test data set at `ACOUSTIC/GRIDDED/pc_{i}`, run the TS(f) estimation and store the output under `ACOUSTIC/GRIDDED/pc_{i}/??.??`.



