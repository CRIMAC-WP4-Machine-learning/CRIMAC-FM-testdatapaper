# CRIMAC-FM-testdatapaper

This repository contains code to download, process and visualize the IMR test data sets.

The code lists avilable data sets, download them, convert from raw to pulse compressed data and plots the channels.

# Preparations

It is recommended to use git for obtaining the latest updates for the code. 

The code use `uv` for managing the environment. Installation instructions for  `uv` is found [`here`](https://docs.astral.sh/uv/getting-started/installation/)

# Scripts

## List test data
List avilable tests data sets. Run
```bash
uv run list -h
```
for instructions.

## Get test data

Download a test data set. Run
```bash
uv run get -h
```
for instructions.

## Validating data downloads

Each data set comes with a text file containing a list of the data set contents with SHA-256 checksums.  E.g., completeness and integrity of the T2019001 dataset can be verified with the command `sha256sum -c T2019001-sha256.txt` on Linux (or Windows with WSL), or with `certUtil -hashfile T2019001-sha256.txt sha256` on Windows.

## Preprocess test data

Preprocesses the test data to pulse compressed and store as netcdf. Run
```bash
uv run raw2pc -h
```
for instructions.

## Visualisation

Plot pulse compressed data. Run
```bash
uv run pc2png -h
```
for instructions.

# Example

Replace the data set id from the list provided by the first step. The example should download the data set to the `/tmp/crimac-scratch' directory, convert to pulse compressed data in NetCDF format, and finally produce a png file of the pulse compressed echogram.

```bash
uv run list
uv run get --dataset-id T2020003 --datadir /tmp/crimac-scratch/
uv run raw2pc --dataset-id T2020003 --datadir /tmp/crimac-scratch/
uv run pc2png --dataset-id T2020003 --datadir /tmp/crimac-scratch/
```
