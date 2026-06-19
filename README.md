# CRIMAC-FM-testdatapaper

This repository contains code to download, process and visualize the IMR test data sets.

The code lists avilable data sets, download them, convert from raw to pulse compressed data and plots the channels.

# Preparations

It is recommended to use git for obtaining the latest updates for the code. 

The code use `uv` for managing the environment. Installation methods for  `uv` is found [`here`](https://docs.astral.sh/uv/getting-started/installation/) 

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

Replace the data set id from the list provided by the first step. The example should download the data set, convert to pulse compressed data in NetCDF format and finally produce a png file of the pulsecompressed echogram.

```bash
uv run list
uv run get --dataset-id T2020003 --datadir /crimac-scratch/tmp/
uv run raw2pc --dataset-id T2020003 --datadir /crimac-scratch/tmp/
uv run pc2png --dataset-id T2020003 --datadir /crimac-scratch/tmp/
```
