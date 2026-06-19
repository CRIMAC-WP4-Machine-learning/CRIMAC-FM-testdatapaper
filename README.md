# CRIMAC-FM-testdatapaper

This repository contains code to download, process and visualize the IMR test data sets.

## Preparations

### Software

It is recommended to use git for obtaining the latest updates for the code. 

The code use `uv` for managing the environment. Installation methods for  `uv` is found [`here`](https://docs.astral.sh/uv/getting-started/installation/) 

### Environmental variables

You need to set the location of the test data as an env variable:

```bash
export CRIMACSCRATCH="/crimac-scratch"`
```

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
