# CRIMAC-FM-testdatapaper

This repository contains code to download, process and visualize the IMR test data sets.

The code lists avilable data sets, download them, convert from raw to pulse compressed data and plots the channels.

# Preparations

It is recommended to use git for obtaining the latest updates for the code. 

The code use `uv` for managing the environment. Installation instructions for  `uv` is found [`here`](https://docs.astral.sh/uv/getting-started/installation/)

On Linux, you need to ensure that you have NetCDF functionality installed:
```bash
sudo apt install netcdf-bin
```

Some of the processing modules depend on the KORONA library, shipped as part of the [LSSS acoustics analysis software](https://marec.no/downloads.htm) package.
This should be handled automatically via the dependency on the [KoronaScript](https://github.com/CRIMAC-WP4-Machine-learning/CRIMAC-KoronaScript) library, but it is also possible to specify a separate LSSS installation to use, via the LSSS environment variable.
E.g., you can do
```bash
export LSSS=/opt/lsss-3.2.0
```
or, in Windows PowerShell, 
```powershell
$env:LSSS = 'C:\Program Files\Marec\LSSS 3.2.0'
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
Note that a bug in the Python interface to the `netCDF4` library, it is necessary to specify the drive letter on Windows.

# Testing 

This section is for the developers.

## Linting and type checking

**Formatting, linting, and type checking** are automated ways to improve code quality at different levels. **Formatting** ensures code has a consistent layout (spacing, line breaks, ordering) so it’s easy to read and review, without changing behavior. **Linting** analyzes code for likely mistakes, bad patterns, or style issues—such as unused variables or error-prone constructs—and often suggests or applies fixes. **Type checking** goes deeper by verifying that values are used consistently with their declared types, helping catch subtle bugs that only appear at runtime in dynamically typed languages like Python.


| Tool                        | Purpose                      |
| --------------------------- | ---------------------------- |
| **Formatter** (ruff format) | Changes how code looks       |
| **Linter** (ruff check)     | Finds suspicious or bad code |
| **Type checker** (ty)     | Checks if types make sense   |

Usual order:

1. **Format**
2. **Lint**
3. **Type-check**

### Reformatting and linting

* `uv run ruff format .` → Reformatting code
* `uv run ruff check .` → linting
* `uv run ruff check . --fix` → linting with **safe** auto-fixes

### Type checking

* `uv run ty check` → check for typing errors


## Testing

Run all tests with:

```bash
uv run pytest
```


### Integration tests

The integration tests provide an end-to-end smoke test of the data processing pipeline. For each dataset currently available from the data repository, the tests:

1. download and unpack the dataset with `get` to the pytest `/tmp` folder,
2. convert the raw EK80 data to pulse-compressed NetCDF with `raw2pc`, and
3. generate echogram images with `pc2png`.

The integration tests also verify that requesting a non-existing dataset returns an error.

Run the integration tests for all test data sets:

```bash
uv run pytest -m integration
```

The integration tests require network access and may approximately 3 hours, as all the test datasets are downloaded and processed. As a firste step, run instead a single integration test on one test data set:

```bash
uv run pytest "tests/integration/test_smoke.py::test_full_pipeline[T2019001]" -vv
```
