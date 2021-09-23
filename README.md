# wavy

[![Tests (conda)](https://github.com/bohlinger/wavy/actions/workflows/python-conda-test.yml/badge.svg)](https://github.com/bohlinger/wavy/actions/workflows/python-conda-test.yml)
[![Tests (python)](https://github.com/bohlinger/wavy/actions/workflows/python.yml/badge.svg)](https://github.com/bohlinger/wavy/actions/workflows/python.yml)
[![Lint](https://github.com/bohlinger/wavy/actions/workflows/lint.yml/badge.svg)](https://github.com/bohlinger/wavy/actions/workflows/lint.yml)
[![Coverage Status](https://coveralls.io/repos/github/bohlinger/wavy/badge.svg?branch=master)](https://coveralls.io/github/bohlinger/wavy?branch=master)
[![Documentation Status](https://readthedocs.org/projects/wavyopen/badge/?version=latest)](https://wavyopen.readthedocs.io/en/latest/?badge=latest)

### Main developer and moderation:
Patrik Bohlinger, Norwegian Meteorological Institute, patrikb@met.no

## Purpose
Package to aid wave model validation using satellite altimetry and other sources. This README only illustrates some usage examples.

This workshop comprises:
1. Downloading satellite data
2. Quicklook examples
3. Usage examples on collocation and validation

## Additional info
The collocation method follows Bohlinger et al. (2019): https://www.sciencedirect.com/science/article/pii/S1463500319300435. The satellite data is obtained from http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001.

## Getting Started
### Installing wavy with conda
1. First download or clone the wavy github repository:
```
git clone --single-branch --branch wavyMini https://github.com/bohlinger/wavy.git
```
Info on how-to clone a repository:
https://help.github.com/en/articles/cloning-a-repository

2. To make it consistent with the description in this README please use as target location your home directory e.g.: ~/wavy.

3. Install wavy using conda using the environment.yml like:
```
cd ~/wavy
conda env create -f environment.yml
conda activate wavy
```
Info on installing conda, e.g.:
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

4. Configuration files are organized under wavy/config and might need adjustments according to your plans. Examples are the locations of your wave model output files or observation data (e.g. satellite altimetry data). What is needed for this workshop is shown below.

5. Prepare access to Copernicus products. Enter your account credentials into the .netrc-file. Your .netrc should look something like:
```
machine nrt.cmems-du.eu    login {USER}  password {PASSWORD}
```

6. Prepare your wavy environment with providing the directories for satellite data and model data. There are multiple config files but we only need to worry about a few for now. Explore the config file for satellites like this:
```
cd ~/wavy/config
vim satellite_specs.yaml
```
Add your path for satellite data here under cmems:
```
    cmems:
        level: 3
        satellite: s3a, s3b, c2, al, j3, h2b, cfo
        local:
            path: /home/patrikb/tmp_altimeter
```
Then download satellite data using the download.py script:
```
cd ~/wavy/apps/standalone
```
To get help check ...
```
./download.py -h
```
... then download some satellite altimeter data:
```
./download.py -sat s3a -sd 2020103000 -ed 2020111000
```

### HELP
Executable files usually have help function which can be read using e.g.:
```
./{YourExecutable}.py -h
```

e.g.:
```
cd ~/wavy/apps/standalone
./download.py -h
```
### Quicklook examples
The script "wavyQuick.py" is designed to provide quick and easy access to information regarding satellite coverage and basic validation.
0. Checkout the help:
```
cd ~/wavy/apps/standalone
./wavyQuick.py -h
```
1. Browse for satellite data and show footprints on map for one time step:
```
./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110112 --show
```
2. Browse for satellite data and show footprints on map for time period:
```
./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 --show
```
3. Same as above but now dump data to netcdf-file for further use and sharing
```
./wavyQuick.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 -dump /home/patrikb/tmp_altimeter/quickdump/
```
4. Browse for satellite data, collocate and show footprints and model for one time step:
```
./wavyQuick.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020110112 -lt 0 -twin 30 --col --show
```
You can also check wind speed instead of significant wave height by using the abbreviation U instead of Hs
```
./wavyQuick.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020110112 -lt 0 -twin 30 -var U --col --show
```
5. Use multiple satellite missions in one line e.g. all of them. First download the data. This time use the key "-nproc" to speed up the downloading.
```
./download.py -sat s3b -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat al -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat h3b -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat cfo -sd 2020103000 -ed 2020111000 -nproc 2
```
Note that the delivery of files from j3 and c2 is temporarily interupted. Now execute wavyQuick.py while altering the leadtime (-lt) and time constraints (-twin):
```
./wavyQuick.py -sat all -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 60 --col --show
```

6. Or list of satellites:
```
./wabyQuick.py -sat multi -l s3a,s3b,al -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 --col --show
```

7. Now concentrate on a subregion of your model domain by specifiying a region different from the model domain. The subregion is defined in region_specs.yaml.
```
./wavyQuick.py -sat s3a -reg NorwegianSea -mod mwam4 -sd 2020110112 -lt 0 -twin 30 --col --show
```
### Setup of operational usage: examples
1. Collocation and systematically dump to netcdf-file:
```
coming soon
```
2. Validation and systematically dump to netcdf-file:
```
coming soon
```
3. Plotting basic validation figures:
```
coming soon
```
