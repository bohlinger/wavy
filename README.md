# wavy

### Main developer and moderation:
Patrik Bohlinger, Norwegian Meteorological Institute, patrikb@met.no

## Purpose
Package to aid wave model validation using satellite altimetry and other sources. This README does not describe all functionalities but is tailored to the Vietnam workshop comprising:
1. downloading satellite data
2. quicklook examples
3. usage examples on collocation and validation
4. simple html setup

## Additional info
The collocation method follows Bohlinger et al. (2019): https://www.sciencedirect.com/science/article/pii/S1463500319300435. The satellite data is obtained from http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001. In the end of this exercise validation figures are produced comparable to https://cmems.met.no/ARC-MFC/Wave3kmValidation/2020-10/index.html.

## Getting Started
### Installing wavy with conda
1. First download or clone the wavy github repository: https://github.com/bohlinger/wavy
Info on how-to clone a repository:
https://help.github.com/en/articles/cloning-a-repository
2. To make it consistent with the description in this README please use as target location your home directory e.g.: ~/wavy.
3. Install wavy using conda using the environment.yml like:
```
conda env create -f environment.yml
```
Additiona functionality can be investigated using:
```
conda create --help 
```
4. Configuration files are organized under wavy/config and might need adjustments according to your plans. Examples are the locations of your wave model output files or observation data (e.g. satellite altimetry data)

### HELP
Executable files usually have help function which can be read using e.g.:
./{YourExecutable}.py -h

### Quicklook examples
1. browse for satellite data and show footprints on map for one time step:
```
./check_sat.py -sat c2 -reg Vietnam -sd 2019100118 --show
```
2. browse for satellite data and show footprints on map for time period:
```
./check_sat.py -sat c2 -reg Vietnam -sd 2019100118 -ed 2019100318 --show
```
3. same as above but now dump data to netcdf-file for further use and sharing
```
./check_sat.py -sat c2 -reg Vietnam -sd 2019100118 -ed 2019100318 -dump /home/vietadm/wavyMini/data/altimetry/quickdump/
```
4. browse for satellite data, collocate and show footprints and model for one time step:
```
./check_sat.py -sat c2 -mod SWAN -reg Vietnam -sd 2019100118 -lt 30 -twin 30 -col --show
```
5. use multiple satellite missions in one line e.g. all of them:
```
./check_sat.py -sat all -mod SWAN -reg Vietnam -sd 2019100118 -lt 30 -twin 30 -col --show
```
6. or list of stellites:
```
./check_sat.py -sat multi -l s3a,c2,al -mod SWAN -reg Vietnam -sd 2019100118 -lt 30 -twin 30 -col --show
```
### Setup of operational usage: examples
1. Collocation and systematically dump to netcdf-file:
```
./collocate.py -sd 2019093012 -ed 2019100318 -sat c2 -mod SWAN -reg Vietnam
```
2. Validation and systematically dump to netcdf-file:
```
./validate.py -sd 2019093012 -ed 2019100318 -sat c2 -mod SWAN
```
3. Plotting basic validation figures:
```
./figures.py -d 201910 -sat c2 -mod SWAN
```
### Creating of simple html file for publishing and/or sharing results
