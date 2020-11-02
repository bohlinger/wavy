# wavy

### Main developer and moderation:
Patrik Bohlinger, Norwegian Meteorological Institute, patrikb@met.no

## Purpose
Package to aid wave model validation using satellite altimetry and other sources. This README does not describe all functionalities but is tailored to the Vietnam workshop. There are still files remaining due to legacy which will not be discussed here. They are still included to make an upgrade to the next version easier for the user.

This workshop comprises:
1. Downloading satellite data
2. Quicklook examples
3. Usage examples on collocation and validation
4. Simple html setup with basic validation figures

## Additional info
The collocation method follows Bohlinger et al. (2019): https://www.sciencedirect.com/science/article/pii/S1463500319300435. The satellite data is obtained from http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001. In the end of this exercise validation figures are produced comparable to https://cmems.met.no/ARC-MFC/Wave3kmValidation/2020-10/index.html.

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
```
Info on installing conda, e.g.:
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

4. Configuration files are organized under wavy/config and might need adjustments according to your plans. Examples are the locations of your wave model output files or observation data (e.g. satellite altimetry data)

5. Prepare your wavy environment with providing the directories for satellite data and model data. There are multiple config files but we only need to worry about a few for now. Explore the config file for satellites like this:
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
cd ~/wavy/wavy
```
To get help check ...
```
./download.py -h
```
... then download some satellite altimeter data:
```
./download.py -sat s3a -sd 2020110100 -ed 2020111000
```

### HELP
Executable files usually have help function which can be read using e.g.:
```
./{YourExecutable}.py -h
```

e.g.:
```
cd ~/wavy/wavy
./download.py -h
```
### Quicklook examples
The script "check_sat.py" is designed to provide quick and easy access to information regarding satellite coverage and basic validation. 
1. Browse for satellite data and show footprints on map for one time step:
```
cd ~/wavy/wavy
./check_sat.py -sat s3a -reg mwam4 -sd 2020110112 --show
```
2. Browse for satellite data and show footprints on map for time period:
```
cd ~/wavy/wavy
./check_sat.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 --show
```
3. Same as above but now dump data to netcdf-file for further use and sharing
```
cd ~/wavy/wavy
./check_sat.py -sat s3a -reg mwam4 -sd 2020110100 -ed 2020110300 -dump /home/patrikb/tmp_altimeter/quickdump/
```
4. Browse for satellite data, collocate and show footprints and model for one time step:
```
cd ~/wavy/wavy
./check_sat.py -sat s3a -reg mwam4 -mod mwam4 -sd 2020110112 -lt 0 -twin 30 --col --show
```
5. Use multiple satellite missions in one line e.g. all of them:
```
cd ~/wavy/wavy
./check_sat.py -sat all -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 --col --show
```
(Caution: As for now, all chosen satellites should have some data, if not there will be an error)
6. Or list of stellites:
```
cd ~/wavy/wavy
./check_sat.py -sat multi -l s3a,c2,al -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 -col --show
```

### Setup of operational usage: examples
1. Collocation and systematically dump to netcdf-file:
```
cd ~/wavy/op/support
./op_collocate.py -sat s3a -mod mwam4 -sd 2020110100 -ed 2020111000 -path ~/tmp_validation/
```
2. Validation and systematically dump to netcdf-file:
```
cd ~/wavy/op/support
./op_validate.py -sat s3a -mod mwam4 -sd 2020110100 -ed 2020111000 -path ~/tmp_validation/
```
3. Plotting basic validation figures:
```
cd ~/wavy/op/support
./op_figures.py -sat s3a -mod mwam4 -d 202011 -path ~/tmp_validation/
```
### Creating of simple html file for publishing and/or sharing results
```
cd ~/wavy/web/op
sh webpage.sh -i ~/tmp_validation -w ~/wavy/web/op -m mwam4 -s s3a -Y 2020 -M 11
```
Now, copy and paste the location of the index.html into your browser to test the html-page.

### Exercises
1. Exercise:
Add your own altimeter directory and wave model to wavy.
2. Exercise:
Try to execute these commands using your own wave model.
3. Exercise:
Make your own html example using wavy.
