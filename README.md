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
cd ~/wavy/wavy
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
cd ~/wavy/wavy
./download.py -h
```
### Quicklook examples
The script "check_sat.py" is designed to provide quick and easy access to information regarding satellite coverage and basic validation. 
0. Checkout the help:
```
cd ~/wavy/wavy
./check_sat.py -h
```
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
5. Use multiple satellite missions in one line e.g. all of them. First download the data. This time use the key "-nproc" to speed up the downloading.
```
cd ~/wavy/wavy
./download.py -sat s3b -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat al -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat h3b -sd 2020103000 -ed 2020111000 -nproc 2
./download.py -sat cfo -sd 2020103000 -ed 2020111000 -nproc 2
```
Note that the delivery of files from j3 and c2 is temporarily interupted. Now execute check_sat.py while altering the leadtime (-lt) and time constraints (-twin):
```
./check_sat.py -sat all -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 60 --col --show
```

6. Or list of satellites:
```
cd ~/wavy/wavy
./check_sat.py -sat multi -l s3a,s3b,al -mod mwam4 -reg mwam4 -sd 2020110112 -lt 30 -twin 30 --col --show
```

7. Now concentrate on a subregion of your model domain by specifiying a region different from the model domain. The subregion is defined in region_specs.yaml.
```
cd ~/wavy/wavy
./check_sat.py -sat s3a -reg NorwegianSea -mod mwam4 -sd 2020110112 -lt 0 -twin 30 --col --show
```
### Customizing quicklook figures and save figure
The quicklook figures created with check_sat.py can be customized to some degree. The range of the colorbar can be adjusted as well as the levels of the contour lines. Point of interests can be introduced to navigate more easily. These settings can be introduced/altered in a config file called quicklook_specs.yaml in the config directory. So let's have a look there first:
```
cd ~/wavy/config
vim quicklook_specs.yaml
```
This file is organized such that all settings are tied to a region which in this example is the model domain swan_vietnam. As can be seen two point of interests (poi) are introduced, different colorbar levels, the index where the contour lines start, and the increment to specify for which values contour lines shall be drawn. The following line illustrates the introduced changes:
```
./check_sat.py -sat all -reg swan_vietnam -mod swan_vietnam -sd 2020110212 -lt 0 -twin 60 -dist 6 --col --show
```
Note that the contour lines start at Hs = 2m because this is the 9th index in the custom colormap levels (cm_levels) that we introduced into the config file.

Now do the same but save the figure to a defined location/directory:
```
./check_sat.py -sat all -reg swan_vietnam -mod swan_vietnam -sd 2020110212 -lt 0 -twin 60 -dist 6 --col --show -save /home/patrikb/tmp_pics
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
Add your own wave model to wavy. I already added the ecwam example, now add your SWAN model. This is done in the config file for models:
```
cd ~/wavy/config
vim model_specs.yaml
```
2. Exercise:
Try to execute the above described steps (quicklook and operational examples) using your own wave model.

3. Exercise:
Make your own html example using wavy. Add some describing text and use figures for your own wave model.

4. Add a preferred region for validation. Here the following config file for region specification is important:
```
cd ~/wavy/config
vim region_specs.yaml
```
In this file, either a regular region (lat/lon) or a polygon can be defined. Caution: when defining a polygon, the order of nodes is important.
