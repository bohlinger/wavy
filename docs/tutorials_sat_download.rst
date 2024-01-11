Download Satellite data
#######################

One strength of **wavy** is to ease obtaining and using satellite altimetry data. To get an overview over the supported missions or data sources go to wavy/apps/standalone and execute:

.. code-block:: bash

   $ ./wavyDownload.py -h


or:
.. code-block:: bash

   $ ./wavyQuick.py -h


The help-message displayed would give you, among other information, the following options:


.. code-block:: bash

    cmems_L3_NRT:            
     s3a - Sentinel-3A            
     s3b - Sentinel-3B            
     j3 - Jason-3 (reference mission)            
     c2 - Cryosat-2            
     al - SARAL/AltiKa            
     cfo - CFOSAT            
     h2b - HaiYang-2B            
                
    cmems_L3_s6a:            
     s6a - Sentinel-6A Michael Freilich            
                
    eumetsat_L2:            
     s3a - Sentinel-3A            
     s3b - Sentinel-3B            
                
    cci_L2P:            
     j1 - Jason-1            
     j2 - Jason-2            
     j3 - Jason-3            
     c2 - Cryosat-2            
     envisat - Envisat            
     ers1 - European Remote-Sensing Satellite-1            
     ers2 - European Remote-Sensing Satellite-2            
     topex - TOPEX/Poseidon            
     al - SARAL/AltiKa            
     gfo - GEOSAT Follow-On            
        
    cci_L3:            
     multi - multimission product 1991-2018 

    cfo_swim_L2P:
     cfo - CFOSAT

This means that for product cmems_L3_NRT you can choose among 7 satellite missions. Unfortunatley, most of the satellite data is not accessible via thredds or similar options but needs to be downloaded from e.g. a VPN server. To do that you would need the credentials for Copernicus CMEMS, CEDA, or for the AVISO cataloque as these are the main sources that **wavy** currently exploits.

**wavy** relies on you to store the respective usernames and passwords in your local .netrc file. This could look like:

.. code::

   machine nrt.cmems-du.eu    login {USER}  password {PASSWORD}
   machine my.cmems-du.eu     login {USER}  password {PASSWORD}
   machine ftp.ceda.ac.uk    login {USER}  password {PASSWORD}
   machine ftp-access.aviso.altimetry.fr    login {USER}  password {PASSWORD}

Now, prepare your **wavy** environment with providing the directories for satellite data and model data. Add your path for satellite data here demonstrated for CMEMS with my user:

.. code-block:: yaml

   cmems_L3_NRT:
      dst:
         path_template: /home/patrikb/tmp_altimeter/L3/mission


There exists also something called strsub which defines strings that are o substituted. In this case some are predefined as:

.. code-block:: yaml

   strsub: ['varalias','mission','region']

The str "mission" in your path_template will be replaced by the satellite mission that you download. So for Sentinel-3a the final path for your downloaded files will be automatically /home/patrikb/tmp_altimeter/L3/s3a with subfolders on year and month.

You can proceed now and download CMEMS NRT L3 data using the wavyDownload.py script:

.. code-block:: bash

   $ cd ~/wavy/apps/standalone

To get help check ...

.. code-block:: bash

   $ ./wavyDownload.py -h

... or download some satellite altimeter data:

.. code-block:: bash

   $ ./wavyDownload.py -sat s3a -sd 2020110100 -ed 2020111000 -product cmems_L3_NRT

You can find the downloaded files in your chosen download directory.

Similarily one can download L2P and L3 multi-mission altimetry data from the CEDA Climate Change Initiative. This spans a long time period from 1991 to 2018 and enables climate related research and wave model hindcast validation.

.. code-block:: bash

   $ ./wavyDownload.py -sat multi -sd 2017112000 -ed 2017112100 -product cci_L3
   
You can also download altimeter data directly from python with the following lines. 

.. code-block:: bash

   >>> from wavy.satellite_module import satellite_class as sc
   >>> nID = 'cmems_L3_NRT'
   >>> name = 's3a'
   >>> sd = '2023-11-10 00'
   >>> ed = '2023-11-10 10'
   >>> # Initialize sc object
   >>> sco = sc(sd=sd,ed=ed,nID=nID,name=name)
   >>> # Download the data to a chosen directory
   >>> path = '/home/patrikb/tmp_altimeter/L3/s3a'
   >>> sco.download(nproc=4, path=path)


