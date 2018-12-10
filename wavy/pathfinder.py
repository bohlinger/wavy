"""
python file including needed paths
the integrated wave verification system should be able to run
including only the paths listed here
"""

stationpath_lustre_om = '/lustre/storeA/project/fou/om/' + \
                        'waveverification/data/'
stationpath_lustre_hi = '/lustre/storeA/project/fou/hi/' + \
                        'waveverification/data/'

satpath_lustre = '/lustre/storeA/project/fou/om/altimeter/'
satpath_copernicus = '/lustre/storeB/project/copernicus/ARC-MFC/sentinel3a/'
satpath_ftp_008_052 = ('/Core/'
                    + 'SEALEVEL_GLO_WAV_L3_NRT_OBSERVATIONS_008_052/'
                    + 'dataset-wav-alti-l3-rt-global-s3a/'
                    ) # from 20170709
satpath_ftp_014_001 = ('/Core/'
                    + 'WAVE_GLO_WAV_L3_SWH_NRT_OBSERVATIONS_014_001/'
                    + 'dataset-wav-alti-l3-swh-rt-global-s3a/'
                    ) # from 20180320
satpath_sentinel3b = ('...') # when it is available
station_d22_starc = ('/lustre/storeB/immutable/archive/projects/'\
                   + 'metproduction/DNMI_OFFSHORE/')
                   # +rig+'/d22/'+YY+'/'+dy+'.d22'
station_d22_opdate = ('/vol/gorgon/offshore/')#+rig+'/d22/'+dy+'.d22'
