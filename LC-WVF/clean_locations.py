import numpy as np
import pandas as pd

df=pd.read_pickle("./location_pickle.pkl")

lats = df['latitude']
lons = df['longitude']
statid = df['station:id']
idx = df['idx']
idy = df['idy']

Data = {'latitude':lats,
        'longitudes':lons,
        'station:id':statid,
        'idx':idx,
        'idy':idy}
    
newdf = pd.DataFrame (Data, columns = ['station:id','longitudes','latitude','idx','idy'])
newdf = newdf.drop_duplicates()
newdf = newdf.dropna()

newdf.to_pickle('locations_mwam4_clean.pkl')
