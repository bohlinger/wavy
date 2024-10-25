import requests
import xarray as xr
import pandas as pd
import os
import json
from wavy.credentials import credentials_from_netrc
from wavy.utils import parse_date


def get_AIS_data(bbox, sd=None, ed=None, minspeed=0.5):
    '''
    Args:
        bbox: list of strings [lon min, lat min,
                               lon max, lat max]
        sd, ed: start and end dates, as strings
                under "yyyymmddHHMM" format
    '''
    sd = parse_date(sd)
    ed = parse_date(ed)
    sd_str = sd.strftime('%Y%m%d%H%M')
    ed_str = ed.strftime('%Y%m%d%H%M')

    # get token
    auth = credentials_from_netrc('kystdatahuset.no')
    url = "https://kystdatahuset.no/ws/api/auth/login"
    data = {"username": auth[0], "password": auth[1]}
    headers = {'accept': '*/*', 'Content-Type': 'application/json'}
    response = requests.post(url, json=data, headers=headers, auth=auth)
    resp_json = response.json()
    token = resp_json['data']['JWT']

    # prepare curl request
    curl_com = "curl -X 'POST' "
    url_req = "'https://kystdatahuset.no/ws/api/ais/positions/within-bbox-time' "
    head_1 = "-H 'accept: text/plain' "
    head_2 = "-H 'Authorization: Bearer " + token + "' "
    head_3 = "-H 'Content-Type: application/json' "

    dict_data = {
                 "bbox": ",".join(bbox),
                 "start": sd_str,
                 "end": ed_str,
                 "minSpeed": minspeed
                }

    data = "-d '" + str(dict_data).replace("'", '"') + "'"

    curl_request = curl_com + url_req + head_1 + head_2 + head_3 + data

    # get result content of curl request
    stream = os.popen(curl_request)
    output = stream.read()
    dict_output = json.loads(output)

    df_request = pd.DataFrame(
        data=dict_output["data"],
        columns=["mmsi", "datetime", "lon",
                 "lat", "x1", "x2", "x3", "x4",
                 "x5", "x6"])

    # Get ship info
    unique_mmsi = [int(mm) for mm in df_request["mmsi"].unique()]
    url_ships_nsr = 'https://kystdatahuset.no/ws/api/ship/data/nsr/for-mmsis-imos'
    json_ships_nsr = {"mmsi": unique_mmsi, "imo": [], "callsign": []}
    headers_ships_nsr = {'accept': 'test/plain',
                         'Content-Type': 'application/json'}
    response_ships_nsr = requests.post(url_ships_nsr,
                                       json=json_ships_nsr,
                                       headers=headers_ships_nsr)
    resp_ships_json_nsr = response_ships_nsr.json()
    ships_df = pd.DataFrame(resp_ships_json_nsr['data'])
    ships_df = ships_df[['mmsino', 'grosstonnage', 'dwt', 'length',
                         'breadth', 'draught', 'shiptypeeng']]

    # Merge the two dataframes
    df_merge = pd.merge(df_request, ships_df, left_on="mmsi",
                        right_on="mmsino", how="left")
    df_merge = df_merge[['mmsi', 'datetime', 'lon', 'lat',
                         'grosstonnage', 'dwt', 'length', 'breadth',
                         'draught', 'shiptypeeng']]
    df_merge = df_merge.rename(columns={'lon': 'lons', 'lat': 'lats',
                                        'datetime': 'time'})

    # Turn the dataframe into a dataset
    df_merge['time'] = pd.to_datetime(df_merge['time'])
    ds = df_merge.set_index(['time'])
    ds = ds.to_xarray().sortby("time")

    return ds
