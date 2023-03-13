#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn
import json
import csv
from datetime import date, datetime, timedelta

import os
import requests
from bs4 import BeautifulSoup


# # Coordinates conversion API

# In[2]:


#https://geodesy.geo.admin.ch/reframe/navref?format=json&easting=2649930&northing=1177380&altitude=NaN&input=lv95&output=etrf93-ed
#https://geodesy.geo.admin.ch/reframe/navref?format=json&easting=649930&northing=177380&altitude=NaN&input=lv03&output=etrf93-ed

def swiss_coordinates_to_wgs84(easting, northing, system="MN95"):
    """Transforms swiss MN03 or MN95 coordinates to standard GPS WGS84,
    using https://geodesy.geo.admin.ch/reframe/navref API.
    (interface: https://www.swisstopo.admin.ch/en/maps-data-online/calculation-services/navref.html).
    Coordinates (MN03 and MN95) satisfy: easting > northing, everywhere in Switzerland."""
    E=str(int(easting))
    N=str(int(northing))
    service_url = "https://geodesy.geo.admin.ch/reframe/navref"
    
    if system=="MN95":
        if(len(E)!=7 or len(N)!=7):
            raise(ValueError("Length of MN95 coordinates must be 7 each."))
        resp = requests.get(service_url+"?format=json&easting="+E+"&northing="+N+"&altitude=NaN&input=lv95&output=etrf93-ed").json()
    
    elif system=="MN95-short":
        if(len(E)!=6 or len(N)!=6):
            raise(ValueError("Length of 'MN95-short' coordinates must be 6 each."))
        resp = requests.get(service_url+"?format=json&easting=2"+E+"&northing=1"+N+"&altitude=NaN&input=lv95&output=etrf93-ed").json()
    
    elif system=="MN03":
        if(len(E)!=6 or len(N)!=6):
            raise(ValueError("Length of MN03 coordinates must be 6 each."))
        resp = requests.get(service_url+"?format=json&easting="+E+"&northing="+N+"&altitude=NaN&input=lv03&output=etrf93-ed").json()
    
    else:
        raise(ValueError("Accepted systems: 'MN03' 'MN95' 'MN95-short'"))
    
    result = {"longitude" : float(resp["easting"]),
             "latitude" : float(resp["northing"])}
    return(result)


# In[3]:


resp=swiss_coordinates_to_wgs84(2649930,1177380)
resp


# # Webscraping stations infos

# In[4]:


#req = requests.get("https://www.hydrodaten.admin.ch/en/current-situation-table-discharge-and-water-levels.html")
#soup = BeautifulSoup(req.text,'html.parser')
#stations_list = soup.find("table", class_="table table-bordered").find("tbody").find_all("tr")


# In[5]:


req = requests.get("https://www.hydrodaten.admin.ch/en/2019.html")
req.encoding = "utf-8"
soup = BeautifulSoup(req.text,'html.parser')
Stations_list_soup = soup.find("div", class_="well search").find("select").find_all("option")


# In[6]:


Stations_names = pd.DataFrame({"BAFU_ID":[int(station_soup["value"]) for station_soup in Stations_list_soup],
             "Name":[station_soup.string for station_soup in Stations_list_soup]})
Stations_names


# In[7]:


def webscrape_station(station_id):
    id_str = str(station_id)
    while(len(id_str)<4):
        id_str = "0"+id_str
    try:
        req = requests.get("https://www.hydrodaten.admin.ch/en/"+id_str+".html")
        req.encoding = "utf-8"
        soup = BeautifulSoup(req.text,'html.parser')
        station_info_soup = soup.find_all(class_="col-md-6")[-1].table.tbody.find_all("tr")
        infos = [info_.td.string for info_ in station_info_soup]
        
        gps03 = swiss_coordinates_to_wgs84(easting=int(infos[4][0:6]),northing=int(infos[4][9:15]),
                                         system="MN03")
        gps95 = swiss_coordinates_to_wgs84(easting=int(infos[4][0:6]),northing=int(infos[4][9:15]),
                                         system="MN95-short")
        data=pd.Series({"Elevation":float(infos[0][:-8]),
                        "Surface_catchment":float(infos[1][:-3]),
                        "Mean_elevation_catchment":float(infos[2][:-8]),
                        "Glaciation_percent":float(infos[3][:-1]),
                        "easting":int(infos[4][0:6]),
                        "northing":int(infos[4][9:15]),
                        "longitude03":gps03["longitude"],
                        "latitude03":gps03["latitude"],
                        "longitude95":gps95["longitude"],
                        "latitude95":gps95["latitude"]})
    except:
        data=pd.Series({"Elevation":np.NaN,
                        "Surface_catchment":np.NaN,
                        "Mean_elevation_catchment":np.NaN,
                        "Glaciation_percent":np.NaN,
                        "easting":np.NaN,
                        "northing":np.NaN,
                        "longitude03":np.NaN,
                        "latitude03":np.NaN,
                        "longitude95":np.NaN,
                        "latitude95":np.NaN})
    return(data)


# In[8]:


Stations_infos = Stations_names.BAFU_ID.apply(webscrape_station)
Stations_BAFU = Stations_names.join(Stations_infos)
Stations_BAFU.to_csv("data_wrangled/Webscraped_BAFU_stations_info.csv",index=False)
Stations_BAFU.info()


# In[9]:


Stations_BAFU = pd.read_csv("data_wrangled/Webscraped_BAFU_stations_info.csv")
Stations_BAFU


# In[10]:


Stations_BAFU[Stations_BAFU.easting.notnull()].BAFU_ID.values


# In[11]:


Stations_BAFU[Stations_BAFU.easting.isna()].BAFU_ID.values


# # Merge test

# In[12]:


Stations_ours = pd.read_csv("data_wrangled/info_discharges_gps.csv", index_col="id")
Stations_ours


# In[13]:


Stations_BAFU.BAFU_ID.values


# In[14]:


Stations_ours.BAFU_ID.values


# In[15]:


Stations_ours.reset_index(inplace=True)
Discharges_info = Stations_ours.merge(Stations_BAFU, how="left", left_on="BAFU_ID", right_on="BAFU_ID",
                                      suffixes=('', '_BAFU'))
#CSV save to verify all the matching info
#Discharges_info.to_csv("data_wrangled/Discharges_info_merge_debug.csv",index=False)
Discharges_info.set_index("id")

