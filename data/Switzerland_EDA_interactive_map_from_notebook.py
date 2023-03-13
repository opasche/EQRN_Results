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
import folium
import branca.colormap as cm


# In[2]:


def show_nan(df, size=(16,20), y_ticks=None, show_axis=False, only_date=False):
    plt.figure(figsize=size)
    plt.imshow(df.isnull().values, aspect='auto', interpolation="none")#"nearest""antialiased"
    plt.xticks(range(df.shape[1]),df.columns.tolist(),rotation=90)
    if y_ticks:
        if only_date:
            plt.yticks(y_ticks,df.index[y_ticks].date.tolist())
        else:
            plt.yticks(y_ticks,df.index[y_ticks].tolist())
    if show_axis:
        ax = plt.gca()
        ax.set_xticks(np.arange(0.5, df.shape[1],1),minor=True)
        plt.grid(which="minor",axis="x",color="w",linewidth=0.5)
    plt.show()


# # Import

# ## River Discharges

# In[3]:


Stations_info = pd.read_csv("data_wrangled/info_discharges_initial.csv", encoding="latin", index_col=0)
Stations_info = Stations_info.rename(columns={"Name.y":"Name","ID":"BAFU_ID"})
Stations_info.index.name = "id"

Discharges = pd.read_csv("data_wrangled/Discharges_all.csv")
Discharges["Date"]=pd.to_datetime(Discharges.Date)
Discharges.set_index("Date",inplace=True)
Discharges.sort_index(inplace=True)

Summer = pd.read_csv("data_wrangled/Discharges_summer.csv")
Summer["Date"]=pd.to_datetime(Summer.Date)
Summer.set_index("Date",inplace=True)
Summer.sort_index(inplace=True)


# In[4]:


Stations_info


# In[5]:


Discharges


# In[6]:


Summer


# ## Precipitations

# In[7]:


Precip = pd.read_csv("Data/precip.csv")
Precip.rename(columns={"Unnamed: 0":"Date"},inplace=True)
Precip["Date"]=pd.to_datetime(Precip.Date)
Precip.set_index("Date",inplace=True)
Precip.sort_index(inplace=True)

PrecipInfo = pd.read_csv("Data/precip_info.csv",encoding="latin")
PrecipInfo.set_index("code",inplace=True)
PrecipInfo.drop(["Unnamed: 0","folder"],axis=1,inplace=True)


# In[8]:


Precip


# In[9]:


PrecipInfo


# In[10]:


Precip21_1 = pd.read_table("Data/Precipitations_112021/order_98339_data.txt",sep=";")
Precip21_2 = pd.read_table("Data/Precipitations_112021/order_98340_data.txt",sep=";")
Precip21_3 = pd.read_table("Data/Precipitations_112021/order_98341_data.txt",sep=";")


# In[11]:


Precip21 = pd.concat([Precip21_1,Precip21_2,Precip21_3])
Precip21["time"]=pd.to_datetime(Precip21.time, format='%Y%m%d', errors="coerce")
Precip21.dropna(inplace=True)
Precip21.reset_index(drop=True, inplace=True)
#Precip21.set_index(["stn","time"],inplace=True)
#Precip21.sort_index(inplace=True)
Precip21


# ## Stats

# In[12]:


Stations_stats = Discharges.describe().transpose()
Stations_stats


# In[13]:


Stations_stats[["count","min"]].min()


# In[14]:


Precip_stats = Precip.describe().transpose()
Precip_stats


# In[15]:


Precip_stats[["count","min"]].min()


# In[16]:


Precip.min().sort_values()


# # Dates Range and Missing Values

# In[17]:


newyears_disch = Discharges.reset_index()[Discharges.reset_index().Date.dt.dayofyear==1].index.tolist()
newyears_summ = Summer.reset_index()[(Summer.reset_index().Date.dt.month==6)&(Summer.reset_index().Date.dt.day==1)].index.tolist()
newyears_precip = Precip.reset_index()[Precip.reset_index().Date.dt.dayofyear==1].index.tolist()


# In[18]:


show_nan(Discharges, y_ticks=newyears_disch, only_date=True)


# In[19]:


show_nan(Summer, y_ticks=newyears_summ, only_date=True)


# In[20]:


show_nan(Precip, y_ticks=newyears_precip, only_date=True)


# # Coordinates conversion

# In[21]:


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


# In[22]:


#mn95 = Stations_info[["X","Y"]].apply(lambda x: pd.Series(
#    swiss_coordinates_to_wgs84(x["X"], x["Y"], system="MN95-short")), axis=1)
#mn03 = Stations_info[["X","Y"]].apply(lambda x: pd.Series(
#    swiss_coordinates_to_wgs84(x["X"], x["Y"], system="MN03")), axis=1)


# In[23]:


#mn95 = mn95.add_suffix("95")
#mn03 = mn03.add_suffix("03")
#gps_coord = mn95.join(mn03)
#Stations_info = Stations_info.join(gps_coord)
#Stations_info.to_csv("data_wrangled/info_discharges_gps.csv")


# In[24]:


Stations_info = pd.read_csv("data_wrangled/info_discharges_gps.csv", index_col="id")


# In[25]:


Stations_info


# # Merging BAFU info

# In[26]:


Stations_BAFU = pd.read_csv("data_wrangled/Webscraped_BAFU_stations_info.csv")
Stations_BAFU


# In[27]:


Discharges_info = Stations_info.reset_index().merge(Stations_BAFU, how="left",
                                                    left_on="BAFU_ID", right_on="BAFU_ID",
                                                    suffixes=('', '_BAFU'))
Discharges_info.set_index("id", inplace=True)
#CSV save to verify all the matching info
Discharges_info.convert_dtypes().to_csv("data_wrangled/Discharges_info_merge_debug.csv")
#To have the right dtypes in base formats
Discharges_info = pd.read_csv("data_wrangled/Discharges_info_merge_debug.csv",index_col=0)
Discharges_info


# In[28]:


Discharges_info.isnull().sum()


# In[29]:


(Stations_BAFU.Mean_elevation_catchment%1).value_counts(dropna=False)


# ### Divergences in information

# In[30]:


Discharges_info[(Discharges_info.X!=Discharges_info.easting)]


# In[31]:


Discharges_info[(Discharges_info.Y!=Discharges_info.northing)]


# In[32]:


Stations_info_keep = Discharges_info[['Name','BAFU_ID','Ave','Elevation','Surface_catchment',
                                      'Mean_elevation_catchment','Glaciation_percent',
                                      'X','Y','longitude95','latitude95','longitude03','latitude03']]
Stations_info_keep.to_csv("data_wrangled/Discharges_info_merged.csv")


# # Interactive Map

# In[33]:


#https://ogre.adc4gis.com/convert
#http://geojson.io/#map=12/46.2122/6.1111
map_ch=folium.Map(location=[46.8,8.2],zoom_start=8,tiles='cartodbpositron')#'Stamen Toner'
#folium.TileLayer('openstreetmap').add_to(map_ge)
#geo_communes = json.load(open(r'SITG/CAD_COMMUNE.geo.json'))#.geo.json

DischLayer = folium.FeatureGroup(name="Discharges", overlay=True, control=True, show=True)
for disch_stati in Stations_info_keep.iterrows():
    folium.Marker(
        location=[disch_stati[1]["latitude95"], disch_stati[1]["longitude95"]],
        icon=folium.Icon(color='darkblue', icon_color='white', icon='tint',
                         angle=0, prefix='fa'),#prefix='glyphicon'
        #tooltip="hey!",
        popup=disch_stati[1][["BAFU_ID","Name","Ave",'Elevation','Surface_catchment',
                              'Mean_elevation_catchment','Glaciation_percent',]].to_frame().to_html(
            classes='table table-striped table-hover table-condensed table-responsive')
    ).add_to(DischLayer)
map_ch.add_child(DischLayer)


PrecipLayer = folium.FeatureGroup(name="Precipitations", overlay=True, control=True, show=True)
for precip_stati in PrecipInfo.iterrows():
    folium.Marker(
        location=[precip_stati[1]["latitude"], precip_stati[1]["longitude"]],
        icon=folium.Icon(color='lightblue', icon_color='white', icon='cloud',
                         angle=0, prefix='fa'),#prefix='glyphicon'
        #tooltip="hey!",
        popup=precip_stati[1][["name","Elevation","maxMean","totMean",
                               "startYear","endYear","pctNA"]].to_frame().to_html(
            classes='table table-striped table-hover table-condensed table-responsive')
    ).add_to(PrecipLayer)
map_ch.add_child(PrecipLayer)

#['red', 'blue', 'green', 'purple', 'orange', 'darkred', 'lightred', 'beige', 'darkblue', 'darkgreen',
#'cadetblue', 'darkpurple', 'white', 'pink', 'lightblue', 'lightgreen', 'gray', 'black', 'lightgray']

map_ch.add_child(folium.LatLngPopup())
folium.LayerControl(collapsed=False).add_to(map_ch)
map_ch.save(os.path.join('Interactive_Swiss_map', 'Switzerland_Discharges_and_Precipitations_map.html'))
map_ch


# # Time Series of Discharges

# In[34]:


Discharges.sort_index(inplace=True)
for sta in Discharges.columns:
    Discharges[sta].dropna().plot(style="b.", figsize=(16,6), title=sta)#"bo"
    plt.show()


# In[35]:


Summer.sort_index(inplace=True)
for sta in Summer.columns:
    Summer[sta].dropna().plot(style="b.", figsize=(16,6), title=sta)#"bo"
    plt.show()


# # Time Series of Precipitations

# In[36]:


Precip.min().sort_values()


# In[37]:


Precip["GSG"].dropna().plot(style="b.", figsize=(16,6), title=sta)#"bo"
plt.show()#negative min station


# In[38]:


Precip.sort_index(inplace=True)
for sta in Precip.columns:
    Precip[sta].dropna().plot(style="b.", figsize=(16,6), title=sta)#"bo"
    plt.show()

