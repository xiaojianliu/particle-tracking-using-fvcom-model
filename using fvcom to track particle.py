# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 15:50:45 2017

@author: bling
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:34:20 2016

@author: hxu
"""

import datetime as dt
import pytz
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from pytz import timezone
import numpy as np
import csv
from scipy import  interpolate
from matplotlib.dates import date2num,num2date
from barycentric_polygonal_interpolation import get_drifter_track,get_fvcom,get_roms,calculate_SD,drifterhr
######## Hard codes ##########
Model='30yr' # JiM changed from array to string
wind_get_type='FVCOM'
wind=0 
st_lat=[]
st_lon=[]
latc=np.linspace(44.5,45,10)
lonc=np.linspace(-66.8,-66,10)
for aa in np.arange(len(lonc)):
    for bb in np.arange(len(latc)):
        st_lat.append(latc[bb])
        st_lon.append(lonc[aa])
'''
# alternative method of creating random set of points
num = 70
st_lat = np.random.uniform(44.5,45,num)[:]
st_lon = np.random.uniform(-66.8,-66,num)[:]
'''
days=15 # number of days to track
jia=14*0 # always zero
end_times=[]
start_time=[]
start_times=[dt.datetime(2015,5,1,0,0,0,0),dt.datetime(2014,5,1,0,0,0,0)]#,dt.datetime(2013,5,1,0,0,0,0),dt.datetime(2012,5,1,0,0,0,0),dt.datetime(2011,5,1,0,0,0,0),dt.datetime(2010,5,1,0,0,0,0),dt.datetime(2009,5,1,0,0,0,0),dt.datetime(2008,5,1,0,0,0,0),dt.datetime(2007,5,1,0,0,0,0)]
###############################################  END of HARDCODES  #########################################

for a in np.arange(len(start_times)):
    start_time.append(start_times[a]+timedelta(hours=jia*24))
m_ps =dict(lon=[],lat=[],time=[])
for a in np.arange(len(start_time)):
    print a
    end_time=start_time[a]+timedelta(hours=days*24)
    
    #i=Model[0] #why do this???
    GRIDS= ['GOM3','massbay','30yr']
    if Model in GRIDS:
        try:
            get_obj =  get_fvcom(i)
            url_fvcom = get_obj.get_url(start_time[a],end_time)                
            b_points = get_obj.get_data(url_fvcom) # boundary points?
            #print '##########################'
            #print drifter_points['lon_hr'][(nday-1)*24]
            #print drifter_points['lat_hr'][(nday-1)*24]
            #print start_time
            model_points =dict(lon=[],lat=[],time=[]) # all point in the track of one start location
            for b in np.arange(len(st_lon)):
                print b
                modelpoints = dict(lon=[],lat=[],time=[]) # dictionary holding  one point for one of the 100 points
                modelpoints,windspeed= get_obj.get_track(st_lon[b],st_lat[b],0,start_time[a],wind,wind_get_type)
                model_points['lon'].append(modelpoints['lon'])
                model_points['lat'].append(modelpoints['lat'])  
                model_points['time'].append(modelpoints['time'])    
        except:
            print 'There is no model-point near the given-point'
            continue
    m_ps['lon'].append(model_points['lon'])
    m_ps['lat'].append(model_points['lat'])
    m_ps['time'].append(model_points['time'])
np.save('m_ps'+str(startime[0].year)+'-'+str(startime[1].year),m_ps)# this 
