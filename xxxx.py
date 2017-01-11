# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:22:54 2016

@author: hxu
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:34:20 2016

@author: xiaojian
"""

import datetime as dt
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
from matplotlib.path import Path
from Back_forecast_function import get_fvcom,sh_bindata1,nearest_point1
######## Hard codes ##########
Model='30yr'
point_data=[]
num=1#how many point do you want forecast
end_times=dt.datetime(2010,5,1)
w=0.0001 #m/s
wind_get_type='FVCOM'#'NCEP'
wind=0
data_type='hourly'
FNCL='necscoast_worldvec.dat'
add_model_w='yes'
data = np.genfromtxt('sea.csv',dtype=None,names=['x','y','h'],delimiter=',')    
x=[]
y=[]
h=[]
x=data['x']
y=data['y']
h=data['h']
xi = np.arange(-78.,-55.00,0.2)
yi = np.arange(34,46,0.2)
xb,yb,hb_mean,hb_median,hb_std,hb_num = sh_bindata1(x, y, h, xi, yi)
hb_median=abs(hb_mean)
xxxb,yyyb = np.meshgrid(xb, yb)
plt.figure()
CS=plt.contour(xxxb, yyyb, -abs(hb_mean.T),levels=[-200,-100])        
plt.clabel(CS, inline=1, fontsize=8,fmt='%4.0f')           
plt.axis([-78,-55,34,46],linewidth=0.5)
CL=np.genfromtxt(FNCL,names=['lon','lat'])
plt.plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
tongji=0
for a in np.arange(312):#313
    print 'a',a 
    point_id=a*10
    try: 
        if data_type=='hourly':
            
            if point_data==[]:
                if '30yr' in Model:
                    point_data='nes_lon_lat.csv'
            
            data = np.genfromtxt(point_data,dtype=None,names=['local','lon','lat'],delimiter=',',skip_header=1) 
            #plt.scatter(data['lon'][0:300],data['lat'][0:300])
            #num=len(data['lon'])
            
            lonp,latp = nearest_point1(data['lon'][point_id], data['lat'][point_id], x, y,0.8)
            for vv in np.arange(len(x)):
                if x[vv]==lonp and y[vv]==latp:
                    end_deep=-h[vv]+0
            if end_deep>=0:
                continue
            start_times =end_times-timedelta(hours=int(np.abs(end_deep)/float(w*3600)))
            
            plt.title('end_time=%s-%s-%s  w=%sm/s'%(end_times.year,end_times.month,end_times.day,w))
            if a==3:
                plt.scatter(data['lon'][point_id],data['lat'][point_id],marker='o',color='red',label='end bottom points',s=10)
            lonn=[]
            latt=[]
            for i in range(num):
                back=dict(lon=[],lat=[],time=[],deep=[]) 
                get_obj =  get_fvcom(Model)
                url_fvcom = get_obj.get_url(start_times,end_times)                
                b_points = get_obj.get_data(url_fvcom)
                try:
                    if add_model_w=='yes': 
                        back,windspeed= get_obj.get_track(data['lon'][point_id+i],data['lat'][point_id+i],end_deep,start_times,wind,wind_get_type,w)
                except:
                    continue
                try:
                    if add_model_w!='yes': 
                        back,windspeed= get_obj.get_track1(data['lon'][point_id+i],data['lat'][point_id+i],end_deep,start_times,wind,wind_get_type,w)
                except:
                    continue
                lonn.append(data['lon'][point_id+i])
                latt.append(data['lat'][point_id+i])
                plt.plot([data['lon'][point_id+i],back['lon'][-1]],[data['lat'][point_id+i],back['lat'][-1]],'r-',linewidth=0.2)
                plt.scatter(data['lon'][point_id+i],data['lat'][point_id+i],marker='o',color='red',s=10)
                plt.scatter(back['lon'][-1],back['lat'][-1],marker='o',color='green',s=10)
                lonn.append(back['lon'][-1])
                latt.append(back['lat'][-1])
                tongji=tongji+1
                '''
                plt.figure()
                plt.title('%s'%(data['local'][point_id+i]))
                plt.plot(back['lon'],back['lat'],'bo-')
                plt.scatter(data['lon'][point_id+i],data['lat'][point_id+i],marker='o',color='red',s=70)
                '''
            if a==3:
                plt.scatter(back['lon'][-1],back['lat'][-1],marker='o',color='green',s=10,label='start surface points')
                plt.legend(loc='best')
            
    except:
        continue
plt.savefig('back forecast333',dpi=400)
plt.show()
print tongji

