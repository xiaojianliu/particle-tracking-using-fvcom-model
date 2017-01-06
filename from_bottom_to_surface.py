# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 16:34:20 2016
@author: xiaojian
Routine to backtrack particles that have landed on the bottom and determine where they came from.
Involves prescribing (ie harcoding) the "descent rate".  
Note: Uses a set of functions in "Back_forecast_function.py". 
"""

import datetime as dt
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import csv
import math
from matplotlib.path import Path
from Back_forecast_function import get_fvcom
def sh_bindata(x, y, z, xbins, ybins):
    """
    Bin irregularly spaced data on a rectangular grid.
    """
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_median=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_std=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_num=np.zeros((len(xbins)-1,len(ybins)-1),dtype=int)    
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
#            k=np.where((ix==iix) and (iy==iiy)) # wrong syntax
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
            zb_median[iix-1,iiy-1]=np.median(z[k])
            zb_std[iix-1,iiy-1]=np.std(z[k])
            zb_num[iix-1,iiy-1]=len(z[k])
            
    return xb,yb,zb_mean,zb_median,zb_std,zb_num
def nearest_point( lon, lat, lons, lats, length):  #0.3/5==0.06
        '''Find the nearest point to (lon,lat) from (lons,lats),
           return the nearest-point (lon,lat)
           author: Bingwei'''
        p = Path.circle((lon,lat),radius=length)
        #numpy.vstack(tup):Stack arrays in sequence vertically
        points = np.vstack((lons.flatten(),lats.flatten())).T  
        
        insidep = []
        #collect the points included in Path.
        for i in xrange(len(points)):
            if p.contains_point(points[i]):# .contains_point return 0 or 1
                insidep.append(points[i])  
        # if insidep is null, there is no point in the path.
        if not insidep:
            print 'There is no model-point near the given-point.'
            raise Exception()
        #calculate the distance of every points in insidep to (lon,lat)
        distancelist = []
        for i in insidep:
            ss=math.sqrt((lon-i[0])**2+(lat-i[1])**2)
            distancelist.append(ss)
        # find index of the min-distance
        mindex = np.argmin(distancelist)
        # location the point
        lonp = insidep[mindex][0]; latp = insidep[mindex][1]
        
        return lonp,latp
######## Hard codes ##########
plt.figure()
tongji=0
for a in np.arange(312):#313
    print 'a',a
    try:
        Model='30yr'
        point_data=[]
        num=1#how many point do you want forecast
        point_id=a*10 # here we are subsampling every tenth particle
        #if have drifter_ID it will be use how many drifter have how many data 
        end_times=dt.datetime(2010,5,1)
        w=0.00005 #descent rate in m/s
        wind_get_type='FVCOM'#'NCEP'
        wind=0#if you don't want have wind correct,wind =0 
        data_type='hourly'#hourly or month
        FNCL='necscoast_worldvec.dat' #coastline data
        data = np.genfromtxt('sea.csv',dtype=None,names=['x','y','h'],delimiter=',')    
        x=[]
        y=[]
        h=[]
        x=data['x']
        y=data['y']
        h=data['h']
        xi = np.arange(-78.,-55.00,0.2)
        yi = np.arange(34,46,0.2)
        xb,yb,hb_mean,hb_median,hb_std,hb_num = sh_bindata(x, y, h, xi, yi)
        hb_median=abs(hb_mean)
        xxxb,yyyb = np.meshgrid(xb, yb)
        CS=plt.contour(xxxb, yyyb, -abs(hb_mean.T),levels=[-200,-100])
        
        plt.clabel(CS, inline=1, fontsize=8,fmt='%4.0f')
            
        plt.axis([-78,-55,34,46],linewidth=0.5)
        
        ########back forecast ########## 
        CL=np.genfromtxt(FNCL,names=['lon','lat'])
        plt.plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
        if data_type=='hourly':
            
            if point_data==[]:
                if '30yr' in Model:
                    point_data='nes_lon_lat.csv'
            
            data = np.genfromtxt(point_data,dtype=None,names=['local','lon','lat'],delimiter=',',skip_header=1) 
            #plt.scatter(data['lon'][0:300],data['lat'][0:300])
            #num=len(data['lon'])
            
            lonp,latp = nearest_point(data['lon'][point_id], data['lat'][point_id], x, y,0.8)
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
                    
                    back,windspeed= get_obj.get_track(data['lon'][point_id+i],data['lat'][point_id+i],end_deep,start_times,wind,wind_get_type)
                
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
plt.savefig('back forecast30',dpi=400)
plt.show()
print tongji

