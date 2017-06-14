# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 09:13:45 2017

@author: bling
"""
import sys
import datetime as dt
from matplotlib.path import Path
import netCDF4
from dateutil.parser import parse
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from datetime import datetime, timedelta
from math import radians, cos, sin, atan, sqrt  
import numpy as np
import sys
import datetime as dt
from matplotlib.path import Path
import netCDF4
from dateutil.parser import parse
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from datetime import datetime, timedelta
from math import radians, cos, sin, atan, sqrt  
from matplotlib.dates import date2num,num2date
def haversine(lon1, lat1, lon2, lat2): 
    """ 
    Calculate the great circle distance between two points  
    on the earth (specified in decimal degrees) 
    """   
    #print 'lon1, lat1, lon2, lat21',lon1, lat1, lon2, lat2
    #print 'lon1, lat1, lon2, lat22',lon1, lat1, lon2, lat2
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])  
    #print 34
    dlon = lon2 - lon1   
    dlat = lat2 - lat1   
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2  
    c = 2 * atan(sqrt(a)/sqrt(1-a))   
    r = 6371 
    d=c * r
    #print 'd',d
    return d
def calculate_SD(dmlon,dmlat):
    '''compare the model_points and drifter point(time same as model point)
    (only can pompare one day)!!!'''
    #print modelpoints,dmlon,dmlat,drtime
    #print len(dmlon)
    
    dd=0
    
    for a in range(len(dmlon)-1):
        #print 12
        #dla=(dmlat[a+1]-dmlat[a])*111
        #dlo=(dmlon[a+1]-dmlon[a])*(111*np.cos(dmlat[a]*np.pi/180))
        #d=sqrt(dla**2+dlo**2)#Calculate the distance between two points 
        #print model_points['lon'][a][j],model_points['lat'][a][j],dmlon[a][j],dmlat[a][j],d           
        #print 'd',d
        d=haversine(dmlon[a+1],dmlat[a+1],dmlon[a],dmlat[a])
        dd=dd+d
    #print 'dd',dd
    return dd
m26=np.load('m_ps2009-2007_830haha.npy')#'m_ps2011-2010_630.npy'
p=m26.tolist()
print len(p['lon'])
FN='necscoast_worldvec.dat'
'''
x=np.load('xx.npy')
m=np.load('mm.npy')
xx=[]
mm=[]
for a in np.arange(len(x)/2):
    xx.append(x[a*2])
    mm.append(m[a*2])
'''
CL=np.genfromtxt(FN,names=['lon','lat'])
fig,axes=plt.subplots(3,3,figsize=(15,15))#sharex=True,sharey=True,dpi=800,figsize=(15,15))
plt.subplots_adjust(wspace=0.1,hspace=0.1)
axes[0,0].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][2])):
    axes[0,0].scatter(p['lon'][2][a][0],p['lat'][2][a][0],color='green')
    if len(p['lon'][2][a])>=361:
        
        dis=calculate_SD(p['lon'][2][a][:360],p['lat'][2][a][:360])
        print 'dis',dis
        axes[0,0].scatter(p['lon'][2][a][360],p['lat'][2][a][360],color='red')
        axes[0,0].plot([p['lon'][2][a][0],p['lon'][2][a][360]],[p['lat'][2][a][0],p['lat'][2][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        #print 'j',j
        cl=plt.Circle((p['lon'][2][a][360],p['lat'][2][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[0,0].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][2][a][:],p['lat'][2][a][:])
        axes[0,0].scatter(p['lon'][2][a][-1],p['lat'][2][a][-1],color='red')
        axes[0,0].plot([p['lon'][2][a][0],p['lon'][2][a][-1]],[p['lat'][2][a][0],p['lat'][2][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][2][a][-1],p['lat'][2][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[0,0].add_patch(cl)

axes[0,0].scatter(p['lon'][2][a][0],p['lat'][2][a][0],label='start',color='green')
axes[0,0].scatter(p['lon'][2][a][-1],p['lat'][2][a][-1],label='end',color='red')
axes[0,0].plot(CL['lon'],CL['lat'],'b-',linewidth=0.5) 
axes[0,0].legend(loc='best') 
axes[0,0].set_xlabel('2007')
axes[0,0].set_xlim([-70,-65])
axes[0,0].set_ylim([42,47]) 
axes[0,0].xaxis.tick_top() 

axes[0,1].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][1])):
    axes[0,1].scatter(p['lon'][1][a][0],p['lat'][1][a][0],color='green')
    if len(p['lon'][1][a])>=361:
        dis=calculate_SD(p['lon'][1][a][:360],p['lat'][1][a][:360])
        axes[0,1].scatter(p['lon'][1][a][360],p['lat'][1][a][360],color='red')
        axes[0,1].plot([p['lon'][1][a][0],p['lon'][1][a][360]],[p['lat'][1][a][0],p['lat'][1][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][360],p['lat'][1][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[0,1].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][1][a][:],p['lat'][1][a][:])
        axes[0,1].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],color='red')
        axes[0,1].plot([p['lon'][1][a][0],p['lon'][1][a][-1]],[p['lat'][1][a][0],p['lat'][1][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][-1],p['lat'][1][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[0,1].add_patch(cl)
axes[0,1].scatter(p['lon'][1][a][0],p['lat'][1][a][0],label='start',color='green')
axes[0,1].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],label='end',color='red')
axes[0,1].plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
axes[0,1].set_xlabel('2008')
axes[0,1].xaxis.tick_top() 
axes[0,1].set_yticklabels([])
axes[0,1].set_xlim([-70,-65])
axes[0,1].set_ylim([42,47]) 
#plt.legend(loc='best') 


axes[0,2].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][0])):
    axes[0,2].scatter(p['lon'][0][a][0],p['lat'][0][a][0],color='green')
    if len(p['lon'][0][a])>=361:
        dis=calculate_SD(p['lon'][0][a][:360],p['lat'][0][a][:360])
        axes[0,2].scatter(p['lon'][0][a][360],p['lat'][0][a][360],color='red')
    #print len(lon[a])
        axes[0,2].plot([p['lon'][0][a][0],p['lon'][0][a][360]],[p['lat'][0][a][0],p['lat'][0][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][360],p['lat'][0][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[0,2].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][0][a][:],p['lat'][0][a][:])
        axes[0,2].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],color='red')
        axes[0,2].plot([p['lon'][0][a][0],p['lon'][0][a][-1]],[p['lat'][0][a][0],p['lat'][0][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][-1],p['lat'][0][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[0,2].add_patch(cl)
axes[0,2].scatter(p['lon'][0][a][0],p['lat'][0][a][0],label='start',color='green')
axes[0,2].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],label='end',color='red')  
axes[0,2].set_xlabel('2009')
axes[0,2].xaxis.tick_top() 

axes[0,2].set_yticklabels([])
axes[0,2].set_xlim([-70,-65])
axes[0,2].set_ylim([42,47]) 

m26=np.load('m_ps2011-2010_830haha.npy')#'m_ps2011-2010_630.npy'
p=m26.tolist()
axes[1,0].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][1])):
    axes[1,0].scatter(p['lon'][1][a][0],p['lat'][1][a][0],color='green')
    if len(p['lon'][1][a])>=361:
        dis=calculate_SD(p['lon'][1][a][:360],p['lat'][1][a][:360])
        axes[1,0].scatter(p['lon'][1][a][360],p['lat'][1][a][360],color='red')
        axes[1,0].plot([p['lon'][1][a][0],p['lon'][1][a][360]],[p['lat'][1][a][0],p['lat'][1][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][360],p['lat'][1][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[1,0].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][1][a][:],p['lat'][1][a][:])
        axes[1,0].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],color='red')
        axes[1,0].plot([p['lon'][1][a][0],p['lon'][1][a][-1]],[p['lat'][1][a][0],p['lat'][1][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][-1],p['lat'][1][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[1,0].add_patch(cl)
axes[1,0].scatter(p['lon'][1][a][0],p['lat'][1][a][0],label='start',color='green')
axes[1,0].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],label='end',color='red')
axes[1,0].plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
axes[1,0].set_xlabel('2010')
axes[1,0].set_xticklabels([])
#axes[0,1].set_xticklabels([])
axes[1,0].set_xlim([-70,-65])
axes[1,0].set_ylim([42,47]) 

axes[1,1].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][0])):
    axes[1,1].scatter(p['lon'][0][a][0],p['lat'][0][a][0],color='green')
    if len(p['lon'][0][a])>=361:
        dis=calculate_SD(p['lon'][0][a][:360],p['lat'][0][a][:360])
        axes[1,1].scatter(p['lon'][0][a][360],p['lat'][0][a][360],color='red')
    #print len(lon[a])
        axes[1,1].plot([p['lon'][0][a][0],p['lon'][0][a][360]],[p['lat'][0][a][0],p['lat'][0][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][360],p['lat'][0][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[1,1].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][0][a][:],p['lat'][0][a][:])
        axes[1,1].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],color='red')
        axes[1,1].plot([p['lon'][0][a][0],p['lon'][0][a][-1]],[p['lat'][0][a][0],p['lat'][0][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][-1],p['lat'][0][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[1,1].add_patch(cl)
axes[1,1].scatter(p['lon'][0][a][0],p['lat'][0][a][0],label='start',color='green')
axes[1,1].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],label='end',color='red')  
axes[1,1].set_xlabel('2011')
axes[1,1].set_yticklabels([])
axes[1,1].set_xticklabels([])
axes[1,1].set_xlim([-70,-65])
axes[1,1].set_ylim([42,47]) 



m26=np.load('m_ps2013-2012_830haha.npy')#'m_ps2011-2010_630.npy'
p=m26.tolist()
axes[1,2].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][1])):
    axes[1,2].scatter(p['lon'][1][a][0],p['lat'][1][a][0],color='green')
    if len(p['lon'][1][a])>=361:
        dis=calculate_SD(p['lon'][1][a][:360],p['lat'][1][a][:360])
        axes[1,2].scatter(p['lon'][1][a][360],p['lat'][1][a][360],color='red')
        axes[1,2].plot([p['lon'][1][a][0],p['lon'][1][a][360]],[p['lat'][1][a][0],p['lat'][1][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][360],p['lat'][1][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[1,2].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][1][a][:],p['lat'][1][a][:])
        axes[1,2].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],color='red')
        axes[1,2].plot([p['lon'][1][a][0],p['lon'][1][a][-1]],[p['lat'][1][a][0],p['lat'][1][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][-1],p['lat'][1][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[1,2].add_patch(cl)
axes[1,2].scatter(p['lon'][1][a][0],p['lat'][1][a][0],label='start',color='green')
axes[1,2].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],label='end',color='red')
axes[1,2].plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
axes[1,2].set_xlabel('2012')
axes[1,2].set_yticklabels([])
axes[1,2].set_xticklabels([])
axes[1,2].set_xlim([-70,-65])
axes[1,2].set_ylim([42,47]) 

axes[2,0].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][0])):
    axes[2,0].scatter(p['lon'][0][a][0],p['lat'][0][a][0],color='green')
    if len(p['lon'][0][a])>=361:
        dis=calculate_SD(p['lon'][0][a][:360],p['lat'][0][a][:360])
        axes[2,0].scatter(p['lon'][0][a][360],p['lat'][0][a][360],color='red')
    #print len(lon[a])
        axes[2,0].plot([p['lon'][0][a][0],p['lon'][0][a][360]],[p['lat'][0][a][0],p['lat'][0][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][360],p['lat'][0][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[2,0].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][0][a][:],p['lat'][0][a][:])
        axes[2,0].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],color='red')
        axes[2,0].plot([p['lon'][0][a][0],p['lon'][0][a][-1]],[p['lat'][0][a][0],p['lat'][0][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][-1],p['lat'][0][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[2,0].add_patch(cl)
axes[2,0].scatter(p['lon'][0][a][0],p['lat'][0][a][0],label='start',color='green')
axes[2,0].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],label='end',color='red')  
axes[2,0].set_xlabel('2013')
axes[2,0].set_xticklabels([])
axes[2,0].set_xlim([-70,-65])
axes[2,0].set_ylim([42,47]) 


m26=np.load('m_ps2015-2014_830haha.npy')#'m_ps2011-2010_630.npy'
p=m26.tolist()
axes[2,1].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][1])):
    axes[2,1].scatter(p['lon'][1][a][0],p['lat'][1][a][0],color='green')
    if len(p['lon'][1][a])>=361:
        dis=calculate_SD(p['lon'][1][a][:360],p['lat'][1][a][:360])
        axes[2,1].scatter(p['lon'][1][a][360],p['lat'][1][a][360],color='red')
        axes[2,1].plot([p['lon'][1][a][0],p['lon'][1][a][360]],[p['lat'][1][a][0],p['lat'][1][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][360],p['lat'][1][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[2,1].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][1][a][:],p['lat'][1][a][:])
        axes[2,1].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],color='red')
        axes[2,1].plot([p['lon'][1][a][0],p['lon'][1][a][-1]],[p['lat'][1][a][0],p['lat'][1][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][1][a][-1],p['lat'][1][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[2,1].add_patch(cl)
axes[2,1].scatter(p['lon'][1][a][0],p['lat'][1][a][0],label='start',color='green')
axes[2,1].scatter(p['lon'][1][a][-1],p['lat'][1][a][-1],label='end',color='red')
axes[2,1].plot(CL['lon'],CL['lat'],'b-',linewidth=0.5)
axes[2,1].set_xlabel('2014')
axes[2,1].set_yticklabels([])
axes[2,1].set_xticklabels([])
axes[2,1].set_xlim([-70,-65])
axes[2,1].set_ylim([42,47]) 

axes[2,2].plot(CL['lon'],CL['lat'],'b-')
for a in np.arange(len(p['lon'][0])):
    axes[2,2].scatter(p['lon'][0][a][0],p['lat'][0][a][0],color='green')
    if len(p['lon'][0][a])>=361:
        dis=calculate_SD(p['lon'][0][a][:360],p['lat'][0][a][:360])
        #axes[2,2].scatter(p['lon'][0][a][360],p['lat'][0][a][360],color='red')
    #print len(lon[a])
        axes[2,2].plot([p['lon'][0][a][0],p['lon'][0][a][360]],[p['lat'][0][a][0],p['lat'][0][a][360]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][360],p['lat'][0][a][360]),70*0.009009009,color='red',alpha=0.2)
        axes[2,2].add_patch(cl)
    else:
        dis=calculate_SD(p['lon'][0][a][:],p['lat'][0][a][:])
        #axes[2,2].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],color='red')
        axes[2,2].plot([p['lon'][0][a][0],p['lon'][0][a][-1]],[p['lat'][0][a][0],p['lat'][0][a][-1]],'y-')#,linewidth=0.5)
        #for j in np.arange(len(xx)):
        cl=plt.Circle((p['lon'][0][a][-1],p['lat'][0][a][-1]),70*0.009009009,color='red',alpha=0.2)
        axes[2,2].add_patch(cl)
axes[2,2].scatter(p['lon'][0][a][0],p['lat'][0][a][0],label='start',color='green')
axes[2,2].scatter(p['lon'][0][a][-1],p['lat'][0][a][-1],label='end',color='red')  
axes[2,2].set_xlabel('2015')
axes[2,2].set_yticklabels([])
axes[2,2].set_xticklabels([])
axes[2,2].set_xlim([-70,-65])
axes[2,2].set_ylim([42,47]) 
plt.savefig('''8-1 to 8-151p''',dpi=200)
