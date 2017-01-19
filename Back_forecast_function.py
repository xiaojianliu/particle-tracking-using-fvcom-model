# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:15:03 2015
the function about how to get and compare the drifter and model data
@author: qianran
"""
import sys
import datetime as dt
from matplotlib.path import Path
import netCDF4
from dateutil.parser import parse
import numpy as np
import math
import pandas as pd
from datetime import datetime, timedelta
from math import radians, cos, sin, atan, sqrt  
from matplotlib.dates import date2num,num2date
######## function ##########
def getrawdrift(did,filename):
   '''
   routine to get raw drifter data from ascii files posted on the web
   '''
   url='http://nefsc.noaa.gov/drifter/'+filename
   df=pd.read_csv(url,header=None, delimiter="\s+")
   # make a datetime
   dtime=[]
   index = np.where(df[0]==int(did))[0]
   newData = df.ix[index]
   for k in newData[0].index:
      #dt1=dt.datetime(int(filename[-10:-6]),df[2][k],df[3][k],df[4][k],df[5][k],0,0,pytz.utc)
      dt1=datetime(2015, newData[2][k],newData[3][k],newData[4][k],newData[5][k],0,0)
      dtime.append(dt1)
   return newData[8],newData[7],dtime,newData[9]
def getdrift(did):
    """
    routine to get drifter data from archive based on drifter id (did)
    -assumes "import pandas as pd" has been issued above
    -get remotely-stored drifter data via ERDDAP
    -input: deployment id ("did") number where "did" is a string
    -output: time(datetime), lat (decimal degrees), lon (decimal degrees), depth (meters)
    
    note: there is another function below called "data_extracted" that does a similar thing returning a dictionary
    
    Jim Manning June 2014
    """
    url = 'http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?time,latitude,longitude,depth&id="'+did+'"&orderBy("time")'
    df=pd.read_csv(url,skiprows=[1]) #returns a dataframe with all that requested
    # generate this datetime 
    for k in range(len(df)):
       df.time[k]=parse(df.time[k]) # note this "parse" routine magically converts ERDDAP time to Python datetime
    return df.latitude.values,df.longitude.values,df.time.values,df.depth.values 
def __cmptime(time, times):
    '''
    return indies of specific or nearest time in times.
    '''
    tdelta = []
    #print len(times)
    for t in times:
        tdelta.append(abs((time-t).total_seconds()))
        
    index = tdelta.index(min(tdelta))
    
    return index
def sh_bindata1(x, y, z, xbins, ybins):
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
def nearest_point1( lon, lat, lons, lats, length):  #0.3/5==0.06
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
def get_drifter_track(method_of_drifter,start_time, days,drifter_ID):  
    dr_points=dict(lon=[],lat=[],time=[]) 
    drpoints=dict(ids=[],lat_hr=[],lon_hr=[],lon=[],lat=[],time=[],distance=[])
    drifter_points = dict(lon=[],lat=[],time=[])
    #print start_time
    if method_of_drifter=='raw':
        ids=drifter_ID
        dr_points = get_drifter_raw(start_time,days,drifter_ID)
        drifter_points['time'].extend(dr_points['time'])
    if method_of_drifter=='csv':
        ids=drifter_ID
        starttime=start_time.strftime("%Y-%m-%d")
        dr_points = get_drifter_csv(starttime,drifter_ID,days)
        drifter_points['time'].extend(dr_points['time'])
    if method_of_drifter=='erddap':
        dr_point=dict(lon=[],lat=[],time=[]) 
        drtime=[]
        id=drifter_ID
        ids=id
        dr_point['time'],dr_point['lat'],dr_point['lon'] =get_drifter_erddap(id,start_time,days+1)
        #print dr_points['time']        
        for w in range(len(dr_point['time'])):       
            times=[]       
            times=dr_point['time'][w].replace(tzinfo=None)
            time=times-timedelta(hours=4)
            #print times
            drtime.append(time)  
        drtimenp = np.array(drtime)
        #print drtimenp

        dst = drtimenp-start_time; dstindex = np.argmin(abs(dst))
        det = drtimenp-(start_time+timedelta(days)); detindex = np.argmin(abs(det))#every compare days drifter end index
        #print dstindex,detindex
        dr_points['lon']=dr_point['lon'][dstindex:detindex];dr_points['lat']=dr_point['lat'][dstindex:detindex]
        drifter_points['time']=drtime[dstindex:detindex]        
        #print drifter_points['time'],dr_points
        dr_points['lon'].tolist;dr_points['lat'].tolist
    if method_of_drifter=='npy':
        ids=drifter_ID
        dr_points = get_drifter_npy(start_time,days,drifter_ID)
        #print dr_points
        drifter_points['time'].extend(dr_points['time'])
    drpoints=dr_points ;drpoints['ids']=ids;drpoints['time']=drifter_points['time']
    #print drpoints
    return drpoints
def get_drifter_npy(starttime, days,drifterID):  
    '''return drifter_points lon lat time
    get from vitalii step 2'''
    drpoints=dict(lon=[],lat=[],time=[])
    source='driftfvcom_data2/'
    filename='ID_%s.npz' %drifterID
    #print filename
    drfile=source+filename
    Z=np.load(drfile)
    sttime=date2num(starttime)
    endtime=date2num(starttime+timedelta(days))
    times=np.array(Z['tdh'])
    #print times,sttime,endtime
    #tdh days since 0001-00-00 00:00:00, GMT 
    lons=np.array(Z['londh'])
    lats=np.array(Z['latdh'])
    sindex = np.argmin(abs(sttime-times))
    eindex = np.argmin(abs(endtime-times))
    #print times[sindex:eindex]
    for i in times[sindex:eindex]:
        a=num2date(i)
        b=a.replace(tzinfo=None)
        drpoints['time'].append(b)
    #print sindex,eindex
    drpoints['lon']=lons[sindex:eindex]
    drpoints['lat']=lats[sindex:eindex]
    #print drpoints
    return drpoints
def get_drifter_raw(starttime, days,drifterID):
    '''
        return drifter nodes
        if starttime is given, return nodes started from starttime
        if both starttime and days are given, return nodes of the specific time period
    '''
    filename = 'drift_X.dat' 
    if filename:
        temp=getrawdrift(drifterID,filename)
    else:
        temp=getdrift(drifterID)  
    nodes = {}
    nodes['lon'] = np.array(temp[1])
    nodes['lat'] = np.array(temp[0])
    nodes['time'] = np.array(temp[2])
        #starttime = np.array(temp[2][0])
    if not starttime:
        starttime = np.array(temp[2][0])
    if days:
        endtime = starttime + timedelta(days=days)
        i = __cmptime(starttime, nodes['time'])
        j = __cmptime(endtime, nodes['time'])
        nodes['lon'] = nodes['lon'][i:j+1]
        nodes['lat'] = nodes['lat'][i:j+1]
        nodes['time'] = nodes['time'][i:j+1]
    else:
        i = __cmptime(starttime, nodes['time'])
        nodes['lon'] = nodes['lon'][i:-1]
        nodes['lat'] = nodes['lat'][i:-1]
        nodes['time'] = nodes['time'][i:-1]
    return nodes
def get_drifter_erddap(id,start_time,days):
    """
     get data from url, return ids latitude,longitude, times
     input_time can either contain two values: start_time & end_time OR one value:interval_days
     and they should be timezone aware
     example: input_time=[dt(2012,1,1,0,0,0,0,pytz.UTC),dt(2012,2,1,0,0,0,0,pytz.UTC)]
     """
    df=dict(id=[],lon=[],lat=[],time=[])
    mintime=start_time.strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')  # change time format
    endtime=start_time+timedelta(days)    
    maxtime=endtime.strftime('%Y-%m-%d'+'T'+'%H:%M:%S'+'Z')    
    # open url to get data
    #url='http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?id,time,latitude,longitude&id=%22100390731%22&time>='+str(mintime)+'&time<='+str(maxtime)+'&orderBy(%22time%22)'
    url='http://comet.nefsc.noaa.gov:8080/erddap/tabledap/drifters.csv?id,time,latitude,longitude&time>='\
    +str(mintime)+'&time<='+str(maxtime)+'&id="'+str(id)+'"&orderBy("time")'
    df=pd.read_csv(url,skiprows=[1])
    for k in range(len(df)):
        #print df.time[k]
        df.time[k]=parse(df.time[k][:-1])
    df=df[df.longitude <=-20]
    #print df.time.values
    return df.time.values,df.latitude.values,df.longitude.values
    

def get_drifter_csv(start_time,drifter_ID,days):
    dt_starttime = datetime.strptime(start_time, "%Y-%m-%d")
    nodes=dict(lon=[],lat=[],time=[])
    FN='ID_%s.csv' %drifter_ID 
    D = np.genfromtxt(FN,dtype=None,names=['ID','TimeRD','TIME_GMT','YRDAY0_GMT','LON_DD','LAT_DD','TEMP','DEPTH_I'],delimiter=',')    

    nodes['lon'] = D['LON_DD']
    nodes['lat'] =  D['LAT_DD']
    for i in range(len(D['YRDAY0_GMT'])):
        a=[]
        a=dt.datetime(2010,01,01,0,0,0,0)+timedelta(D['YRDAY0_GMT'][i])
        nodes['time'].append(a)
    #print nodes['time']
        #starttime = np.array(temp[2][0])
    if days:
        #print start_time
        endtime = dt_starttime+timedelta(days)
        i = __cmptime(dt_starttime, nodes['time'])
        j = __cmptime(endtime, nodes['time'])
        nodes['lon'] = nodes['lon'][i:j+1]
        nodes['lat'] = nodes['lat'][i:j+1]
        nodes['time'] = nodes['time'][i:j+1]
    else:
        print 'I need the days'
    return nodes
    
def drifterhr(dr_points,days):
    '''get drifter data on the hour'''
    #print dr_points
    drifter_points_hr = dict(ids=[],lat_hr=[],h_hr=[],lon_hr=[],lon=[],lat=[],time=[],distance=[])
    cprtime=[];
    #print dr_points['time']
    cprstime=dr_points['time'][0].replace(minute=0)+timedelta(hours=1)
    cpretime=cprstime+timedelta(days=days)
    for i in range((cpretime-cprstime).days*24):
        cprtime.append(cprstime+timedelta(hours=i))
    drifter_points_hr['h_hr'] = np.array(cprtime)
    '''npdrlons = np.array(dr_points['lon'])
    npdrlats = np.array(dr_points['lat'])'''
    npdrtimes = np.array(dr_points['time'])
    for i in drifter_points_hr['h_hr']:
        td=npdrtimes-i
        index = np.argmin(abs(td))
        #print npdrtimes[index],i
        if npdrtimes[index]>i:
            erindex=index-1
            laindex=index
        else:
            erindex=index
            laindex=index+1

        if (npdrtimes[laindex]-npdrtimes[erindex]).seconds==0:
            #print i,npdrtimes[erindex],npdrtimes[laindex+1],(npdrtimes[laindex+1]-npdrtimes[erindex]).seconds
            hrlon=dr_points['lon'][erindex]+float((i-npdrtimes[erindex]).seconds)/(npdrtimes[laindex+1]-npdrtimes[erindex]).seconds*(dr_points['lon'][laindex+1]-dr_points['lon'][erindex])
            hrlat=dr_points['lat'][erindex]+float((i-npdrtimes[erindex]).seconds)/(npdrtimes[laindex+1]-npdrtimes[erindex]).seconds*(dr_points['lat'][laindex+1]-dr_points['lat'][erindex])
        else:
            hrlon=dr_points['lon'][erindex]+float((i-npdrtimes[erindex]).seconds)/(npdrtimes[laindex]-npdrtimes[erindex]).seconds*(dr_points['lon'][laindex]-dr_points['lon'][erindex])
            hrlat=dr_points['lat'][erindex]+float((i-npdrtimes[erindex]).seconds)/(npdrtimes[laindex]-npdrtimes[erindex]).seconds*(dr_points['lat'][laindex]-dr_points['lat'][erindex])
        drifter_points_hr['lon_hr'].append(hrlon);
        drifter_points_hr['lat_hr'].append(hrlat)
        #print dr_points['lon'][erindex],dr_points['lon'][laindex],hrlon


    drifter_points_hr['lat']=dr_points['lat']
    drifter_points_hr['lon']=dr_points['lon']
    drifter_points_hr['time']=dr_points['time']
    drifter_points_hr['ids']=dr_points['ids']

    
    '''npdrtime=npdrtimes[index]
    npdrlon=npdrlons[index]    
    npdrlat=npdrlats[index]
    npdrdellon.append(npdrlon)
    npdrdellat.append(npdrlat)
    npdrdeltime.append(npdrtime)'''

    #print drifter_points_hr,len(drifter_points_hr['time']),len(drifter_points_hr['lon']),len(drifter_points_hr['lon_hr'])
    return drifter_points_hr
    
#######################model track###
def dm2dd(lat,lon):
    """
    convert lat, lon from decimal degrees,minutes to decimal degrees
    """
    (a,b)=divmod(float(lat),100.)   
    aa=int(a)
    bb=float(b)
    lat_value=aa+bb/60.

    if float(lon)<0:
        (c,d)=divmod(abs(float(lon)),100.)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60.)
        lon_value=-lon_value
    else:
        (c,d)=divmod(float(lon),100.)
        cc=int(c)
        dd=float(d)
        lon_value=cc+(dd/60.)
    return lat_value, -lon_value
def get_nc_data(url, *args):
    '''
    get specific dataset from url

    *args: dataset name, composed by strings
    ----------------------------------------
    example:
        url = 'http://www.nefsc.noaa.gov/drifter/drift_tcs_2013_1.dat'
        data = get_url_data(url, 'u', 'v')
    '''
    nc = netCDF4.Dataset(url)
    data = {}
    for arg in args:
        try:
            data[arg] = nc.variables[arg]
        except (IndexError, NameError, KeyError):
            print 'Dataset {0} is not found'.format(arg)
    #print data
    return data
class get_fvcom():
    def __init__(self, mod):
        self.modelname = mod
            
    def nearest_point(self, lon, lat, lons, lats, length):  #0.3/5==0.06
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
        
        
    def get_url(self, starttime, endtime):
        '''
        get different url according to starttime and endtime.
        urls are monthly.
        '''
        self.hours = int(round((endtime-starttime).total_seconds()/60/60))
        
                
        if self.modelname == "GOM3":
            timeurl = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?Times[0:1:144]'''
            url = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?
            lon[0:1:51215],lat[0:1:51215],lonc[0:1:95721],latc[0:1:95721],siglay[0:1:39][0:1:51215],h[0:1:51215],nbe[0:1:2][0:1:95721],
            u[{0}:1:{1}][0:1:39][0:1:95721],v[{0}:1:{1}][0:1:39][0:1:95721],zeta[{0}:1:{1}][0:1:51215]'''
            '''urll = http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc?
            u[{0}:1:{1}][0:1:39][0:1:95721],v[{0}:1:{1}][0:1:39][0:1:95721],zeta[{0}:1:{1}][0:1:51215]'''
            try:
                mTime = netCDF4.Dataset(timeurl).variables['Times'][:]              
            except:
                print '"GOM3" database is unavailable!'
                raise Exception()
            Times = []
            for i in mTime:
                strt = '201'+i[3]+'-'+i[5]+i[6]+'-'+i[8]+i[9]+' '+i[11]+i[12]+':'+i[14]+i[15]
                Times.append(datetime.strptime(strt,'%Y-%m-%d %H:%M'))
            fmodtime = Times[0]; emodtime = Times[-1]         
            if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
                print 'Time: Error! Model(GOM3) only works between %s with %s(UTC).'%(fmodtime,emodtime)
                raise Exception()
            npTimes = np.array(Times)
            tm1 = npTimes-starttime; #tm2 = mtime-t2
            index1 = np.argmin(abs(tm1))
            index2 = index1 + self.hours#'''
            url = url.format(index1, index2)
            self.mTime = Times[index1:index2+1]
            
            self.url = url
            
        elif self.modelname == "massbay":
            timeurl = '''http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?Times[0:1:144]'''
            url = """http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?
            lon[0:1:98431],lat[0:1:98431],lonc[0:1:165094],latc[0:1:165094],siglay[0:1:9][0:1:98431],h[0:1:98431],
            nbe[0:1:2][0:1:165094],u[{0}:1:{1}][0:1:9][0:1:165094],v[{0}:1:{1}][0:1:9][0:1:165094],zeta[{0}:1:{1}][0:1:98431]"""
            
            try:
                mTime = netCDF4.Dataset(timeurl).variables['Times'][:]              
            except:
                print '"massbay" database is unavailable!'
                raise Exception()
            Times = []
            for i in mTime:
                strt = '201'+i[3]+'-'+i[5]+i[6]+'-'+i[8]+i[9]+' '+i[11]+i[12]+':'+i[14]+i[15]
                Times.append(datetime.strptime(strt,'%Y-%m-%d %H:%M'))
            fmodtime = Times[0]; emodtime = Times[-1]         
            if starttime<fmodtime or starttime>emodtime or endtime<fmodtime or endtime>emodtime:
                print 'Time: Error! Model(massbay) only works between %s with %s(UTC).'%(fmodtime,emodtime)
                raise Exception()
            npTimes = np.array(Times)
            tm1 = npTimes-starttime; #tm2 = mtime-t2
            index1 = np.argmin(abs(tm1))
            index2 = index1 + self.hours#'''
            url = url.format(index1, index2)
            self.mTime = Times[index1:index2+1]
            
            self.url = url

        elif self.modelname == "30yr": #start at 1977/12/31 23:00, end at 2014/1/1 0:0, time units:hours
            timeurl = """http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?time[0:1:316008]"""
            url = '''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?h[0:1:48450],
            lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],nbe[0:1:2][0:1:90414],siglay[0:1:44][0:1:48450],
            u[{0}:1:{1}][0:1:44][0:1:90414],v[{0}:1:{1}][0:1:44][0:1:90414],zeta[{0}:1:{1}][0:1:48450],
            ww[{0}:1:{1}][0:1:44][0:1:48450]'''
            
            try:
                mtime = netCDF4.Dataset(timeurl).variables['time'][:]
                
            except:
                print '"30yr" database is unavailable!'
                raise Exception
            # get model's time horizon(UTC).
            '''fmodtime = datetime(1858,11,17) + timedelta(float(mtime[0]))
            emodtime = datetime(1858,11,17) + timedelta(float(mtime[-1]))
            mstt = fmodtime.strftime('%m/%d/%Y %H:%M')
            mett = emodtime.strftime('%m/%d/%Y %H:%M') #'''
            # get number of days from 11/17/1858
            #print starttime
            t1 = (starttime - datetime(1858,11,17)).total_seconds()/86400 
            t2 = (endtime - datetime(1858,11,17)).total_seconds()/86400
            if not mtime[0]<t1<mtime[-1] or not mtime[0]<t2<mtime[-1]:
                #print 'Time: Error! Model(massbay) only works between %s with %s(UTC).'%(mstt,mett)
                print 'Time: Error! Model(30yr) only works between 1978-1-1 with 2014-1-1(UTC).'
                raise Exception()
            
            tm1 = mtime-t1; #tm2 = mtime-t2
            #print mtime,tm1
            index1 = np.argmin(abs(tm1)); #index2 = np.argmin(abs(tm2)); print index1,index2
            index2 = index1 + self.hours
            url = url.format(index1, index2)
            Times = []
            for i in range(self.hours+1):
                Times.append(starttime+timedelta(hours=i))
            #print 'Times', Times
            self.mTime = Times
            self.url = url
        
        #print url
        return url
    def get_data(self,url):
        '''
        "get_data" not only returns boundary points but defines global attributes to the object
        '''
        self.data = get_nc_data(url,'lat','lon','latc','lonc','siglay','h','nbe','u','v','zeta','ww')#,'nv'
        self.lonc, self.latc = self.data['lonc'][:], self.data['latc'][:]  #quantity:165095
        self.lons, self.lats = self.data['lon'][:], self.data['lat'][:]
        self.h = self.data['h'][:]; self.siglay = self.data['siglay'][:]; #nv = self.data['nv'][:]
        self.u = self.data['u']; self.v = self.data['v']; self.zeta = self.data['zeta']
        self.ww = self.data['ww']
        nbe1=self.data['nbe'][0];nbe2=self.data['nbe'][1];
        nbe3=self.data['nbe'][2]
        pointt = np.vstack((nbe1,nbe2,nbe3)).T; self.pointt = pointt
        wl=[]
        for i in pointt:
            if 0 in i: 
                wl.append(1)
            else:
                wl.append(0)
        self.wl = wl
        tf = np.array(wl)
        inde = np.where(tf==True)
        #b_index = inde[0]
        lonb = self.lonc[inde]; latb = self.latc[inde]        
        self.b_points = np.vstack((lonb,latb)).T#'''
        #self.b_points = b_points
        return self.b_points #,nv lons,lats,lonc,latc,,h,siglay
        
    def shrink_data(self,lon,lat,lons,lats,rad):
        lont = []; latt = []
        p = Path.circle((lon,lat),radius=rad)
        pints = np.vstack((lons,lats)).T
        for i in range(len(pints)):
            if p.contains_point(pints[i]):
                lont.append(pints[i][0])
                latt.append(pints[i][1])
        lonl=np.array(lont); latl=np.array(latt)#'''
        if not lont:
            print 'point position error! shrink_data'
            #sys.exit()
        return lonl,latl
    
    def boundary_path(self,lon,lat):
        p = Path.circle((lon,lat),radius=0.03)
        dis = []
        for i in self.b_points:
            if p.contains_point(i):
                d = math.sqrt((lon-i[0])**2+(lat-i[1])**2)
                dis.append(d)
        if dis:
            md = min(dis)
            pa = Path.circle((lon,lat),radius=md+0.005)
            return pa
        else: return None
        
    
    def eline_path(self,lon,lat):
        '''
        When drifter close to boundary(less than 0.1),find one nearest point to drifter from boundary points, 
        then find two nearest boundary points to previous boundary point, create a boundary path using that 
        three boundary points.
        '''
        def boundary_location(locindex,pointt,wl):
            '''
            Return the index of boundary points nearest to 'locindex'.
            '''
            loca = []
            dx = pointt[locindex]; #print 'func',dx 
            for i in dx: # i is a number.
                #print i  
                if i ==0 :
                    continue
                dx1 = pointt[i-1]; #print dx1
                if 0 in dx1:
                    loca.append(i-1)
                else:
                    for j in dx1:
                        if j != locindex+1:
                            if wl[j-1] == 1:
                                loca.append(j-1)
            return loca
        
        p = Path.circle((lon,lat),radius=0.02) #0.06
        dis = []; bps = []; pa = []
        tlons = []; tlats = []; loca = []
        for i in self.b_points:
            if p.contains_point(i):
                bps.append((i[0],i[1]))
                d = math.sqrt((lon-i[0])**2+(lat-i[1])**2)
                dis.append(d)
        bps = np.array(bps)
        if not dis:
            return None
        else:
            print "Close to boundary."
            dnp = np.array(dis)
            dmin = np.argmin(dnp)
            lonp = bps[dmin][0]; latp = bps[dmin][1]
            index1 = np.where(self.lonc==lonp)
            index2 = np.where(self.latc==latp)
            elementindex = np.intersect1d(index1,index2)[0] # location 753'''
            #print index1,index2,elementindex  
            loc1 = boundary_location(elementindex,self.pointt,self.wl) ; #print 'loc1',loc1
            loca.extend(loc1)
            loca.insert(1,elementindex)               
            for i in range(len(loc1)):
                loc2 = boundary_location(loc1[i],self.pointt,self.wl); #print 'loc2',loc2
                if len(loc2)==1:
                    continue
                for j in loc2:
                    if j != elementindex:
                        if i ==0:
                            loca.insert(0,j)
                        else:
                            loca.append(j)
            
            for i in loca:
                tlons.append(self.lonc[i]); tlats.append(self.latc[i])
                       
            for i in xrange(len(tlons)):
                pa.append((tlons[i],tlats[i]))
            path = Path(pa)#,codes
            return path
            

        
    def uvt(self,u1,v1,u2,v2):
        t = 2
        a=0; b=0
        if u1==u2:
            a = u1
        else:
            ut = np.arange(u1,u2,float(u2-u1)/t)
            for i in ut:
                a += i
            a = a/t  
        
        if v1==v2:
            b = v1
        else:
            c = float(v2-v1)/t
            vt = np.arange(v1,v2,c)
            for i in vt:
                b += i
            b = b/t
               
        return a, b


    def nearpoint_index(self,lat,lon):
        dists=[];points=[]
        url_wind="""http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/meteorology_v1/m_forcing_2011.nc?XLAT[0:1:140][0:1:182],XLONG[0:1:140][0:1:182]"""  #%(year,ftime,stime,lats,lats,lons,lons)
        vdata=get_nc_data(url_wind,'XLAT','XLONG')
        LAT=vdata['XLAT'][:];LON=vdata['XLONG'][:]
        for a in range(140):#have 141 group data
            #print max(LAT[a]), min(LAT[a]) ,max(LON[a]) ,min(LON[a])
            if lat<max(LAT[a]) and lat>min(LAT[a]) and lon<max(LON[a]) and lon>min(LON[a]):
                #print a
                for o in range(182):#every group have 182 data
                    if lon-0.1<LON[a][o]<lon+0.1 and lat-0.1<LAT[a][o]<lat+0.1:
                        #print a,o,LON[a][o],LAT[a][o]
                        points.append((a,o))
                        dist=haversine(LON[a][o],LAT[a][o],lon,lat)
                        dists.append(dist)
        #print dists,points
        if len(dists)>1:
            for i in range(len(dists)):
                if dists[i]-min(dists)==0:
                    index=i
        else :
            index=0
        #print  points[i]
        return points[index]
            
    def get_wind_fvcom(self,starttime,lat,lon):
        #print starttime
        year=starttime.year
        index=self.nearpoint_index(lat,lon)

        cptime="%i,01,01,00,00"  %year
        cptimes=datetime.strptime(cptime, '%Y,%m,%d,%H,%M')
        time_s=(starttime-cptimes).total_seconds()
        timeindex=int(time_s/60/60)


        url_wind="""http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/meteorology_v1/m_forcing_%s.nc?U10[%i:1:%i][%i:1:%i][%i:1:%i],V10[%i:1:%i][%i:1:%i][%i:1:%i]"""  %(year,timeindex,timeindex,index[0],index[0],index[1],index[1],timeindex,timeindex,index[0],index[0],index[1],index[1])
        data=get_nc_data(url_wind,'V10','U10')
        self.v_wind=data['V10'][:];self.u_wind=data['U10'][:]
        #print self.v_wind,self.u_wind
       
        #print vwind
        #print self.v_wind[0][0][0],self.u_wind[0][0][0]
        return self.v_wind,self.u_wind
        
    def get_wind_ncep(self,starttime,lat,lon):
        year=starttime.year
        #print starttime
        lats=(lat-18)/0.1875
        if lats-int(lats)>0.5:
            lats=int(lats)+1
        else:
            lats=int(lats)
        lons=(lon+91.9375)/0.1875 ;
        if lons-int(lons)>0.5:
            lons=int(lons)+1
        else:
            lons=int(lons)
        cptime="%i,01,01,00,00"  %year

        cptimes=datetime.strptime(cptime, '%Y,%m,%d,%H,%M')
        #print cptimes
        time_d=(starttime-cptimes).days
        time_s=(starttime-cptimes).seconds
        
        stime=time_d*8+int(time_s/60/60/3)+1
        ftime=time_d*8+int(time_s/60/60/3)
        

        url_vwind="""http://tds.marine.rutgers.edu:8080/thredds/dodsC/met/narr/Vwind_narr_NENA_%s.nc?Vwind[%i:1:%i][%i:1:%i][%i:1:%i]"""  %(year,ftime,stime,lats,lats,lons,lons)
        vdata=get_nc_data(url_vwind,'Vwind')
        self.v_wind=vdata['Vwind'][:]
        
        url_uwind="""http://tds.marine.rutgers.edu:8080/thredds/dodsC/met/narr/Uwind_narr_NENA_%s.nc?Uwind[%i:1:%i][%i:1:%i][%i:1:%i]"""  %(year,ftime,stime,lats,lats,lons,lons)    
        udata=get_nc_data(url_uwind,'Uwind')
        self.u_wind=udata['Uwind'][:]
        #print (time_s/60.0/60.0/3.0),int(time_s/60/60/3),self.v_wind[0][0][0],self.v_wind[1][0][0]
        vwind=self.v_wind[0][0][0]+((time_s/60.0/60.0/3.0)-int(time_s/60/60/3))*(self.v_wind[1][0][0]-self.v_wind[0][0][0])
        uwind=self.u_wind[0][0][0]+((time_s/60.0/60.0/3.0)-int(time_s/60/60/3))*(self.u_wind[1][0][0]-self.u_wind[0][0][0])
        #print vwind
        #print self.v_wind[0][0][0],self.u_wind[0][0][0]
        return vwind,uwind
        
    def get_track(self,lon,lat,depth,starttime,wind,wind_get_type,w1): #,b_index,nvdepth, 
        '''
        Get forecast points start at lon,lat
        '''
        windspeed=[]
        modpts = dict(lon=[lon], lat=[lat], layer=[], time=[],deep=[depth]) #model forecast points
        #uvz = netCDF4.Dataset(self.url)
        #u = uvz.variables['u']; v = uvz.variables['v']; zeta = uvz.variables['zeta']
        #print 'len u',len(u)
        #print modpts['lon']
        if lon>90:
            lon, lat = dm2dd(lon, lat)
        
        lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc,1)#1 day elements(blue)
        lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats,1)#1 day nodes(red)
        print 'in'
        try:
            if self.modelname == "GOM3" or self.modelname == "30yr":
                lonp,latp = self.nearest_point(lon, lat, lonl, latl,1)#elements(blue)
                
                lonn,latn = self.nearest_point(lon,lat,lonk,latk,1.2)#nodes(red)
                
            if self.modelname == "massbay":
                lonp,latp = self.nearest_point(lon, lat, lonl, latl,0.06)
                lonn,latn = self.nearest_point(lon,lat,lonk,latk,0.1)    
            index1 = np.where(self.lonc==lonp)
            index2 = np.where(self.latc==latp)
            elementindex = np.intersect1d(index1,index2)#nearest index elements(blue)
            index3 = np.where(self.lons==lonn)
            index4 = np.where(self.lats==latn)
            nodeindex = np.intersect1d(index3,index4)#nearest index nodes(red)
            
            
            pa = self.eline_path(lon,lat)# boundary 

            if self.modelname == "30yr":
                waterdepth = self.h[nodeindex]
            else:
                waterdepth = self.h[nodeindex]+self.zeta[0,nodeindex]
            
            modpts['time'].append(self.mTime[0])
            depth_total = self.siglay[:,nodeindex]*waterdepth; #print 'Here one' 
            for xx in np.arange(600):
                    if waterdepth<(abs(depth)):
                        depth=depth+5
                    if waterdepth<(abs(depth)):
                        continue
                    else:
                        break
            layer = np.argmin(abs(depth_total+depth)); #print 'layer',layer
            modpts['layer'].append(layer); 
            
            
            if waterdepth<(abs(depth)):
                print 'This point is too shallow.Less than %d meter.'%abs(depth)
                raise Exception()
        except:
            #print 12345
            return modpts
        sum_deep=0    
        t = abs(self.hours)
        print 'hour',t         
        for ii in xrange(t):  
            i=t-ii
            
                #print 'layer,lon,lat,i',layer,lon,lat,i
            lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc,1)
            lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats,1)
            u_t1 = self.u[i,layer,elementindex][0]; v_t1 = self.v[i,layer,elementindex][0]
            u_t2 = self.u[i-1,layer,elementindex][0]; v_t2 = self.v[i-1,layer,elementindex][0]
            u_t = -(u_t1+u_t2)/2; v_t = -(v_t1+v_t2)/2
            w_t1=self.ww[i,layer,nodeindex][0]
            w_t2=self.ww[i,layer,nodeindex][0]
            #u_t,v_t = self.uvt(u_t1,v_t1,u_t2,v_t2)
            w=-(w_t1+w_t2)/2
            starttimes=starttime+timedelta(hours=i)
            #print 'starttime',starttimes
            if wind==0:

                dx = 60*60*u_t; dy = 60*60*v_t;dz=60*60*(w+w1)
            else:
                if wind_get_type=='NCEP':
                    v_wind,u_wind=self.get_wind_ncep(starttimes,lat,lon)
                    meanwindspeed=0
                if wind_get_type=='FVCOM':
                    v_wind,u_wind=self.get_wind_fvcom(starttimes,lat,lon)

                dx = 60*60*u_t+60*60*u_wind[0][0][0]*wind; dy = 60*60*v_t+60*60*v_wind[0][0][0]*wind
                windspeed.append(u_wind[0][0][0]*u_wind[0][0][0]+v_wind[0][0][0]*v_wind[0][0][0])
                
                #windspeed.append(v_wind[0][0][0]*v_wind[0][0][0]+u_wind[0][0][0]*u_wind[0][0][0])
            #u_t = (u_t1+u_t2)/2; v_t = (v_t1+v_t2)/2
            '''if u_t==0 and v_t==0: #There is no water
                print 'Sorry, hits the land,u,v==0'
                return modpts,1 #'''
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #temlon,temlat = mapx(x+dx,y+dy,inverse=True)            
            
            temlon = lon + (dx/(111111*np.cos(lat*np.pi/180)))#convert decimal to degree
            temlat = lat + dy/111111 #'''
            
            sum_deep=sum_deep+dz
            #print '%d,lat,lon,layer'%(i+1),temlat,temlon,layer
            #########case for boundary 1 #############
            if pa:
                teml = [(lon,lat),(temlon,temlat)]
                tempa = Path(teml)
                if pa.intersects_path(tempa): 
                    print 'Sorry, point hits land here.path'
                    if len(windspeed)<2:
                        print 12
                        meanwindspeed=windspeed
                    else:
                        print 34
                        meanwindspeed=np.mean(windspeed)                 
                    return modpts,meanwindspeed


            lon = temlon; lat = temlat
            #if i!=(t-1):                
            try:
                #print 12345
                if self.modelname == "GOM3" or self.modelname == "30yr":
                    lonp,latp = self.nearest_point(lon, lat, lonl, latl,1)
                    lonn,latn = self.nearest_point(lon,lat,lonk,latk,1.2)
                    
                if self.modelname == "massbay":
                    lonp,latp = self.nearest_point(lon, lat, lonl, latl,0.06)
                    lonn,latn = self.nearest_point(lon,lat,lonk,latk,0.1)
                index1 = np.where(self.lonc==lonp)
                index2 = np.where(self.latc==latp)
                elementindex = np.intersect1d(index1,index2); #print 'elementindex',elementindex
                index3 = np.where(self.lons==lonn)
                index4 = np.where(self.lats==latn)
                nodeindex = np.intersect1d(index3,index4)
                
            
                pa = self.eline_path(lon,lat)#boundary
                if self.modelname == "30yr":            
                    waterdepth = self.h[nodeindex]
                else:
                    waterdepth = self.h[nodeindex]+self.zeta[i+1,nodeindex]
                modpts['time'].append(self.mTime[i-1])            
                
                depth_total = self.siglay[:,nodeindex]*waterdepth 
                for xx in np.arange(600):
                    if waterdepth<(abs(depth)):
                        depth=depth+5
                    if waterdepth<(abs(depth)):
                        continue
                    else:
                        break
                layer = np.argmin(abs(depth_total+depth+sum_deep)) 
                modpts['lon'].append(lon); modpts['lat'].append(lat); modpts['layer'].append(layer); 
                modpts['layer'].append(layer)
                modpts['deep'].append(depth+sum_deep)
                if depth+sum_deep>=0:
                 
                    return modpts,meanwindspeed
                    break
                
                if waterdepth<(abs(depth)):
                    print 'This point hits the land here.Less than %d meter.'%abs(depth)
                    raise Exception()
            except:
                return modpts 
        if len(windspeed)<2:
            print 12
            meanwindspeed=windspeed
        else:
            print 34
            meanwindspeed=np.mean(windspeed)                 
        return modpts,meanwindspeed
    def get_track1(self,lon,lat,depth,starttime,wind,wind_get_type,w1): #,b_index,nvdepth, 
        '''
        Get forecast points start at lon,lat
        '''
        windspeed=[]
        modpts = dict(lon=[lon], lat=[lat], layer=[], time=[],deep=[depth]) #model forecast points
        #uvz = netCDF4.Dataset(self.url)
        #u = uvz.variables['u']; v = uvz.variables['v']; zeta = uvz.variables['zeta']
        #print 'len u',len(u)
        #print modpts['lon']
        if lon>90:
            lon, lat = dm2dd(lon, lat)
        
        lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc,1)#1 day elements(blue)
        lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats,1)#1 day nodes(red)
        print 'in'
        try:
            if self.modelname == "GOM3" or self.modelname == "30yr":
                lonp,latp = self.nearest_point(lon, lat, lonl, latl,1)#elements(blue)
                
                lonn,latn = self.nearest_point(lon,lat,lonk,latk,1.2)#nodes(red)
                
            if self.modelname == "massbay":
                lonp,latp = self.nearest_point(lon, lat, lonl, latl,0.06)
                lonn,latn = self.nearest_point(lon,lat,lonk,latk,0.1)    
            index1 = np.where(self.lonc==lonp)
            index2 = np.where(self.latc==latp)
            elementindex = np.intersect1d(index1,index2)#nearest index elements(blue)
            index3 = np.where(self.lons==lonn)
            index4 = np.where(self.lats==latn)
            nodeindex = np.intersect1d(index3,index4)#nearest index nodes(red)
            
            
            pa = self.eline_path(lon,lat)# boundary 

            if self.modelname == "30yr":
                waterdepth = self.h[nodeindex]
            else:
                waterdepth = self.h[nodeindex]+self.zeta[0,nodeindex]
            
            modpts['time'].append(self.mTime[0])
            depth_total = self.siglay[:,nodeindex]*waterdepth; #print 'Here one' 
            for xx in np.arange(600):
                    if waterdepth<(abs(depth)):
                        depth=depth+5
                    if waterdepth<(abs(depth)):
                        continue
                    else:
                        break
            layer = np.argmin(abs(depth_total+depth)); #print 'layer',layer
            modpts['layer'].append(layer); 
            
            
            if waterdepth<(abs(depth)):
                print 'This point is too shallow.Less than %d meter.'%abs(depth)
                raise Exception()
        except:
            #print 12345
            return modpts
        sum_deep=0    
        t = abs(self.hours)
        print 'hour',t         
        for ii in xrange(t):  
            i=t-ii
            
                #print 'layer,lon,lat,i',layer,lon,lat,i
            lonl,latl = self.shrink_data(lon,lat,self.lonc,self.latc,1)
            lonk,latk = self.shrink_data(lon,lat,self.lons,self.lats,1)
            u_t1 = self.u[i,layer,elementindex][0]; v_t1 = self.v[i,layer,elementindex][0]
            u_t2 = self.u[i-1,layer,elementindex][0]; v_t2 = self.v[i-1,layer,elementindex][0]
            u_t = -(u_t1+u_t2)/2; v_t = -(v_t1+v_t2)/2
            #w_t1=self.ww[i,layer,nodeindex][0]
            #w_t2=self.ww[i,layer,nodeindex][0]
            #u_t,v_t = self.uvt(u_t1,v_t1,u_t2,v_t2)
            #w=(w_t1+w_t2)/2
            starttimes=starttime+timedelta(hours=i)
            #print 'starttime',starttimes
            if wind==0:

                dx = 60*60*u_t; dy = 60*60*v_t;dz=60*60*(w1)
            else:
                if wind_get_type=='NCEP':
                    v_wind,u_wind=self.get_wind_ncep(starttimes,lat,lon)
                    meanwindspeed=0
                if wind_get_type=='FVCOM':
                    v_wind,u_wind=self.get_wind_fvcom(starttimes,lat,lon)

                dx = 60*60*u_t+60*60*u_wind[0][0][0]*wind; dy = 60*60*v_t+60*60*v_wind[0][0][0]*wind
                windspeed.append(u_wind[0][0][0]*u_wind[0][0][0]+v_wind[0][0][0]*v_wind[0][0][0])
                
                #windspeed.append(v_wind[0][0][0]*v_wind[0][0][0]+u_wind[0][0][0]*u_wind[0][0][0])
            #u_t = (u_t1+u_t2)/2; v_t = (v_t1+v_t2)/2
            '''if u_t==0 and v_t==0: #There is no water
                print 'Sorry, hits the land,u,v==0'
                return modpts,1 #'''
            #mapx = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')                        
            #x,y = mapx(lon,lat)
            #temlon,temlat = mapx(x+dx,y+dy,inverse=True)            
            
            temlon = lon + (dx/(111111*np.cos(lat*np.pi/180)))#convert decimal to degree
            temlat = lat + dy/111111 #'''
            
            sum_deep=sum_deep+dz
            #print '%d,lat,lon,layer'%(i+1),temlat,temlon,layer
            #########case for boundary 1 #############
            if pa:
                teml = [(lon,lat),(temlon,temlat)]
                tempa = Path(teml)
                if pa.intersects_path(tempa): 
                    print 'Sorry, point hits land here.path'
                    if len(windspeed)<2:
                        print 12
                        meanwindspeed=windspeed
                    else:
                        print 34
                        meanwindspeed=np.mean(windspeed)                 
                    return modpts,meanwindspeed


            lon = temlon; lat = temlat
            #if i!=(t-1):                
            try:
                #print 12345
                if self.modelname == "GOM3" or self.modelname == "30yr":
                    lonp,latp = self.nearest_point(lon, lat, lonl, latl,1)
                    lonn,latn = self.nearest_point(lon,lat,lonk,latk,1.2)
                    
                if self.modelname == "massbay":
                    lonp,latp = self.nearest_point(lon, lat, lonl, latl,0.06)
                    lonn,latn = self.nearest_point(lon,lat,lonk,latk,0.1)
                index1 = np.where(self.lonc==lonp)
                index2 = np.where(self.latc==latp)
                elementindex = np.intersect1d(index1,index2); #print 'elementindex',elementindex
                index3 = np.where(self.lons==lonn)
                index4 = np.where(self.lats==latn)
                nodeindex = np.intersect1d(index3,index4)
                
            
                pa = self.eline_path(lon,lat)#boundary
                if self.modelname == "30yr":            
                    waterdepth = self.h[nodeindex]
                else:
                    waterdepth = self.h[nodeindex]+self.zeta[i+1,nodeindex]
                modpts['time'].append(self.mTime[i-1])            
                
                depth_total = self.siglay[:,nodeindex]*waterdepth 
                for xx in np.arange(600):
                    if waterdepth<(abs(depth)):
                        depth=depth+5
                    if waterdepth<(abs(depth)):
                        continue
                    else:
                        break
                layer = np.argmin(abs(depth_total+depth+sum_deep)) 
                modpts['lon'].append(lon); modpts['lat'].append(lat); modpts['layer'].append(layer); 
                modpts['layer'].append(layer)
                modpts['deep'].append(depth+sum_deep)
                if depth+sum_deep>=0:
                 
                    return modpts,meanwindspeed
                    break
                
                if waterdepth<(abs(depth)):
                    print 'This point hits the land here.Less than %d meter.'%abs(depth)
                    raise Exception()
            except:
                return modpts 
        if len(windspeed)<2:
            print 12
            meanwindspeed=windspeed
        else:
            print 34
            meanwindspeed=np.mean(windspeed)                 
        return modpts,meanwindspeed
