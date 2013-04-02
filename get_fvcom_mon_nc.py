# -*- coding: utf-8 -*-
"""
access FVCOM data (individual monthly nc files) via pydap or netCDF4 and save to local disk as
YYYY/variable/YYYYMMDDHRMISE_variable.npy

@author: Vitalii Sheremet, FATE Project, 2012, vsheremet@whoi.edu
"""
"""
Accessing FVCOM GOM3 30yr hindcast data
DDS:

Dataset {
    Float32 a1u[four = 4][nele = 90415];
    Float32 a2u[four = 4][nele = 90415];
    Float32 art1[node = 48451];
    Float32 art2[node = 48451];
    Float32 aw0[three = 3][nele = 90415];
    Float32 awx[three = 3][nele = 90415];
    Float32 awy[three = 3][nele = 90415];
    Float32 cc_hvc[nele = 90415];
    Float32 h[node = 48451];
    Float32 lat[node = 48451];
    Float32 latc[nele = 90415];
    Float32 lon[node = 48451];
    Float32 lonc[nele = 90415];
    Int32 nbe[three = 3][nele = 90415];
    Int32 nbsn[maxnode = 11][node = 48451];
    Int32 nbve[maxelem = 9][node = 48451];
    Float32 nn_hvc[node = 48451];
    Int32 nprocs;
    Int32 ntsn[node = 48451];
    Int32 ntve[node = 48451];
    Int32 nv[three = 3][nele = 90415];
    Int32 partition[nele = 90415];
    Float32 siglay[siglay = 45][node = 48451];
    Float32 siglev[siglev = 46][node = 48451];
    Float32 x[node = 48451];
    Float32 xc[nele = 90415];
    Float32 y[node = 48451];
    Float32 yc[nele = 90415];
    Float32 z0b[nele = 90415];
    Int32 Itime[time = 171883];
    Int32 Itime2[time = 171883];
    String Times[time = 171883];
    String file_date[time = 171883];
    Int32 iint[time = 171883];
    Float32 kh[time = 171883][siglev = 46][node = 48451];
    Float32 km[time = 171883][siglev = 46][node = 48451];
    Float32 kq[time = 171883][siglev = 46][node = 48451];
    Float32 l[time = 171883][siglev = 46][node = 48451];
    Float32 net_heat_flux[time = 171883][node = 48451];
    Float32 omega[time = 171883][siglev = 46][node = 48451];
    Float32 q2[time = 171883][siglev = 46][node = 48451];
    Float32 q2l[time = 171883][siglev = 46][node = 48451];
    Float32 salinity[time = 171883][siglay = 45][node = 48451];
    Float32 short_wave[time = 171883][node = 48451];
    Float32 temp[time = 171883][siglay = 45][node = 48451];
    Float32 time[time = 171883];
    Float32 u[time = 171883][siglay = 45][nele = 90415];
    Float32 ua[time = 171883][nele = 90415];
    Float32 uwind_stress[time = 171883][nele = 90415];
    Float32 v[time = 171883][siglay = 45][nele = 90415];
    Float32 va[time = 171883][nele = 90415];
    Float32 vwind_stress[time = 171883][nele = 90415];
    Float32 ww[time = 171883][siglay = 45][nele = 90415];
    Float32 zeta[time = 171883][node = 48451];
} fvcom%2fhindcasts%2f30yr_gom3;
"""
import numpy as np
from SeaHorseLib import *
from datetime import *
import os
#from pydap.client import open_url  # pydap not yet supported by epd
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
#import matplotlib.mlab as ml

#http://www.smast.umassd.edu:8080/thredds/catalog/models/fvcom/NECOFS/Archive/catalog.html
#http://www.smast.umassd.edu:8080/thredds/catalog.html
#http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_gom3.html
#URL='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3.html'
#URL='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
#http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?
#a1u[0:1:3][0:1:90414],a2u[0:1:3][0:1:90414],art1[0:1:48450],art2[0:1:48450],
#aw0[0:1:2][0:1:90414],awx[0:1:2][0:1:90414],awy[0:1:2][0:1:90414],cc_hvc[0:1:90414],
#h[0:1:48450],lat[0:1:48450],latc[0:1:90414],lon[0:1:48450],lonc[0:1:90414],
#nbe[0:1:2][0:1:90414],nbsn[0:1:10][0:1:48450],nbve[0:1:8][0:1:48450],
#nn_hvc[0:1:48450],nprocs,ntsn[0:1:48450],ntve[0:1:48450],nv[0:1:2][0:1:90414],
#partition[0:1:90414],siglay[0:1:44][0:1:48450],siglev[0:1:45][0:1:48450],
#x[0:1:48450],xc[0:1:90414],y[0:1:48450],yc[0:1:90414],z0b[0:1:90414],
#Itime[0:1:171882],Itime2[0:1:171882],Times[0:1:171882],file_date[0:1:171882],
#iint[0:1:171882],kh[0:1:171882][0:1:45][0:1:48450],
#km[0:1:171882][0:1:45][0:1:48450],kq[0:1:171882][0:1:45][0:1:48450],
#l[0:1:171882][0:1:45][0:1:48450],net_heat_flux[0:1:171882][0:1:48450],
#omega[0:1:171882][0:1:45][0:1:48450],q2[0:1:171882][0:1:45][0:1:48450],
#q2l[0:1:171882][0:1:45][0:1:48450],salinity[0:1:171882][0:1:44][0:1:48450],
#short_wave[0:1:171882][0:1:48450],temp[0:1:171882][0:1:44][0:1:48450],
#time[0:1:171882],u[0:1:171882][0:1:44][0:1:90414],ua[0:1:171882][0:1:90414],
#uwind_stress[0:1:171882][0:1:90414],v[0:1:171882][0:1:44][0:1:90414],
#va[0:1:171882][0:1:90414],vwind_stress[0:1:171882][0:1:90414],
#ww[0:1:171882][0:1:44][0:1:90414],zeta[0:1:171882][0:1:48450]

#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/gom3_199701.nc?time[0:1:744],zeta[0:1:744][0:1:48450],u[0:1:744][0][0:1:90414],v[0:1:744][0][0:1:90414],ua[0:1:744][0:1:90414],va[0:1:744][0:1:90414],temp[0:1:744][0:1:44][0:1:48450],uwind_stress[0:1:744][0:1:90414],vwind_stress[0:1:744][0:1:90414]
#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/gom3_199702.nc?time[0:1:672],zeta[0:1:672][0:1:48450],u[0:1:672][0][0:1:90414],v[0:1:672][0][0:1:90414],ua[0:1:672][0:1:90414],va[0:1:672][0:1:90414],temp[0:1:672][0:1:44][0:1:48450],uwind_stress[0:1:672][0:1:90414],vwind_stress[0:1:672][0:1:90414]
#ds=open_url(URL)

#xxx=ds['xxx']; np.save('gom3.xxx.npy',np.array(xxx))

#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/
#gom3_199701.nc?
#time[0:1:744],zeta[0:1:744][0:1:48450],
#u[0:1:744][0][0:1:90414],v[0:1:744][0][0:1:90414],
#ua[0:1:744][0:1:90414],va[0:1:744][0:1:90414],
#temp[0:1:744][0:1:44][0:1:48450],
#uwind_stress[0:1:744][0:1:90414],vwind_stress[0:1:744][0:1:90414]

URL0='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/'

#Available range 1978-2012
#for yr in range(1978,2013):
for yr in range(2001,2002):
    YYYY=str(yr)

    DN='GOM3_'+YYYY
    if os.path.exists(DN)==False:
        os.mkdir(DN)
    DN='GOM3_'+YYYY+'/'+'zeta'
    if os.path.exists(DN)==False:
        os.mkdir(DN)
    DN='GOM3_'+YYYY+'/'+'u44'
    if os.path.exists(DN)==False:
        os.mkdir(DN)
    DN='GOM3_'+YYYY+'/'+'v44'
    if os.path.exists(DN)==False:
        os.mkdir(DN)
    DN='GOM3_'+YYYY+'/'+'temp44'
    if os.path.exists(DN)==False:
        os.mkdir(DN)    

#    for mo in range(1,13):
    for mo in range(1,2):
        FN0='gom3_'+str(yr)+str(mo).zfill(2)
        FN1=FN0+'.nc'
        print FN1
        
#
        URL=URL0+FN1
#        ds=open_url(URL)                 # pydap version 
        ds = Dataset(URL,'r').variables   # netCDF4 version

# 
        Times=ds['Times'][:]
        for i in range(len(Times)):
            Time=Times[i]
            T=Time.tostring()
#            TS=T[0:4]+T[5:7]+T[8:10]+T[11:13]+T[14:16]+T[17:19]
            yr,mo,da,hr,mi,se=sh_parse_timestamp(T)
            tRD=RataDie(yr,mo,da)+(hr*60*60+mi*60+se)/86400.

            if True: #da>26:

    # sometimes fails 199407005959
    # need to round to 1 hour
                tt=np.round(tRD*24.)/24.
                ti=datetime.fromordinal(int(tt))
                YEAR=str(ti.year)
                MO=str(ti.month).zfill(2)
                DA=str(ti.day).zfill(2)
                hr=(tt-int(tt))*24
                HR=str(int(np.round(hr))).zfill(2)            
                TS=YEAR+MO+DA+HR+'0000'
                FN2=TS
    
                tchk=RataDie(int(YEAR),int(MO),int(DA))+int(HR)/24.
                if (tRD-tchk)*24. > 0.5 :
                    print 'warning tchk ', kdr,TS
            
    #            time=ds['time'][i];    np.save('GOM3_'+YYYY+'/'+'time/'+FN2+'_time.npy',np.array(time))
                zeta=ds['zeta'][i,:]; np.save('GOM3_'+YYYY+'/'+'zeta/'+FN2+'_zeta.npy',np.array(zeta))
    #            u0=ds['u'][i,0,:];   np.save('GOM3_'+YYYY+'/'+'u0/'+FN2+'_u0.npy',np.array(u0))
    #            v0=ds['v'][i,0,:];   np.save('GOM3_'+YYYY+'/'+'v0/'+FN2+'_v0.npy',np.array(v0))
    #            u44=ds['u'][i,44,:];   np.save('GOM3_'+YYYY+'/'+'u44/'+FN2+'_u44.npy',np.array(u44))
    #            v44=ds['v'][i,44,:];   np.save('GOM3_'+YYYY+'/'+'v44/'+FN2+'_v44.npy',np.array(v44))
        # averaged velocities
    #            ua=ds['ua'][i,:]; np.save('GOM3_'+YYYY+'/'+'ua/'+FN2+'_ua.npy',np.array(ua))
    #            va=ds['va'][i,:]; np.save('GOM3_'+YYYY+'/'+'va/'+FN2+'_va.npy',np.array(va))
        # surface temperature
    #            temp0=ds['temp'][i,0,:]; np.save('GOM3_'+YYYY+'/'+'temp0/'+FN2+'_temp0.npy',np.array(temp0))
    #            temp44=ds['temp'][i,44,:]; np.save('GOM3_'+YYYY+'/'+'temp44/'+FN2+'_temp44.npy',np.array(temp44))
    #            uwind_stress=ds['uwind_stress'][i,:]; np.save('GOM3_'+YYYY+'/'+'uwind_stress/'+FN2+'_uwind_stress.npy',np.array(uwind_stress))
    #            vwind_stress=ds['vwind_stress'][i,:]; np.save('GOM3_'+YYYY+'/'+'vwind_stress/'+FN2+'_vwind_stress.npy',np.array(vwind_stress))
    #           
                print TS
#
