# -*- coding: utf-8 -*-
"""
@author: Vitalii Sheremet, SeaHorse Project
"""
import numpy as np

def sh_rmtide(f,dt=1.,ends=0.):
    """
    removes solar and lunar tidal signal by sequentially averaging 
    over their periods: of 24h and 24h50.47m. This is equivalent to
    applying convolution with a trapezoidal shaped window (almost triagular).
    
    f - input tidal sequence
    dt - uniform time step of input series, in hours, dt=1. hr default
    ends - fill value for the ends: 0. or nan
    
    fm = sh_rmtide(f,dt=1.,ends=np.nan)
    """
    TS=24. # Solar hour angle period 24h
    TL=24.+50.47/60. # Lunar hour angle period
    N=len(f)
    fm=np.zeros(N)*ends

# remove solar period    
    T=TS/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
# all weights =1 except at ends      
#    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
        
# remove lunar period
    f=fm*1.0  # deep copy!  
    T=TL/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
# all weights =1 except at ends      
#    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
    
    return fm

def sh_interp3(ti,t,f):
    """
    interpolates irragularly spaced 1D data to uniform grid
    using cubuc Hermite spline.
    http://en.wikipedia.org/wiki/Cubic_Hermite_spline
    tangents are defined as a finite difference approximation of derivative
    interpolation formula is expressed in terms of 
    Bernstein polynomials of 3rd order: B_0,B_1,B_2,B_3.
    """
    m=f*0. # derivative
    N=len(f)
    for i in range(1,N-1):
# In this version the tangent m is a plain average of differences from two sides
# which is not consistent with the finite difference
# approximation of a derivative on a nonuniform grid
#        m[i]=((f[i+1]-f[i])/(t[i+1]-t[i])+(f[i]-f[i-1])/(t[i]-t[i-1]))*0.5
# This version of the tangent m is the finite difference
# approximation of a derivative on a nonuniform grid 
        ha=t[i+1]-t[i]
        hb=t[i]-t[i-1]        
        m[i]=((f[i+1]-f[i])*hb/ha+(f[i]-f[i-1])*ha/hb)/(ha+hb)
    # ends
    m[0]  =2.*(f[1]  -f[0])  /(t[1]  -t[0])  -m[1]
    m[N-1]=2.*(f[N-1]-f[N-2])/(t[N-1]-t[N-2])-m[N-2]    
        
    # 
    fi=ti*0.
    NN=len(ti)
    for ii in range(NN):
        i=max(np.argwhere(t<=ti[ii])) # left end of segment
        dt=t[i+1]-t[i]
        a=(ti[ii]-t[i])/dt  # 0<a<1 normalized coordinate within a segment 
        b=1.-a
        # multiply f and m*dt must have the same dimensions
        fi[ii]=(b*b*b+3.*a*b*b)*f[i]+(a*a*a+3.*b*a*a)*f[i+1]+(a*b*b*m[i]-b*a*a*m[i+1])*dt
    
    return fi