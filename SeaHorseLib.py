# SeaHorse Project Python Library
# Vitalii Sheremet, 2008-2012
# report bugs to vsheremet@whoi.edu
import numpy as np
from scipy import interpolate
import re

def RataDie(Y,M,D):
    """

RataDie(Y,M,D) -> RD    

returns serial day number in the (proleptic) Gregorian calendar
starting with January 1, year 1, as day 1 (Monday).

Gregorian calendar has a 400 year cycle (of which 97 are leap years)
of 146097 days, 20871 seven-day weeks. 
RataDie(1,1,1) -> 1  (Monday)
RataDie(2001,1,1) -> 146097*5+1 = 730486 (Monday)

RataDie is same as Python date2num or toordinal.
MATLAB datenum has January 1, year 0 as day 1, therefore
MATLAB datenum = RataDie + 366

M>0,D may be outside valid calendar range, e.g.,
2001,1,121 is a yearday 121 of 2001;
2001,1,0 is 2000-12-31.

D may have fractional part.
If TIMESTAMP='2001-01-121' is a string, use
(YYYY,MM,DD)=TIMESTAMP.split('-')
Y,M,D=int(YYYY),int(MM),int(DD)
or
sh_parse_timestamp(TIMESTAMP) -> yr,mo,da,hr,mi,se

Vitalii Sheremet, SeaHorse Project, 2008-2012
    """
#
# all divisions are integer, floor is not required
    RD=367*Y-(7*(Y+((M+9)/12))/4)-(3*(((Y+(M-9)/7)/100)+1)/4)+(275*M/9)+D-396;
    return RD
    
def sh_parse_timestamp(TIMESTAMP):
    """

sh_parse_timestamp(TIMESTAMP) -> yr,mo,da,hr,mi,se

parse TIMESTAMP string and convert to yr,mo,da,hr,mi,se

Acceptable formats:
'YYYY-MM-DD HR:MI:SE'
'YYYY-MM-DDTHR:MI:SE'
'YYYY-MM-DDTHR:MI:SEZ'
'YYYY/MM/DD HR:MI:SE'
'YYYY/MM/DDTHR:MI:SE'    
'YYYY/MM/DDTHR:MI:SEZ'

If MM,DD,HR,MI,SE are omited, 
then default values are assumed 
   01,01,00,00,00, respectively.
DD,HR,MI,SE may be fractional and outside formal ranges, e.g., 
TIMESTAMP='2001-01-121.25 24.5:120.4:90.1234'
 - 2001, yearday 121 00:00:00 + 0.25d + 24.5h + 120.4m +90.1234s
TIMESTAMP='2001-01 24.5'
 - 2001-01-01 00:00:00 + 24.5h

Vitalii Sheremet, SeaHorse Project, 2008-2012
    """
    YR='0000';MO='01';DA='01';HR='00';MI='00';SE='00';TIME='00:00:00'
    TIMESTAMP=TIMESTAMP.strip()
    if TIMESTAMP[-1]=='Z':
        TIMESTAMP=TIMESTAMP[0:-1]
        
    if TIMESTAMP.find(' ')>-1:
        DATE,TIME=TIMESTAMP.split(' ')
    elif TIMESTAMP.find('T')>-1:
        DATE,TIME=TIMESTAMP.split('T')
    else:
        DATE=TIMESTAMP
            
    if DATE.find('/')>-1:    
        CS='/'
    else:
        CS='-'
        
    DATE=DATE.split(CS)
    if len(DATE)==3:
        YR=DATE[0];MO=DATE[1];DA=DATE[2]
    elif len(DATE)==2:
        YR=DATE[0];MO=DATE[1]
    elif len(DATE)==1:
        YR=DATE[0]
    else:
        print('sh_parse_timestamp: error: unknown date format')
    
    TIME=TIME.split(':')
    if len(TIME)==3:
        HR=TIME[0];MI=TIME[1];SE=TIME[2]
    elif len(TIME)==2:
        HR=TIME[0];MI=TIME[1]
    elif len(TIME)==1:
        HR=TIME[0]
    else:
        print('sh_parse_timestamp: error: unknown time format')
             
    yr=int(YR);mo=int(MO);da=float(DA);hr=float(HR);mi=float(MI);se=float(SE)
    return yr,mo,da,hr,mi,se    


def sh_readhobo(FN):
#"Plot Title: 2039067_100212_WB_1-50_LRBY "
#"#"	"Date"	"Time, GMT+00:00"	"X Accel, g(LGR S/N: 2039067)"	"Y Accel, g(LGR S/N: 2039067)"	"Z Accel, g(LGR S/N: 2039067)"
#1	2/12/2010	23:00:00	-0.999	-0.325	0.125
# open file
    f = open(FN,'r')
# read header and 1st data line    
    nhdr=0;hdr='';
    line = f.readline()
    while line[0]=='"':
        nhdr=nhdr+1
        hdr=hdr+line
#        linevarnames=line
        line = f.readline()
    f.close()    
#    f.seek(0)


#'(LGR S/N: 9953408)'
    S=re.search('(LGR S/N: '+'\d+'+')',hdr)
    if S:
        S=S.group()
        S=re.search('\d+',S)
        SN=S.group()
    else:
        print(FN+' has no S/N in header')
        SN=''

# get Time Zone Info to do

    if hdr.find('X Accel, g')==-1:
        print(FN + ' has no X Accel, g column')
    if hdr.find('Y Accel, g')==-1:
        print(FN + ' has no Y Accel, g column')
    if hdr.find('Z Accel, g')==-1:
        print(FN + ' has no Z Accel, g column')


# sniff 1st data line for the date format: mm/dd/yyyy or yyyy-mm-dd
#1	2012-02-09	17:55:00	-0.850	0.600	0.300
#1	02/09/2012	17:55:00	-0.850	0.600	0.300
#    print(line)    
#    CS='\t'
    CS=None    
    if line.find(',')>-1:
        CS=','
    elif line.find(';')>-1:
        CS=';'
        
    S = line.split(CS)
    scan=S[0];date=S[1];time=S[2];xac=S[3];yac=S[4];zac=S[5]
    
    DS='-'
    if date.find('/')>-1:
        DS='/'
    
    s0,s1,s2=date.split(DS)
    if len(s0)==4:
        iyr=0;imo=1;ida=2 # standard format yyyy-mm-dd
    elif len(s2)==4:
        iyr=2;imo=0;ida=1 
    else:
        print(FN + ' unrecognized Date format')

    
# open file again and read data    
    f = open(FN,'r')
    hdr=[];sc=[];yr=[];mo=[];da=[];hr=[];mi=[];se=[];X=[];Y=[];Z=[];
    for line in f:
#  for line in f.readlines(): # equivalent, read all lines      
        if line[0] == '"':
            hdr.append(line)
        else:
            S = line.split(CS)
            scan=S[0];date=S[1];time=S[2];xac=S[3];yac=S[4];zac=S[5]
            sc.append(int(scan))
            yr.append(int(date.split(DS)[iyr]))
            mo.append(int(date.split(DS)[imo]))
            da.append(int(date.split(DS)[ida]))        
            hr.append(int(time.split(':')[0]))
            mi.append(int(time.split(':')[1]))        
            se.append(float(time.split(':')[2])) # seconds may have fractional part
            X.append(float(xac))
            Y.append(float(yac))
            Z.append(float(zac))
# Close the file
    f.close()
    sc=np.array(sc);
    yr=np.array(yr);mo=np.array(mo);da=np.array(da)
    hr=np.array(hr);mi=np.array(mi);se=np.array(se)
    X=np.array(X);Y=np.array(Y);Z=np.array(Z);
    return SN,sc,yr,mo,da,hr,mi,se,X,Y,Z
    
def sh_readcsv(FN):
    """

reads generic csv file with a header containing variable names    
tripple quotes comment text fragments
pound sign comments a line
# Vars: Field0, Field1, Field3, ... - specify variable names   

# Vars: DataFile, SN, ValidT1, ValidT2, AzYAx, MAvg, MODEL, SiteCode, SiteLatitude, SiteLongitude, SiteDepth
2039068_111201_1p50_ChRBr.txt,,,,,,,,,,
2039068_111201_1p50_ChRBr.txt,2039068,2011-12-01 18:00,2011/12/31T20:54:00,,,1p50,ChRBr,,,
2039068_111201_1p50_ChRBr.txt,2039068,2011-12-01 18:00,2011/12/31T20:54:00,2039068_111201_1p50_ChRBr_BaseOrient.txt,,1p50,ChRBr,,,
2039068_111201_1p50_ChRBr.txt,,,,,,,,,,
    """
    
# open file
    f = open(FN,'r')

    Items=[]
    skipflag=0  
    for fline in f:
        fline=fline.translate(None,'\n\r').strip()
#        print fline
        if fline[0:3] == '"""':
            if skipflag:
                skipflag=0 
            else:
                skipflag=1
            continue
        elif not skipflag:
            if fline.startswith('# Vars:'): 
# column names 
                fline=fline[len('# Vars:'):]
                CS=None
                Fields=fline.split(CS)                
                if fline.find(',')>-1:
                    CS=','
                    Fields=fline.translate(None,' ').split(CS)
                Items=[] # clear, use only the last list
            elif fline.startswith('#'):
                pass  
            else:
# read info for each datafile
                s=fline.split(CS)
# append if nonempty
#                print s
                if not ((s == ['']) or (s==[])): 
                    Items.append(s)
#
    f.close()
    return Items,Fields

def sh_getcfg(SN):
    FN='SH_'+SN+'.cfg'
#    try:
    Items,Fields=sh_readcsv(FN)
#    except:    
#        Items,Fields=sh_readcsv('SH_catalog.cfg')
    
    iSN=Fields.index('SN')
    iQX0=Fields.index('QX0')
    iQX1=Fields.index('QX1')
    iQY0=Fields.index('QY0')
    iQY1=Fields.index('QY1')
    iQZ0=Fields.index('QZ0')
    iQZ1=Fields.index('QZ1')
    iV0X=Fields.index('V0X')
    iV0Y=Fields.index('V0Y')
    iV0Z=Fields.index('V0Z')
    iMODEL=Fields.index('MODEL')

    for k in range(len(Items)):
        Item=Items[k]

        while len(Item)<len(Fields):
            Item.append('') # append optional fields

        SN1=Item[iSN]
        if SN1==SN:            
            QX0=float(Item[iQX0])
            QX1=float(Item[iQX1])
            QY0=float(Item[iQY0])
            QY1=float(Item[iQY1])
            QZ0=float(Item[iQZ0])
            QZ1=float(Item[iQZ1])
            V0X=float(Item[iV0X])
            V0Y=float(Item[iV0Y])
            V0Z=float(Item[iV0Z])
            MODEL=Item[iMODEL]
            return QX0,QX1,QY0,QY1,QZ0,QZ1,V0X,V0Y,V0Z,MODEL
        
def sh_tilt2vel(Tilt,MODEL):
    """function Vel=sh_tilt2vel(Tilt,MODEL)
% model drag curve parameters file
FN=['SHTiltCM' MODEL '.cfg'];
[Items,Fields]=sh_readcsv(FN);
iTilt=strmatch('Tilt',Fields,'exact');
iVel=strmatch('Vel',Fields,'exact');
N=length(Items);
TiltTab=zeros(N,1);
VelTab=zeros(N,1);
for k=1:N
    Item=Items{k};
    TiltTab(k)=str2double(Item{iTilt});
    VelTab(k)=str2double(Item{iVel});
end
Vel=interp1(TiltTab,VelTab,Tilt,'cubic','extrap');
end
    """
    FN='SHTiltCM'+MODEL+'.cfg'
    Items,Fields=sh_readcsv(FN)
    
    iTilt=Fields.index('Tilt')
#    print iTilt
    iVel=Fields.index('Vel')
#    print iVel
    N=len(Items)
#    print N
    TiltTab=np.zeros(N)
    VelTab=np.zeros(N)
    for k in range(N):
        Item=Items[k]
        TiltTab[k]=float(Item[iTilt])
        VelTab[k]=float(Item[iVel])

    tilt2vel=interpolate.interp1d(TiltTab,VelTab,kind='cubic',bounds_error=False) #returns function not array
#    tilt2vel=interpolate.interp1d(TiltTab,VelTab,kind='linear',bounds_error=False) #returns function not array
    Vel=tilt2vel(Tilt)
#    print TiltTab
#    print VelTab
#    print Tilt[0]
#    print Vel[0]
    return Vel

def sh_moments(X,Y):
    """
function [sx,sy,sxx,sxy,syy]=moments(X,Y)
j=find(~isnan(X) & ~isnan(Y));X=X(j);Y=Y(j);
N=length(X);
sx=mean(X);
sy=mean(Y);
X=X-sx;
Y=Y-sy;
sxx=mean(X.*X);
sxy=mean(X.*Y);
syy=mean(Y.*Y);
    """
    j=np.argwhere((X != np.NaN)&(Y != np.NaN))
    X=X[j].flatten();Y=Y[j].flatten();
    sx=np.mean(X);
    sy=np.mean(Y);
    X=X-sx;
    Y=Y-sy;
    sxx=np.mean(X*X);
    sxy=np.mean(X*Y);
    syy=np.mean(Y*Y);
    return sx,sy,sxx,sxy,syy
