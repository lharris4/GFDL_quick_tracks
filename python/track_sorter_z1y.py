#!/usr/bin/env python

## Version 0.1 17oct13 initial version
## Version 1.0b 30oct13: Many changes: new output formatting,
##     output files for each basin, option to not output
##     storm data before TS strength, option to not write figs
##     or do timing
## Version 1.0c 31oct13: Criteria as the existed before were those of
##     Zhao et al 2009, not Chen and Lin 2013. Altered so that 
##     only the genesis region need be equatorward of 40 latitude
## Version 1.1 6feb14: Changed hurricane criteria; must also satisfy
##     conditions for a TS. Also changed noPrecursor to False.
## Version 1.1b 20feb14: Removed bug preventing the 'World' files from
##     being written to. Also added the number of records (ie. number
##     of times a storm is identified) to header line of ASCII output.
## Version 1.1c 19mar14: Output text files now give four-digit-years
## Version 1.1d 23jul14: Added '-h' or '--hurThresh' option to select
##     hurricane threshold at run time; creates a file with a 'hThresh'
##     in its filename.
##     NOTE: need for a hurricane to also be a TS has been removed at
##     some point
## Version 1.1e 23jul14: Optimized the to_string call to greatly
##     speed up this bottleneck. Main loop now twice as fast!!
##     to_string still takes up nearly half of the main loop time, tho
## Version 1.1f 14aug14: added dates to output files, and '-t' or
##     '--TSThresh' option
## Version 1.1g 3sep14: replaced contains_point by a pre-computed
##      grid of basin locations, speeding up the sorter by 30%.
##      Also now have option to compute land mask.
## BRANCH 10sep14: Following Zhitao Yu's rules as closely as possible
##      modulo some bug fixes and some rules to permit general
##      world-wide year-round tracking.
## Version 1.1i 14apr16: Fixed bugs with new libraries, including a
##      bug that caused problems plotting hurricane tracks

import pandas as pd
import numpy as np
import itertools
from datetime import datetime, timedelta
from dateutil.parser import parse
import getopt, sys
import re
import matplotlib.path as path
import string as s
import time
import matplotlib
import math

#Allows script to run without a $DISPLAY set
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

filename='filename has not been set'
startdate=None
enddate=None
getHeader=False
doWarmCore=True
noPrecursor=True
doFigs=True
doTiming=True
plotStormNumber=False
debug=False
#writeAllStorms includes complete track for each storm,
#  including precursor, in the *.allstorms.* file, 
#  even if noPrecursor is True
#  The resulting file may get very large and performance
#   will be much slower
writeAllStorms=True # False

pd.set_option('chained_assignment', 'raise')

hurThresh = 32.5
TSThresh = 17.5

total_lifetime_thresh    =  timedelta(hours=72) # or attains hurricane strength while WC
WC_cumul_thresh          =  timedelta(hours=48)
WC_wind_contig_thresh    =  timedelta(hours=36)
wind_contig_thresh       =  timedelta(hours=36)

equatorward_of_lat = 40 #for warm-core storms

cumulativeTimers = {}
runningTimers = {}

init_grid_called = False
maskres = 0.25 # degrees
#Text input created using ncks:
#ncks -C -h -s '%3.2f\n' -H -v land_mask /archive/jhc/Hiram_C384L63/nudge_nh_0827_2/pp/atmos/ts/month/atmos.20001201.nc >! /work/lmh/research/seasonal/quick_tracks/land_mask_c384_1440x720.dat
#ncks                       -v land_mask /archive/jhc/Hiram_C384L63/nudge_nh_0827_2/pp/atmos/ts/month/atmos.20001201.nc >! /work/lmh/research/seasonal/quick_tracks/land_mask_c384_1440x720.text
landmaskfile = '/work/lmh/research/seasonal/quick_tracks/land_mask_c384_1440x720.dat'
global lon_grid
global lat_grid
global land_mask
global basinpoints
global nx_grid
global ny_grid

#http://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array
def run_of_true(bits):
    for bit, group in itertools.groupby(bits):
        if bit:
            yield sum(group) #returns an iterator; use 'return' to return a value
        else:
            yield 0

def max_run_of_true(bits):
    return max(run_of_true(bits))

def sum_run_of_true(bits):
    return sum(run_of_true(bits))

def date_parser(ss):
    year, month, day, hour = ss.split()
    year = int(year)
    #Ugly hack
    if year < 1900:
        year = 1900 + np.mod(year,300)
    st = '%04d-%02d-%02d %02d:00:00' % (int(year), int(month), int(day), int(hour))
    return pd.Timestamp(st) #, min, sec))    
    #return datetime(year, month, day, hour, 0, 0)
    

def plot_the_track(m, lon, lat, TCmask, hurmask, startind):
    x, y = m(lon[startind:], lat[startind:])

    xm = np.ma.MaskedArray(x,mask=~TCmask[startind:])
    ym = np.ma.MaskedArray(y,mask=~TCmask[startind:])
    xh = np.ma.MaskedArray(x,mask=~hurmask[startind:])
    yh = np.ma.MaskedArray(y,mask=~hurmask[startind:])
    
    #Contiguity requirement for TCs largely obviates
    # the need to draw lines for isolated TC points
    #xm, ym = half_interp_track(xm, ym, TCmask)
    xh, yh = half_interp_track(xh.data, yh.data, hurmask[startind:])

    m.plot(x,y,color='gray',linewidth=1.0)
    m.plot(xm,ym,color='k',linewidth=1.25)
    m.plot(xh,yh,color='r',linewidth=1.5)

    if ~any(TCmask) and plotStormNumber:
        yoffset = 0.022*(m.ymax-m.ymin)
        xl = plt.xlim()
        yl = plt.ylim()
        if m.xmin <= x[0] and x[0] <= m.xmax and m.ymin <= y[0] and y[0] <= m.ymax:
            plt.text(x[0],y[0],'o',fontsize=9,
                     ha='center',va='center')
            plt.text(x[0],y[0]-yoffset,str(number),fontsize=9,
                     ha='center',va='top')
        if m.xmin <= x[-1] and x[-1] <= m.xmax and m.ymin <= y[-1] and y[-1] <= m.ymax:
            plt.text(x[-1],y[-1],'X',fontsize=9,
                     ha='center',va='center')
            plt.text(x[-1],y[-1]-yoffset,str(number),fontsize=9,
                     ha='center',va='top')

#NOTE: with more recent versions of numpy library
#  weird things can happen if masked data is passed
# to this routine; pass the unmasked data instead.
def half_interp_track(x, y, mask):
    x2 = np.zeros(2*len(x) - 1)
    y2 = np.zeros(2*len(y) - 1)

    x2[0::2] = x
    x2[1::2] = 0.5*( x[0:-1]+x[1:] )

    y2[0::2] = y
    y2[1::2] = 0.5*( y[0:-1]+y[1:] )

    mask2 =  np.zeros(2*len(x) - 1,dtype=np.bool)
    mask2[0::2] = mask
    mask2[1::2] = mask[0:-1] & mask[1:]

    xm = np.ma.MaskedArray(x2,mask=~mask2)
    ym = np.ma.MaskedArray(y2,mask=~mask2)
    #print x.data
    #print xm[0::2]
    #print xm[1::2]

    return (xm, ym)

def plot_density(H, xedges, yedges, cmap):
    plt.figure().set_size_inches((8,3))
    m =  Basemap(llcrnrlon=0., llcrnrlat=-60.,
            urcrnrlon=360., urcrnrlat=60)
    x, y = m(xedges,yedges)
    m.pcolormesh(x,y,H.T,cmap=plt.get_cmap(cmap))
    m.drawcoastlines(color='0.75')
    m.drawparallels(np.linspace(-60,60,7),labels=[1,0,0,0])
    m.drawmeridians(np.linspace(0,360,7),labels=[0,0,0,1])
    plt.colorbar()

def start_clock(clockname):
    if doTiming:
        runningTimers[clockname] = time.time()
        if clockname not in cumulativeTimers:
            cumulativeTimers[clockname] = 0.

def stop_clock(clockname):
    if doTiming:
        cumulativeTimers[clockname] += time.time() - runningTimers[clockname]
        del runningTimers[clockname]

def print_clocks():
    if doTiming:
        for key  in cumulativeTimers.keys():
            print '%(key)30s: %(time)5.2f seconds' % {"key":key, "time": cumulativeTimers[key]}

def init_grid():
    global init_grid_called, lon_grid, lat_grid, nx_grid, ny_grid

    if init_grid_called:
        return

    #Define a grid by referencing midpoints
    lx  = np.arange(  0.+maskres*0.5, 360.+maskres, maskres)
    ly  = np.arange(-90.+maskres*0.5, 90.+maskres , maskres)
    nx_grid = lx.size
    ny_grid = ly.size
    lon_grid, lat_grid = np.meshgrid(lx,ly)

    init_grid_called = True

def init_land_mask():
    global land_mask, lon_grid

    init_grid()

    #lon_grid re-defined to cover one extra point in each direction
    #  to handle edge cases (lon 360 and lat 90)
    land_mask = np.zeros(np.array(lon_grid.shape))
    land_mask[:-1,:-1] = np.loadtxt(landmaskfile).reshape(np.array(lon_grid.shape) - 1) > 0.5
    land_mask[-1,:] = land_mask[0,:]
    land_mask[:,-1] = land_mask[:,0]

def over_land(lons, lats):
    global land_mask

    #np.floor no longer returns ints, just rounded floats
    iis = map(int, np.floor(lons/maskres).tolist())
    jjs = map(int, np.floor((lats + 90.)/maskres).tolist())

    return land_mask[jjs,iis]

def init_basin_points():
    global basinpoints, lon_grid, lat_grid

    init_grid()

    p = zip(lon_grid.ravel(),lat_grid.ravel())
    qs = []
    for n, zoom in enumerate(zooms):
        if zoom['domain'] is None:
            qs += [np.zeros(lon_grid.shape) < 1000.]
        else:
            qs += [zoom['domain'].contains_points(p).reshape(lon_grid.shape)] 
    basinpoints = qs # np.array(zip( *qs ) ).reshape((len(zooms),) + lon_grid.shape)

def inbasin(point):
    global nx_grid, ny_grid
    #Assumes lons start at 0
    # and lats start at -90
    #NOTE: what to do for a lon of 360.00 or a lat of 90.00?
    i = int(math.floor(point[0]/maskres))
    j = int(math.floor((point[1] + 90.)/maskres))

    return [basinpoints[n][j,i] for n, zoom in enumerate(zooms)]
        
#Begin main

start_clock('Init')

try:
    opts, args = getopt.getopt(sys.argv[1:], "s:e:gwch:t:d", ["startdate=", "enddate=", "getheader", "getHeader","WC","wc","doWarmCore","dowarmcore","NWC","nwc","noWarmCore","nowarmcore","hurThresh=","TSThresh="])
except getopt.GetoptError as err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit(2)
for o, a in opts:
    if o in ("-s", "--startdate"):
        startdate=a
    elif o in ("-e", "--enddate"):
        enddate=a
    elif o in ("-g", "--getHeader", "--getheader"):
        getHeader = True
    elif o in ("-w", "--WC", "--wc", "doWarmCore", "dowarmcore"):
        doWarmCore=True
    elif o in ("-c", "--NWC", "--nwc", "noWarmCore", "nowarmcore"):
        doWarmCore=False
    elif o in ("--hurThresh", "-h"):
        hurThresh = float(a)
    elif o in ("--TSThresh", "-t"):
        TSThresh = float(a)
    elif o in ("-d"):
        debug = True
    else:
        assert False, "unhandled option"

if args:
    filename = args[0]
    if len(args) >= 2:
        outfileprefix = args[1]
    else:
        outfileprefix = filename
        
        
#Can modify this so that we can determine column data from a file header.
#From that we can then import almost anything!
unames = ['key', 'year', 'month', 'day', 'hour', 'stormID', 'lon', 'lat', 'minSLP', 'maxWind', 'maxPrec', 'maxRain', 'maxSnow', 'maxVort', 'numBlz', 'WC', 'Area']
dtype  = {'key': np.int32,
          'year': 'S4',
          'month': 'S2',
          'day': 'S2',
          'hour': 'S2',
          'stormID': np.int32,
          'lon': np.float64,
          'lat': np.float64,
          'minSLP': np.float32,
          'maxWind': np.float64,
          'maxPrec': np.float64,
          'maxSnow': np.float64,
          'maxVort': np.float64,
          'numBlz': np.int32,
          'WC': np.float64,
          'Area': np.float64 }

start_clock('Get data')
if getHeader:
    trackerData = pd.read_table(filename,sep='\s+',lineterminator='\n',
                            index_col='date',comment='#',
                            date_parser=date_parser,
                            parse_dates={'date' : ['year','month','day','hour']},
                            keep_date_col=True)
else:
    trackerData = pd.read_table(filename,sep='\s+',lineterminator='\n',
                                names=unames,index_col='date',comment='#',
                                date_parser=date_parser,
                                parse_dates={'date' : ['year','month','day','hour']},
                                keep_date_col=True)

#timepart = trackerData[['year','month','day','hour']].astype('str')
#times = pd.to_datetime((timepart['year'] + ' ' + timepart['month'] + ' ' + timepart['day'] + ' ' + timepart['hour']).values, format="%Y %m %d %H")
    #times = pd.to_datetime(str(trackerData.year.values) + ' ' + str(trackerData.month.values) + ' ' + str(trackerData.day.values) + ' ' + str(trackerData.hour.values),"%Y %m %d %H")
#trackerData.set_index(times, inplace=True)
#trackerData = trackerData.drop(['year','month','day','hour'], axis=1)
    
zooms = ( {'Name': 'World', 'ShortName': 'world', 'corners': (0, -80, 360, 80),
           'projection': 'mill','parallels': np.linspace(-90,90,5),
           'meridians': np.linspace(0,360,7),
           'domain': None}, 
           {'Name': 'Northern Hemisphere', 'ShortName': 'NH', 'corners': (0, 0, 360, 60),
           'projection': 'mill','parallels': np.linspace(0,60,4),
           'meridians': np.linspace(0,360,7),
           'domain': path.Path(( (0,0), (0, 90), (360, 90), (360, 0) ))},
           {'Name': 'Western Pacific', 'ShortName': 'WPac', 'corners': (100, 0, 180, 60),
           'projection': 'mill','parallels': np.linspace(0,60,4),
           'meridians': np.linspace(100,180,5),
           'domain': path.Path(( (100, 0), (100, 90), (180, 90), (180, 0) ))},
           #{'Name': 'Western Pacific Nest', 'ShortName': 'WPacNest', 'corners': (100, 0, 150, 60),
           #'projection': 'mill','parallels': np.linspace(0,60,4),
           #'meridians': np.linspace(100,180,5),
           #'domain': path.Path(( (100, 0), (100, 90), (150, 90), (150, 0) ))},
           {'Name': 'Eastern Pacific', 'ShortName': 'EPac', 'corners': (220, 0, 270, 40),
           'projection': 'mill','parallels': np.linspace(0,40,5),
           'meridians': np.linspace(220,270,6),
           'domain': path.Path(( (295, 0), (260,20), (260,90), (220,90), (220, 0) ))},
           {'Name': 'North Atlantic', 'ShortName': 'NAtl', 'corners': (250,0,370,60),
           'projection': 'mill','parallels': np.linspace(0,60,4),
           'meridians': np.linspace(250,370,4),
           'domain': path.Path(( (295, 0), (260, 20), (260, 90), (360, 90), (360, 0) ))},
           #{'Name': 'North Atlantic Nest', 'ShortName': 'NAtlNest', 'corners': (250,0,310,60),
           #'projection': 'mill','parallels': np.linspace(0,60,4),
           #'meridians': np.linspace(250,370,4),
           #'domain': path.Path(( (295, 0), (260, 20), (260, 90), (310, 90), (310, 0) ))},
          #  {'Name': 'E North America', 'ShortName': 'ENAm', 'corners': (250,20,300,60),
          #  'projection': 'mill','parallels': np.linspace(0,60,4),
          #  'meridians': np.linspace(250,370,4),
          #  'domain': path.Path(( (250, 20), (300, 20), (300, 60), (250,60)   ))},
          #  {'Name': 'W North America', 'ShortName': 'WNAm', 'corners': (220,0,250,60),
          #  'projection': 'mill','parallels': np.linspace(0,60,4),
          #  'meridians': np.linspace(250,370,4),
          #  'domain': path.Path(( (250, 20), (220, 20), (220, 60), (250,60)   ))},
          # #Region defined by Chang 2013 JClim
          #  {'Name': 'North America', 'ShortName': 'NAm', 'corners': (260,10,330,70),
          #  'projection': 'mill','parallels': np.linspace(10,70,4),
          #  'meridians': np.linspace(270,330,4),
          #  'domain': path.Path(( (260, 10), (330, 10), (330,70), (260,70)   ))},
        )

init_basin_points()
init_land_mask() # not yet used

stop_clock('Get data')

stormcount = {'cy':  np.zeros(len(zooms),dtype=np.int32),
              'TC':  np.zeros(len(zooms),dtype=np.int32),
              'hur': np.zeros(len(zooms),dtype=np.int32)}
ShortNames = [ zoom['ShortName'] for zoom in zooms]

if startdate is None:
    sd = trackerData[0:1].index[0].to_datetime()
else:
    sd = parse(startdate)
    
if enddate is None:
    ed = trackerData[-1:].index[0].to_datetime()
else:
    ed = parse(enddate)

sdDispStr = sd.strftime('%Y-%m-%d')
edDispStr = ed.strftime('%Y-%m-%d')
sdShortStr = sd.strftime('%Y%m%d')
edShortStr = ed.strftime('%Y%m%d')

start_clock('Sort data')
# Remove un-grouped storms, if present
trackerData = trackerData[trackerData['stormID'] > 0]
times = trackerData.index
#trackerData = trackerData[sd:ed] #unreliable?
trackerData = trackerData[(times >= sd) & (times <= ed)]
stormGroups = trackerData.groupby('stormID')
stop_clock('Sort data')

figs = []
ms = []
if doFigs:

    for zoom in zooms:
        figs += [plt.figure()]
        m = Basemap(llcrnrlon=zoom['corners'][0], llcrnrlat=zoom['corners'][1],
                    urcrnrlon=zoom['corners'][2], urcrnrlat=zoom['corners'][3])
        ms += [m]
        m.drawcoastlines(color='0.75')
        #m.drawmapboundary(fill_color='#99ffff')
        #m.drawmapboundary(fill_color='aqua')
        m.fillcontinents(color='0.75')#,lake_color='#99ffff')
        m.drawparallels(zoom['parallels'],labels=[1,1,0,0])
        m.drawmeridians(zoom['meridians'],labels=[0,0,0,1])
        plt.title(zoom['Name'] + ' cyclone tracks: ' + sdDispStr + ' to ' + edDispStr)
    
if not doWarmCore:
    outfileprefix  += '.z1y-noWC'
else:
    outfileprefix += '.z1y-warm'

if hurThresh != 32.5:
    outfileprefix += ('.h' + str(hurThresh).replace('.','_'))
if TSThresh != 17.5:
    outfileprefix += ('.t' + str(TSThresh).replace('.','_'))

datesuffix =  sdShortStr + '-' + edShortStr

sumfmt = '%05d++++++++++++++++++++%05d   %06d    %12.4f  %4d %4d\n'
recfmt = '%05 %02d%02d%02d%02d  %03d %04d %4d %5.2f  %1d\n'
  
#Files for each basin
fstorms={}
fTCsum={}
fHTsum={}
nstorm = {}
nTC = {}
nhur = {}
for zoom in zooms:
    basin = zoom['ShortName']
    outfile = outfileprefix  + '.allstorms.' + basin + '.' + datesuffix + '.txt'
    TCsumfile = outfileprefix  + '.TS.' + basin + '.' + datesuffix + '.txt'
    HTsumfile = outfileprefix  + '.C15w.' + basin + '.' + datesuffix + '.txt'

    fstorms[basin] = open(outfile  ,'w')
    fTCsum[basin]  = open(TCsumfile,'w')
    fHTsum[basin]  = open(HTsumfile,'w')

    nstorm[basin] = 0
    nTC[basin] = 0
    nhur[basin] = 0

#Looks for date format output from to_string and then fixes it to conform to earlier trackplot output
#datefmt = re.compile('([0-9][0-9][0-9][0-9])-([0-9][0-9])-([0-9][0-9])\s+([0-9][0-9]):([0-9][0-9]):([0-9][0-9])')
datefmt = re.compile('([0-9][0-9][0-9][0-9])-([0-9][0-9])-([0-9][0-9])\s+([0-9][0-9]):([0-9][0-9]):([0-9][0-9])')

hist_nx = 72
hist_ny = 45
Hcy  = np.zeros((hist_nx,hist_ny))
Htc  = np.zeros((hist_nx,hist_ny))
Hhur = np.zeros((hist_nx,hist_ny))
genesis = np.zeros((hist_nx,hist_ny))
#Set up edges in case histograms are not made
dum, xedges, yedges = np.histogram2d([],[],bins=[hist_nx,hist_ny],range=[[0,360],[-90,90]])
del dum

stop_clock('Init')

start_clock('Main Loop')

start_clock('Loop head')

#Get storms
#Assumes hours interval is difference between the first two values
#(ie. cannot accept tracker data at varying intervals)
for stormID, thisstorm in stormGroups:

    #Eliminate identified centers never made into a track
    if stormID <= 0:
        continue
    if len(thisstorm) < 2:
        continue

    stop_clock('Loop head')
    start_clock('Storm proc')


    times = thisstorm.index
    #Some IBTrACS storms use a 12-hour interval instead
    # of a 6-hour interval. This re-computes dt for each storm.
    dt = (times[1].to_datetime() - times[0].to_datetime())
    dtdays = dt.seconds/86400.

    #This code is duplicated from below. Can be put into a function?
    if writeAllStorms:

        start_clock('Storm count')
        start_clock('contains point')
        indomain = inbasin( thisstorm[['lon','lat']][0:1].values[0] )
        stop_clock('contains point')
        stormcount['cy'] += indomain
        inbasins = [ShortName for ind, ShortName in zip(indomain, ShortNames)
                    if ind > 0 ]
        stop_clock('Storm count')

        start_clock('to_string 1')
        ar = thisstorm[['lon','lat','minSLP','maxWind','WC']]
        #st = ar.to_string(header=False,index_names=False) + '\n' #+ '\n'
        #st = '\n'.join(ar.to_string().split('\n')[2:]) + '\n' 
        st = '\n'.join(ar.to_string().split('\n')[2:]) + '\n' 
        stop_clock('to_string 1')
        st = datefmt.sub(r'    \1\2\3\4',st)

        for basin in inbasins:
            nstorm[basin] += 1
            sumstr = sumfmt % (nstorm[basin], nstorm[basin], stormID,
                               np.max(thisstorm['maxWind'].values), 
                               np.min(thisstorm['minSLP'].values), len(thisstorm)) 


            fstorms[basin].write(sumstr);
            fstorms[basin].write(st)
        

    #Performing Z1Y checks.

    # Check I: 36 consecutive hours of max winds >= TSThresh
    TSwind = thisstorm['maxWind'] >= TSThresh
    wind_contig = max_run_of_true(TSwind)*dt
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print thisstorm
    #     print 'CHECK I:', wind_contig, WC_wind_contig_thresh+dt
    # ### END DEBUG CODE
    if wind_contig < WC_wind_contig_thresh+dt:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue

    # Check II: 36 consecutive hours of warm core
    WC = thisstorm['WC'] > 0
    WC_contig = max_run_of_true(WC)*dt
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print 'CHECK II: ', WC_contig, WC_wind_contig_thresh+dt
    # ### END DEBUG CODE
    if WC_contig < WC_wind_contig_thresh+dt:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue

    # Check III: First time of winds over threshold must be over water (misses storms?)
    # This step actually commented out in Z1Y tracker

    is_land = over_land(thisstorm['lon'], thisstorm['lat'])
    if is_land[0]:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue        

    # Check IV: Genesis Equatorward of 40 latitude
    # Here we deviate from Z1Y, which just uses latitude of first detection.
    #  Doing this will count a lot of polar lows in the winter.
    TS       = WC & TSwind
    if not any(TS):
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue
        
    if noPrecursor:
        firstTS  = list(TS).index(True)
    else:
        firstTS  = 0
    latitude = np.abs(thisstorm['lat']) <= equatorward_of_lat
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print 'CHECK IV: ', latitude[firstTS], thisstorm['lat'], firstTS
    # ### END DEBUG CODE
    if not latitude[firstTS]:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue

    # Check V: Total wind >= TCthresh lifetime 
    # NOTE: The Z1Y sorter will NOT count a storm
    #  which exceeds the wind threshold for 12
    #  consecutive (6-hourly) times, but
    #  WILL count a storm which exceeds the wind threshold
    #  for 12 NON-consecutive times. This sorter counts neither.
    # See storms 4340 (2005 06 23 12 00 00   274.35    17.30    18.94    1  1004.26)
    #        and 6247 (2005 08 19 06 00 00   307.67    13.67    15.59    1  1009.11)
    #  in tracker test L3
    #  /home/Lucas.Harris/research/seasonal/quick_tracks/trunk/little_tests.sh
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print 'CHECK V: ', sum_run_of_true(TSwind)
    # ### END DEBUG CODE
    #if len(thisstorm)*dt < total_lifetime_thresh+dt:
    if sum_run_of_true(TSwind)*dt < total_lifetime_thresh+dt:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue
        
    # Check VI: Total WC time
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print 'CHECK VI: ', sum_run_of_true(WC)
    # ### END DEBUG CODE
    if sum_run_of_true(WC)*dt < WC_cumul_thresh+dt:
        stop_clock('Storm proc')
        start_clock('Loop head')
        continue        

    start_clock('Storm info')

    lon =  thisstorm['lon'].values
    lon[lon < 0] = lon[lon < 0] + 360.
#    thisstorm['lon'] = lon  #no longer assign back
    lat =  thisstorm['lat'].values

    start_clock('Storm count')
    indomain = inbasin ( (lon[TS][0], lat[TS][0] ))
#    indomain = inbasin ( thisstorm[['lon','lat']][TS].values[0] )
    stormcount['TC'] += indomain
    inbasins = [ShortName for ind, ShortName in zip(indomain, ShortNames)
                if ind > 0 ]
    basinstring = s.join(inbasins, ' ')
    
    #latitude is still TRUE for values equatorward of 40
    hur = (thisstorm['maxWind'] >= hurThresh) & WC
    # ### DEBUG CODE
    # if stormID == 6247:
    #     print 'Hurricane check: ', hur & latitude
    # ### END DEBUG CODE
    if any(hur & latitude):
        stormcount['hur'] += inbasin((lon[hur][0],lat[hur][0]))
#        stormcount['hur'] += inbasin(thisstorm[['lon','lat']][hur].values[0])

    #Set up wrap-around
    df = np.append(0,np.diff(lon))
    side = np.cumprod(np.sign(180-abs(df))) < 0
    if any(side):
        sg = -np.sign(df[side][0])
        lon2 = lon.copy()
        lon3 = lon.copy()
        lon3[side] += 360*sg
        lon2[~side] -= 360*sg

    stop_clock('Storm count')
    stop_clock('Storm info')

    start_clock('Storm I/O')
    if doFigs:
        for m, fig in zip(ms, figs):
            plt.figure(fig.number)
            if any(side):
                plot_the_track(m, lon2, thisstorm['lat'].values,
                          TS, hur, firstTS)
                plot_the_track(m, lon3, thisstorm['lat'].values,
                           TS, hur, firstTS)
            else:
                plot_the_track(m, lon, thisstorm['lat'].values,
                               TS, hur, firstTS)


    # Storm I/0
    start_clock('to_string 2')
    ar = thisstorm[firstTS:][['lon','lat','minSLP','maxWind','WC']]
    st = '\n'.join(ar.to_string().split('\n')[2:]) + '\n' #+ '\n'
    stop_clock('to_string 2')
    st = datefmt.sub(r'    \1\2\3\4',st)

    #Is nstorm the same as stormcount?
    for basin in inbasins:
        if not writeAllStorms:
            nstorm[basin] += 1
            sumstr = sumfmt % (nstorm[basin], nstorm[basin], stormID,
                               np.max(thisstorm[TS]['maxWind'].values), 
                               np.min(thisstorm[TS]['minSLP'].values), len(thisstorm[firstTS:])) 

            fstorms[basin].write(sumstr);
            fstorms[basin].write(st)

        nTC[basin] += 1
        sumstr = sumfmt % (nTC[basin], nTC[basin], stormID,
                           np.max(thisstorm[TS]['maxWind'].values), 
                           np.min(thisstorm[TS]['minSLP'].values), len(thisstorm[firstTS:])) 
        fTCsum[basin].write( sumstr);
        fTCsum[basin].write(st)

        if any(hur & latitude):
            nhur[basin] += 1
            sumstr = sumfmt %  (nhur[basin], nhur[basin], stormID,
                           np.max(thisstorm[TS]['maxWind'].values), 
                           np.min(thisstorm[TS]['minSLP'].values), len(thisstorm[firstTS:])) 
            fHTsum[basin].write(sumstr);
            fHTsum[basin].write(st)

    
    # Histogram
    if doFigs:
        Hd, xedges, yedges= np.histogram2d(lon, lat,
                                       bins=[hist_nx,hist_ny],range=[[0,360],[-90,90]])
        Hcy += Hd*dtdays
        Hd, xedges, yedges= np.histogram2d(
            lon, lat,
            bins=[hist_nx,hist_ny],range=[[0,360],[-90,90]])
        Htc += Hd*dtdays
        if any(hur & latitude):
            Hd, xedges, yedges= np.histogram2d(
                lon, lat,
                bins=[hist_nx,hist_ny],range=[[0,360],[-90,90]])
            Hhur += Hd*dtdays
        Hg, xedges, yedges = np.histogram2d(lon[TS][0:1], lat[TS][0:1],
            bins=[hist_nx,hist_ny],range=[[0,360],[-90,90]])
        genesis += Hg

    stop_clock('Storm I/O')
    stop_clock('Storm proc')
    start_clock('Loop head')

stop_clock('Loop head')
stop_clock('Main Loop')

start_clock('Finish')

if doFigs:

    for m, fig, zoom in zip(ms, figs, zooms):
        plt.figure(fig.number)
        plt.savefig(outfileprefix  + '.'+ zoom['ShortName'] + '.' + datesuffix + '.png',
                    bbox_inches='tight')

    # Plot histograms
    plot_density(Hcy, xedges, yedges, 'Greens')
    plt.title('Cyclone storm-days:\n ' + sdDispStr + ' to ' + edDispStr)
    plt.savefig(outfileprefix  + '.cy.' + datesuffix + '.png',
                bbox_inches='tight')

    plot_density(Htc, xedges, yedges, 'Greys')
    plt.title('Tropical cyclone storm-days:\n ' + sdDispStr + ' to ' + edDispStr)
    plt.savefig(outfileprefix  + '.TC.' + datesuffix + '.png',
                bbox_inches='tight')

    plot_density(Hhur, xedges, yedges, 'Reds')
    plt.title('Hurricane/typhoon storm-days:\n ' + sdDispStr + ' to ' + edDispStr)
    plt.savefig(outfileprefix  + '.hur.' + datesuffix + '.png',
                bbox_inches='tight')

    plot_density(genesis, xedges, yedges, 'Reds')
    plt.title('Tropical cyclone genesis events:\n ' + sdDispStr + ' to ' + edDispStr)
    plt.savefig(outfileprefix  + '.gen.' + datesuffix + '.png',
                bbox_inches='tight')

    #Save histograms
    np.savez(outfileprefix  + '.hist' + datesuffix,
             Hcy=Hcy, Htc=Htc, Hhur=Hhur, genesis=genesis,
             xedges=xedges,yedges=yedges,zooms=zooms,stormcount=stormcount,
             nstorm=nstorm,nTC=nTC,nhur=nhur)

stop_clock('Finish')
print_clocks()

#plt.show()

for f in fstorms.itervalues():
    f.close()

for f in fTCsum.itervalues():
    f.close()

for f in fHTsum.itervalues():
    f.close()
