# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 11:32:59 2016

@author: Matt
"""
# Special thanks to Greg Reda for his website explaning pandas
# Point matching, high efficiency, variable with height generator

# Select folder of Level 3 Data
import numpy as np
from metpy.cbook import get_test_data
from metpy.io.nexrad import Level3File
import os
import numpy.ma as ma
from pyproj import Geod #Assign lat, Lon to each point using WGS 84 ellipsoid
import pandas as pd

g = Geod(ellps = 'WGS84')
RE = 6371; # radius of Earth
Ke = 4/3; # Equivalent Earth radius

# Bring in the file
dataPath = 'C:\\Users\\Matt\\Desktop\\KumjianCode\\PyTest\\'
file = 'C:\\Users\\Matt\\Desktop\\KumjianCode\\2015_Concatenated_FieldData.xlsx'
df = pd.read_excel(file)
groups = df.groupby('ID')
latGR = groups['LAT'].first()
lonGR = groups['LON'].first()
dataArray = []

for l in range(len(latGR)):
    print(l,latGR[l],lonGR[l])
    for filename in os.listdir(dataPath):
        print(filename)
        # Grab every file in the folder
        a = get_test_data(dataPath+filename)
     
        #Extract Level 3 data and place in an array format 
        b = Level3File(a)
        datadict = b.sym_block[0][0]
        data = b.map_data(datadict['data'])
        time = b.metadata.get('vol_time')
     
         # Turn into an array, then mask
        data[np.isnan(data)] = ma.masked
    
        # Grab azimuths and calculate a range based on number of gates
        az = np.array(datadict['start_az'] + [datadict['end_az'][-1]])
        az = np.delete(az,0,0) #Clear duplicate azimuth
        rng = np.linspace(0, b.max_range, data.shape[-1])
        theta = b.metadata['el_angle']
        beam = np.transpose(np.repeat(rng[:,np.newaxis],360,1))
 
        center_lat = np.ones([len(az),len(rng)])*b.lat    
        center_lon = np.ones([len(az),len(rng)])*b.lon
        print(center_lat[0][0],center_lon[0][0])
        az2D = np.ones_like(center_lat)*az[:,None]
        rng2D = np.ones_like(center_lat)*np.transpose(rng[:,None])*1000
        lon,lat,back=g.fwd(center_lon,center_lat,az2D,rng2D)
        
        # Use pyproj's Geod function to use it to find the distance between
        # the selected location where the hail was found and the rest of the 
        # points. We will use this to select what points we want.
        latsRe = np.broadcast_to(latGR[l],(lat.shape))
        lonsRe = np.broadcast_to(lonGR[l],(lon.shape))
        a1,a2,dist = g.inv(lonsRe,latsRe,lon,lat)
        dist = dist/1000 # Turn dist array to km from meters
        goodX,goodY = np.where(dist <= 2)
        if all(goodX == 0):
            print('ERROR: NO POINTS FOUND')
        else:
            dataPts = np.ones((len(goodX),2))
            timePts = []
            for d in range(len(goodX)):
                dataPts[d,0] = data[goodX[d],goodY[d]]
                timePts.append(time)
                dataPts[d,1] = dist[goodX[d],goodY[d]]
            dfdataCols = ['Data','Distance (km)']
            dfPts = pd.DataFrame(dataPts, columns = dfdataCols)
            timeDF = pd.DataFrame(timePts,columns=['Time'])
            finalDF = pd.concat([dfPts,timeDF],axis=1)
            dataArray.append(finalDF)
            
                
                
            
        
        
        # Now we have the final points that we want for our plots. Let's store
        # the values in a list and then we can go back to make the box plot.
        
        
        
            
            
        
         
        
            
            
