# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 23:32:53 2023

@author: ALEXIS
"""
#Source: https://www.mathworks.com/matlabcentral/fileexchange/43435-rotated-grid-transform
#https://gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
from math import *

def rotated_grid_transform(grid_in, option, SP_coor):
    lon = grid_in[0]
    lat = grid_in[1];

    lon = (lon*pi)/180; # Convert degrees to radians
    lat = (lat*pi)/180;

    SP_lon = SP_coor[0];
    SP_lat = SP_coor[1];

    theta = 90+SP_lat; # Rotation around y-axis
    phi = SP_lon; # Rotation around z-axis

    theta = (theta*pi)/180;
    phi = (phi*pi)/180; # Convert degrees to radians

    x = cos(lon)*cos(lat); # Convert from spherical to cartesian coordinates
    y = sin(lon)*cos(lat);
    z = sin(lat);

    if option == 1: # Regular -> Rotated

        x_new = cos(theta)*cos(phi)*x + cos(theta)*sin(phi)*y + sin(theta)*z;
        y_new = -sin(phi)*x + cos(phi)*y;
        z_new = -sin(theta)*cos(phi)*x - sin(theta)*sin(phi)*y + cos(theta)*z;

    else:  # Rotated -> Regular

        phi = -phi;
        theta = -theta;

        x_new = cos(theta)*cos(phi)*x + sin(phi)*y + sin(theta)*cos(phi)*z;
        y_new = -cos(theta)*sin(phi)*x + cos(phi)*y - sin(theta)*sin(phi)*z;
        z_new = -sin(theta)*x + cos(theta)*z;



    lon_new = atan2(y_new,x_new); # Convert cartesian back to spherical coordinates
    lat_new = asin(z_new);

    lon_new = (lon_new*180)/pi; # Convert radians back to degrees
    lat_new = (lat_new*180)/pi;

    #print (lon_new,lat_new)
    return lon_new,lat_new

##FOR SP_COOR : rotated _latitude_longitude  - verify the variable in the domain
##FOR SAM domain: ## define or change according to domain
grid_north_pole_latitude = 70.60
grid_north_pole_longitude = -56.06
## order of input: Longitude,Latitude, option 1, north pole longitude, north pole latitude
rlon=180+grid_north_pole_longitude #0 is south pole
rlat=0-grid_north_pole_latitude #negative because is ecuator and minus cause is Southern hemisphere
test=rotated_grid_transform((-63.5258,-20.0067), 1, (rlon,rlat))
print(test)

#%%
import os
import pandas as pd
directory=os.getcwd()
print(directory)
pathcsv=(r"D:\UPWORK_FREELANCE\BOENCO\INPUT-LL") ## Change directory of csv file
os.chdir(pathcsv)
Coord_stations=pd.read_csv('Estaciones_coord.csv')

ss=Coord_stations
ss['rlon']=""
ss['rlat']=""
for i in range(len(ss)):
    lt_i=ss.loc[i,'Lat']
    ln_i=ss.loc[i,'Lon']
    var=rotated_grid_transform((ln_i,lt_i), 1, (rlon,rlat))
    ss.loc[i,'rlon']=round(var[0],2)
    ss.loc[i,'rlat']=round(var[1],2)

convert_dict = {'rlon': float,
                'rlat': float
                }
ss=ss.astype(convert_dict)
print(ss)
ss.to_csv(pathcsv +'\ '+'01_Rotated_Estaciones_Coords.csv',encoding = 'utf-8-sig')
