# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 23:26:30 2023

@author: ALEXIS
"""

###OPEN NETCDF FILES FOR THE REGIONAL CLIMATIC MODELS####
### EXTRACT INFORMATION FROM RCM SAM Rotated coordinates###

import xarray as xr
import os
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
directory=os.getcwd()
print(directory)

## open netcdf file, define directories
path=(r"F:\RCMs\Merged_SA\RCMs_Bolivia\RCP85")
pathcsv=(r"D:\UPWORK_FREELANCE\BOENCO\INPUT-LL")
Pathout=(r"D:\UPWORK_FREELANCE\BOENCO\INPUT-LL\0_Extracted\RCP85")
##Here change the file name
os.chdir(path)
#label name of the RCM
RCM_name0='1_pr_rcp85_CCCma-CanESM2_1'
RCM_name1='2_pr_rcp85_CSIRO-QCCCE-CSIRO-Mk3-6-0_2'
RCM_name2='3_pr_rcp85_ICHEC-EC-EARTH_3'
RCM_name3='4_pr_rcp85_IPSL-IPSL-CM5A-MR_4'
RCM_name4='5_pr_rcp85_MIROC-MIROC5_5'
RCM_name5='6_pr_rcp85_MPI-M-MPI-ESM-LR_6'
RCM_name6='7_pr_rcp85_NCC-NorESM1-M_7'
RCM_name7='8_pr_rcp85_NOAA-GFDL-GFDL-ESM2M_8'
#models
rcm0='7_pr_SAM-44_CCCma-CanESM2_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm1='5_pr_SAM-44_CSIRO-QCCCE-CSIRO-Mk3-6-0_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm2='4_pr_SAM-44_ICHEC-EC-EARTH_rcp85_r12i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm3='3_pr_SAM-44_IPSL-IPSL-CM5A-MR_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm4='2_pr_SAM-44_MIROC-MIROC5_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm5='15_pr_SAM-44_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm6='14_pr_SAM-44_NCC-NorESM1-M_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
rcm7='16_pr_SAM-44_NOAA-GFDL-GFDL-ESM2M_rcp85_r1i1p1_SMHI-RCA4_v3_day_Combined.nc'
##open and select models
###################################
data= xr.open_dataset(rcm6) ###Chose which RCM to process
RCM_name=RCM_name6
###################################
#print (data)
print(data.data_vars)
#print(data.variables.keys())
# Storing the lat and lon data into the variables 
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]

### Storing the lat and lon of station into variables 
os.chdir(pathcsv)
Coord_stations=pd.read_csv(' 01_Rotated_Estaciones_Coords.csv')  ##Change the name of the coordinate file
cds=Coord_stations

##Rotated coordinates
lat_station = -0.76
lon_station =  172.99

##Variables plot ncfile with stations coordinates
pcp=data.pr.isel(time=9150) ## change according the day you want to evaluate, according to the lenght spell of data of the model
pcp.plot()
plt.scatter(cds['rlon'],cds['rlat'],color='yellow',s=5, label='Estaciones Hidrologicas')
#plt.xlim([168,169]) #change according to coordinates limit
#plt.ylim([2,3]) #change according to coordinates limit
plt.title('Ubicacion de Estaciones Hidrologicas')
plt.legend( loc='upper left',labelcolor='w')
plt.show()
plt.clf()
#%%

date_range=data['time'].to_index()
df=pd.DataFrame(index=date_range)
df1=pd.DataFrame(index=date_range)
os.chdir(Pathout)
for i in range(len(cds)):
    lt_i=cds.loc[i,'rlat']
    ln_i=cds.loc[i,'rlon']
    name=cds.loc[i,'ID']
    da2=da.sel(rlat=lt_i,rlon=ln_i,method='nearest') ###Extractor of information
    da2=da2.to_dataframe()
    #df[name]=da2['pr']  #In case precipitation is in mm/d
    df[name]=da2['pr']*86400  #In case the precipitation is NOT in mm/d
    print('Process finished for: ' +str(i)+'_'+ str(name))

df.to_excel(str(RCM_name)+'_estaciones'+'.xlsx',encoding = 'utf-8-sig')    
print('Process finished for RCM Model: '+RCM_name)



