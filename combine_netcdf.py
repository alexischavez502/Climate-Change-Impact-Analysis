# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 13:11:13 2023

@author: ALEXIS
"""

### JOIN NETCDF FILES FOR THE REGIONAL CLIMATIC MODELS
import netCDF4
import xarray
import os

directory=os.getcwd()
print(directory)

##Files directories
ncpath=(r"F:\RCMs\South America")  ## Change Directory ncfiles
outputpath=(r"F:\RCMs\Merged_SA\RCMs_Bolivia\RCP85")  ## Change output directory

##Open ncfiles
##the * representing that whatever comes there is a don’t-care string
os.chdir(ncpath)
directory1=os.getcwd()
print(directory1)
### RCMs
rcm='pr_SAM-44_MPI-M-MPI-ESM-LR_historical_r1i1p1_MPI-CSC-REMO2009_v1_day_*.nc'
rcm1='pr_SAM-44_MOHC-HadGEM2-ES_historical_r1i1p1_SMHI-RCA4_v3_day_*.nc'
rcm2='pr_SAM-44_MIROC-MIROC5_rcp85_r1i1p1_SMHI-RCA4_v3_day_*.nc'
rcm3='pr_SAM-44_IPSL-IPSL-CM5A-MR_rcp85_r1i1p1_SMHI-RCA4_v3_day_*.nc'
rcm4='pr_SAM-44_ICHEC-EC-EARTH_rcp85_r12i1p1_SMHI-RCA4_v3_day_*.nc'


lt=[rcm,rcm1,rcm2,rcm3,rcm4]
listRCM=[rcm[:-5],rcm1[:-5],rcm2[:-5],rcm3[:-5],rcm4[:-5]]

##Function to combine files

def Combine(x):
    ds=xarray.open_mfdataset(x,
                             combine='nested', concat_dim="time")
    return ds

#Id codes:1,2,3,4,5
##change the number
id_codes=[1,2,3,4,5]
## Apply function 
dsf=Combine(rcm14) ## Here combines the RCM selected, so change the RCM variable according to need

os.chdir(outputpath) ## change directory to output
##export into a combined netCDF
i=1  ## to set the id code from the list
dsf.to_netcdf(str(id_codes[i])+'_'+listRCM[i]+'_Combined'+'.nc')
print('Processed finished for: '+str(id_codes[i])+'_'+listRCM[i])



