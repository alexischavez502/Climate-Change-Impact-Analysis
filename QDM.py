# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 19:02:48 2023

@author: ALEXIS
"""

## DOWNSCALING QUANTILE DELTA MAPPING 
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymannkendall as pmk 
import seaborn as sns


#read file
## CHANGE DIRECTORIES
path_observations=r'D:\UPWORK_FREELANCE\BOENCO\INPUT-LL'
path_historical=r'D:\UPWORK_FREELANCE\BOENCO\INPUT-LL\0_Extracted\Historical'
path_future=r'D:\UPWORK_FREELANCE\BOENCO\INPUT-LL\0_Extracted\RCP45'
salida=r'D:\UPWORK_FREELANCE\BOENCO\INPUT-LL\02_Perturbed_Series\RCP45'


#QUANTILE DELTA MAPPING
##Input
i=1 ## change i to evaluate another station
##CHANGE NAME OF MODELS
model_name='1_pr_rcp45_CCCma-CanESM2_1'
os.chdir(path_observations)
obs=pd.read_excel('LaPaz_Pr_Obs.xlsx',usecols=[0,i],index_col=0)
os.chdir(path_historical)
curr_clim=pd.read_csv('1_pr_historical_CCCma-CanESM2_1_estaciones.csv',usecols=[0,i],index_col=0)
curr_clim.index=pd.to_datetime(curr_clim.index,format='%Y-%m-%d', errors='coerce')
os.chdir(path_future)
fut_clim=pd.read_csv('1_pr_rcp45_CCCma-CanESM2_1_estaciones.csv',usecols=[0,i],index_col=0)
fut_clim.index=pd.to_datetime(fut_clim.index,format='%Y-%m-%d', errors='coerce')


Nyears=round(len(fut_clim)/365.5,0) 
Obsyears=round(len(obs)/365.5,0)
#create output dataframe
df_salida=pd.DataFrame()
########
##Filtering Time periods

def filter_date(df):
    start_date='1980-01-01'
    end_date='2018-12-31'
    mask = (df.index >= start_date) & (df.index <= end_date)
    df=df.loc[mask]
    return df

def filter_date_HIST(df):
    start_date='1965-01-01'
    end_date='2005-01-01'
    mask = (df.index >= start_date) & (df.index <= end_date)
    df=df.loc[mask]
    return df

obs=filter_date(obs)
curr_clim=filter_date_HIST(curr_clim)

def filter_date1(df):
    start_date='2030-01-01'
    end_date='2070-01-01'
    mask = (df.index >= start_date) & (df.index <= end_date)
    df=df.loc[mask]
    return df

fut_clim=filter_date1(fut_clim)

##get the name of the column
name_station=obs.columns.values
name_station=name_station[0]
####mankendal test
df_yearly_Obs=obs.copy()
df_yearly_Obs['year']=df_yearly_Obs.index.year
df_yearly_Obs=df_yearly_Obs.groupby('year')[name_station].sum()
df_yearly_Obs=df_yearly_Obs.to_frame()

plt.scatter(df_yearly_Obs.index,df_yearly_Obs[name_station])
plt.show()
plt.clf()

test = pmk.original_test(df_yearly_Obs[name_station], alpha=0.05)

print("\n"+'Precipitation Observed at: '+str(name_station), test)


def POT(df,nam):
    POT=pd.DataFrame()
    POT[nam]=df[nam]
    POT.sort_values(by=nam,ascending=False, inplace=True)
    POT['Ranked']=range(1,len(POT)+1) 
    POT= POT.astype({'Ranked':'float'})
    n=len(POT) #sample size
    years=round(n/365) #number of years
    def ERP(x):
        return (years/(x))
    POT['Return_Period']=POT['Ranked'].apply(ERP)
    return POT

C_climate=POT(curr_clim,name_station)
Obs_climate=POT(obs,name_station)

plt.scatter(C_climate['Return_Period'],C_climate[name_station],s=10,label='Model_Historical')
plt.scatter(Obs_climate['Return_Period'],Obs_climate[name_station],s=10,label='Observations')
plt.xscale('log')
plt.legend()
plt.title(str(name_station))
plt.show()
plt.clf()

#%%
dailyO=pd.DataFrame()
dailyO['DailyO']=obs[name_station]
dailyMC=curr_clim.copy()
dailyMT=fut_clim.copy()
datesO=['1980-01-01','2005-12-31'] ##Change Time period according to historical observations
datesO = [pd.to_datetime(date) for date in datesO]

dailyO=dailyO['DailyO'].tolist()
dailyMC=dailyMC[name_station].to_list()
dailyMT=dailyMT[name_station].to_list()

datesO = obs.index.tolist()
datesC = curr_clim.index.tolist()
datesT = fut_clim.index.tolist()

def QDM(var, dailyO, dailyMC, dailyMT, datesO, datesC, datesT):
    
    # Input
    # dailyO: serie obs
    # dailyMC: serie sim historica
    # dailyMT: serie sim futura
    #datesO: lista fechas serie obs
    #datesC: lista fecha clima historico modelo
    #datexT: lista fechas clima futuro modelo
    # dailyO: [1970,2005]
    # dailyMC: [1970,2005]
    # dailyMT: [2036,2065]
    

    # Written by Martin Gomez Garcia Alvestegui(ゴメス　マルティン)
               
    startYear, endYear = datesO[0].year, datesO[-1].year
    CDFO               = np.zeros((12,31*3*(endYear-startYear+1)),dtype = np.float64) #*3 window 3 months   
    startYear, endYear = datesC[0].year, datesC[-1].year
    CDFMC              = np.zeros((12,31*3*(endYear-startYear+1)),dtype = np.float64)
    startYear, endYear = datesT[0].year, datesT[-1].year
    CDFMT              = np.zeros((12,31*3*(endYear-startYear+1)),dtype = np.float64)

    arrLenCO,arrLenC,arrLenT = np.zeros((3,12),dtype = np.int32)

    
    # For precipitation  
    if var == 'pr':

        for im in range(1,13):
                
            # Define the months for calibration
            monBef = 12 if im == 1  else im-1
            monAct = im
            monAft = 1  if im == 12 else im+1    

            indxsCO = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesO])
            indxsC  = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesC])
            arrLenCO[im-1] = len(indxsCO[0])
            arrLenC[im-1]  = len(indxsC[0])

            indxsT  = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesT])
            arrLenT[im-1] = len(indxsT[0])

            # Obs
            tempDaily = np.sort(np.array(dailyO)[indxsCO])
            numZeros = len(np.where(tempDaily<0.05)[0])
            randNum  = np.random.uniform(0.00001,0.05,size=numZeros)
            tempDaily[np.where(tempDaily<0.05)] = np.sort(randNum)
            CDFO[im-1][:arrLenCO[im-1]] = tempDaily    

            # Mod Calib
            tempDaily = np.sort(np.array(dailyMC)[indxsC])
            numZeros = len(np.where(tempDaily<0.05)[0])
            randNum = np.random.uniform(0.00001,0.05,size=numZeros)
            tempDaily[np.where(tempDaily<0.05)] = np.sort(randNum)
            CDFMC[im-1][:arrLenC[im-1]] = tempDaily

            # Mod Target
            tempDaily = np.sort(np.array(dailyMT)[indxsT])
            numZeros = len(np.where(tempDaily<0.05)[0])
            randNum  = np.random.uniform(0.00001,0.05,size = numZeros)
            tempDaily[np.where(tempDaily<0.05)] = np.sort(randNum)
            CDFMT[im-1][:arrLenT[im-1]] = tempDaily

        #tildeDaily = np.zeros_like(dailyMT)   #Focus on the Future Climatic model
        tildeDaily=np.zeros_like(dailyO)
        #for ii,iDay in enumerate(datesT):   #dates length of the Future Climatic Model
        for ii,iDay in enumerate(datesO):
                
            # Loop over each day in the target period : datesT . Determine this month.
            im = iDay.month
            # Copy the CDFs
            ecdfO  =  CDFO[im-1][:arrLenCO[im-1]]  #Obs
            ecdfMC = CDFMC[im-1][:arrLenC[im-1]]
            ecdfMT = CDFMT[im-1][:arrLenT[im-1]]

            # Create a probability array for the calibration period (pCO;pC) and target period (pT)
            pCO = (np.arange(arrLenCO[im-1])-0.5)/arrLenCO[im-1] #Obs
            pC  = (np.arange(arrLenC[im-1])-0.5)/arrLenC[im-1]  #Current Climate
            pT  = (np.arange(arrLenT[im-1])-0.5)/arrLenT[im-1]  #Fut climate

            # Calculate the absolute change of the quantile (qDelta)
            thisQuantile = np.interp(dailyMT[ii],ecdfMT,pT)
            daily_MC     = np.interp(thisQuantile,pC,ecdfMC)
            #plt.plot(thisQuantile)
            #plt.show()
            #plt.clf()

            qDelta       = dailyMT[ii] / daily_MC    
            print(round(qDelta,3))
            #print(thisQuantile)
            #print(pCO)
            #print(ecdfO)
            # Calculate the corrected value
            tildeDaily[ii] = np.interp(thisQuantile,pCO,ecdfO) * qDelta 
            #print(tildeDaily[ii])
            tildeDaily[np.where(tildeDaily<0.05)] = 0                
    
    
    # For temperature CAL
    else:
        for im in range(1,13):
                
            # Define the months for calibration
            monBef = 12 if im == 1  else im-1
            monAct = im
            monAft = 1  if im == 12 else im+1    

            indxsCO = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesO])
            indxsC  = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesC])
            arrLenCO[im-1] = len(indxsCO[0])
            arrLenC[im-1]  = len(indxsC[0])

            indxsT  = np.where([(dt.month==monBef or dt.month==monAct or dt.month==monAft)for dt in datesT])
            arrLenT[im-1] = len(indxsT[0])

            # Obs
            tempDaily = np.sort(dailyO[indxsCO])
            CDFO[im-1][:arrLenCO[im-1]] = tempDaily    

            # Mod Calib
            tempDaily = np.sort(dailyMC[indxsC])
            CDFMC[im-1][:arrLenC[im-1]] = tempDaily

            # Mod Target
            tempDaily = np.sort(dailyMT[indxsT])
            CDFMT[im-1][:arrLenT[im-1]] = tempDaily

        tildeDaily = np.zeros_like(dailyMT)
        for ii,iDay in enumerate(datesT):
                
            # Loop over each day in the target period : datesT . Determine this month.
            im = iDay.month
            # Copy the CDFs
            ecdfO  =  CDFO[im-1][:arrLenCO[im-1]]
            ecdfMC = CDFMC[im-1][:arrLenC[im-1]]
            ecdfMT = CDFMT[im-1][:arrLenT[im-1]]

            # Create a probability array for the calibration period (pCO;pC) and target period (pT)
            pCO = (np.arange(arrLenCO[im-1])-0.5)/arrLenCO[im-1]
            pC  = (np.arange(arrLenC[im-1])-0.5)/arrLenC[im-1]
            pT  = (np.arange(arrLenT[im-1])-0.5)/arrLenT[im-1]

            # Calculate the absolute change of the quantile (qDelta)
            thisQuantile = np.interp(dailyMT[ii],ecdfMT,pT)
            daily_MC     = np.interp(thisQuantile,pC,ecdfMC)

            qDelta       = dailyMT[ii] - daily_MC           
            # Calculate the corrected value
            tildeDaily[ii] = np.interp(thisQuantile,pCO,ecdfO) + qDelta         

    return tildeDaily
     

#%%
Aldunate_series=QDM('pr',dailyO,dailyMC,dailyMT,datesO,datesC,datesT)
Aldunate_series[4]
obs['QDM_perturbed']=Aldunate_series


plt.plot(obs.QDM_perturbed,label='QDM')
plt.plot(obs[name_station],label='historic obs')
plt.legend()
plt.show()
plt.clf()


#%%
#Seasonality Functions for Precipitation

def Monthly(df):
    mon_dict={1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',
              6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
    df.index=pd.to_datetime(df.index,errors='coerce')
    df_mon=df.groupby([df.index.month,df.index.year],sort=False).sum()
    df_mon.index.rename(['Month','Year'],inplace=True)
    df_mon.reset_index(level=0,inplace=True)
    df_mon.replace({'Month':mon_dict},inplace=True)
    return df_mon

def Season(df):
    mon_dict={1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',
              6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
    df_season=df.groupby([df['Month']],sort=False).mean()
    #df_season.reset_index(level=0,inplace=True)
    df_season=df_season.rename(index=mon_dict)
    return df_season

QDM=pd.DataFrame(Aldunate_series,columns=['QDM'])
QDM.index=obs.index.copy()



QDM_POT=POT(QDM,'QDM')


plt.scatter(Obs_climate['Return_Period'],Obs_climate[name_station],s=10,label='Observations')
plt.scatter(QDM_POT['Return_Period'],QDM_POT['QDM'],s=10,label='QDM')
plt.xscale('log')
plt.title('Quantiles Delta Mapping'+' '+name_station)
plt.legend()
plt.show()
plt.clf()

mon_obs=Monthly(obs)
season_obs=Season(mon_obs)

plt.plot(season_obs.QDM_perturbed,label='QDM')
ax=sns.lineplot(data=mon_obs,x='Month',y=str(name_station),
                legend=False, lw=1.5,color='black',label='Historical obs')
plt.plot(season_obs[name_station],label='historic obs')
plt.legend()
plt.title(name_station)
plt.show()
plt.clf()



##Output file to excel Matrix

name=str(i)+'_'+str(name_station)
df_salida[name]=obs['QDM_perturbed']
#df_salida=obs.QDM_perturbed.copy().to_frame()
#df_salida=df_salida.rename(columns={'QDM_perturbed':'QDM_'+str(name_station)})

print('station number is: '+ str(i))

os.chdir(salida)
#df_salida.to_excel(model_name+'_Perturbed.xlsx')
print('station number is: '+ str(i)+' and output success!')
