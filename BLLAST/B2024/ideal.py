import netCDF4 as nc

import numpy as np
from scipy.signal import savgol_filter

from metpy.calc import mixing_ratio_from_relative_humidity, relative_humidity_from_mixing_ratio
from metpy.units import units
import matplotlib.pyplot as plt

import dephycf.thermo as thermo
import dephycf.constants as CC

window = 51
order = 2

def ideal_theta(z):
    #dataI['zh_theta'] =    np.array([0,  65, 85,   150, 500,  1500, 2500, 4000, 10000])
    #dataI['theta']    =    np.array([20, 21, 23.5,  24,  25.5,  32,   35,   42,    62])

    theta = z*0. + np.nan
    theta = np.where(z < 17500., 66   + (z-11000)/(17500-11000)*(155-66), theta)
    theta = np.where(z < 11000., 42   + (z-4000)/(11000-4000)*(66-42), theta)
    theta = np.where(z <  4000., 35   + (z-2500)/(4000-2500)*(42-35), theta)
    theta = np.where(z <  2500., 32   + (z-1500)/(2500-1500)*(35-32), theta)
    theta = np.where(z <  1500., 25.5 + (z-500)/(1500-500)*(32-25.5), theta)
    theta = np.where(z <   500., 24   + (z-150)/(500-150)*(25.5-24), theta)
    theta = np.where(z <   150., 23.5 + (z-85)/(150-85)*(24-23.5), theta)
    theta = np.where(z <    85., 21   + (z-65.)/(85.-65.)*(23.5-21.), theta)
    theta = np.where(z <    65., 20   + (z-0.)/(65.-0.)*(21.-20.), theta)
    
    return theta

def ideal_hur(z):
    #dataI['zh_hur'] = np.array([ 0, 45, 60, 150, 550, 6000, 12000])
    #dataI['hur']    = np.array([78, 72, 64,  58,  50,   50,     0])

    hur = z*0. + np.nan
    hur = np.where(z < 17500.,  0, hur)
    hur = np.where(z < 12000., 50 + (z-6000)/(12000-6000)*(0-50), hur)
    hur = np.where(z <  6000., 50, hur)
    hur = np.where(z <   550., 58 + (z-150)/(550-150)*(50-58), hur)
    #hur = np.where(z <   150., 64 + (z-60)/(150-60)*(58-64), hur)
    hur = np.where(z <   150., 62 + (z-85)/(150-85)*(58-62), hur)
    #hur = np.where(z <    60., 72 + (z-45.)/(60.-45.)*(64-72.), hur)
    hur = np.where(z <    85., 72 + (z-45.)/(85.-45.)*(62-72.), hur)
    hur = np.where(z <    45., 78 + (z-0.)/(45.-0.)*(72.-78.), hur)
    
    return hur

def ideal_ua(z):

    ua = z*0. + np.nan
    ua = np.where(z < 17500., 18 + (z-14000)/(17500-14000)*(5-18.), ua)
    ua = np.where(z < 14000., 18, ua)
    ua = np.where(z <  5000., 15 + (z-3500)/(5000-3500)*(18-15.), ua)
    ua = np.where(z <  3500., 0 + (z-0)/(3500-0)*(15-0.), ua)
    
    return ua

def ideal_va(z):

    va = z*0. + np.nan
    va = np.where(z < 17500.,  0, va)
    va = np.where(z <  5000., -4 + (z-4000)/(5000-4000)*(0+4), va)
    va = np.where(z <  4000., 0 + (z-400)/(4000-400)*(-4-0), va)
    va = np.where(z <   400., 0.5 + (z-0)/(400-0)*(0-0.5), va)
    
    return va


data = {}

data['zh'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=1)
data['ua'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=4)
data['va'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=5)
data['ta'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=11)
data['hur'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=14)
data['pa'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.cor',dtype=None,skip_header=1,skip_footer=1,usecols=15)

data['theta'] = (data['ta']+273.15)*(1000/data['pa'])**(2./7.)-273.15

p_loc = data['pa']*units.hPa
ta_loc = data['ta']*units.degC
hur_loc = data['hur']/100.
nlev, = p_loc.shape
data['rv'] = np.array([mixing_ratio_from_relative_humidity(p_loc[i], ta_loc[i], hur_loc[i]) for i in range(nlev)])

data['qv'] = data['rv']/(1+data['rv'])

data['rv'] *= 1000.
data['qv'] *= 1000.

data2 = {}
for var in ['ta', 'theta', 'qv', 'rv', 'hur', 'ua', 'va']: 
    data2[var] = savgol_filter(data[var], window, order)

data_fil = {}

data_fil['zh'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=1)
data_fil['ua'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=4)
data_fil['va'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=5)
data_fil['ta'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=10)
data_fil['hur'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=13)
data_fil['pa'] = np.genfromtxt('data/RS_20110620_0515_site1_MODEM_CRA.fil',dtype=None,skip_header=1,skip_footer=1,usecols=14)

data_fil['theta'] = (data_fil['ta']+273.15)*(1000/data_fil['pa'])**(2./7.)-273.15

p_loc = data_fil['pa']*units.hPa
ta_loc = data_fil['ta']*units.degC
hur_loc = data_fil['hur']/100.
nlev, = p_loc.shape
data_fil['rv'] = np.array([mixing_ratio_from_relative_humidity(p_loc[i], ta_loc[i], hur_loc[i]) for i in range(nlev)])

data_fil['qv'] = data_fil['rv']/(1+data_fil['rv'])

data_fil['rv'] *= 1000.
data_fil['qv'] *= 1000.



data3 = {}
with nc.Dataset('../NOADV/BLLAST_NOADV_SCM_driver.nc') as ds:
    for var in ['zh','ua','va','ta','theta','rv','pa']:
        data3[var] = np.squeeze(ds[var][:])

data3['qv'] = data3['rv']/(1+data3['rv'])

p_loc = data3['pa']*units.Pa
ta_loc = data3['ta']*units.K
rv_loc = data3['rv']*1.
nlev, = p_loc.shape
data3['hur'] = np.array([relative_humidity_from_mixing_ratio(p_loc[i], ta_loc[i], rv_loc[i]) for i in range(nlev)])

data3['ta'] -= 273.15
data3['theta'] -= 273.15
data3['rv'] *= 1000.
data3['qv'] *= 1000.
data3['hur'] *= 100.

for var in ['ta','theta','rv','qv','hur','ua','va']:
    data3[var] = np.where(data3['zh'] > 3000., np.nan, data3[var])


data4 = {}
with nc.Dataset('../NOADV_MOIST/BLLAST_NOADV_MOIST_SCM_driver.nc') as ds:
    for var in ['zh','ua','va','ta','theta','rv','pa']:
        data4[var] = np.squeeze(ds[var][:])

data4['qv'] = data4['rv']/(1+data4['rv'])

p_loc = data4['pa']*units.Pa
ta_loc = data4['ta']*units.K
rv_loc = data4['rv']*1.
nlev, = p_loc.shape
data4['hur'] = np.array([relative_humidity_from_mixing_ratio(p_loc[i], ta_loc[i], rv_loc[i]) for i in range(nlev)])

data4['ta'] -= 273.15
data4['theta'] -= 273.15
data4['rv'] *= 1000.
data4['qv'] *= 1000.
data4['hur'] *= 100.

for var in ['ta','theta','rv','qv','hur','ua','va']:
    data4[var] = np.where(data4['zh'] > 3000., np.nan, data4[var])

data5 = {}
with nc.Dataset('BLLAST_B2024_SCM_driver.nc') as ds:
    for var in ['zh','ua','va','ta','theta','rv','pa']:
        data5[var] = np.squeeze(ds[var][:])

data5['qv'] = data5['rv']/(1+data5['rv'])

p_loc = data5['pa']*units.Pa
ta_loc = data5['ta']*units.K
rv_loc = data5['rv']*1.
nlev, = p_loc.shape
data5['hur'] = np.array([relative_humidity_from_mixing_ratio(p_loc[i], ta_loc[i], rv_loc[i]) for i in range(nlev)])

data5['ta'] -= 273.15
data5['theta'] -= 273.15
data5['rv'] *= 1000.
data5['qv'] *= 1000.
data5['hur'] *= 100.

header = ['year' ,'month' ,'day' ,'h30' ,'mm30' ,'offsetc' ,'gainc' ,'offsetq' ,'gainq' ,'hveg' ,'d' ,'zo' ,'temp2' ,'huc2' ,'temp15' ,'huc15' ,'temp30' ,'huc30' ,'temp45' ,'huc45' ,'temp60' ,'huc60' ,'r2' ,'r15' ,'r30' ,'r45' ,'r60' ,'q' ,'lic' ,'kh20' ,'oz' ,'CO2' ,'Ps' ,'Plic' ,'rho' ,'Rn' ,'offseto' ,'gaino' ,'FF15' ,'FF30' ,'FFG45' ,'FF45' ,'FF60' ,'DD15' ,'DD30' ,'DDG45' ,'DD45' ,'DD60' ,'G1' ,'G2' ,'G3' ,'rain']
tmp = np.array([2011,6,20,5,15,0,1,0,1,0.05,np.nan,np.nan,16.273,75.025,15.685,76.133,15.788,73.885,16.197,72.264,17.197,64.721,9.2226,9.0252,8.8503,8.8999,8.5017,np.nan,512.02,-0.042409,15.038,15.984,950.36,945,1.1374,-6.989,0,1,0.26057,1.1965,2.2675,2.2586,2.9016,184.37,190.65,176.96,177.88,190.63,-29.937,-30.233,-32.414,0])
dataM = {}
dataM['zh'] = [2, 15, 30, 45, 60]
for var in ['temp','huc','r','FF','DD']:
    dataM[var] = []
    for height in dataM['zh']:
        var_loc = f'{var}{height}'
        ind = [s == var_loc for s in header]
        if np.sum(ind) == 1:
            dataM[var].append(float(tmp[ind]))
        else:
            dataM[var].append(np.nan)

    dataM[var] = np.array(dataM[var])

dataM['ta'] = dataM['temp']
dataM['hur'] = dataM['huc']
dataM['rv'] = dataM['r']
dataM['qv'] = dataM['rv']/(1+dataM['rv']/1000.)

dataM['ua'] = -1*dataM['FF']*np.sin(np.deg2rad(dataM['DD']))
dataM['va'] = -1*dataM['FF']*np.cos(np.deg2rad(dataM['DD']))

sites = ['grass','wheat','corn','forest']
dataS = {}

site = 'grass'
dataS[site] = {}
dataS[site]['ta'] = 288.07000732421875
dataS[site]['qv'] = 0.009695500135421753
dataS[site]['ua'] = -0.21229732036590576
dataS[site]['va'] = 0.5204781293869019
dataS[site]['rho'] = 1.145300030708313


site = 'wheat'
dataS[site] = {}
dataS[site]['ta'] = 288.9800109863281
dataS[site]['qv'] = 0.010250000283122063
dataS[site]['ua'] = -0.33480384945869446
dataS[site]['va'] = 0.4236309826374054
dataS[site]['rho'] = 1.1401000022888184


site = 'corn'
dataS[site] = {}
dataS[site]['ta'] = 291.6400146484375
dataS[site]['qv'] = 0.009568300098180771
dataS[site]['ua'] = -0.1659073680639267
dataS[site]['va'] = 2.1270394325256348
dataS[site]['rho'] = 1.1233999729156494


site = 'forest'
dataS[site] = {}
dataS[site]['ta'] = 288.6199951171875
dataS[site]['qv'] = 0.009139000438153744
dataS[site]['ua'] = 0.07069340348243713
dataS[site]['va'] = 2.3138201236724854
dataS[site]['rho'] = 1.135200023651123

for site in sites:
    dataS[site]['rv'] = dataS[site]['qv']/(1-dataS[site]['qv'])
    dataS[site]['pa'] = dataS[site]['rho'] * (dataS[site]['qv']*CC.Rv + (1-dataS[site]['qv'])*CC.Rd)*dataS[site]['ta']
    p_loc = dataS[site]['pa']*units.Pa
    ta_loc = dataS[site]['ta']*units.K
    rv_loc = dataS[site]['rv']*1.
    dataS[site]['hur'] = float(relative_humidity_from_mixing_ratio(p_loc, ta_loc, rv_loc))*100.
    dataS[site]['rv'] *= 1000.
    dataS[site]['qv'] *= 1000.
    dataS[site]['theta'] = dataS[site]['ta']*(100000./dataS[site]['pa'])**(2./7.)-273.15
    dataS[site]['ta'] -= 273.15


colors = {}
colors['grass'] = 'lawngreen'
colors['wheat'] = 'peru'
colors['corn'] = 'gold'
colors['forest'] = 'forestgreen'

dataI = {}
dataI['zh'] = np.linspace(0,17500,17501)
dataI['theta'] = ideal_theta(dataI['zh'])+273.15
dataI['hur'] = ideal_hur(dataI['zh'])
dataI['ua'] = ideal_ua(dataI['zh'])
dataI['va'] = ideal_va(dataI['zh'])

#dataI['zh_theta'] =    np.array([0,  65, 85,   150, 500,  1500, 2500, 4000, 10000])
#dataI['theta']    =    np.array([20, 21, 23.5,  24,  25.5,  32,   35,   42,    62])

#dataI['zh_hur'] = np.array([ 0, 45, 60, 150, 550, 6000, 12000])
#dataI['hur']    = np.array([78, 72, 64,  58,  50,   50,     0])

dataI['qv'] = dataI['theta']*0.
niter = 20
for i in range(niter):
    print('iteration #',i)
    dataI['pa'] = thermo.z2p(theta=dataI['theta'], z=dataI['zh'],ps=95000, qv=dataI['qv'])
    dataI['ta'] = thermo.theta2t(theta=dataI['theta'],p=dataI['pa'])

    p_loc = dataI['pa']*units.Pa
    ta_loc = dataI['ta']*units.K
    hur_loc = dataI['hur']/100.
    nlev, = p_loc.shape
    dataI['rv'] = np.array([mixing_ratio_from_relative_humidity(p_loc[i], ta_loc[i], hur_loc[i]) for i in range(nlev)])

    dataI['qv'] = dataI['rv']/(1+dataI['rv'])

dataI['rv'] *= 1000.
dataI['qv'] *= 1000.
dataI['ta'] -= 273.15
dataI['theta'] -= 273.15

dataERA5 = {}
with nc.Dataset('../aux/ERA5/ERA5_P2OA_20110620000000-20110620230000.nc') as ds:
    dates = nc.num2date(ds['time'][:], ds['time'].units, calendar=ds['time'].calendar)
    ind = [d.hour == 5 for d in dates]
    dataERA5['pa'] = ds['plev'][:]
    for var in ['zg','ta','qv','ua','va']:
        dataERA5[var] = np.squeeze(ds[var][ind])

dataERA5['rv'] = dataERA5['qv']/(1-dataERA5['qv'])
p_loc = dataERA5['pa']*units.Pa
ta_loc = dataERA5['ta']*units.K
rv_loc = dataERA5['rv']*1.
nlev, = p_loc.shape
dataERA5['hur'] = [float(relative_humidity_from_mixing_ratio(p_loc[ilev], ta_loc[ilev], rv_loc[ilev]))*100. for ilev in range(nlev)]
dataERA5['theta'] = (dataERA5['ta'])*(100000/dataERA5['pa'])**(2./7.)

dataERA5['rv'] *= 1000.
dataERA5['qv'] *= 1000.
dataERA5['ta'] -= 273.15
dataERA5['theta'] -= 273.15
dataERA5['zh'] = dataERA5['zg']

for config in ['low','low2','bl','mid','all']:

    print('config:', config)

    if config == 'low':
        hmax = 200.
        ta_range = (10,20)
        theta_range = (15,25)
        rv_range = (7,11)
        hur_range = (50,90)
        ua_range = (-10,10)
        va_range = (-10,10)
    if config == 'low2':
        hmax = 1000.
        ta_range = (10,20)
        theta_range = (15,30)
        rv_range = (5,11)
        hur_range = (0,110)
        ua_range = (-10,10)
        va_range = (-10,10)
    elif config == 'bl':
        hmax = 2500.
        ta_range = (5,20)
        theta_range = (15,40)
        rv_range = (3,11)
        hur_range = (0,110)
        ua_range = (-10,20)
        va_range = (-10,10)
    elif config == 'mid':
        hmax = 6000.
        ta_range = (-20,20)
        theta_range = (15,60)
        rv_range = (0,11)
        hur_range = (0,110)
        ua_range = (-10,30)
        va_range = (-10,10)
    elif config == 'all':
        hmax = 20000.
        ta_range = (-80,20)
        theta_range = (15,160)
        rv_range = (0,11)
        hur_range = (0,110)
        ua_range = (-10,30)
        va_range = (-10,10)
    elif config == 'bernard':
        hmax = 4000.
        ta_range = (268-273.15,293-273.15)
        theta_range = (15,50)
        rv_range = (0,10.5)
        hur_range = (0,110)
        ua_range = (-10,25)
        va_range = (-10,10)


    fig, axs = plt.subplots(nrows=1, ncols=6, sharey=True, figsize=(20,4))

    for ivar,var in enumerate(['ta', 'theta', 'rv', 'hur', 'ua', 'va']):
        axs[ivar].plot(data[var], data['zh']-data['zh'][0], color='lightgrey', label='raw')
        axs[ivar].plot(data_fil[var], data_fil['zh']-data_fil['zh'][0], color='k', label='filtered')
        axs[ivar].plot(dataERA5[var], dataERA5['zh']-data_fil['zh'][0], color='k', linestyle='--', label='ERA5')
        #axs[ivar].plot(data2[var], data['zh']-data['zh'][0], color='b', label='filtered')
        
        if var in ['ta','rv','hur','ua','va']:
            axs[ivar].plot(dataM[var], dataM['zh']-data4['zh'][0], color='salmon', label='60-m Mast')
        for site in sites:
            axs[ivar].plot(dataS[site][var], 2, color = colors[site], marker='o', linestyle='None', label=site)

        axs[ivar].plot(data3[var], data3['zh']-data3['zh'][0], color='r', label='Darbieu et al. (2015)')
        #axs[ivar].plot(data4[var], data4['zh']-data4['zh'][0], color='cyan', label='NOADV_MOIST')
        axs[ivar].plot(dataI[var], dataI['zh']-dataI['zh'][0], color='blue', label='IDEAL')
        axs[ivar].plot(data5[var], data5['zh']-data5['zh'][0], color='cyan', label='B2024')
            
                
        axs[ivar].grid(color='grey', linestyle='--')

    axs[0].set_xlabel('Temperature [°C]')
    axs[0].set_xlim(*ta_range)
    axs[1].set_xlabel(r'$\theta$ [°C]')
    axs[1].set_xlim(*theta_range)
    axs[2].set_xlabel(r'$r_v$ [g kg$^{-1}$]')
    axs[2].set_xlim(*rv_range)
    axs[3].set_xlabel('Relative humidity [%]')
    axs[3].set_xlim(*hur_range)
    if hur_range[1] > 100:
        axs[3].axvline(100, color='grey', linestyle='--')
    axs[4].set_xlabel(r'Zonal wind [m s$^{-1}$]')
    axs[4].set_xlim(*ua_range)
    axs[5].set_xlabel(r'Meridional wind [m s$^{-1}$]')
    axs[5].set_xlim(*va_range)


    axs[0].set_ylabel('Altitude above ground [m]')
    axs[0].set_ylim(0,hmax)
    axs[5].legend(bbox_to_anchor=(1,0), loc="lower left")

    fig.suptitle(f'Radiosounding for 2011-06-20 05:15:00')


    plt.savefig(f'images/radiosounding_iter{niter:0>2}_{config}.png')
    plt.close()
