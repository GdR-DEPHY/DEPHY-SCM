#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 9 November 2022

@author: Catherine Rio

Modification
  2022/11/08, C. Rio: EUROCS case
  2024/02/05, F. Couvreux: adaptation au dernier format
  2026/06/02, N. Villefranque: clean for publication
"""

## EUROCS original case definition
## From  https://rmets.onlinelibrary.wiley.com/doi/abs/10.1256/qj.03.145

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=96
tmin = datetime(1997, 6, 27, 11, 30)
tmax = tmin + timedelta(hours=duration) #  "19970701113000"

case = Case('EUROCS/REF',
        lat=36.61,
        lon=-97.49,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=360.)

case.set_title("Forcing and initial conditions for EUROCS case - Original definition")
case.set_reference("Guichard et al. (2004, QJRMS)")
case.set_author("C. Rio")
case.set_script("DEPHY-SCM/EUROCS/REF/driver_DEF.py")
case.set_comment("https://rmets.onlinelibrary.wiley.com/doi/abs/10.1256/qj.03.145")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97285.89
case.add_init_ps(ps)

# Pressure, potential temperature, specific humidity profiles from aux file
P,th,rv = np.genfromtxt('../aux/init_thrv.txt', dtype=float, skip_header=0, usecols=[0,1,2]).transpose()
case.add_init_theta(th, lev=P, levtype='pressure')
case.add_init_rv(rv, lev=P, levtype='pressure')

# Pressure, horizontal winds profiles from aux file
Pu,u,v = np.genfromtxt('../aux/init_uv.txt', dtype=None, skip_header=0, usecols=[0,1,2]).transpose()
case.add_init_wind(u=u,v=v, lev=Pu, levtype='pressure')

################################################
# 3. Forcing
################################################

# Advection forcing (+ radiative tendency)
#       t (s), A_theta (K hour-1) R_theta (K hour-1) A_rt (g kg-1 hour-1)
nk=20
yyyy,mm,dd,timess = np.genfromtxt('../aux/Forc_Tq_sfc.txt', dtype=float, skip_header=0, usecols=[0,1,2,3]).transpose()
nt = len(timess)
Timeforc=np.zeros(nt)
for it, (y,m,d,t) in enumerate(zip(yyyy,mm,dd,timess)):
    dateforc = datetime(int(y),int(m),int(d))+timedelta(hours=t/3600.)
    if it==0:
        dateinit=dateforc
    difftime=dateforc-dateinit
    Timeforc[it]=difftime.total_seconds()
# !!!!!!!!!!!!!!!!!!!!
# !!!! EN COURS !!!!!!
# !!!!!!!!!!!!!!!!!!!!
print(Timeforc)
exit()
alt=np.genfromtxt('Forc_Tq_sfc.txt',dtype=None,skip_header=0, usecols=4)
psforc=np.genfromtxt('Forc_Tq_sfc.txt',dtype=None,skip_header=0, usecols=5)
thsforc=np.genfromtxt('Forc_Tq_sfc.txt',dtype=None,skip_header=0, usecols=6)
rvsforc=np.genfromtxt('Forc_Tq_sfc.txt',dtype=None,skip_header=0, usecols=7)

ntimes=len(Timeforc)

Pfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=0)
Pfrck = np.array(Pfrck,dtype=np.float64)
Pfrck = np.reshape(Pfrck, (ntimes,nk))
Pfrck1D=Pfrck[0,:]
print('Pfrc 1D=',Pfrck1D.shape,'nk=',nk)
ufrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=1)
ufrck = np.array(ufrck,dtype=np.float64)
ufrck = np.reshape(ufrck, (ntimes,nk))
vfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=2)
vfrck = np.array(vfrck,dtype=np.float64)
vfrck = np.reshape(vfrck, (ntimes,nk))
thfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=3)
thfrck = np.array(thfrck,dtype=np.float64)
thfrck = np.reshape(thfrck, (ntimes,nk))
rvfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=4)
rvfrck = np.array(rvfrck,dtype=np.float64)
rvfrck = np.reshape(rvfrck, (ntimes,nk))
dthdtfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=6)
dthdtfrck = np.array(dthdtfrck,dtype=np.float64)
dthdtfrck = np.reshape(dthdtfrck, (ntimes,nk))
drvdtfrck=np.genfromtxt('Forc_Tq.txt',dtype=None,skip_header=0, usecols=7)
drvdtfrck = np.array(drvdtfrck,dtype=np.float64)
drvdtfrck = np.reshape(drvdtfrck, (ntimes,nk))

##################################################################################
#ATTENTION POUR L INSTANT POUR ARCHIVE ON A GARDE l INFO SUR LA tENDANCE RADIATIVE
#ISSU DE DEFINITION DU CAS EUROCS DANS DE VIEILLES VERSION DE LMDZ
# ON NE L UTILISE PAS DANS LA DEFINITION DU CAS CAR LE CAS TOURNE AVEC RAYONNEMENT INTERACTIF
# PBM aussi que rad_T est défini sur 18 niveaux alors que les forcages précédents
#sont fournis sur 20 niveaux
###############################################################################"
rad_T = [-.3867370E+01,  -.9437038E+00,   .3279568E+01,   .3491005E+01,   .3179960E+00,
  -.3563539E+01,  -.4864929E+01,  -.4350163E+01,  -.3867370E+01,  -.9437038E+00,
   .3279568E+01,   .3491005E+01,   .3179960E+00,  -.3563539E+01,  -.4864929E+01,
  -.4350163E+01,  -.3867370E+01, -.9437038E+00,   .3279568E+01,   .3491005E+01,   .3179960E+00,
  -.3563539E+01,  -.4864929E+01,  -.4350163E+01,  -.3867370E+01,  -.9437038E+00,
   .3279568E+01,   .3491005E+01,   .3179960E+00,  -.3563539E+01,  -.4864929E+01,
  -.4350163E+01,  -.3867370E+01,
  -.3087036E+01,  -.8644809E+00,   .1899312E+01,   .1892331E+01,  -.8682765E-01,
  -.1842751E+01,  -.3281245E+01,  -.3491231E+01,  -.3087036E+01,  -.8644809E+00,
   .1899312E+01,   .1892331E+01,  -.8682765E-01,  -.1842751E+01,  -.3281245E+01,
  -.3491231E+01,  -.3087036E+01, -.8644809E+00,   .1899312E+01,   .1892331E+01,  -.8682765E-01,
  -.1842751E+01,  -.3281245E+01,  -.3491231E+01,  -.3087036E+01,  -.8644809E+00,
   .1899312E+01,   .1892331E+01,  -.8682765E-01,  -.1842751E+01,  -.3281245E+01,
  -.3491231E+01,  -.3087036E+01,
  -.2838080E+01,  -.9134414E+00,   .1061002E+01,   .1546524E+01,  -.1895735E+00,
  -.1981932E+01,  -.2895933E+01,  -.3198640E+01,  -.2838080E+01,  -.9134414E+00,
   .1061002E+01,   .1546524E+01,  -.1895735E+00,  -.1981932E+01,  -.2895933E+01,
  -.3198640E+01,  -.2838080E+01, -.9134414E+00,   .1061002E+01,   .1546524E+01,  -.1895735E+00,
  -.1981932E+01,  -.2895933E+01,  -.3198640E+01,  -.2838080E+01,  -.9134414E+00,
   .1061002E+01,   .1546524E+01,  -.1895735E+00,  -.1981932E+01,  -.2895933E+01,
  -.3198640E+01,  -.2838080E+01,
  -.2796604E+01,  -.6982704E+00,   .1151711E+01,   .8755094E+00,  -.3523057E+00,
  -.1913864E+01,  -.2814059E+01,  -.3134714E+01,  -.2796604E+01,  -.6982704E+00,
   .1151711E+01,   .8755094E+00,  -.3523057E+00,  -.1913864E+01,  -.2814059E+01,
  -.3134714E+01,  -.2796604E+01, -.6982704E+00,   .1151711E+01,   .8755094E+00,  -.3523057E+00,
  -.1913864E+01,  -.2814059E+01,  -.3134714E+01,  -.2796604E+01,  -.6982704E+00,
   .1151711E+01,   .8755094E+00,  -.3523057E+00,  -.1913864E+01,  -.2814059E+01,
  -.3134714E+01,  -.2796604E+01,
  -.2680587E+01,  -.5366148E+00,   .1186146E+01,   .6456953E+00,  -.1644712E+01,
  -.3054163E+01,  -.3096808E+01,  -.3154263E+01,  -.2680587E+01,  -.5366148E+00,
   .1186146E+01,   .6456953E+00,  -.1644712E+01,  -.3054163E+01,  -.3096808E+01,
  -.3154263E+01,  -.2680587E+01, -.5366148E+00,   .1186146E+01,   .6456953E+00,  -.1644712E+01,
  -.3054163E+01,  -.3096808E+01,  -.3154263E+01,  -.2680587E+01,  -.5366148E+00,
   .1186146E+01,   .6456953E+00,  -.1644712E+01,  -.3054163E+01,  -.3096808E+01,
  -.3154263E+01,  -.2680587E+01,
  -.2442485E+01,  -.2916573E+00,   .1181932E+01,   .7098126E+00,  -.2031150E+01,
  -.3241811E+01,  -.3113050E+01,  -.2991056E+01,  -.2442485E+01,  -.2916573E+00,
   .1181932E+01,   .7098126E+00,  -.2031150E+01,  -.3241811E+01,  -.3113050E+01,
  -.2991056E+01,  -.2442485E+01, -.2916573E+00,   .1181932E+01,   .7098126E+00,  -.2031150E+01,
  -.3241811E+01,  -.3113050E+01,  -.2991056E+01,  -.2442485E+01,  -.2916573E+00,
   .1181932E+01,   .7098126E+00,  -.2031150E+01,  -.3241811E+01,  -.3113050E+01,
  -.2991056E+01,  -.2442485E+01,
  -.2167791E+01,  -.2415075E+00,   .9281422E+00,   .5617926E+00,  -.1065317E+01,
  -.2690809E+01,  -.3053746E+01,  -.2645833E+01,  -.2167791E+01,  -.2415075E+00,
   .9281422E+00,   .5617926E+00,  -.1065317E+01,  -.2690809E+01,  -.3053746E+01,
  -.2645833E+01,  -.2167791E+01, -.2415075E+00,   .9281422E+00,   .5617926E+00,  -.1065317E+01,
  -.2690809E+01,  -.3053746E+01,  -.2645833E+01,  -.2167791E+01,  -.2415075E+00,
   .9281422E+00,   .5617926E+00,  -.1065317E+01,  -.2690809E+01,  -.3053746E+01,
  -.2645833E+01,  -.2167791E+01,
  -.1712053E+01,  -.2189050E+00,   .7260165E+00,   .4900903E+00,  -.5656150E+00,
  -.1544243E+01,  -.2284874E+01,  -.2072463E+01,  -.1712053E+01,  -.2189050E+00,
   .7260165E+00,   .4900903E+00,  -.5656150E+00,  -.1544243E+01,  -.2284874E+01,
  -.2072463E+01,  -.1712053E+01, -.2189050E+00,   .7260165E+00,   .4900903E+00,  -.5656150E+00,
  -.1544243E+01,  -.2284874E+01,  -.2072463E+01,  -.1712053E+01,  -.2189050E+00,
   .7260165E+00,   .4900903E+00,  -.5656150E+00,  -.1544243E+01,  -.2284874E+01,
  -.2072463E+01,  -.1712053E+01,
  -.1406971E+01,  -.1938429E+00,   .7387693E+00,   .5495740E+00,  -.5620997E+00,
  -.1545807E+01,  -.1935528E+01,  -.1930948E+01,  -.1406971E+01,  -.1938429E+00,
   .7387693E+00,   .5495740E+00,  -.5620997E+00,  -.1545807E+01,  -.1935528E+01,
  -.1930948E+01,  -.1406971E+01, -.1938429E+00,   .7387693E+00,   .5495740E+00,  -.5620997E+00,
  -.1545807E+01,  -.1935528E+01,  -.1930948E+01,  -.1406971E+01,  -.1938429E+00,
   .7387693E+00,   .5495740E+00,  -.5620997E+00,  -.1545807E+01,  -.1935528E+01,
  -.1930948E+01,  -.1406971E+01,
  -.1569944E+01,  -.3956247E+00,   .5938859E+00,   .4954197E+00,  -.6265966E+00,
  -.1628802E+01,  -.1902504E+01,  -.1967809E+01,  -.1569944E+01,  -.3956247E+00,
   .5938859E+00,   .4954197E+00,  -.6265966E+00,  -.1628802E+01,  -.1902504E+01,
  -.1967809E+01,  -.1569944E+01, -.3956247E+00,   .5938859E+00,   .4954197E+00,  -.6265966E+00,
  -.1628802E+01,  -.1902504E+01,  -.1967809E+01,  -.1569944E+01,  -.3956247E+00,
   .5938859E+00,   .4954197E+00,  -.6265966E+00,  -.1628802E+01,  -.1902504E+01,
  -.1967809E+01,  -.1569944E+01,
  -.1612240E+01,  -.4640218E+00,   .3844041E+00,   .5113403E+00,  -.3453308E+00,
  -.1097269E+01,  -.1659311E+01,  -.1792680E+01,  -.1612240E+01,  -.4640218E+00,
   .3844041E+00,   .5113403E+00,  -.3453308E+00,  -.1097269E+01,  -.1659311E+01,
  -.1792680E+01,  -.1612240E+01, -.4640218E+00,   .3844041E+00,   .5113403E+00,  -.3453308E+00,
  -.1097269E+01,  -.1659311E+01,  -.1792680E+01,  -.1612240E+01,  -.4640218E+00,
   .3844041E+00,   .5113403E+00,  -.3453308E+00,  -.1097269E+01,  -.1659311E+01,
  -.1792680E+01,  -.1612240E+01,
  -.1374078E+01,  -.3583995E+00,   .4245075E+00,   .7145232E+00,  -.1200377E+00,
  -.1306179E+01,  -.1692640E+01,  -.1740910E+01,  -.1374078E+01,  -.3583995E+00,
   .4245075E+00,   .7145232E+00,  -.1200377E+00,  -.1306179E+01,  -.1692640E+01,
  -.1740910E+01,  -.1374078E+01, -.3583995E+00,   .4245075E+00,   .7145232E+00,  -.1200377E+00,
  -.1306179E+01,  -.1692640E+01,  -.1740910E+01,  -.1374078E+01,  -.3583995E+00,
   .4245075E+00,   .7145232E+00,  -.1200377E+00,  -.1306179E+01,  -.1692640E+01,
  -.1740910E+01,  -.1374078E+01,
  -.1371049E+01,  -.3670744E+00,   .4365064E+00,   .5257573E+00,  -.2752405E+00,
  -.1400615E+01,  -.2153996E+01,  -.2349175E+01,  -.1371049E+01,  -.3670744E+00,
   .4365064E+00,   .5257573E+00,  -.2752405E+00,  -.1400615E+01,  -.2153996E+01,
  -.2349175E+01,  -.1371049E+01, -.3670744E+00,   .4365064E+00,   .5257573E+00,  -.2752405E+00,
  -.1400615E+01,  -.2153996E+01,  -.2349175E+01,  -.1371049E+01,  -.3670744E+00,
   .4365064E+00,   .5257573E+00,  -.2752405E+00,  -.1400615E+01,  -.2153996E+01,
  -.2349175E+01,  -.1371049E+01,
  -.1419775E+01,  -.4842038E+00,   .3604439E+00,   .5668442E+00,  -.4210227E+00,
  -.1250045E+01,  -.2228679E+01,  -.2824064E+01,  -.1419775E+01,  -.4842038E+00,
   .3604439E+00,   .5668442E+00,  -.4210227E+00,  -.1250045E+01,  -.2228679E+01,
  -.2824064E+01,  -.1419775E+01, -.4842038E+00,   .3604439E+00,   .5668442E+00,  -.4210227E+00,
  -.1250045E+01,  -.2228679E+01,  -.2824064E+01,  -.1419775E+01,  -.4842038E+00,
   .3604439E+00,   .5668442E+00,  -.4210227E+00,  -.1250045E+01,  -.2228679E+01,
  -.2824064E+01,  -.1419775E+01,
  -.7939980E+00,  -.2438151E+00,   .8957900E+00,   .1890437E+01,   .2338486E+00,
  -.7301186E+00,  -.1994794E+01,  -.1858190E+01,  -.7939980E+00,  -.2438151E+00,
   .8957900E+00,   .1890437E+01,   .2338486E+00,  -.7301186E+00,  -.1994794E+01,
  -.1858190E+01,  -.7939980E+00, -.2438151E+00,   .8957900E+00,   .1890437E+01,   .2338486E+00,
  -.7301186E+00,  -.1994794E+01,  -.1858190E+01,  -.7939980E+00,  -.2438151E+00,
   .8957900E+00,   .1890437E+01,   .2338486E+00,  -.7301186E+00,  -.1994794E+01,
  -.1858190E+01,  -.7939980E+00,
  -.1698492E-01,   .3160324E+00,   .1086991E+01,   .1564187E+01,   .5056278E+00,
  -.5262386E+00,  -.1222305E+01,  -.9789044E+00,  -.1698492E-01,   .3160324E+00,
   .1086991E+01,   .1564187E+01,   .5056278E+00,  -.5262386E+00,  -.1222305E+01,
  -.9789044E+00,  -.1698492E-01,  .3160324E+00,   .1086991E+01,   .1564187E+01,   .5056278E+00,
  -.5262386E+00,  -.1222305E+01,  -.9789044E+00,  -.1698492E-01,   .3160324E+00,
   .1086991E+01,   .1564187E+01,   .5056278E+00,  -.5262386E+00,  -.1222305E+01,
  -.9789044E+00,  -.1698492E-01,
   .6604483E+00,   .6242185E+00,   .1274897E+01,   .1056411E+01,   .8297247E+00,
   .4579909E+00,  -.2120029E+00,  -.3305153E+00,   .6604483E+00,   .6242185E+00,
   .1274897E+01,   .1056411E+01,   .8297247E+00,   .4579909E+00,  -.2120029E+00,
  -.3305153E+00,   .6604483E+00,  .6242185E+00,   .1274897E+01,   .1056411E+01,   .8297247E+00,
   .4579909E+00,  -.2120029E+00,  -.3305153E+00,   .6604483E+00,   .6242185E+00,
   .1274897E+01,   .1056411E+01,   .8297247E+00,   .4579909E+00,  -.2120029E+00,
  -.3305153E+00,   .6604483E+00,
   .1250477E+01,   .1237773E+01,   .1704626E+01,   .1607207E+01,   .1545491E+01,
   .1177068E+01,   .1890299E+00,   .6328749E-01,   .1250477E+01,   .1237773E+01,
   .1704626E+01,   .1607207E+01,   .1545491E+01,   .1177068E+01,   .1890299E+00,
   .6328749E-01,   .1250477E+01, .1237773E+01,   .1704626E+01,   .1607207E+01,   .1545491E+01,
   .1177068E+01,   .1890299E+00,   .6328749E-01,   .1250477E+01,   .1237773E+01,
   .1704626E+01,   .1607207E+01,   .1545491E+01,   .1177068E+01,   .1890299E+00,
   .6328749E-01,   .1250477E+01]
#rad_T = np.array(rad_T,dtype=np.float64)
#rad_T = np.reshape(rad_T, (18,33))
#case.add_temp_radiation_tendency(rad_T/3600./24.,time=timeF,lev=p,levtype='pressure',levid='lev') # converted in K s-1


case.add_theta_advection(dthdtfrck,time=Timeforc,lev=Pfrck1D,levtype='pressure',include_rad=False) # converted in K s-1
case.add_rv_advection(drvdtfrck,time=Timeforc,lev=Pfrck1D,levtype='pressure') # converted in K s-1


#Wind nudging
#unudg = u_wind
#vnudg = v_wind

#unudg = extend('ua',unudg,init=False,time=timeForc)
#vnudg = extend('va',vnudg,init=False,time=timeForc)

case.add_wind_nudging(unudg=ufrck,vnudg=vfrck,timescale=3600.*2.,p_nudging=110000.,time=Timeforc,timeid='time',lev=Pfrck1D,levtype='pressure',levid='lev')


timeSfc = np.genfromtxt('sfceforcing_time.txt',dtype=None,skip_header=0,usecols=0)
sensSfc = np.genfromtxt('sfceforcing_fluxsens.txt',dtype=None,skip_header=0,usecols=0)
lateSfc = np.genfromtxt('sfceforcing_fluxlat.txt',dtype=None,skip_header=0,usecols=0)
tsSfc = np.genfromtxt('sfceforcing_ts.txt',dtype=None,skip_header=0,usecols=0)

case.add_surface_fluxes(sens=sensSfc,lat=lateSfc,time=timeSfc,forc_wind='z0',z0=0.15)
#case.add_forcing_ts(tsSfc,time=timeSfc,z0=0.15)
case.add_rad_ts(tsSfc,time=timeSfc)

#besoin de rajouter definition de Emissivite=0.994 et Albedo=0.17
alb=0.17
emis=0.994
case.add_forcing_variable('alb',alb)
case.add_forcing_variable('emis',emis)
################################################
# 4. Writing file
################################################

case.write('EUROCS_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
