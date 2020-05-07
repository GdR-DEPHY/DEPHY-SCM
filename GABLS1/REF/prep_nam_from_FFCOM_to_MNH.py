#-*- coding:UTF-8 -*-
# utilisation du module netCDF4 pour la lecture et l'ecriture de fichiers netcdf3
import matplotlib as mpl
# "backend" pour sortie fichier
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from itertools import izip,count
#import cdms2
#import cdtime

# pour enlever les dims = 1 dans un tableau
#new_array=np.squeeze(array)

casname='ARMCUB'
casname='SANDUFAST'
casname='SANDUSLOW'
casname='ARMCUSCM'
casname='RICOSCM'
casname='SANDUREF'
if casname =='ARMCUDEF':
    file='ARMCU_REF_DEF_driver.nc'
    DX=100.
    DY=100.
    DZ=100.
    NAM1='ARMCU'
    NAM2='DEFT1'
    CLOUDSC='ICE3'
    NI=10
    NJ=10
    NZ=40
    keyword='ZUVTHDMR'
    zmax=5000.
    
elif casname =='ARMCUSCM':
    file='ARMCU_REF_SCM_driver.nc'
    DX=100.
    DY=100.
    DZ=100.
    NAM1='ARMCU'
    NAM2='SCMT1'
    CLOUDSC='ICE3'
    NI=10
    NJ=10
    NZ=40
    keyword='ZUVTHDMR'
    zmax=5000.
     
elif casname =='SANDUREF':
    file='SANDU_REF_SCM_driver.nc'
    DX=35.
    DY=35.
    DZ=100.
    NAM1='SANDU'
    NAM2='REFT1'
    CLOUDSC='KHKO'
    NI=256
    NJ=256
    NZ=110
    keyword='ZUVTHLMR'
    zmax=3300.

elif casname =='SANDUFAST':
    file='Composite_FAST_driver_RR.nc'
    DX=35.
    DY=35.
    DZ=100.
    NAM1='SANDU'
    NAM2='FAST1'
    CLOUDSC='KHKO'
    NI=256
    NJ=256
    NZ=110
    keyword='ZUVTHLMR'
    zmax=5000.

elif casname =='SANDUSLOW':
    file='Composite_SLOW_driver_RR.nc'
    DX=35.
    DY=35.
    DZ=100.
    NAM1='SANDU'
    NAM2='SLOW1'
    CLOUDSC='KHKO'
    NI=256
    NJ=256
    NZ=110
    keyword='ZUVTHLMR'
    zmax=5000.

elif casname =='RICOSCM':
    file='RICO_SHORT_SCM_driver.nc'
    DX=25.
    DY=25.
    DZ=25.
    NAM1='RICON'
    NAM2='SCMT1'
    CLOUDSC='ICE3'
    NI=100
    NJ=100
    NZ=130
    keyword='ZUVTHDMR'
    zmax=5000.

f = netCDF4.Dataset(file,'r')
liste_var=f.variables
print('liste des variables du fichier driver:')
for nomvar in liste_var:
    print(nomvar)
startDate = f.startDate
year = int(startDate[0:4])
month = int(startDate[4:6])
day = int(startDate[6:8])
hour = int(startDate[8:10])
minute = int(startDate[10:12])
second = int(startDate[12:14])
timedeb=hour*3600.+minute*60.+second
print('timedeb',timedeb)
endDate = f.endDate
tunits0 = 'seconds since ' + str(year) + '-' + str(month) + '-' + str(day) + ' ' + str(hour) + ':' + str(minute) + ':0.0'

seconds = minute*60 + hour*3600

yearf = int(endDate[0:4])
monthf = int(endDate[4:6])
dayf = int(endDate[6:8])
hourf = int(endDate[8:10])
minutef = int(endDate[10:12])
#besoin de faire facilement la diff des 2 temps

attributes = {}
for att in ['adv_temp','adv_qt','adv_thetal','adv_rv','adv_rt','adv_theta','adv_qv','rad_temp','rad_theta','rad_thetal','forc_omega','forc_w','forc_geo','nudging_temp','nudging_qv','nudging_u','nudging_v','nudging_theta','nudging_thetal','nudging_qt','nudging_rv','nudging_rt','zorog','z0','surfaceType','surfaceForcing','surfaceForcingWind']:
  attributes[att] = 0
for att in f.ncattrs():
  attributes[att] = getattr(f,att)  
  print(att,attributes[att])
# champs qui ne varie pas en fonction du temps
time=f.variables['time'][:]
nt=time.shape[0]
lat = f.variables['lat'][0]
lon =f.variables['lon'][0]
theta=f.variables['theta'][:]
pressure_f=f.variables['pressure'][:]
height=f.variables['height'][:]
u=f.variables['u'][:]
v=f.variables['v'][:]
qv=f.variables['qv'][:]
qt=f.variables['qt'][:]
rv=f.variables['rv'][:]
rt=f.variables['rt'][:]
temp=f.variables['temp'][:]
ql=f.variables['ql'][:]
qi=f.variables['qi'][:]
ps=f.variables['ps'][:]
if (keyword=='ZUVTHLMR'):
    thetal=theta-(100000./pressure_f)**0.286*2.5E6/1004.25*ql
#### A METTRE EN COMMENTAIRE
###omega=f.variables['omega'][:]
###w=f.variables['w'][:]
############################"
print('lat',str(lat))
variables3D = []
#champs qui varie en fonction du temps
#for nomvar in liste_var:
#    variables3D = f.variables[nomvar][:]
ps_forc=f.variables['ps_forc'][:]
height_forc=f.variables['height_forc'][:]
pressure_forc=f.variables['pressure_forc'][:]
ts=f.variables['ts'][:]
print('lecture ug vg',attributes['forc_geo'])
if attributes['forc_geo']==1:
    ug=f.variables['ug'][:]
    vg=f.variables['vg'][:]
print('lecture u_nudg',attributes['nudging_u'])   
if attributes['nudging_u']>0:
    unudg=f.variables['u_nudging'][:]
print('lecture v_nudg',attributes['nudging_v'] )  
if attributes['nudging_v']>0:
    vnudg=f.variables['v_nudging'][:]
print('lecture t_nudg',attributes['nudging_theta'])   
if attributes['nudging_theta']>0:
    thnudg=f.variables['theta_nudging'][:]
print('lecture rv_nudg', attributes['nudging_rv']  ) 
if attributes['nudging_rv']>0:
    qnudg=f.variables['rv_nudging'][:]    
print('lecture theta_adv', attributes['adv_theta']  ) 
if attributes['adv_theta']==1:
    thadv=f.variables['theta_adv'][:]
print('lecture rad_theta',attributes['rad_theta'])
if attributes['rad_theta']==1:
    thadv=thadv+f.variables['theta_rad'][:]
if attributes['rad_theta']=='adv':
    print('le terme de rayonnement est déjà pris en compte dans l advection')
print('lecture rv_adv', attributes['adv_rv'])   
if attributes['adv_rv']==1:
    rvadv=f.variables['rv_adv'][:]
    print('lecture rv_adv')
print('lecture forc_w', attributes['forc_w'])   
if attributes['forc_w']==1:
    w=f.variables['w'][:]
    print('lecture w')
print('lecture forc_omega', attributes['forc_omega'])   
if attributes['forc_omega']==1:
    omega=f.variables['omega'][:]
    print('lecture omega')
startDate = str(f.startDate)
endDate = str(f.endDate)
print('lecture surface_forc', attributes['surfaceForcing'])   
print('lecture surface_forcwind', attributes['surfaceForcingWind'])   
if attributes['surfaceForcing'] == 'ts':
    sst = f.variables['ts'][:]
    print('lecture sst')
elif attributes['surfaceForcing'] == 'surfaceFlux':
  hfls = f.variables['sfc_lat_flx'][:]
  hfss = f.variables['sfc_sens_flx'][:]
  print('lecture flux')
if attributes['surfaceForcingWind'] == 'ustar':
  ustar = f.variables['ustar'][:]
  print('lecture ustar')
f.close()


print(time)
print('time shape=',time.shape)
nt = time.shape[0]
nk = temp.shape[1]
if attributes['surfaceForcing'] == 'ts':
    print(' sst shape=',sst.shape)
print('height=',height.shape)
print('height_forc=',height_forc.shape)
print('nk=',nk)
print('u=',u.shape)
print('qv',qv.shape)

#ii1 = 0
#ii2 = nt-1
#for it in range(0,nt):
#  tt = cdtime.reltime(time[it],time.units)
#  if tt.torel(time.units).value == t0.torel(time.units).value:
#    ii1 = it
#  if tt.torel(time.units).value == tend0.torel(time.units).value:
#    ii2 = it

#nt0 = nt-ii1-(nt-1-ii2)
nt0=nt

g = open('PRE_IDEA1.nam_'+casname,'w')

g.write('&NAM_DIMn_PRE NIMAX='+str(NI)+', NJMAX='+str(NJ)+' /\n')
g.write('&NAM_CONF_PRE LCARTESIAN=.TRUE., NVERB=10,'+'\n')
g.write(' CIDEAL="RSOU",  CZS="FLAT", LFORCING=.TRUE.,'+'\n')
g.write(' NHALO=3,JPHEXT=3,\n')
g.write(' LBOUSS=.FALSE., CEQNSYS="DUR", LPERTURB=.TRUE. /\n')
g.write('&NAM_PERT_PRE CPERT_KIND="WH",XAMPLIWH=0.1 /\n')
g.write('&NAM_CONFn LUSERV=.TRUE.,NSV_USER=3,\n')
if (keyword=='ZUVTHLMR'):
    g.write('LUSERC=.TRUE. /\n')
else :
    g.write('/\n')
g.write('&NAM_GRIDH_PRE XDELTAX='+str(DX)+', XDELTAY='+str(DY)+' /\n')
g.write('&NAM_VER_GRID  LTHINSHELL=.TRUE., NKMAX='+str(NZ)+',\n')
if casname =='SANDUREF' or casname =='SANDUFAST' or casname =='SANDUSLOW':
    g.write(' YZGRID_TYPE="MANUAL"/\n')
else:
    g.write(' ZDZGRD='+str(DZ)+', ZDZTOP='+str(DZ)+',\n')
    g.write(' ZZMAX_STRGRD=1000. , ZSTRGRD=0., ZSTRTOP=0. /\n')
g.write('&NAM_LUNITn CINIFILE="'+casname+'_3D_REF",\n')
g.write('            CINIFILEPGD="'+casname+'_PGD" /\n')
g.write('&NAM_CONFZ MPI_BUFFER_SIZE=400/\n')
g.write('&NAM_LBCn_PRE CLBCX=2*"CYCL", CLBCY=2*"CYCL" /\n') 
g.write('&NAM_GRn_PRE CSURF="EXTE"/\n')

g.write('&NAM_PGD_SCHEMES\n')
daymonth=[31,28,31,30,31,30,31,31,30,31,30,31]
day0=int(day)
month0=int(month)
print('day0=',day0,'month0=',month0,daymonth[month0-1])
timeref=timedeb
if attributes['surfaceType'] == 'ocean' and attributes['surfaceForcing'] == 'ts':
  g.write('  CNATURE = "NONE"  ,\n')
  g.write('  CSEA    = "SEAFLX",\n')
  g.write('  CWATER  = "NONE"  ,\n')
  g.write('  CTOWN   = "NONE"  ,\n')
  g.write('/\n')
  g.write('&NAM_SEABATHY XUNIF_SEABATHY=5. /\n')
  g.write('&NAM_PREP_SEAFLUX XSST_UNIF = '+str(float(ts[0]))+' /\n')
  g.write('&NAM_DATA_SEAFLUX  LSST_DATA=T, \n')
  g.write('NTIME_SST= '+str(nt)+' ,\n')
  for it in range(0,nt):
    if (it >=1):
        timeref=timeref+int(time[it])-int(time[it-1])
    if (timeref >= 86400):
        timeref=timeref-86400
        day0=day0+1
        if (day0 >daymonth[month0-1]):
            day0=day0-daymonth[month0-1]
            month0=month0+1    
    g.write('XUNIF_SST('+str(it+1)+')='+str(float(ts[it])) +',\n')
    g.write('NYEAR_SST('+str(it+1)+')='+str(int(year)) +',\n')
    g.write('NMONTH_SST('+str(it+1)+')='+str(month0) +',\n')
    g.write('NDAY_SST('+str(it+1)+')='+str(day0) +',\n')
    g.write('XTIME_SST('+str(it+1)+')='+str(timeref) +',\n')
  g.write('/\n')
if attributes['surfaceType'] == 'ocean' and attributes['surfaceForcing'] == 'surfaceFlux':
  g.write('  CNATURE = "NONE"  ,\n')
  g.write('  CSEA    = "FLUX" ,\n')
  g.write('  CWATER  = "NONE"  ,\n')
  g.write('  CTOWN   = "NONE"  ,\n')
  g.write('/\n')
if attributes['surfaceType'] == 'land' and attributes['surfaceForcing'] == 'surfaceFlux':
  g.write('  CNATURE = "NONE" ,\n')
  g.write('  CSEA    = "FLUX"  ,\n')
  g.write('  CWATER  = "NONE"  ,\n')
  g.write('  CTOWN   = "NONE"  ,\n')
  g.write('/\n')



g.write('&NAM_GRID_PRE'+' \n')
g.write('  XLAT0 = ' + str(lat) + ',\n')
g.write('  XLON0 = ' + str(lon) + '/\n')

#g.write('&NAM_FRAC')
#g.write('  LECOCLIMAP = F,')
#if surfaceType == 'ocean' or surfaceType == 'land':
#  g.write('  XUNIF_SEA    = 1.,')
#else:
#  g.write('  XUNIF_SEA    = 0.,')
#g.write('  XUNIF_WATER  = 0.,')
#g.write('  XUNIF_TOWN   = 0.,')
#g.write('  XUNIF_NATURE = 0.,')
#if surfaceType == 'land':
#  g.write('  XUNIF_NATURE = 1.,")
#else:
#  g.write('  XUNIF_NATURE = 0.,")
#g.write('/')

g.write('&NAM_COVER'+' \n')
g.write('  XUNIF_COVER(1) = 1./\n')

g.write('&NAM_ZS\n')
g.write('  XUNIF_ZS = ' + str(attributes['zorog']) + '/\n')

#g.write('&NAM_PREP_SURF_ATM')
#g.write('  NYEAR= ' + str(int(year)) + ',')
#g.write('  NMONTH=' + str(int(month)) + ',')
#g.write('  NDAY=' + str(int(day)) + ',')
#g.write('  XTIME=' + str(int(seconds)) + '/')

#if surfaceType == 'ocean':
#  g.write('&NAM_PREP_SEAFLUX')
#  if surfaceForcing == 'ts':
#    g.write('  XSST_UNIF=%(sst)6.2f,'%{'sst':sst[0]})
#  else:
#    g.write('  XSST_UNIF=300.,')
#  g.write('  NYEAR=' + str(int(year)) + ',')
#  g.write('  NMONTH=' + str(int(month)) + ',')
#  g.write('  NDAY=' + str(int(day)) + ',')
#  g.write('  XTIME=' + str(int(seconds)) + '/')

#if surfaceForcing == 'ts':
#  g.write('&NAM_DATA_SEAFLUX')
#  g.write('  LSST_DATA = T,')
#  g.write('  NTIME_SST = ' + str(nt0) + ',')
#
#  for it in range(0,nt0):
#    g.write('  XUNIF_SST(%(ii)3.i) = %(sst)6.2f,'%{'ii': it+1, 'sst': sst[it]})
#
#  for it in range(0,nt0):
#    g.write('  CFTYP_SST(%(ii)3.i) = "DIRECT",'%{"ii": it+1})
#
#  for it in range(0,nt0):
#    g.write('  NYEAR_SST(%(ii)3.i)=%(year)4.4i,  NMONTH_SST(%(ii)3.i)=%(month)2.2i,  NDAY_SST(%(ii)3.i)=%(day)2.2i , XTIME_SST(%(ii)3.i)=%(seconds)7.1f,'%{'ii': it+1, 'year': year, 'month': month , 'day': day ,'seconds': time[it]})
#
#  g.write('/')
if casname =='SANDUREF' or casname =='SANDUFAST' or casname =='SANDUSLOW':
    g.write(' ZHAT \n')
    g.write('0.0\n')
    g.write('25.0\n')
    g.write('50.0\n')
    g.write('75.0\n')
    g.write('100.0\n')
    g.write('125.0\n')
    g.write('150.0\n')
    g.write('175.0\n')
    g.write('200.0\n')
    g.write('225.0\n')
    g.write('250.0\n')
    g.write('275.0\n')
    g.write('300.0\n')
    g.write('325.0\n')
    g.write('350.0\n')
    g.write('375.0\n')
    g.write('400.0\n')
    g.write('425.0\n')
    g.write('450.0\n')
    g.write('475.0\n')
    g.write('500.0\n')
    g.write('525.0\n')
    g.write('550.0\n')
    g.write('575.0\n')
    g.write('600.0\n')
    g.write('625.0\n')
    g.write('650.0\n')
    g.write('675.0\n')
    g.write('700.0\n')
    g.write('725.0\n')
    g.write('750.0\n')
    g.write('775.0\n')
    g.write('800.0\n')
    g.write('825.0\n')
    g.write('850.0\n')
    g.write('875.0\n')
    g.write('900.0\n')
    g.write('925.0\n')
    g.write('950.0\n')
    g.write('975.0\n')
    g.write('1000.0\n')
    g.write('1025.0\n')
    g.write('1050.0\n')
    g.write('1075.0\n')
    g.write('1100.0\n')
    g.write('1125.0\n')
    g.write('1150.0\n')
    g.write('1175.0\n')
    g.write('1200.0\n')
    g.write('1225.0\n')
    g.write('1250.0\n')
    g.write('1275.0\n')
    g.write('1300.0\n')
    g.write('1325.0\n')
    g.write('1350.0\n')
    g.write('1375.0\n')
    g.write('1400.0\n')
    g.write('1425.0\n')
    g.write('1450.0\n')
    g.write('1475.0\n')
    g.write('1500.0\n')
    g.write('1525.0\n')
    g.write('1550.0\n')
    g.write('1575.0\n')
    g.write('1600.0\n')
    g.write('1625.0\n')
    g.write('1650.0\n')
    g.write('1675.0\n')
    g.write('1700.0\n')
    g.write('1725.0\n')
    g.write('1750.0\n')
    g.write('1775.0\n')
    g.write('1800.0\n')
    g.write('1825.0\n')
    g.write('1850.0\n')
    g.write('1875.0\n')
    g.write('1900.0\n')
    g.write('1925.0\n')
    g.write('1950.0\n')
    g.write('1975.0\n')
    g.write('2000.0\n')
    g.write('2025.0\n')
    g.write('2050.0\n')
    g.write('2075.0\n')
    g.write('2100.0\n')
    g.write('2125.0\n')
    g.write('2150.0\n')
    g.write('2175.0\n')
    g.write('2200.0\n')
    g.write('2225.0\n')
    g.write('2250.0\n')   
    g.write('2275.0\n')
    g.write('2300.0\n')
    g.write('2325.0\n')
    g.write('2350.0\n')
    g.write('2375.0\n')
    g.write('2400.0\n')
    g.write('2425.0\n')
    g.write('2450.0\n')
    g.write('2475.0\n')
    g.write('2500.0\n')
    g.write('2530.0\n')
    g.write('2566.0\n')
    g.write('2609.2\n')
    g.write('2661.0\n')
    g.write('2723.2\n')
    g.write('2797.9\n')
    g.write('2887.5\n')
    g.write('2995.0\n')
    g.write('3124.0\n')
    g.write('3278.8\n')
#ecriture du fichier initial
g.write('RSOU\n')
g.write(str(int(year)) + " "+ str(int(month)) + " " + str(int(day)) + " "+ str(int(timedeb))+'\n')
g.write(keyword+'\n')
g.write('%.2f' % attributes['zorog']+'\n')
g.write('%.2f' % ps[0]+'\n')
if (keyword=='ZUVTHDMR'):
    g.write('%.2f' % theta[0,0,0,0]+'\n')
elif (keyword=='ZUVTHLMR'):
    g.write('%.2f' % thetal[0,0,0,0]+'\n')
if (keyword=='ZUVTHDMR'):
    g.write('%.8f' % rv[0,0,0,0]+'\n')
elif (keyword=='ZUVTHLMR'):
    g.write('%.8f' % rt[0,0,0,0]+'\n')
nknew=np.where(height[0,:,0,0]<=zmax)[0].shape[0]
g.write('%d' % nknew+'\n')
for ik in range(0,nknew):
    g.write('%.1f %.2f %.2f' % (float(height[0,ik,0,0]),float(u[0,ik,0,0]),float(v[0,ik,0,0]))+'\n')
    #g.write('%.1f %.2f %.2f' % (height[ik][0],u[0,ik,0,0][0],v[0,ik,0,0][0])
g.write('%d' % nknew+'\n')
for ik in range(1,nknew):
    if (keyword=='ZUVTHDMR'):
        g.write('%.1f %.2f %.8f' % (height[0,ik,0,0],theta[0,ik,0,0],rv[0,ik,0,0])+'\n')
    elif (keyword=='ZUVTHLMR'):
        g.write('%.1f %.2f %.8f' % (height[0,ik,0,0],thetal[0,ik,0,0],rt[0,ik,0,0])+'\n')


#ecriture des forcages
g.write('ZFRC \n')
g.write('%d' %  nt+'\n')
# attention ne va marcher que pour 1 seul jour sinon besoin d'updater le fichier pour tenir compte de la date variable
# Attention pour l'instant pas de traitement particulier pour le trad
tadvtemp=np.zeros((nt,nk,1,1))
qadvtemp=np.zeros((nt,nk,1,1))
day0=int(day)
month0=int(month)
timeref=timedeb
nknew=np.where(height_forc[0,:,0,0]<=zmax)[0].shape[0]
print('nknew=',np.where(height_forc[0,:,0,0]<=zmax)[0].shape[0])
print('height_forc=',height_forc.shape)
for it in range(0,nt):
    if (it >=1):
        timeref=timeref+time[it]-time[it-1]
    if (timeref >= 86400):
        timeref=timeref-86400
        day0=day0+1
        if (day0 >daymonth[month0-1]):
            day0=day0-daymonth[month0-1]
            month0=month0+1
    print('timeref=',timeref,'day0=',day0,'month0=',month0)
    g.write(str(int(year)) + " "+ str(month0) + " " + str(day0) + " "+ str(timeref)+'\n')
    g.write('%.2f' % attributes['zorog']+'\n')
    g.write('%.2f' % ps_forc[it]+'\n')
    g.write('%.2f' % theta[0,0,0,0]+'\n')
    g.write('%.8f' % rv[0,0,0,0]+'\n')
    g.write('%d' % nknew+'\n')
    if ((attributes['adv_theta']==1) or (attributes['rad_theta']==1)):
        tadvtemp=thadv
    else:
        tadvtemp=np.zeros((nt,nk,1,1))
    if (attributes['adv_rv']==1):
        qadvtemp=rvadv
    else:
        qadvtemp=np.zeros((nt,nk,1,1))    
    if (attributes['forc_w']==1):
        wtemp=w
    elif (attributes['forc_omega']==1):
        wtemp=omega/(-1.*288*temp*9.81)
    else:
        wtemp=tadvtemp*0.
    if (attributes['forc_geo']==1):
        utemp=ug
        vtemp=vg
    elif (attributes['nudging_u']>0):
        utemp=unudg
        vtemp=vnudg
    else:
        utemp=np.zeros((nt,nk,1,1))
        vtemp=np.zeros((nt,nk,1,1))
    if (attributes['nudging_theta']>0):
        thtemp=thnudg
        qtemp=rvnudg
    else:
        thtemp=np.zeros((nt,nk,1,1))
        qtemp=np.zeros((nt,nk,1,1))
    for ik in range(0,nknew):
        g.write('%.1f %.2f %.2f %.2f %.6f %.7f %.7f %.10f %.2f %.2f' % (height_forc[it,ik,0,0],utemp[it,ik,0,0],vtemp[it,ik,0,0],thtemp[it,ik,0,0],qtemp[it,ik,0,0], wtemp[it,ik,0,0],tadvtemp[it,ik,0,0],qadvtemp[it,ik,0,0],0.,0.)+'\n')

g.close()
g = open('EXSEG1.nam_'+casname,'w')
g.write('&NAM_LUNITn CINIFILE ="'+casname+'_3D_REF",\n')
g.write('		 CINIFILEPGD="'+casname+'_PGD" /\n')
g.write('&NAM_CONFn LUSERV=T,NSV_USER=3/\n')
g.write('&NAM_DYNn XTSTEP=1.,XT4DIFU=300.,CPRESOPT="ZRESI" /\n')
g.write('&NAM_ADVn  CUVW_ADV_SCHEME = "CEN4TH",CTEMP_SCHEME="RKC4",\n')
g.write('           CMET_ADV_SCHEME = "PPM_01", CSV_ADV_SCHEME = "PPM_01",/\n')
g.write('&NAM_PARAMn  CTURB="TKEL", CRAD="NONE", CCLOUD="'+CLOUDSC+'", CSCONV="NONE",\n')
g.write('             CDCONV="NONE"  /\n')
g.write('&NAM_LBCn    CLBCX = 2*"CYCL", CLBCY = 2*"CYCL"/\n'		)
g.write('&NAM_TURBn XIMPL=1., CTURBLEN="DEAR", CTURBDIM="3DIM",\n')
g.write('           LTURB_FLX=T, LTURB_DIAG=T, LSUBG_COND=F,\n')
g.write('           XKEMIN=1E-10,\n')
g.write('           LSIGMAS=F, LSIG_CONV=F, LRMC01=T /\n')
if (CLOUDSC=='KHKO'):
    g.write('&NAM_PARAM_C2R2 HPARAM_CCN="CPB", HINI_CCN="CCN",\n')
    g.write('XCHEN=0.173E+09, XKHEN=1.403, XMUHEN=0.834,\n')		    
    g.write('XBETAHEN=25.499, LRAIN= F,LSEDC= F/\n')
g.write('&NAM_CONF CCONF="START", CEQNSYS ="DUR", LFLAT=T,\n')
g.write('          NMODEL=1, NVERB=6, CEXP="'+NAM1+'", CSEG="'+NAM2+'", LFORCING=T/\n')
# besoin de modifier pour avoir la bonne durée
g.write('&NAM_DYN XSEGLEN=25200., LCORIO=T,\n'	 	)
g.write('             LNUMDIFU=T,\n'   	)
g.write('             LNUMDIFTH = .F.,\n')
# ca aussi ca doit pouvoir etre paramétré => pas vraiment l'info dans les cas 1D
g.write('         XALKTOP=0.001, XALZBOT=3000 /\n')
#sorties 3D toutes les h
g.write('&NAM_BACKUP XBAK_TIME_FREQ(1) = 3600.0/\n')
# a modifier pour tenir compte des clés
g.write('&NAM_FRC \n')
if attributes['forc_geo']==1:
    g.write('         LGEOST_UV_FRC=.TRUE.,\n')
else:
    g.write('         LGEOST_UV_FRC=.FALSE.,\n')
if (attributes['adv_theta']==1) or (attributes['adv_rv']==1) or (attributes['rad_theta']==1):
    g.write('         LTEND_THRV_FRC=.TRUE.,\n')
else:
    g.write('         LTEND_THRV_FRC=.FALSE.,\n')
if (attributes['forc_w']==1) or (attributes['forc_omega']==1):
    g.write('         LVERT_MOTION_FRC=.TRUE.,\n')
else:
    g.write('         LVERT_MOTION_FRC=.FALSE.,\n')    
g.write('         LGEOST_TH_FRC=.FALSE.,\n')
if (attributes['nudging_u']>0) or (attributes['nudging_v']>0):
    g.write('         LRELAX_UV_FRC=.TRUE.,XRELAX_TIME_FRC='+str(attributes["nudging_u"])+',\n')
else:
    g.write('         LRELAX_UV_FRC=.FALSE.,\n')
if (attributes['nudging_theta']>0) or (attributes['nudging_rv']>0):
    g.write('         LRELAX_THRV_FRC=.TRUE.,XRELAX_TIME_FRC='+str(attributes["nudging_theta"])+',\n')
else:
    g.write('         LRELAX_THRV_FRC=.FALSE.,\n'    )
g.write('/\n')
g.write('&NAM_BLANK /\n')
g.write('&NAM_LES LLES_MEAN=.TRUE., LLES_SUBGRID=.TRUE.,LLES_RESOLVED=.TRUE.,\n')
g.write('         LLES_NEB_MASK = .TRUE.,\n')
g.write('            LLES_CORE_MASK = .TRUE.,\n')
g.write('            LLES_CS_MASK = .TRUE.,\n')
g.write('         XLES_TEMP_SAMPLING=60.,XLES_TEMP_MEAN_START=1.,\n' )
g.write('         XLES_TEMP_MEAN_END=25200., XLES_TEMP_MEAN_STEP=3600.  /\n')
g.write('&NAM_CONDSAMP LCONDSAMP=T, NCONDSAMP=3 /\n')
if attributes['surfaceForcing'] == 'surfaceFlux':
    g.write('&NAM_IDEAL_FLUX\n' )
    g.write('NFORCT = '+str(nt)+',\n')
    g.write('NFORCF = '+str(nt)+',\n')
    for it in range(0,nt):
        g.write('XTIMET('+str(it+1)+') = '+str(time[it])+',\n')
    for it in range(0,nt):
        g.write('XTIMEF('+str(it+1)+') = '+str(time[it])+',\n')
    for it in range(0,nt):
        g.write('XSFTH('+str(it+1)+') = '+str(float(hfss[it]))+',\n')
    g.write('CSFTQ="W/m2",\n')
    for it in range(0,nt):
        g.write('XSFTQ('+str(it+1)+') = '+str(float(hfls[it]))+',\n')
    for it in range(0,nt):
        g.write('XSFCO2('+str(it+1)+')=0.\n')
    if attributes['surfaceForcingWind'] == 'z0': 
        g.write('CUSTARTYPE = "Z0   ",\n')
        g.write('XZ0='+str(attributes["z0"])+',\n')
    if attributes['surfaceForcingWind'] == 'ustar':
        g.write('CUSTARTYPE = "USTAR",\n')
        for it in range(0,nt):
            g.write('USTAR('+str(it+1)+')= '+str(float(ustar[it]))+',\n')  
    g.write('XALB   = 0.,\n')
    g.write('XEMIS  = 1.,\n')
    for it in range(0,nt):
        g.write('XTSRAD('+str(it+1)+') = '+str(float(ts[it]))+',\n')
    g.write('/\n')
g.close()  
