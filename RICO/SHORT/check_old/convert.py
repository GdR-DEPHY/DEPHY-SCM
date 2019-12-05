import os

import netCDF4 as nc

os.system('rm rico_driver_RR_new3_converted.nc')
os.system('ncrename -v tadvh,tadv -v qvadvh,qvadv -v ps,ps_forc rico_driver_RR_new3.nc rico_driver_RR_new3_converted.nc')

f = nc.Dataset('rico_driver_RR_new3_converted.nc','a')
z = f['height'][0,:,0,0]
f['lev'][:] = z[:]
f['lev'].setncattr('units','m')

f.close()

