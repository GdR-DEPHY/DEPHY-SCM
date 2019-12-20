import os

import netCDF4 as nc

os.system('rm Composite_REF_driver_RR_new3_converted.nc')
os.system('ncrename -v temp_nudg,temp_nudging -v qv_nudg,qv_nudging -v u_nudg,u_nudging -v v_nudg,v_nudging -v ps,ps_forc Composite_REF_driver_RR_new3.nc Composite_REF_driver_RR_new3_converted.nc')

f = nc.Dataset('Composite_REF_driver_RR_new3_converted.nc','a')
z = f['height'][0,:,0,0]
f['lev'][:] = z[:]
f['lev'].setncattr('units','m')

f.close()

