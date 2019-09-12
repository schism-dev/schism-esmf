from pylab import *
import netCDF4

nc = netCDF4.Dataset('outputs/1_hvel.nc')
ncv=nc.variables

ul = ncv['u'][:,-1,22]
ur = ncv['u'][:,-1,42]

time = ncv['time'][:]/86400.

figure(figsize=(10,3))

plot(time,ul,'k-',lw=2.0,label='left')
plot(time,ur,'-',lw=2.0,color='orange',label='right')
ylabel('[m/s]')
xlabel('days')
legend()

show()
nc.close()
