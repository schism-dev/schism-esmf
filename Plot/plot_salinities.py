from pylab import *
import netCDF4

runs = {}

runs['upwind'] = 'o001'
runs['tvd'] = 'o002'
runs['tvd 0.1*dt'] = 'o003'
runs['tvd 10*dt'] = 'outputs'

lss = ['-','--',':','-.']

fig = figure(figsize=(8,5))

for label,folder in runs.iteritems():

  nc = netCDF4.Dataset('%s/1_salt.nc'%folder)
  ncv=nc.variables

  t = ncv['time'][:]
  y = ncv['y'][:]
  yidx = where(y==0.0)
  s = ncv['salt'][-1,-1][yidx]-35.
  x = ncv['x'][yidx]/1000.

  
  if label=='upwind':
    col = 'r'
    ls = '-'
    s0 = ncv['salt'][0,-1][yidx]-35.
    plot(x,s0,'-',lw=3.0,color=(0.6,0.6,0.6),label='initial condition')
  else:
    col = 'k'
    ls = lss.pop()
    

  plot(x,s,ls=ls,color=col,lw=2.0,label=label)

  nc.close()

ylim(-0.1,1.1)
ylabel('Tracer [1/1]')
xlabel('km')
legend(loc='upper left',frameon=False)

savefig('transport_comparison.pdf')
show()
