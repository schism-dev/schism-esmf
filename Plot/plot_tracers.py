from pylab import *
import netCDF4
import sys

runs = {}
if len(sys.argv)>1:
  name = sys.argv[1]
else:
  name='test'


runs[1000] = name+'_1000m'
runs[500] = name+'_500m'
runs[200] = name+'_200m'
runs[100] = name+'_100m'
runs[50] = name+'_50m'
runs[20] = name+'_20m'

lss = ['-','--',':','-.','-..']
tidx = 477

fig = figure(figsize=(8,5))

keys = sorted(runs)

for dx in keys:
  folder = runs[dx]
  label = '%d m'%dx

  print(label)
  nc = netCDF4.Dataset('%s/outputs/1_salt.nc'%folder)
  ncv=nc.variables

  t = ncv['time'][:]
  y = ncv['y'][:]
  yidx = where(y==0.0)
  s = ncv['salt'][tidx,-1][yidx]-35.
  x = ncv['x'][yidx]/1000.

  
  if label=='20 m':
    col = 'r'
    ls = '-'
    s0 = ncv['salt'][0,-1][yidx]-35.
    plot(x,s0,'-',lw=3.0,color='orange',label='initial condition')
  else:
    col = (log10(int(label[:-2]))/3.-0.2)*ones((3))
    ls = '-'#lss.pop()
    

  plot(x,s,ls=ls,color=col,lw=2.0,label=label)

  nc.close()

ylim(-0.1,1.1)
ylabel('Tracer [1/1]')
xlabel('km')
legend(loc='upper left',frameon=False)

savefig('transport_comparison_%s.pdf'%name)
show()
