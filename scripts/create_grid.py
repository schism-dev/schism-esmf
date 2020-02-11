from pylab import *
from matplotlib.tri import Triangulation
import os,sys

def write_gr3(tri,fname,const=0.1):
  fh = open(fname,'w')
  fh.write('%s description\n'%fname)
  fh.write('%d %d\n'%(len(tri.triangles),len(tri.x)))
  for n,(x,y) in enumerate(zip(tri.x,tri.y)):
    fh.write('%d %0.2f %0.2f %f\n'%(n+1,x,y,const))
  fh.close()

def write_bctides(tri,amp=0.1):
  fh = open('bctides.in','w')
  fh.write("""01/01/2016 00:00:00 PST
0 40. ntip
1  nbfr
M2
 1.405189e-04  0.98160 98.26966
2 nope
3 3 0 0 2 !left
M2
%0.2f 32.25
%0.2f 32.25
%0.2f 32.25
35.
1.
3 3 0 0 2 !right
M2
%0.2f 0.0
%0.2f 0.0
%0.2f 0.0
36.
1.
  """%(amp,amp,amp,amp,amp,amp))
  fh.close()

f = open('hgrid.gr3','w')

length=20000. # 20km
if len(sys.argv)>1:
  dx=float(sys.argv[1])
else:
  dx = 500. # m
depth=10. # m


nodes = []
nv = []
pointsx=[0.0,dx]
pointsx.extend(list(arange(dx+2*dx,length,2*dx)))
pointsx.append(length)
pointsy = list(0.0*ones(len(pointsx)))

upperx = list(arange(0.0,length,2*dx))
upperx.append(length)
pointsx.extend(upperx)
pointsx.extend(upperx)

pointsy.extend(list(dx*ones(len(upperx))))
pointsy.extend(list(-dx*ones(len(upperx))))

tri = Triangulation(pointsx,pointsy)
triplot(tri)


f.write('hgrid.gr3\n')
f.write('%d %d\n'%(len(tri.triangles),len(tri.x)))

for n,(x,y) in enumerate(zip(tri.x,tri.y)):
  f.write('%d %0.2f %0.2f %0.2f\n'%(n+1,x,y,depth))

for n,(i,j,k) in enumerate(tri.triangles):
  f.write('%d 3 %d %d %d\n'%(n+1,i+1,j+1,k+1))

# number of open boundaries
f.write('2\n')
# number of open boundary nodes
f.write('6\n')
# number of left boundary nodes and nodes
f.write('3\n')
for i in where(tri.x==0.0)[0]:
  f.write('%d\n'%(i+1))

# number of left boundary nodes and nodes
f.write('3\n')
for i in where(tri.x==length)[0]:
  f.write('%d\n'%(i+1))


nlower = len(where(tri.y==-dx)[0])
nupper = len(where(tri.y==dx)[0])
# number of land boundaries
f.write('2\n')
# number of land boundary nodes
f.write('%d\n'%(nupper+nlower-4))
# number of lower land boundary nodes and nodes
f.write('%d\n'%(nlower-2))
print(nupper)
for i in where(tri.y==-dx)[0][1:-1]:
  f.write('%d\n'%(i+1))

# number of upper land boundary nodes and nodes
f.write('%d\n'%(nupper-2))
for i in where(tri.y==dx)[0][1:-1]:
  f.write('%d\n'%(i+1))

# number of island boundaries
f.write('0\n')

f.close()

# write salt initial condition
fs = open('salt.ic','w')
fs.write('salinity gradient !description\n')
fs.write('%d %d\n'%(len(tri.triangles),len(tri.x)))

saltv=[]
for n,(x,y) in enumerate(zip(tri.x,tri.y)):
  salt = interp(x,[0.0,0.45*length,0.55*length,length],[35.0,35.0,36.0,36.0])
  saltv.append(salt)
  fs.write('%d %0.2f %0.2f %0.2f\n'%(n+1,x,y,salt))
fs.close()


write_gr3(tri,'diffmax.gr3',const=0.1)
write_gr3(tri,'diffmin.gr3',const=0.000001)
write_gr3(tri,'rough.gr3',const=0.001)
write_gr3(tri,'xlsc.gr3',const=0.5)
write_gr3(tri,'temp.ic',const=10.)
os.system('ln -sf hgrid.gr3 hgrid.ll')
write_bctides(tri,amp=0.1)

# write tvd.prop
fh = open('tvd.prop','w')
for i in range(len(tri.triangles)):
  fh.write('%d 1\n'%(i+1))
fh.close()

# plot grid
fig = figure(figsize=(10,2))
tripcolor(tri,saltv,shading='flat',cmap=cm.RdYlGn_r)
triplot(tri,'k-',lw=1.0)
xlim(-1000,length+1000)
ylim(-500,dx+500)
ax=gca()
ax.set_yticks([0,500])
ax.set_aspect('equal')
xlabel('x [m]')
ylabel('y [m]')
savefig('hgrid.pdf')
