import os,sys

if len(sys.argv)>1:
  name=sys.argv[1]
else:
  name='dt'

dx = 200.
cfls=[0.25,0.5,0.75,1.0,2.0,4.0,8.0,16.]
#cfls=[4.0,8.0]
#dxs=[1000.,500.]

for cfl in cfls:
  dt = cfl*dx/10. # 10 m/s group velocity
  dtmin = dt/5
  nspool = int(1800./dt)
  stack = int(864000/(nspool*dt))*nspool
  print('  run %0.0f s timestep'%dt)
  os.system('Scripts/make_setup.scr %s_%ds %0.0f %0.1f %0.2f %d %d'%(name,int(dt),dx,dt,dtmin,nspool,stack))



