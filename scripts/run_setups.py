import os,sys

if len(sys.argv)>1:
  name=sys.argv[1]
else:
  name='test'

dxs=[1000.,500.,200.,100.,50.,20.]
#dxs=[1000.,500.]

for dx in dxs:
  dt = 2*dx/10. # 10 m/s group velocity
  dtmin = dt/5
  nspool = int(1800./dt)
  stack = int(864000/dt)
  print('  run %0.0f m resolution'%dx)
  os.system('Scripts/make_setup.scr %s_%dm %0.0f %0.1f %0.1f %d %d'%(name,int(dx),dx,dt,dtmin,nspool,stack))



