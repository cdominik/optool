import optool
import os.path

if (os.path.exists("grid")):
    grid   = optool.particle(dir='grid')
else:
    grid   = optool.particle('./optool -nl 30 -s -d',       'grid')
if (os.path.exists("one")):
    one    = optool.particle(dir='one')
else:
    one    = optool.particle('./optool -nl 30 -s   ',       'one')
if (os.path.exists("grid_r")):
    grid_r = optool.particle(dir='grid-r')
else:
    grid_r = optool.particle('./optool -nl 30 -s -d -radmc','grid-r')
if (os.path.exists("grid_r")):
    one_r  = optool.particle(dir='one-r')
else:
    one_r  = optool.particle('./optool -nl 30 -s    -radmc','one-r')


nn = grid.a1**(-2.5)/1000
p = grid.sizedist(nn)
one.plot()
p.plot()

nn = grid_r.a1**(-2.5)/1000
q = grid_r.sizedist(nn)
one_r.plot()
q.plot()
