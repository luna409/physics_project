import math
import numpy as np
from namelist import*
region=np.zeros(2)
region=(1./math.sqrt(density))*initucell  
nMol=initucell[0]*initucell[1]  
velMag=math.sqrt(NDIM*(1.-1./nMol)*temperature)  
rcut=2**(1./6.)  
rrcut=rcut**2
r=rcut/2
class Prop():
 def __init__(self, val,sum,sum2):
     self.val = val
     self.sum=sum
     self.sum2=sum2
kinEnergy=Prop(0,0,0)
totEnergy=Prop(0,0,0)
pressure=Prop(0,0,0)

cells=(region/rcut).astype(int)
prod_cells=cells[0]*cells[1]
xx=region[0]/2
yy=region[1]/2
xylim=max(xx+r,yy+r)
timeList=[]
hFuncList=[]
histVelList=[]
t=['']
step=['']
vx=['']
vy=['']
tot=['']
tot2=['']
kin=['']
kin2=['']
p=['']
p2=['']
