import math
import numpy as np
from random import random
from namelist import*
from setparams import*
#allocarray
#
#initcoords
gap=region/initucell
n=0
molr=np.zeros((nMol,2))
molrv=np.zeros((nMol,2))
molra=np.zeros((nMol,2))
molru=np.zeros(nMol)
c=np.zeros(2)
for ny in range(initucell[1]):
    for nx in range(initucell[0]):
        c[0]=nx+0.5
        c[1]=ny+0.5
        c=c*gap
        c+=-0.5*region
        molr[n]=c  #n up to nMol
        n+=1     
#initVels
np.random.seed(17)        
vsum=np.zeros(2)
randomang=2*math.pi*np.random.rand(nMol)
molrv[:,0]=np.cos(randomang)
molrv[:,1]=np.sin(randomang)
molrv=molrv*velMag
vsum=molrv.sum(axis=0)
molrv-=1./nMol*vsum
#initaccels
  #molra is 0 already
#accumprops(0)
totEnergy.sum=0
kinEnergy.sum=0
pressure.sum=0
totEnergy.sum2=0
kinEnergy.sum2=0
pressure.sum2=0

def VLinear(p,s):
    return (p[:,1])*s[0]+p[:,0]
N_OFFSET=5
vOff=[[0,0],[1,0],[1,1],[0,1],[-1,1]]
m1v=np.zeros((prod_cells,NDIM))
m2v=np.zeros((N_OFFSET,prod_cells,NDIM))
shift=np.zeros((N_OFFSET,prod_cells,NDIM))
m2=np.zeros((N_OFFSET,prod_cells))
i=0
for m1y in range(cells[1]):
    for m1x in range(cells[0]):
        m1v[i]=np.array([m1x,m1y])
        i+=1
m1=(VLinear(m1v,cells)+nMol).astype(int)
for offset in range(N_OFFSET):
    m2v[offset]=m1v+vOff[offset]
    shift[offset]=(m2v[offset]>=cells)*region-(m2v[offset]<0)*region
    m2v[offset]=(m2v[offset]>=0)*(m2v[offset]<cells)*m2v[offset]+(m2v[offset]<0)*(cells-1)
    m2[offset]=(VLinear(m2v[offset],cells)+nMol).astype(int)
    m2=(m2).astype(int)  

# for EvalVelDist
deltaP=0
histSumPos=0
rangePos=region[0]
countPos=0


deltaV=0
histSum=0
rangeVel=6.
hFunction=0
countVel=0


morecycle=1
stepcount=0
timenow=0
