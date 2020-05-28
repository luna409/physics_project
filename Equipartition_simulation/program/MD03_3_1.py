"""
========
MD03_3.1
========
Do 2D MD simulation and plot the position of moleculars.
"""

import math
import numpy as np
import csv
from random import random
from namelist import*
from setparams import*
from setupjob import*
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle
from matplotlib.collections import PatchCollection

# define functions

#use leapforg integration to compute velocity and position
def leapforg(part):
    global molrv,molr,dt
    if part==1:
        molrv+=0.5*dt*molra
        molr+=dt*molrv
    if part==2:
        molrv+=0.5*dt*molra
        
#use the boundary condition
def applyboundary():
    global molr
    molr=molr-((molr+0.5*region)/region).astype(int)*region*(molr>=0.5*region)-((molr-0.5*region)/region).astype(int)*region*(molr<-0.5*region)

#part 0: set all prop to 0
#part 1: add the prop's value to sum
#part 2: averge the prop of a period
def accumprops(part):
    global totEnergy,kinEnergy,pressure
    if part==0:
        totEnergy.sum,kinEnergy.sum,pressure.sum=0,0,0
        totEnergy.sum2,kinEnergy.sum2,pressure.sum2=0,0,0
    if part==1:
        totEnergy.sum+=totEnergy.val
        totEnergy.sum2+=(totEnergy.val)**2
        kinEnergy.sum+=kinEnergy.val
        kinEnergy.sum2+=(kinEnergy.val)**2
        pressure.sum+=pressure.val
        pressure.sum2+=(pressure.val)**2
    if part==2:
        totEnergy.sum/=stepAvg
        totEnergy.sum2=math.sqrt(max(totEnergy.sum2/stepAvg-(totEnergy.sum)**2,0))
        kinEnergy.sum/=stepAvg
        kinEnergy.sum2=math.sqrt(max(kinEnergy.sum2/stepAvg-(kinEnergy.sum)**2,0))
        pressure.sum/=stepAvg
        pressure.sum2=math.sqrt(max(pressure.sum2/stepAvg-(pressure.sum)**2,0))

#use Lennard-Jones potential to compute the force        
def computeforce():
    global celllist,invwid,cc,c,molra,usum,virsum,rrcut,rs,www,molr
    rrcut=(rcut)**2#
    invwid=cells/region#
    celllist=np.arange(prod_cells+nMol)#
    celllist=-1*(celllist>=nMol)#
    rs=molr+0.5*region
    cc=(rs*invwid).astype(int)
    c=VLinear(cc,cells)+nMol
    for n in range(nMol):
        celllist[n]=celllist[c[n]]
        celllist[c[n]]=n
    molra=molra*0
    usum=0
    virsum=0
    for i in range(prod_cells):
        for offset in range(N_OFFSET):
            j1=celllist[m1[i]]
            while j1>=0 :
                j2=celllist[m2[offset,i]]
                while j2>=0:
                    if m1[i]!=m2[offset,i] or j2<j1:
                        dr=molr[j1]-molr[j2]
                        dr=dr-shift[offset,i]
                        rr=(dr*dr).sum()
                        if rr<rrcut and rr>0 :
                            rri=1./rr
                            rri3=rri**3
                            fcval=48*rri3*(rri3-0.5)*rri
                            uval=4.*rri3*(rri3-1)+1
                            molra[j1]=molra[j1]+fcval*dr
                            molra[j2]=molra[j2]-fcval*dr
                            usum+=uval
                            virsum+=fcval*rr
                    j2=celllist[j2]
                j1=celllist[j1]

#use external potential to compute the force
def computeforcePotential(method):
    if method ==0:
        molra[:,0] += -k0/(molr[:,0]*molr[:,0])*np.sign(molr[:,0])*(abs(molr[:,0])>rcut*0.5)    
    if method ==1:
        molra[:,0] += -k1 * (molr[:,0])    
    if method ==2:
        force=k2 * (molr[:,0]**2 + molr[:,1]**2)**1/2
        molra[:,0] += -(molr[:,0]/((molr[:,0]**2 + molr[:,1]**2)**1/2))*force
        molra[:,1] += -(molr[:,1]/((molr[:,0]**2 + molr[:,1]**2)**1/2))*force   
    if method ==3:
        force=(np.sqrt(molr[:,0]**2 + molr[:,1]**2)>rcut*0.5)*k3 / (molr[:,0]**2 + molr[:,1]**2)
        molra[:,0] += -(molr[:,0]/((molr[:,0]**2 + molr[:,1]**2)**1/2))*force
        molra[:,1] += -(molr[:,1]/((molr[:,0]**2 + molr[:,1]**2)**1/2))*force

#evaluate the prop's value                
def evalprop():
    global totEnergy,kinEnergy,pressure,vsum
    vv=np.zeros(NDIM)
    vvsum=0.
    vsum=molrv.sum(axis=0)
    vv=molrv*molrv
    vvsum=vv.sum()
    kinEnergy.val=0.5*vvsum/nMol
    totEnergy.val=kinEnergy.val+usum/nMol
    #print(usum/nMol,'       ',kinEnergy.val)
    pressure.val=density*(vvsum+virsum)/(nMol*NDIM)

#print the result
def printsum():
    print ('stepcount= ',stepcount)
    print ('timenow= ',timenow)
    print (vsum/nMol)
    print ('tot= ',totEnergy.sum)
    print ('tot2= ',totEnergy.sum2)
    print ('kin= ',kinEnergy.sum)
    print ('kin2= ',kinEnergy.sum2)
    print ('p= ',pressure.sum)
    print ('p2= ',pressure.sum2)
    print
    t.append(timenow)
    step.append(stepcount)
    vx.append(vsum[0]/nMol)
    vy.append(vsum[1]/nMol)
    tot.append(totEnergy.sum)
    tot2.append(totEnergy.sum2)
    kin.append(kinEnergy.sum)
    kin2.append(kinEnergy.sum2)
    p.append(pressure.sum)
    p2.append(pressure.sum2)
    
#present the result of velcity's distribution
def printVelDist():
    global vBin,timenow,hFunction
    plt.figure('velocity_distribution')
    plt.plot(histVel)
    plt.savefig('velocity_distribution')
    histVelList.append(histVel)

#evaluate the velcity's distribution
def evalVelDist():
    global histVel, rangeVel, countVel, limitVel, sizeHistVel, stepVel
    global deltaV, histSum,j1,hFunction      
    if countVel == 0:
        histVel=histVel*0   
    deltaV = rangeVel / sizeHistVel
    j1 = np.sqrt((molrv*molrv).sum(axis=1)) / deltaV
    j1=j1.astype(int)
    j1=j1*(j1<sizeHistVel)+(sizeHistVel-1)*(j1>=sizeHistVel)
    for i in range(sizeHistVel):
        histVel[i]+=(j1==i).sum()
    countVel+=1
    if countVel == limitVel :
        histSum=histVel.sum()
        histVel/=histSum
        j2=np.arange(sizeHistVel)
        hFunction = (histVel * np.log((histVel/ (j2+0.5) * deltaV)+(histVel==0))).sum()
        timeList.append(timenow)
        hFuncList.append(hFunction)
        printVelDist()
        countVel=0
        
#output the result
def outputresult():
    t.append('')
    step.append('')
    vx.append('')
    vy.append('')
    tot.append('')
    tot2.append('')
    kin.append('')
    kin2.append('')
    p.append('')
    p2.append('')
    with open("RESULT123.csv", "w") as output:
        output.write('step, ')
        output.write(str(step)+'\n')
        output.write('time, ')
        output.write(str(t)+'\n')
        output.write('vx,')
        output.write(str(vx)+'\n')
        output.write('vy, ')
        output.write(str(vy)+'\n')
        output.write('vz, ')
        output.write(str(tot)+'\n')
        output.write('tot2, ')
        output.write(str(tot2)+'\n')
        output.write('kin, ')
        output.write(str(kin)+'\n')
        output.write('kin2, ')
        output.write(str(kin2)+'\n')
        output.write('p, ')
        output.write(str(p)+'\n')
        output.write('p2, ')
        output.write(str(p2)+'\n')

def plotnow():
    fig=plt.figure('pos_now',figsize=(6,6),clear=True)
    fig.suptitle("--- %s step ---" % (stepcount))
    plt.ion()
    plt.plot([-xx,-xx,xx,xx,-xx],[-yy,yy,yy,-yy,-yy],'--',linewidth=0.5)
    ax=fig.add_subplot(111)   
    patches = []
    for x1, y1 in zip(molr[:,0], molr[:,1]):
        circle = Circle((x1, y1), r)
        patches.append(circle)
    colors = 20*(np.zeros(len(patches)))
    pa = PatchCollection(patches, alpha=0.4)
    pa.set_array(np.array(colors))
    ax.add_collection(pa)
    ax.set_xlim([-xylim,xylim])#
    ax.set_ylim([-xylim,xylim])#
    plt.pause(0.001)


while morecycle==1:
    stepcount+=1
    timenow=stepcount*dt
    leapforg(1)
    applyboundary()
    computeforce()
    computeforcePotential(method)
    leapforg(2)    
    evalprop()            
    accumprops(1)
    if stepcount%5==0:
        plotnow()
    if stepcount%stepAvg==0:
        accumprops(2)
        printsum()
        accumprops(0)          
    if stepcount >= stepEquil and (stepcount-stepEquil)% stepVel == 0 :
        evalVelDist()
    if stepcount>=steplimit:
        morecycle=0
        print('We output velocity_distribution.png, hfunction.png and RESULT.csv ')

outputresult()
print ('end')






