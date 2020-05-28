"""
========
MD03_3.2
========
Do 2D MD simulation. Create a new folder named form the time you run
the program e.g. 2018_828_143933. Save molecular's position, velocity,
accleration and the Lennard-Jones potential in binary files(npz) every
few steps in the new folder.
"""

VelPosDistimport math
import numpy as np
import csv
from random import random
from namelist import*
from setparams import*
from setupjob import*
import pylab as plt
import matplotlib.pyplot as plt

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
    global celllist,invwid,cc,c,molra,usum,virsum,rrcut,rs,www,molr,molru
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
    molru=molru*0
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
                            molru[j1]+=uval
                            molru[j2]+=uval
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
    global vBin,timenow,hFunction,final_histVel
    plt.figure('velocity_distribution')
    plt.plot(histVel,':',linewidth=1)
    plt.savefig(r'{0}_velocity_distribution'.format(foldername))
    final_histVel=histVel
    histVelList.append(histVel)   

#evaluate the velocity's distribution
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
    global t,step,vx,vy,tot,tot2,kin,kin2,p,p2
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
    with open("{0}_RESULT.csv".format(foldername), "a") as output:
        output.write('step, ')
        output.write(str(step)+'\n')
        output.write('time, ')
        output.write(str(t)+'\n')
        output.write('vx,')
        output.write(str(vx)+'\n')
        output.write('vy, ')
        output.write(str(vy)+'\n')
        output.write('tot, ')
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

#Create a folder to store files
import datetime
x = datetime.datetime.now()
foldername='{0}_{1}{2}_{3}{4}{5}'.format(x.year,x.month,x.day,x.hour,x.minute,x.second)
import pathlib
pathlib.Path('{1}'.format(foldername,foldername)).mkdir(parents=True, exist_ok=False)

#Save the position, velocity, accleration and the Lennard-Jones potential 
def savemolr_molrv():
    np.savez_compressed(r'{0}/{1}step.npz'.format(foldername,stepcount),molr=molr,molrv=molrv,molra=molra,molru=molru)
    ## If you want to save as csv files you can uncomment the following codes
    #np.savetxt(r'{0}\{1}step_r.csv'.format(foldername,stepcount),molr,delimiter=",")
    #np.savetxt(r'{0}\{1}step_ra.csv'.format(foldername,stepcount),molra,delimiter=",")
    #np.savetxt(r'{0}\{1}step_rv.csv'.format(foldername,stepcount),molrv,delimiter=",")
    #np.savetxt(r'{0}\{1}step_ru.csv'.format(foldername,stepcount),molru,delimiter=",")
    
# main program start
print (foldername)
print(steplimit)
print(method)
print(k0,k1,k2,k3)

#Save the namelist used.
file_w=open('{0}_namelist.py'.format(foldername),'a')
file_r=open('namelist.py'.format(foldername),'r')
content=file_r.read()
file_w.write('{0}\n'.format(content))
file_w.close()
file_r.close()

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
    if stepcount%stepsave==0:
        savemolr_molrv()
    if stepcount%stepAvg==0:
        accumprops(2)
        printsum()
        accumprops(0)          
    if stepcount >= stepEquil and (stepcount-stepEquil)% stepVel == 0 :
        evalVelDist()
    if stepcount%stepoutputresult==0:
        outputresult()
        h_w=open('{0}_hF.txt'.format(foldername),'a')
        h_w.write('{0}\n'.format(hFuncList))
        h_w.close()
    if stepcount>=steplimit:
        morecycle=0
        plt.figure('hFunctions')
        plt.plot(timeList,hFuncList)
        plt.savefig(r'{0}_hfunction'.format(foldername))
        print('We output velocity_distribution.png, hfunction.png and RESULT.csv ')

outputresult()
print ('end')




