"""
==========
VelPosDist
==========
Load the files(npz) saved from MD03_3.2, and compute the position
distribution in X or Y direction or the velocity distribution.
The figure will be shown and save a csv file. 

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
    
#present the result of velcity's distribution
def printVelDist():
    global vBin,timenow,hFunction,final_histVel
    plt.figure('velocity_distribution')
    plt.plot(np.arange(0,rangeVel,deltaV),histVel,'-',linewidth=1)
    plt.savefig(r'{0}\velocity_distribution2'.format(foldername))
    final_histVel=histVel
    histVelList.append(histVel)
    file=np.savetxt('{0}VelDist_{1}-{2}.csv'.format(foldername,stepcount,stepVel*limitVel),(np.arange(0,rangeVel,deltaV),histVel),delimiter=',',
                    header='stepcount={0},stepPos={1},limitPos={2},Row2 vel,Row3 propotion'.format(stepcount,stepPos,limitPos))
    
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

#present the result of velcity's distribution
def printPosDist(xy):
    global vBin,timenow,hFunction,final_histVel
    plt.figure('Pos_distribution')
    if xy=='x':
        plt.plot(np.arange(0,rangePos,deltaP)-region[0]*0.5,histPos,'-',linewidth=1)
        plt.xlabel('pos_x')
    if xy=='y':
        plt.plot(np.arange(0,rangePos,deltaP)-region[1]*0.5,histPos,'-',linewidth=1)
        plt.xlabel('pos_y')
    if xy=='r':
        plt.plot(np.arange(0,rangePos,deltaP)[0:50],histPos[0:50],'-',linewidth=1)
        plt.xlabel('pos_r')
    plt.plot(0,0)
    plt.plot(0,0.1)
    file=np.savetxt('{0}PosDist{2}_{1}-{3}.csv'.format(foldername,stepcount,xy,stepPos*limitPos),(np.arange(0,rangePos,deltaP)-region[0]*0.5,histPos),delimiter=',',
                    header='stepcount={0},stepPos={1},limitPos={2},Row2 pos,Row3 propotion'.format(stepcount,stepPos,limitPos))
    histVelList.append(histVel)   

#evaluate the velcity's distribution
def evalPosDist(xy):
    global histPos, rangePos, countPos, limitPos, sizeHistPos, stepPos
    global deltaP, histSumPos,j1     
    if countPos == 0:
        histPos=histPos*0   
    deltaP = rangePos / sizeHistPos
    if xy=='x':
        j1 = (molr[:,0]+region[0]*0.5) / deltaP
    if xy=='y':
        j1 = (molr[:,1]+region[1]*0.5) / deltaP
    if xy=='r':
        j1 = np.sqrt((molr*molr).sum(axis=1)) / deltaP
    j1=j1.astype(int)
    j1=j1*(j1<sizeHistPos)+(sizeHistPos-1)*(j1>=sizeHistPos)
    for i in range(sizeHistPos):
        histPos[i]+=(j1==i).sum()
    countPos+=1
    if countPos == limitPos :
        histSumPos=histPos.sum()
        print(histSumPos)
        histPos/=histSumPos
        j2=np.arange(sizeHistPos)
        printPosDist(xy)
        countPos=0


foldername=  '2018_825_63949' ## Choose the folder you want to use.
stepcount=    490000          ## Choose steps you want to start.
steplimit=    500000          ## Choose steps you want to end.

rangeVel=     6               ## Choose the range of velocity. 
sizeHistVel=  50              ## Division of the range of velocity.
stepVel=      200             ## Calculate velocity distribution once how many steps.
limitVel=     500             ## Average velocity distribution once how many times.
histVel=np.zeros(sizeHistVel)

sizeHistPos=  50              ## Division of region. 
stepPos=      200             ## Calculate position distribution once how many steps.
limitPos=     500             ## Average position distribution once how many times.
histPos=np.zeros(sizeHistPos)

while morecycle==1:
    stepcount+=1
    timenow=stepcount*dt
    if stepcount >= stepEquil and (stepcount-stepEquil)% stepPos == 0 :
        data=np.load(r'{0}/{1}step.npz'.format(foldername,stepcount))
        molrv=data['molrv']
        molr=data['molr']
        evalVelDist()         ## Compute velocity distribution.
        evalPosDist('x')      ## Compute position distribution in the direction of 'x','y' or 'r'.
    if stepcount>=steplimit:
        morecycle=0
        plt.show()
        
print ('end')






