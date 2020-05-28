"""
========
MD03_3.3
========
Load the files(npz) saved from MD03_3_2, and plot the position of
moleculars.
You can change the following four as you want:
stepcount=stepcount           ##
steplimit=steplimit           ##
stepplot=10                   ##
foldername='2018_828_135857'  ## 
"""

import math
import numpy as np
import csv
from namelist import*
from setparams import*
from setupjob import*
import pylab as plt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle
from matplotlib.collections import PatchCollection

# define functions
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
    ax.set_xlim([-xylim,xylim])
    ax.set_ylim([-xylim,xylim])
    plt.pause(0.01)

foldername='2018_828_135857'  ## Choose the folder you want to use.
stepcount=stepcount           ## Choose steps you want to start.
steplimit=steplimit           ## Choose steps you want to end.
stepplot=10                   ## Plot the figure once how many steps.
  
while morecycle==1:
    stepcount+=1
    timenow=stepcount*dt
    if stepcount%stepplot==0:
        data=np.load('{0}\{1}step.npz'.format(foldername,stepcount))
        molr=data['molr']
        #molr=np.loadtxt('{0}\{1}step_r.csv'.format(foldername,stepcount),delimiter=",")
        plotnow()
    if stepcount>=steplimit:
        morecycle=0
print ('end')






