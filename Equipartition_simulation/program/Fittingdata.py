"""
===========
Fittingdata
===========
Load the csv files saved from VelPosDist, and fitting the data
with some functions, the fitting parameters will print in the shell.
Note that csv file's name must be right.

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
from scipy import optimize

foldername='2018_820_1902'           ## Choose the folder you want to use.

stepcount=80000
xy='x'
dataP=np.loadtxt('{0}PosDist{2}_{1}-100000.csv'.format(foldername,stepcount,xy),delimiter=',',skiprows=1)
dataV=np.loadtxt('{0}VelDist_{1}-100000.csv'.format(foldername,stepcount),delimiter=',',skiprows=1)

data=dataP                           ## Choose dataP for position, dataV for velocity. 
x_data=data[0]
y_data=data[1]

def test_func(x, a, b, c):           ## Use the function of your choice with parameter being a, b, c.
    return a*np.exp(-b * ((x)**2))
    #return a*np.exp(b/abs(x))
    #return -(a*x)**2+b
    #return a*(x)*np.exp(-b*(x**2))  ##Distribution for the speed

params_initial=[0.001, 2*10**-5, 0]  ##Initial value of parameters

params, params_covariance = optimize.curve_fit(test_func, x_data, y_data, params_initial)

print(params)
print(params_covariance)

## Plot distribution
plt.figure(figsize=(6, 4))
plt.scatter(x_data, y_data,linewidths=0.1)

## Plot distribution with fitting curve
plt.figure(figsize=(6, 4))
plt.scatter(x_data, y_data, label='Data',linewidths=0.1)
plt.plot(x_data, test_func(x_data, params[0], params[1], params[2]),
         label='Fitted function')
plt.legend(loc='best')
plt.show()
print ('end')





