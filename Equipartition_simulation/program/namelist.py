import math
import numpy as np
from random import random
##Namelist
initucell=np.zeros(2)
initucell=initucell.astype(int)

dt=            0.005    #Time for each step 
density=       0.05     #Density of the molecular
initucell[0]=  20       #Number of molecular arrange at each side
initucell[1]=  20
stepAvg=       1000     #The number of step to average the prop
stepEquil=     0        #After the step we start to evaluate the disturbution of velocity
steplimit=     500000   #Total step we want to do
temperature=   1.       #Initial temperature
NDIM=          2        #Dimensions
sizeHistVel=   50       
stepVel=       10       #After stepEquil, interval step to average velocity disturbution
limitVel=      500      

stepoutputresult= 50000 #Steps to saving data of propertity, prevent from accidently shut down of program   
stepsave=      5        #Steps to save the data of position, velocity, accerlation and LJ potential.   
method = 0         #Different kind of external potential, (0): ~1/(distance)^2 in x-axis, (1):~distance in x-axis, (2):~distance to the center, (3): ~1/(distance)^2
k0=10**0           # the const of external potential 0
k1=10**-3          # the const of external potential 1
k2=0.5*10**-3      # the const of external potential 2
k3=10**0           # the const of external potential 3
