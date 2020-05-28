This is a 2D molecular dynamics (MD) simulation program. We use Lennard-Jones potential to approximate the interaction between moleculars and with 4 kinds of external potential can be added. We use leapfrog method to do the intergration and use cell subdivision to help us to have a more efficient way to calculate the force between moleculars. 

Requirements of program language :
Python3
Requirements of Python3 modules :
Numpy, Matplotlib, Scipy

There are some python files.
namelist¡GValues of simulation's initial condition .
setparams¡GSet parameters based on initial condition.
setupjob¡GSet the initial position, velocity and acceleration.
MD03_3.1 : Do 2D MD simulation and display the instant position of molecular. But the history data will not be save.
MD03_3.2 : Do 2D MD simulation and save molecular's position, velocity, accleration and the Lennard-Jones potential in the form of npz files in a new folder every few step.
MD03_3.3 : Load the npz files saved from MD03_3.2, and display the position of molecular.
VelPosDist¡GLoad the npz files saved from MD03_3.2, and compute the position distribution in X or Y direction or the velocity distribution. Plot the figure and save a csv file. 
Fittingdata¡GLoad the csv files saved from VelPosDist, and fitting the data with functions.



You can change the value in namelist to the condiction you want : 

dt=            0.005        #Time for each step 
density=       0.05        #Density of the molecular
initucell[0]=  20           #Number of molecular arrange at each side
initucell[1]=  20
stepAvg=       1000      #The number of step to average the prop
stepEquil=     0             #After the step we start to evaluate the disturbution of velocity
steplimit=     500000    #Total step we want to do
temperature=   1.         #Initial temperature
NDIM=          2             #Dimensions
sizeHistVel=   50       
stepVel=       10       #After stepEquil, interval step to average velocity disturbution
limitVel=      500      

stepoutputresult= 50000 #Steps to saving data of propertity, prevent from accidently shut down of program   
stepsave=      5       #Steps to save the data of position, velocity, accerlation and LJ potential.   
method = 0             #Different kind of external potential, (0): ~1/(distance)^2 in x-axis, (1):~distance in x-axis, (2):~distance to the center, (3): ~1/(distance)^2
k0=10**0               # the const of external potential 0
k1=10**-3             # the const of external potential 1
k2=0.5*10**-3      # the const of external potential 2
k3=10**0              # the const of external potential 3


