from vpython import *
import numpy as np
from histogram import *
import math
import csv
k=1.38E-23
N,T=50,300
m,size =39.948/6E26,188E-12
L = 0.5*((24.4E-3/(6E23))*N)**(1/3.0)/2 + size
t,dt =0,0.5*3E-13
kk=25
average_distance=0
decrease_rate = 1  

dt0 = dt /kk
temperature=T

stage="No Changing"
sigma,epsilon =3.40E-10 ,997/6E23



# histogram setting
dv = 10.
deltav = 25. # slotwidth for v histogram
vdist = graph(x=800, y=0, ymax = N*deltav/100.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left')
theory_T = gcurve(color=color.cyan,graph=vdist) # for plot of the curve for the atom speed distribution

    
def potential_function(r):
    return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

observation = ghistogram(graph = vdist, bins=arange(0.,1500.,deltav), color=color.red) # for the simulation speed distribution

vrms = (2*k*1.5*T/m)**0.5
atoms=[]
plots=[]
potentials=[]
A = -24*epsilon*(sigma**6)
B = sigma**6

#potential graph
dx =0.01*sigma
potential_graph = graph(x=8000,y=0,ymax=3*epsilon,width=500,height=300,xtitle='distance',ytitle='potential_energy',align='left')
potential_pic =gcurve(color=color.orange,graph=potential_graph,width=4)
for x in arange(dx,10*sigma,dx):
    potential_pic.plot(pos=(x,potential_function(x)))

potentials.append(potential_pic)

scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )
p_a ,v_a,a_a = np.zeros((N,3)),np.zeros((N,3)),np.zeros((N,3))



for v in arange(0.,1500.+dv,dv): # theoretical speed distribution
    theory_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*temperature))**1.5)*exp((-0.5*m*v**2)/(k*temperature))*(v**2)*dv))
plots.append(theory_T)

def LJ_force_on_particle(a1p,a2p):   ##force on partircle_1
    vector = a1p-a2p
    distance = sum(vector*vector)**0.5

    if distance==0.0:
        return np.array([0,0,0])
    else:
        return A*( (distance**-7) - 2*B*(distance**-13) ) *vector /distance


    

def potential_on_particle(a1p,a2p):
    vector = a1p-a2p
    distance = sum(vector*vector)**0.5
    if distance==0.0:
        return np.array([0,0,0])
    else:
        return 4*epsilon*((sigma/distance)**12-(sigma/distance)**6)

    

for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N-1: # the last atom is with yellow color and leaves a trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50)
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random()))
    ra = pi*random()
    rb = 2*pi*random()
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms.append(atom)


def keyinput(evt):
    keyword = {'up':0.01,'down':-0.01}
    stop    ={' ':1}
    s=evt.key
    global stage ,decrease_rate
    if s in keyword:
        decrease_rate = decrease_rate + keyword[s]
        if decrease_rate <1:
            stage ="Now Decreasing!"    
        elif decrease_rate ==1:
            stage ="No Changing"
        elif decrease_rate >1:
            stage = "Now Increasing"
        print (' ')
        print ("you are changeing it")  
        print ("now stage=",stage)
        print ('now rate =',decrease_rate )
        print (' ')
    if s in stop:
        decrease_rate = 1
        stage ="No Changing"
        print (' ')
        print ("you are changeing it")  
        print ("now stage=",stage)
        print ('now rate =',decrease_rate )
        print (' ')
scene.bind('keydown', keyinput) # setting for the binding function



def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime

def outputresult():
    
    with open("{0}_RESULT.csv".format(foldername), "a") as output:

        output.write('time, ')
        output.write(str(t)+',')
        output.write('temperature, ')
        output.write(str(temperature)+',')
        output.write('pressure, ')
        output.write(str(pressure)+',')
        output.write('volume, ')
        output.write(str(volume)+',')
        output.write('vrms, ')
        output.write(str(vrms)+',')
        output.write('U, ')
        output.write(str(potential)+'\n')

import datetime
x = datetime.datetime.now()
foldername='{0}_{1}{2}_{3}{4}{5}'.format(x.year,x.month,x.day,x.hour,x.minute,x.second)


total_kinetic_energy=0
t0 = 500*dt # I changed this to remove 'time'
pressure = 0
volume=0
gamma=(5/3)
forcex=0
forcey=0
forcez=0
time = 0
term = 0
potential = 0
count = 0
while True:
    t += dt
    rate(10000)
    a_a = np.zeros((N,3))
    #term = 0
    length,height,width = container.length , container.height , container.width
    
    for i in range(N):
        atoms[i].pos = vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]) # to display atoms at new positions
        if abs(p_a[i][0]) >= (length/2) - size and p_a[i][0]*v_a[i][0] > 0 :  ##x
            forcex += abs(m*(1+decrease_rate)*v_a[i][0])
            v_a[i][0] =  -1*decrease_rate*v_a[i][0]
        if abs(p_a[i][1]) >= (height/2) - size and p_a[i][1]*v_a[i][1] > 0 :  ##y
            v_a[i][1] =  -1*decrease_rate*v_a[i][1]
            forcey += abs(m*(1+decrease_rate)*v_a[i][1])
        if abs(p_a[i][2]) >= (width/2) - size and p_a[i][2]*v_a[i][2] > 0 :  ##z
            v_a[i][2] =  -1*decrease_rate*v_a[i][2]
            forcez += abs(m*(1+decrease_rate)*v_a[i][2])
        if sum(v_a[i]*v_a[i]) > 30*vrms**2:
            ra = pi*random()
            rb = 2*pi*random()
            v_a[i] =0.1*(6*random()+7)*np.array([vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]  )
            p_a[i]=np.array([2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L])
            print (i, " OVER~~~~~~~~~~~~~~")
        if sum(p_a[i]*p_a[i]) > 3*(L**2):
            ra = pi*random()
            rb = 2*pi*random()
            v_a[i] =0.1*(6*random()+7)*np.array([vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]  )
            p_a[i]=np.array([2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L])
            print (i," GO AWAY~~~~~~~~~~~")
    observation.plot(data = np.sqrt(np.sum(np.square(v_a),-1)))
    
    for i in range(N-1):
        a_a[i] += np.array([0,-9.8,0])            #consider gravity
        for j in range(i+1,N):
            distance = sum((p_a[i]-p_a[j])**2)**0.5
            if distance > 10*sigma:
                continue
            if distance < 3*sigma:
                continue
            force = (LJ_force_on_particle(p_a[i],p_a[j]))
            a_a[i] += force/m
            a_a[j] += -1*force/m
            rij = p_a[i]-p_a[j]
            term += np.dot(rij,force) # this is the extra term according to eq.(2.12)
    v_a += a_a*dt
    p_a += v_a*dt # calculate new positions for all atoms
    for i in range(N-1):
        for j in range(i+1,N):
            distance = sum((p_a[i]-p_a[j])**2)**0.5
            if distance <= (2**(1/6))*sigma:
                v_a[i],v_a[j] = vcollision(p_a[i],p_a[j],v_a[i],v_a[j])
            elif (2**(1/6))*sigma<distance <= 3*sigma   :
                com_v  = 0.5*(v_a[i]+v_a[j])
                com_p  = 0.5*(p_a[i]+p_a[j])
                com_a  = 0.5*(a_a[i]+a_a[j])
                u1=v_a[i]-com_v
                u2=v_a[j]-com_v
                p1=p_a[i]-com_p
                p2=p_a[j]-com_p
                a1=a_a[i]-com_a
                a2=a_a[j]-com_a
                for l in range(kk):
                    force = (LJ_force_on_particle(p1,p2))
                    u1    += dt0 * force / m
                    u2    -= dt0 * force / m
                    p1    += dt0*(u1)
                    p2    += dt0*(u2)
                    rij = p1- p2
                    term += np.dot(rij,force)/kk # divided by kk to take average
                a1 = force/m
                a2 = -1*force/m
                p_a[i]= com_p + p1
                p_a[j]= com_p + p2
                v_a[i]= com_v + u1
                v_a[j]= com_v + u2
                a_a[i]= com_a + a1
                a_a[j]= com_a + a2
    if t>=time :
        for i in range(N):
            total_kinetic_energy += 0.5 * m *sum(v_a[i]*v_a[i])
        temperature= (total_kinetic_energy/ (1.5*N*k) )
        vrms = (3*k*temperature/m)**0.5 # current vrms
        volume = length*height*width
        pressure = ( (forcex)/(t0*width*height) + (forcey)/(t0*length*width) + (forcez)/(t0*length*height))/6
        for i in range(N-1):
            for j in range(i+1,N):
                potential += potential_on_particle(p_a[i],p_a[j])
                average_distance+=  sum((p_a[i]-p_a[j])**2)**0.5
                count+=1
        average_distance = average_distance/count
        print(' ')
        print ("T=",temperature," P=",pressure," V=",volume)
        print ("PV=",(pressure*volume)," NKT+term=",(N*k*temperature+term/3/300)," NKT=",(N*k*temperature)) # count 300 times take average
        print("P(V^5/3)=",(pressure* ((volume)**gamma) ))
        print ("Vrms=",vrms)
        print ("stage=",stage," ,rate=",decrease_rate)
        print ("term= ",term)
        print ("potential= ", potential)
        print ("average distance is ",average_distance/sigma," times of sigma")
        print ("this is L.J potential")
        print (' ')
        outputresult()
        potential = 0
        total_kinetic_energy=0
        forcex,forcey,forcez = 0,0,0
        time = time+t0
        term = 0 ## use the average of the term
        
        theory_T = gcurve(color=color.orange,graph=vdist) # for plot of the curve for the atom speed distribution
        for v in arange(0.,1500.+dv,dv): # theoretical speed distribution
            theory_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*temperature))**1.5)*exp((-0.5*m*v**2)/(k*temperature))*(v**2)*dv))
        plots.append(theory_T)
        plots[-2].delete()

        potential_pic =gcurve(color=color.orange,graph=potential_graph,width=4)
        for x in arange(dx,average_distance,dx):
            potential_pic.plot(pos=(x,    potential_function(x)  +  potential_function(-(x-average_distance))))
        potentials.append(potential_pic)
        potentials[-2].delete()
        count =0
        average_distance = 0

