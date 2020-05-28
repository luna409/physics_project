from vpython import *
import numpy as np
from histogram import *

k=1.38E-23
N,T=50,298
m,size =39.948/6E26,188E-12
L = 0.5*((24.4E-3/(6E23))*N)**(1/3.0)/2 + size
t,dt =0,0.5*3E-13
kk=360               #球碰到之後多算360次

decrease_rate = 1   #碰到牆壁速度變幾分之幾

dt0 = dt /kk        #分割成20次變更細
temperature=T       #溫度先令成T 等下溫度會改變

stage="No Changing" #等下decrease rate 會改變，按上下鍵
sigma,epsilon =3.40E-10 ,997/6E23 #refer to wiki #交點等於零sigma eps位能最低點


# histogram setting
dv = 10.
deltav = 25. # slotwidth for v histogram
vdist =     graph(x=800, y=0, ymax = N*deltav/100.,width=500, height=300, xtitle='v', ytitle='dN', align = 'left')
theory_T = gcurve(color=color.cyan) # for plot of the curve for the atom speed distribution


observation = ghistogram(graph = vdist, bins=arange(0.,1500.,deltav), color=color.red) # for the simulation speed distribution

vrms = (2*k*1.5*T/m)**0.5 #初始Ar氣的方均根速度
atoms=[]            #畫好一顆球就加進去
plots=[]            #畫理論曲線麥克斯威爾分佈

A = -24*epsilon*sigma**6    #紙筆計算出的
B = sigma**6                #同上

scene = canvas(width=500, height=500, background=vector(0.2,0.2,0), align = 'left')
container = box(length = 2*L, height = 2*L, width = 2*L, opacity=0.2, color = color.yellow )
p_a ,v_a,a_a = np.zeros((N,3)),np.zeros((N,3)),np.zeros((N,3))


#畫理論曲線圖
for v in arange(0.,1500.+dv,dv): # theoretical speed distribution
    theory_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*temperature))**1.5)*exp((-0.5*m*v**2)/(k*temperature))*(v**2)*dv))
plots.append(theory_T)

c1=4*epsilon*12*sigma**12
c2=4*epsilon*6*sigma**6
##V_JL = 4*epsilon*((sigma/distance)**12-(sigma/distance)**6)
##V_JL = 4*epsilon*(sigma**12/distance**12-sigma**6/distance**6)
##F_JL = 4*epsilon*((12*sigma**12)/distance**13-(6*sigma**6)/distance**7)
##F_JL = (c1/distance**14 + c2/distance**8 ) * vector
#粒子受到的力
def force_on_particle(a1p,a2p):   ##for on partircle_2
    vector = a1p-a2p              #一號球跟二號球的相對位置
    distance = sum(vector*vector)**0.5  #距離大小
    if distance==0.0:               # 距離為零會有錯另外定義 讓他重置成0
        return np.array([0,0,0])
    else:
        return (c1/distance**14 + c2/distance**8 ) * vector #手動算LJ 對r微分
    
#粒子potential不會用到
def potential_on_particle(a1p,a2p):
    vector = a1p-a2p            #2u
    distance = sum(vector*vector)**0.5
    if distance==0.0:
        return np.array([0,0,0])
    else:
        return 4*epsilon*((sigma/distance)**12-(sigma/distance)**6) #LJ式子
#粒子初始條件 隨機位置 V
for i in range(N):
    p_a[i] = [2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L] # particle is initially random positioned in container
    if i== N-1: # the last atom is with yellow color and leaves a trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=color.yellow, make_trail = True, retain = 50) #retain50個迴圈 軌跡才不會越來越長
    else: # other atoms are with random color and leaves no trail
        atom = sphere(pos=vector(p_a[i, 0], p_a[i, 1], p_a[i, 2]), radius = size, color=vector(random(), random(), random()))
    #角度令成均勻的的機率 隨機出現
    ra = pi*random() #anale theta
    rb = 2*pi*random() # angle phi
    v_a[i] = [vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)] # particle initially same speed but random direction
    atoms.append(atom)  #每做好一個球ㄧ就加到list裡面

#鍵盤輸入
def keyinput(evt):
    keyword = {'up':0.01,'down':-0.01}  #鍵盤按上加0.01
    stop    ={' ':1}
    s=evt.key
    global stage ,decrease_rate         #
    if s in keyword:
        decrease_rate = decrease_rate + keyword[s] #根據他小於或大於1判斷印出什麼東西
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
scene.bind('keydown', keyinput) # setting for the binding function #壓下去的時候反應


#定義碰撞
def vcollision(a1p, a2p, a1v,a2v): # the function for handling velocity after collisions between two atoms
    v1prime = a1v - (a1p - a2p) * sum((a1v-a2v)*(a1p-a2p)) / sum((a1p-a2p)**2)
    v2prime = a2v - (a2p - a1p) * sum((a2v-a1v)*(a2p-a1p)) / sum((a2p-a1p)**2)
    return v1prime, v2prime

#給他初始值
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

while True:
    t += dt
    rate(10000)
    a_a = np.zeros((N,3))
    #term = 0
    length,height,width = container.length , container.height , container.width
   
   
    #粒子與六面不同牆壁碰撞
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
        
        #速度爆掉的時候把跑出去的東西抓回來  令成隨機速度 位置
        if sum(v_a[i]*v_a[i]) > 30*vrms**2:
            ra = pi*random()
            rb = 2*pi*random()
            v_a[i] =0.1*(6*random()+7)*np.array([vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]  )
            p_a[i]=np.array([2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L])
            print (i, " OVER~~~~~~~~~~~~~~")
    
        #位置超過三倍牆壁 就抓回箱子 令成隨機速度 位置
        if sum(p_a[i]*p_a[i]) > 3*(L**2):
            ra = pi*random()
            rb = 2*pi*random()
            v_a[i] =0.1*(6*random()+7)*np.array([vrms*sin(ra)*cos(rb), vrms*sin(ra)*sin(rb), vrms*cos(ra)]  )
            p_a[i]=np.array([2 * L*random() - L, 2 * L*random() - L, 2 * L*random() - L])
            print (i," GO AWAY~~~~~~~~~~~")

    observation.plot(data = np.sqrt(np.sum(np.square(v_a),-1))) #紅色長條圖隨時更新

    a_a[N-1] += np.array([0,-9.8,0])
    for i in range(N-1):            # 0~N-2
        a_a[i] += np.array([0,-9.8,0])            #consider gravity #每顆都受到y方向重力加速度
        for j in range(i+1,N):      # 算每顆分子間距離
            distance = sum((p_a[i]-p_a[j])**2)**0.5
            #只算一定距離間受力
            if distance > 5*sigma:
                continue
            if distance < 2*sigma:
                continue
            force = (force_on_particle(p_a[i],p_a[j])) #粒子間受力
            a_a[i] += force/m  #第i顆加速度 受到其他顆作用力
            a_a[j] += -1*force/m #反作用力
            rij = p_a[i]-p_a[j] #另外定義距離
            term += np.dot(rij,force) # this is the extra term according to eq.(2.12)
    v_a += a_a*dt
    p_a += v_a*dt # calculate new positions for all atoms
   
    #算kk次
    for i in range(N-1):
        for j in range(i+1,N):
            distance = sum((p_a[i]-p_a[j])**2)**0.5
            if distance <= sigma: ##(2**(1/6))*
                pass
                ##v_a[i],v_a[j] = vcollision(p_a[i],p_a[j],v_a[i],v_a[j])
            elif distance <= 2*sigma: ##(2**(1/6))*
                com_v  = 0.5*(v_a[i]+v_a[j]) #質心速度
                com_p  = 0.5*(p_a[i]+p_a[j]) #質心位置
                com_a  = 0.5*(a_a[i]+a_a[j]) #質心加速度
                u1=v_a[i]-com_v
                u2=v_a[j]-com_v
                p1=p_a[i]-com_p
                p2=p_a[j]-com_p
                a1=a_a[i]-com_a
                a2=a_a[j]-com_a
                for l in range(kk):
                    force = (force_on_particle(p1,p2))
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
        print(' ')
        print ("T=",temperature," P=",pressure," V=",volume)
        #使用LJ potential 會符合PV=NkT+1/3term
        print ("PV=",(pressure*volume)," NKT+term=",(N*k*temperature+term/3/(t0/dt))," NKT=",(N*k*temperature)) # count (t0/dt) times and then take average
        print("P(V^5/3)=",(pressure* ((volume)**gamma) ))
        print ("Vrms=",vrms)
        print ("stage=",stage," ,rate=",decrease_rate)
        print ("term= ",term/(t0/dt))
        print ("potential= ", potential)
        print (' ')
        #歸零
        potential = 0
        total_kinetic_energy=0
        forcex,forcey,forcez = 0,0,0
        time = time+t0 #每跑過一次就會加t0秒
        term = 0 ## use the average of the term
        theory_T = gcurve(color=color.orange) # for plot of the curve for the atom speed distribution
        #重新畫圖
        for v in arange(0.,1500.+dv,dv): # theoretical speed distribution
            theory_T.plot(pos=(v,(deltav/dv)*N*4.*pi*((m/(2.*pi*k*temperature))**1.5)*exp((-0.5*m*v**2)/(k*temperature))*(v**2)*dv))
        plots.append(theory_T)
        plots[-2].delete()          #倒數第二個理論曲線刪掉(舊的刪掉留新的)
        
        

