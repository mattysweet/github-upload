# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""


import numpy as np
from matplotlib import pyplot as plt
import time
from realmachinery import F,G,Jminus1
from numpy import pi,cos,sin


tick=time.perf_counter()

theta1_init=90
theta2_init=180
thetadot1_init=0
thetadot2_init=0
l1=1
l2=1
m1=3
m2=3
g=9.80665
h=0.001
duration=15

theta1_init=(theta1_init/360)*2*pi
theta2_init=(theta2_init/360)*2*pi


starter=np.array([[theta1_init],[theta2_init]])

p_0=np.array([[(m1+m2)*(l1**2)*thetadot1_init+m2*l1*l2*thetadot2_init*cos(theta1_init - theta2_init)],[
              0.5*m2*(l2**2)+m2*l1*l2*thetadot1_init*cos(theta1_init-theta2_init)]])


b=np.arange(1000*duration)
angles=np.array(starter)
lastangleset=starter
lastlastangleset=starter
p=p_0
x1s=np.array([])
y1s=np.array([])
x2s=np.array([])
y2s=np.array([])
#Es=np.array([])
count=np.array([])
for i in b:
    #find the next set of angles
    a=np.arange(15)
    bestguess=2*lastangleset-lastlastangleset
    
    for j in a:
        old=bestguess
        bestguess=bestguess-np.matmul(Jminus1(lastangleset[0,0],lastangleset[1,0],bestguess[0,0],bestguess[1,0],l1,l2,m1,m2,g,h),F(lastangleset[0,0],lastangleset[1,0],bestguess[0,0],bestguess[1,0],p[0,0],p[1,0],l1,l2,m1,m2,g,h))
        error=((bestguess[0,0]-old[0,0])**2 + (bestguess[1,0]-old[1,0])**2)**(1/2)
        
        if error<0.001:
            break
   
    #update p
    p=G(lastangleset[0,0],lastangleset[1,0],bestguess[0,0],bestguess[1,0],l1,l2,m1,m2,g,h)
    #store angle set
    x_1=l1*sin(bestguess[0,0])
    y_1=-l1*cos(bestguess[0,0])
    x_2=x_1+l2*sin(bestguess[1,0])
    y_2=y_1-l2*cos(bestguess[1,0])
    
    x1s=np.append(x1s,x_1)
    y1s=np.append(y1s,y_1)
    x2s=np.append(x2s,x_2)
    y2s=np.append(y2s,y_2)
    
    thetadot1=(bestguess[0,0]-lastangleset[0,0])/h
    thetadot2=(bestguess[1,0]-lastangleset[1,0])/h
    
    #E=0.5 * l1**2 *thetadot1**2 * (m1+m2) + m2*l1*l2*thetadot1*thetadot2*cos(bestguess[0,0]-bestguess[1,0]) +0.5* m2 *l2**2 *thetadot2**2 - (m1+m2)*g*l1*cos(bestguess[0,0]) - m2*g*l2*cos(bestguess[1,0])
    #Es=np.append(Es,E)
    #count=np.append(count,i)
    #plt.plot(x_1,y_1,'ko')
    
    #update angles
    lastlastangleset=lastangleset
    lastangleset=bestguess
    
tock=time.perf_counter()

crunching_time=tock-tick

print("Time taken to crunch angles: %1.4f" %crunching_time)  


#plt.plot(x1s,y1s,'ko',markersize=3)
#plt.plot(x2s,y2s,'ro',markersize=3)
    
    
    
#plt.plot(count,Es,'ko',markersize=2)
tick=time.perf_counter()
x1s_e=np.array([])
for i in range(int(x1s.size/40)):
    x1s_e=.append(x1s_e,x1s[i*40])
y1s_e=np.array([])
for i in range(int(y1s.size/40)):
    y1s_e=np.append(y1s_e,y1s[i*40])

x2s_e=np.array([])
for i in range(int(x2s.size/40)):
    x2s_e=np.append(x2s_e,x2s[i*40])
y2s_e=np.array([])
for i in range(int(y2s.size/40)):
    y2s_e=np.append(y2s_e,y2s[i*40])


#print(x2s_e)
import matplotlib.animation as animation


fig, ax = plt.subplots()



ax = plt.axis([-3,3,-3,3])

xs=[x2s_e[0]]
ys=[y2s_e[0]]
plt.gca().set_aspect('equal', adjustable='box')
redDot, = plt.plot([0,x1s_e[0],x2s_e[0]], [0,y1s_e[0],y2s_e[0]], 'ro-',)
greyLine, = plt.plot(xs, ys,'0.5',linewidth=1)

def animate(i):
    redDot.set_data([0,x1s_e[i],x2s_e[i]], [0,y1s_e[i],y2s_e[i]])
    x=x2s_e[i]
    y=y2s_e[i]
    xs.append(x)
    ys.append(y)
    greyLine.set_data(xs,ys)
    return redDot, greyLine,

# create animation using the animate() function
myAnimation = animation.FuncAnimation(fig, animate, frames=x2s_e.size,
                                      interval=40, blit=False, repeat=True)
    
plt.show()
#myAnimation.save('CLICKME.gif',writer='imagemagick') 



tock=time.perf_counter()

animation_time=tock-tick

print("Time taken to animate: %1.4f"%animation_time)




    



