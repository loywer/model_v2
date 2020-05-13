import Engine as e
import SystemAutomaticControl_release as sys
import aerodym as aero
import atmosphere as atm
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
class Aircraft:
    def __init__(self,H,V,angle,dt):
        self.engine = e.Engine()
        self.atmos = atm.Atmosphere()
        self.atmos.set_H(H)
        self.aerodynamic = aero.Aerodunamic()
        self.control = sys.Control()
        
        
    def run(self,v_zad,gama_zad,theta_zad,dt):
       X,n,w,V,V_abs,gamma,theta,psi,alpha,betta =  self.aerodynamic.get_data()
       self.atmos.set_H(X[1])
       ro = self.atmos.get_density()
       g = self.atmos.get_accel_of_gravity()
       self.control.set_data(theta_zad,n[1],gama_zad,gamma,w[0])
       elevator,aileron = self.control.get_data()
       self.engine.Set_data(V_abs,v_zad,ro)
       P,M,omega = self.engine.Get_data()
       #P = np.array([5000,0,0])
       #aileron = 0.01
       #elevator = -0.07
       self.aerodynamic.set_data(elevator,aileron,0,P,M,ro,g)
       return n[1]

H=2000
angle=np.array([0,0,0])
V = np.array([50,0.0,0])
plane = Aircraft(H,V,angle,0.02)
T=100
t=0
X=[]
TT=[]
while(T>t):
    x=plane.run(55,0.0,1,0.002)
    t+=0.002

    X.append(x)
    TT.append(t)
#X=np.array(X)
plt.plot(TT,X)
plt.grid()
plt.show()
