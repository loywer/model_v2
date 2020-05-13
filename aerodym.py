import numpy as np
import c172 
from math import sin
from math import cos, sqrt
import matplotlib.pyplot as plt 
class Aerodunamic():
    def __init__(self):
        self.alpha = 0
        self.betta = 0 
        self.alpha_dot = 0
        self.V_abs = 50
        self.elevator = -0.0
        self.rudder = 0
        self.w= np.array([0,0,0])
        self.NSK_SSK=np.array([[0,0,0],[0,0,0],[0,0,0]])
        self.VSK_SSK=np.array([[0,0,0],[0,0,0],[0,0,0]])
        self.V=np.array([50,0,0])
        self.g = 9.81
        self.G=np.array([0,-self.g,0])
        self.X = np.array([0,500,0])
        self.gamma = 0
        self.theta = 0
        self.psi = 0
        self.P=np.array([00000,0,0])
        self.ro = 1.25
        self.aileron = 0.0
    def get_acceleration(self):
        Cx = c172.get_Cx(self.alpha,self.betta)
        Cy = c172.get_Cy(self.alpha,self.alpha_dot,self.elevator,self.w,self.V_abs)
        Cz = c172.get_Cz(self.alpha,self.betta,self.w,self.rudder,self.aileron,self.V_abs)
        q = self.ro*self.V_abs**2/(2*c172.mass)*c172.Sw
        F_aero = np.array([-Cx,Cy,Cz])*q
        buff =  np.dot(self.VSK_SSK,F_aero)

        F = F_aero-np.cross(self.w,self.V)+np.dot(self.NSK_SSK,self.G) +self.P/c172.mass
        return F
    def get_acceleration_angle(self):
        J = c172.inertia
        mx = c172.get_mx(self.alpha,self.betta,self.w,self.V_abs,self.aileron,self.rudder)
        my = c172.get_my(self.betta,self.w,self.V_abs,self.rudder,self.aileron)
        mz = c172.get_mz(self.alpha,self.alpha_dot,self.w,self.V_abs,self.elevator)
        q = self.ro*self.V_abs**2/(2.0)*c172.Sw
        M=np.array([mx*c172.b,my*c172.b,mz*c172.c])*q
        buff = np.cross(self.w,np.dot(J,self.w))
        e = np.dot(np.linalg.inv(J),M-buff)
        return e
    def get_NSK_SSK(self):
        result = np.array([[cos(self.theta)*cos(self.psi),sin(self.theta),-cos(self.theta)*sin(self.psi)],
                           [sin(self.gamma)*sin(self.psi)-cos(self.gamma)*sin(self.theta)*cos(self.psi),cos(self.gamma)*cos(self.theta),sin(self.gamma)*cos(self.psi)+cos(self.gamma)*sin(self.psi)*sin(self.theta)],
                           [cos(self.gamma)*sin(self.psi)+sin(self.gamma)*sin(self.theta)*cos(self.psi),-sin(self.gamma)*cos(self.psi),cos(self.gamma)*cos(self.psi)-sin(self.gamma)*sin(self.theta)*sin(self.psi)]])
        return result
    def get_VSK_SSK(self):
        result = np.array([[cos(self.betta)*cos(self.alpha),sin(self.alpha),-sin(self.betta)*cos(self.alpha)],
                           [-cos(self.betta)*sin(self.alpha), cos(self.alpha),sin(self.betta)*sin(self.alpha)],
                           [sin(self.betta),0,cos(self.betta)]])
        return result
    def AngleSpeed_Ailer (self ):
        
        Speed_teta = self.w[1]*sin(self.gamma) + self.w[2]*cos(self.gamma)
        Speed_gamma = self.w[0] - (self.w[1]*cos(self.gamma) - self.w[2]*sin(self.gamma)) * np.tan(self.theta)
        Speed_psi = 1/cos(self.theta) * (self.w[1]*cos(self.gamma) - self.w[2]*sin(self.gamma))

        return Speed_teta, Speed_gamma, Speed_psi
    def get_V_abs(self):
        return sqrt(self.V[0]**2+self.V[1]**2+self.V[2]**2)
    def get_alpha(self):
        self.alpha = -np.arctan2(self.V[1],self.V[0])
    def get_betta(self):
        self.betta = np.arctan2(self.V[2],self.V[0])
    def get_alpha_dot(self):
        self.apha_dot = (self.alpha-self.alpha_last)/0.002
    def Integrator(self,left,right,dt):
        return left+right*dt
    def set_data(self,elevator,aileron,rudder,P,M_d,ro,g):
        self.elevator = elevator
        self.aileron = aileron
        self.rudder = rudder
        self.P = P
        self.ro = ro
        self.g = g
        self.M_d = M_d
        self.G = np.array([0,-self.g,0])
    def get_data(self):
        self.alpha_last = self.alpha
        self.get_alpha()
        self.get_alpha_dot()
        self.get_betta()
        self.betta =-self.betta
        self.NSK_SSK=self.get_NSK_SSK()
        self.VSK_SSK=self.get_VSK_SSK()
        self.V_abs = self.get_V_abs()
        a = self.get_acceleration()
        self.n = a - np.dot(self.NSK_SSK,self.G) 
        self.n = a/self.g
        e = self.get_acceleration_angle()
        s_theta,s_gamma,s_psi = self.AngleSpeed_Ailer()
        self.V = self.Integrator(self.V,a,0.002)
        self.w = self.Integrator(self.w,e,0.002)
        self.X = self.Integrator(self.X,np.dot(self.NSK_SSK.T,self.V),0.002)
        self.theta = self.Integrator(self.theta,s_theta,0.002)
        self.gamma = self.Integrator(self.gamma,s_gamma,0.002)
        self.psi   = self.Integrator(self.psi,s_psi,0.002)
        return self.X,self.n,self.w,self.V,self.V_abs,self.gamma,self.theta,self.psi,self.alpha,self.betta

