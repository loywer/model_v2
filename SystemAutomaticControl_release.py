import numpy as np
from matplotlib import pyplot as plt

#eps_theta = 0.01       # |(theta_now - theta_spec)| > eps_theta
dt = 0.002

kp_elev = 14.0
ki_elev = 2.5
kd_elev = 1.0

kp_eleron = 1.8
kd_eleron = 0.005



t = 0
T_elev = 5
T_eleron = 2.1

I = 0
elev_spec = 0
last = 0
class Control:

    def set_data (self, theta_spec, theta_now, gamma_spec, gamma_now, w0) :
        self.theta_spec = theta_spec
        self.theta_now = theta_now
        self.gamma_spec = gamma_spec
        self.gamma_now = gamma_now
        self.w0 = w0

    def get_data (self) :
        self.get_elev_and_theta_new(self.theta_spec,self.theta_now)
        self.get_GammaAngle_and_eleron_now(self.gamma_spec,self.gamma_now,self.w0)
        return self.elev_new, self.eleron_now


    def aperiodic_link (self, T) :
        dt = 0.002
        return (1 - np.exp(-dt/T))


    # Функция получения положения руля высоты и текущей перегрузки 
    #   theta_spec - заданное значение перегрузки 
    def get_elev_and_theta_new (self, theta_spec, theta_now) :
        global last
        global I
        global elev_spec
        
        self.transition_function_elev = self.aperiodic_link(T_elev)
        self.delta_theta = theta_now - theta_spec
        dd_theta = (self.delta_theta - last) / dt
        I = I + self.delta_theta * (ki_elev*dt) 
        last = self.delta_theta
        elev_spec =  I + dd_theta*kd_elev + self.delta_theta * kp_elev
        #delta_elev = elev_spec - elev_new
        self.elev_new = elev_spec 

        if (self.elev_new * 180.0/np.pi >= 26) :
             self.elev_new = 26/180.0*np.pi

        if (self.elev_new * 180.0/np.pi <= -28) :
            self.elev_new = -28/180.0*np.pi
        
        return self.elev_new

    # Функция получения положения элерона и крена самолета 
    #   gamma_spec - заданное значение крена самолета (желаемое)
    def get_GammaAngle_and_eleron_now (self, gamma_spec, gamma_now, w0) :
           
        self.transition_function_eleron = self.aperiodic_link(T_eleron)
        self.delta_gamma = gamma_spec - gamma_now
        #self.eleron_spec = self.delta_gamma * kp_eleron - 3*w0
        eleron_spec = self.delta_gamma * kp_eleron - kd_eleron * w0
        self.eleron_now = eleron_spec 

        if (self.eleron_now * 180.0/np.pi >= 20) :
            self.eleron_now = 20/180.0*np.pi

        if (self.eleron_now * 180.0/np.pi <= -15) :
            self.eleron_now = -15/180.0*np.pi

        return self.eleron_now


"""
control = Control()

theta = 1
theta_s = 0.5
a1 = control.get_elev_and_theta_new(theta, theta_s)
print(a1)

gamma1 = 0.52
gamma_2 = 0.65
w_0 = 0.01
a2 = control.get_GammaAngle_and_eleron_now(gamma1, gamma_2, w_0)
print(a2)
"""
