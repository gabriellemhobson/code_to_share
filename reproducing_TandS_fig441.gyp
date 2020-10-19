# A script to reproduce figure 4.41 from Geodynamics: Third Edition, Turcotte and Schubert, 2014 

# Figure 4.41 shows temperature profiles at the center of a 2 m wide intrusion at specific times
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

animate = 0    # set animate=1 if you want to plot the solution over time and save frames as jpg to make a video

n = 256
y = np.linspace(0,20,n) # y is the horizontal axis, in m
Q = 8.8e9 # J m^-2
rho = 2900 # kg m^-3
c = 1.2 * 1000 # kJ kg^-1 K^-1 converted to J kg^-1 K^-1
kappa = 0.5 * 1e-6 # mm^2 s^-1 converted to m^2 s^-1

t1 = 100*24*60*60
t2 = 1000*24*60*60

def temp(Q,rho,c,kappa,y,t):
    T = (Q / (2*rho*c*np.sqrt(np.pi*kappa*t))) * (np.exp(-y**2 /(4*kappa*t)))
    return T


plt.plot(y,temp(Q,rho,c,kappa,y,t1))
plt.plot(y,temp(Q,rho,c,kappa,y,t2))

fs = 14
plt.legend(['t = 100 days','t = 1000 days'])
plt.xlabel(r'$y (m)$',fontsize=fs)
plt.ylabel(r'$(T - T_0) (K)$',fontsize=fs)
plt.xlim(0,20.0)
plt.ylim(0,400)
plt.xticks([0,5,10,15,20],fontsize=fs)
plt.yticks([0,100,200,300,400],fontsize=fs)
plt.show()

###################################################

# making a video of the solution changing with time

if animate == 1:
    
    def time(t):
        # convert time from days to seconds
        ts = t*24*60*60;
        return ts

    fs = 14
    t = np.arange(50,1001,1)

    frame = 0
    for k in range(951):
        frame += 1
        plt.cla()

        plt.title('t = {}'.format(t[k]))
        plt.plot(y,temp(Q,rho,c,kappa,y,time(t[k])))
        plt.xlim(0,20.0)
        plt.ylim(0,500)
        plt.xlabel(r'$y (m)$',fontsize=fs)
        plt.ylabel(r'$(T - T_0) (K)$',fontsize=fs)
        plt.xticks([0,5,10,15,20],fontsize=fs)
        plt.yticks([0,100,200,300,400,500],fontsize=fs)

        plt.savefig('/Users/ghobson/Documents/Homework/234/images_fig441_extra/fig441_plot{}.jpg'.format(frame), dpi = 300)