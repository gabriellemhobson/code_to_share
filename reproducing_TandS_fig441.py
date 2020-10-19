'''
A script to reproduce figure 4.41 from Geodynamics: Third Edition, Turcotte and Schubert, 2014 

Previously solved for the temperature distribution up until solidification of the dike. 
Now, solving 1D unsteady heat conduction problem in an infinite region with specified initial temperature
distribution, to get subsequent thermal history. 

The heat conduction equation is nondimensionalized and solved to get the temperature distribution T-T0,
and this temperature distribution is plotted at times t = 100 days and t = 1000 days in figure 4.41. 
The intrusive body under consideration in this example is a 2 m wide intrusion, and the temperature 
profiles at the center of the intrusion are plotted. That intrusion had a solidification time of 10.9 days, 
so our timescale must be much longer than that. 
'''

# importing required libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

animate = 0    # set animate=1 if you want to plot the solution over time and save frames as jpg to make a video

n = 256 # the number of points on the 1D grid
y = np.linspace(0,20,n) # y is the horizontal axis, in m
Q = 8.8e9 # J m^-2, initial heat content of the dike
rho = 2900 # kg m^-3
c = 1.2 * 1000 # kJ kg^-1 K^-1 converted to J kg^-1 K^-1
kappa = 0.5 * 1e-6 # mm^2 s^-1 converted to m^2 s^-1

def temp(Q,rho,c,kappa,y,t): # the function we wish to plot
    T = (Q / (2*rho*c*np.sqrt(np.pi*kappa*t))) * (np.exp(-y**2 /(4*kappa*t)))
    return T

def time(t): # convert time from days to seconds
    ts = t*24*60*60;
    return ts

plt.plot(y,temp(Q,rho,c,kappa,y,time(100)))
plt.plot(y,temp(Q,rho,c,kappa,y,time(1000)))

fs = 14 # set the font size 
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

    t = np.arange(50,1001,1) # timesteps we want to plot

    frame = 0
    for k in range(951):
        frame += 1
        plt.cla() # clear the plot

        plt.title('t = {}'.format(t[k]))
        plt.plot(y,temp(Q,rho,c,kappa,y,time(t[k])))
        plt.xlim(0,20.0)
        plt.ylim(0,500)
        plt.xlabel(r'$y (m)$',fontsize=fs)
        plt.ylabel(r'$(T - T_0) (K)$',fontsize=fs)
        plt.xticks([0,5,10,15,20],fontsize=fs)
        plt.yticks([0,100,200,300,400,500],fontsize=fs)

        plt.savefig('/Users/ghobson/Documents/Homework/234/images_fig441_extra/fig441_plot{}.jpg'.format(frame), dpi = 300)
