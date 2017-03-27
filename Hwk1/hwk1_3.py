import matplotlib.pyplot as plt
import rebound as rb
from math import pi, acos, log10
import numpy as np



aa = 1.0
G = 1.0
mJ = 1.0e-3
msun = 1.0

sim = rb.Simulation()
sim.add(m=mJ)
sim.add(m=msun, a=aa, e=0.1)
sim.move_to_com()

times = np.linspace(0.1, 2.0*pi, 1000)
x_j = np.zeros(1000)
y_j = np.zeros(1000)
f_j = np.zeros(1000)
r_j = np.zeros(1000)
x_s = np.zeros(1000)
y_s = np.zeros(1000)
f_s = np.zeros(1000)
r_s = np.zeros(1000)
vx_j = np.zeros(1000)
vx_s = np.zeros(1000)
particles = sim.particles
for i,time in enumerate(times):
    for p in sim.particles:
		sim.integrate(time, exact_finish_time=1)
		x_j[i] = particles[1].x
		y_j[i] = particles[1].y
		f_j[i] = particles[1].f
		vx_j[i] = particles[1].vx
		r_j[i] = ((x_j[i]**2) + (y_j[i]**2))**(0.5)
		x_s[i] = particles[0].x
		y_s[i] = particles[0].y
		vx_s[i] = particles[0].vx
		r_s[i] = ((x_s[i]**2) + (y_s[i]**2))**(0.5)


#plot radial distance
fig = plt.figure(figsize = (10,10))
fig.subplots_adjust(hspace=0.25)
ax1 = plt.subplot(211)
plt.plot(times, r_j*(1.5e11), color='red', label='Jupiter to bary') #planet to bary
plt.title('Barycentric Distance of Jupiter/Sun')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(times, r_s*(1.5e11), color='blue', label='Sun to bary') #star to bary
plt.xlabel('Time (2*pi=1 period)')
plt.ylabel('R')
plt.legend()
plt.show()


#plot x(t)
fig = plt.figure(figsize = (10,10))
fig.subplots_adjust(hspace=0.25)
ax1 = plt.subplot(211)
plt.plot(times, x_j*(1.5e11), color='red', label='Jupiter') #planet to bary
plt.title('Barycentric Distance of Jupiter/Sun')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(times, x_s*(1.5e11), color='blue', label='Sun') #star to bary
plt.xlabel('Time (2*pi=1 period)')
plt.ylabel('x')
plt.legend()
plt.show()


#plot
fig = plt.figure(figsize = (10,10))
fig.subplots_adjust(hspace=0.25)
ax1 = plt.subplot(211)
plt.plot(times, vx_j*(1.5e11), color='red', label='Jupiter to bary') #planet to bary
plt.title('Barycentric Distance of Jupiter/Sun')
plt.legend()
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(times, vx_s*(1.5e11), color='blue', label='Sun to bary') #star to bary
plt.xlabel('Time (2*pi=1 period)')
plt.ylabel('vx')
plt.legend()
plt.show()