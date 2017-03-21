from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from math import *


def derivs(k, time, params):
	x,  y, vx, vy = k
	GM = params
	R = ((x**2) + (y**2))**0.5
	dvx_dt = -GM/(R**3)*x
	dvy_dt = -GM/(R**3)*y
	dx_dt = vx
	dy_dt = vy
	d_vect = [dx_dt, dy_dt, dvx_dt, dvy_dt]
	return d_vect


#3b
x = -5.0
y = -1.0
GM = 1.0
vy = 0.0
b = 1.0
v_rel_un = sqrt(2.0*GM/b) #order unity perturber
v_rel_sm = 5.0 #small perturber
v_rel_bi = 0.5 #big perturber

#order unity pertubation
params = GM
init_cond_un = (x, y, v_rel_un, vy)
time_space_un = np.linspace(0, 10.0*2.0*b/v_rel_un, 1000)
sol_un = odeint(derivs, init_cond_un, time_space_un, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_un[:,0], sol_un[:,1])
plt.plot([0], [0], '+', label = 'mass')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Perturbed Mass (order unity)')

#small pertubation
params = GM
init_cond_sm = (x, y, v_rel_sm, vy)
time_space_sm = np.linspace(0, 10.0*2.0*b/v_rel_sm, 1000)
sol_sm = odeint(derivs, init_cond_sm, time_space_sm, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_sm[:,0], sol_sm[:,1])
plt.plot([0], [0], '+', label = 'mass')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Perturbed Mass (small)')

#big pertubation
params = GM
init_cond_bi = (x, y, v_rel_bi, vy)
time_space_bi = np.linspace(0, 10.0*2.0*b/v_rel_bi, 1000)
sol_bi = odeint(derivs, init_cond_bi, time_space_bi, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_bi[:,0], sol_bi[:,1])
plt.plot([0], [0], '+', label = 'mass')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Perturbed Mass (big)')


#3d
#Pericenter orbit added kick in +y
e = 0.2
R0 = 1.0
a = R0 / (1.0 + e)
b = sqrt((a**2)*(1-(e**2)))
params = GM
vperi = (9.0*GM/(5.0*R0))**0.5
Tperi = 1.23 * R0 * ((R0/GM)**0.5)
n = 2.0 * pi / Tperi
v_kick = n * a * e
init_cond_peri1 = (1, 0, 0, v_kick + vperi) #x,y,vx,vy
time_space_per = np.linspace(0, 80.0*Tperi, 1000)
sol_per1 = odeint(derivs, init_cond_peri1, time_space_per, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_per1[:,0], sol_per1[:,1], color = 'green', label = 'pericenter')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Added Guiding Center Kick (+y)')
plt.legend(loc = 'lower left')

#Pericenter orbit added kick in +x
e = 0.2
R0 = 1.0
a = R0 / (1.0 + e)
b = sqrt((a**2)*(1-(e**2)))
params = GM
vperi = (9.0*GM/(5.0*R0))**0.5
Tperi = 1.23 * R0 * ((R0/GM)**0.5)
n = 2.0 * pi / Tperi
v_kick = n * a * e
init_cond_peri2 = (1, 0, v_kick, vperi) #x,y,vx,vy
time_space_per = np.linspace(0, 80.0*Tperi, 1000)
sol_per2 = odeint(derivs, init_cond_peri2, time_space_per, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_per2[:,0], sol_per2[:,1], color = 'green', label = 'pericenter')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Added Guiding Center Kick (+x)')
plt.legend(loc = 'lower left')

#Pericenter orbit added kick in -y
e = 0.2
R0 = 1.0
a = R0 / (1.0 + e)
b = sqrt((a**2)*(1-(e**2)))
params = GM
vperi = (9.0*GM/(5.0*R0))**0.5
Tperi = 1.23 * R0 * ((R0/GM)**0.5)
n = 2.0 * pi / Tperi
v_kick = n * a * e
init_cond_peri3 = (1, 0, 0.0, vperi - v_kick) #x,y,vx,vy
time_space_per = np.linspace(0, 80.0*Tperi, 1000)
sol_per3 = odeint(derivs, init_cond_peri3, time_space_per, args=(params,))

#plot
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_per3[:,0], sol_per3[:,1], color = 'green', label = 'pericenter')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Added Guiding Center Kick (-y)')
plt.legend(loc = 'lower left')
plt.show()