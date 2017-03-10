import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi, sqrt, acos
from scipy.integrate import odeint
import matplotlib.patches as patches

#ode functions
def derivs(k, time, params):
    x, y, vx, vy = k
    GM = params
    R = ((x**2) + (y**2))**0.5
    dvx_dt = -GM/(R**3)*x
    dvy_dt = -GM/(R**3)*y
    dx_dt = vx
    dy_dt = vy
    d_vect = [dx_dt, dy_dt, dvx_dt, dvy_dt]
    return d_vect



#1f Plot integrated orbits


#free fall
x0 = 1.0
y0 = 0
vx0 = 0
init_cond_ff = (x0, y0, vx0, 0) #x,y,vx,vy
GM = 1.0
params = GM
time_space_ff = np.linspace(0, 2.0*np.pi, 1000)
sol_ff = odeint(derivs, init_cond_ff, time_space_ff, args=(params,))
print sol_ff


#circular orbit
R0 = ((x0**2) + (y0**2))**0.5
vcir = (GM/R0)**0.5
Tcir = 2.0*np.pi*R0*((R0/GM)**(0.5))
init_cond_cir = (1, 0, 0, vcir) #x,y,vx,vy
time_space_cir = np.linspace(0, Tcir, 1000)
sol_cir = odeint(derivs, init_cond_cir, time_space_cir, args=(params,))

#Pericenter
e = 0.2
a = R0 / (1.0 + e)
b = sqrt((a**2)*(1-(e**2)))
params = GM
vperi = (9.0*GM/(5.0*R0))**0.5
Tperi = 1.23 * R0 * ((R0/GM)**0.5)
init_cond_peri = (1, 0, 0, vperi) #x,y,vx,vy
time_space_per = np.linspace(0, 80.0*Tperi, 1000)
sol_per = odeint(derivs, init_cond_peri, time_space_per, args=(params,))

#Apocenter
e = 0.2
params = GM
vapo = (4.0*GM/(5.0*R0))**0.5
Tapo = 1.85 * R0 * ((R0/GM)**0.5)
init_cond_apo = (1, 0, 0, vapo) #x,y,vx,vy
time_space_apo = np.linspace(0, 4.0*Tapo, 1000)
sol_apo = odeint(derivs, init_cond_apo, time_space_apo, args=(params,))

#Escape Vel
params = GM
vesc = (2.0*GM/R0)**0.5
Tesc = R0 * ((R0/(2.0*GM))**0.5)
init_cond_esc = (1, 0, 0, vesc) #x,y,vx,vy
time_space_esc = np.linspace(0, 20.0*Tesc, 1000)
sol_esc = odeint(derivs, init_cond_esc, time_space_esc, args=(params,))

#Twice Escape
params = GM
v2esc = (2.0*GM/R0)**0.5
T2esc = 0.5 * R0 * ((R0/(2.0*GM))**0.5)
init_cond_2esc = (1, 0, 0, v2esc) #x,y,vx,vy
time_space_2esc = np.linspace(0, 20.0*T2esc, 1000)
sol_2esc = odeint(derivs, init_cond_2esc, time_space_2esc, args=(params,))

#Plot (1f)
plt.figure(figsize=(10,10))
plt.axis('equal')
plt.plot(sol_ff[:,0],sol_ff[:,1], color = 'blue', label = 'free fall')
plt.plot(sol_cir[:,0], sol_cir[:,1], color = 'red', label = 'circular')
plt.plot(sol_per[:,0], sol_per[:,1], color = 'green', label = 'pericenter')
plt.plot(sol_apo[:,0], sol_apo[:,1], color = 'black', label = 'apocenter')
plt.plot(sol_esc[:,0], sol_esc[:,1], color = 'cyan', label = 'escape')
plt.plot(sol_2esc[:,0], sol_esc[:,1], color = 'pink', label = '2 X escape')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Integrated Orbits')
plt.legend(loc = 'lower left')

#1g Overplot analytic expression

#circular
fig1, ax= plt.subplots(figsize=(10,10))
c = patches.Circle(xy=(0,0), radius = 1.0, alpha = 0.2)
c.set_facecolor('blue' )
c.set_edgecolor('blue')
ax.add_artist(c)
ax.plot(sol_cir[:,0], sol_cir[:,1], color = 'black', label = 'numerical')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
leg0 = plt.legend([c], ['analytical'], loc = 'lower right')
leg = ax.legend(loc = 'lower left')
plt.title('Circular Orbit')
plt.gca().add_artist(leg0)

#Apocenter
a = R0 / (1.0 + e)
b = sqrt((a**2)*(1-(e**2)))
fig2, ax= plt.subplots(figsize=(10,10))
e1 = patches.Ellipse(xy=(1.0-a,0), width=2.0*a, height=2.0*b, alpha = 0.2)
e1.set_facecolor('blue' )
e1.set_edgecolor('blue')
ax.add_artist(e1)
ax.plot(sol_apo[:,0], sol_apo[:,1], color = 'black', label = 'numerical')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
leg0 = plt.legend([e1], ['analytical'], loc = 'lower right')
leg = ax.legend(loc = 'lower left')
plt.title('Apocenter Orbit')
plt.gca().add_artist(leg0)

#Pericenter
a = R0 / (1.0 - e)
b = sqrt((a**2)*(1-(e**2)))
fig2, ax= plt.subplots(figsize=(10,10))
e2 = patches.Ellipse(xy=(1.0-a,0), width=2.0*a, height=2.0*b, alpha = 0.2)
e2.set_facecolor('blue' )
e2.set_edgecolor('blue')
ax.add_artist(e2)
ax.plot(sol_per[:,0], sol_per[:,1], color = 'black', label = 'numerical')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
leg0 = plt.legend([e2], ['analytical'], loc = 'lower right')
leg = ax.legend(loc = 'lower left')
plt.title('Pericenter Orbit')
plt.gca().add_artist(leg0)


#1h plot r(t) and f(t)

#Free Fall
r_ff = np.zeros([1000])
f_ff =np.zeros([1000])
for i in range(len(sol_ff)):
    r_ff[i] = ((sol_ff[:,0][i])**2.0 + (sol_ff[:,1][i])**2.0)**(0.5)
    f_ff[i] = acos(sol_ff[:,0][i] / r_ff[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_ff, r_ff)
plt.xlabel('time')
plt.ylabel('r(t)')
plt.title('Free Fall r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_ff, f_ff)
plt.title('Free Fall f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')

#Circular
r_cir = np.zeros([1000])
f_cir =np.zeros([1000])
for i in range(len(sol_cir)):
    r_cir[i] = ((sol_cir[:,0][i])**2.0 + (sol_cir[:,1][i])**2.0)**(0.5)
    f_cir[i] = acos(sol_cir[:,0][i] / r_cir[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_cir, r_cir)
plt.xlabel('time')
plt.ylabel('r(t)')
plt.ylim(0.0,1.1)
plt.title('Circular r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_cir, f_cir)
plt.title('Circular f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')

#Apocenter
r_apo = np.zeros([1000])
f_apo =np.zeros([1000])
for i in range(len(sol_apo)):
    r_apo[i] = ((sol_apo[:,0][i])**2.0 + (sol_apo[:,1][i])**2.0)**(0.5)
    f_apo[i] = acos(sol_apo[:,0][i] / r_apo[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_apo, r_apo)
plt.xlabel('time')
plt.ylabel('r(t)')
#plt.ylim(0.0,1.1)
plt.title('Apocenter r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_apo, f_apo)
plt.title('Apocenter f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')

#Pericenter
r_peri = np.zeros([1000])
f_peri =np.zeros([1000])
for i in range(len(sol_per)):
    r_peri[i] = ((sol_per[:,0][i])**2.0 + (sol_per[:,1][i])**2.0)**(0.5)
    f_peri[i] = acos(sol_per[:,0][i] / r_peri[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_per, r_peri)
plt.xlabel('time')
plt.ylabel('r(t)')
plt.title('Pericenter r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_per, f_apo)
plt.title('Pericenter f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')

#Escape
r_esc = np.zeros([1000])
f_esc =np.zeros([1000])
for i in range(len(sol_esc)):
    r_esc[i] = ((sol_esc[:,0][i])**2.0 + (sol_esc[:,1][i])**2.0)**(0.5)
    f_esc[i] = acos(sol_esc[:,0][i] / r_esc[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_esc, r_esc)
plt.xlabel('time')
plt.ylabel('r(t)')
plt.title('Escape r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_esc, f_esc)
plt.title('Escape f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')

#Twice Escape
r_2esc = np.zeros([1000])
f_2esc =np.zeros([1000])
for i in range(len(sol_2esc)):
    r_2esc[i] = ((sol_2esc[:,0][i])**2.0 + (sol_2esc[:,1][i])**2.0)**(0.5)
    f_2esc[i] = acos(sol_2esc[:,0][i] / r_2esc[i])
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot(211)
plt.plot(time_space_2esc, r_2esc)
plt.xlabel('time')
plt.ylabel('r(t)')
plt.title('Twice Escape r(t)')
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(time_space_2esc, f_esc)
plt.title('Twice Escape f(t)')
plt.ylabel('f(t)')
plt.xlabel('time')
plt.show()