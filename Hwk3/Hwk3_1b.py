import numpy as np
from math import *
import matplotlib.pyplot as plt
import rebound as rb
import random

def orb_fun(times, x_prime, y_prime, vx_prime, vy_prime):
	sim = rb.Simulation()
	sim.add(m = m1)
	sim.add(m = m2, a = 1)
	sim.move_to_com()
	sim.add(x=x_prime, y=y_prime, vx=vx_prime, vy=vy_prime)
	particles = sim.particles
    

	x_parray = np.zeros(len(times))
	y_parray = np.zeros(len(times))
	vx_parray = np.zeros(len(times))
	vy_parray = np.zeros(len(times))
	for i,time in enumerate(times):
		for p in sim.particles:
			sim.integrate(time, exact_finish_time=1)
			x_parray[i] = particles[2].x
			y_parray[i] = particles[2].y
			vx_parray[i] = particles[2].vx
			vy_parray[i] = particles[2].vy
	return x_parray, y_parray, vx_parray, vy_parray

m_sun = 1.0
m1 = 1.0*m_sun
m2 = 1.0e-3 * m_sun
G = 1.0
a = 1.0
mu = G * (m1 + m2)
n = sqrt(mu)
mu1 = G*m1
mu2 = G*m2
U_L1 = 1.52061452
U_L2 = 1.51967042
U_L3 = 1.50099899
U_L4 = 1.50000150148

#find x_m1 and x_m2 from barycenter
x_m1 = -(m2 / (m1+m2)) * a
x_m2 = (m1 / (m1+m2)) * a


#i
times = np.linspace(0.0, 2.0*pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = x_m2 + 0.1*RH
vx_prime = 0
y_prime = 0
vy_prime = sqrt(G*m2/(0.1*RH)) + (n * x_prime)
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle Inside Hill Sphere (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle Inside Hill Sphere (Rotating RF)')

#ii
times = np.linspace(0.0, 2.0*pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = x_m2 + 0.1*RH
vx_prime = 0
y_prime = 0
vy_prime = 1.2*sqrt(G*m2/(0.1*RH)) + (n * x_prime)
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle Inside Hill Sphere 1.2*v_rot (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle Inside Hill Sphere 1.2*v_rot (Rotating RF)')


#iii
times = np.linspace(0.0, 2.0*pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = x_m2 + 0.1*RH
y_prime = 0
r1 = sqrt(x_prime + mu2)**2.0
r2 = sqrt(x_prime - mu1)**2.0
U = (n**2.0)*(x_prime**2.0 + y_prime**2.0)/2.0 + G*m1/r1 + G*m2/r2
vx_prime = 0
vy_prime = 0.99*sqrt(2.0*(-U_L1 + U))
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('H = -U_L1 (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('H = -U_L1 (Rotating RF)')

#iv
times = np.linspace(0.0,  100.0 * 2 *pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = x_m2 + 0.1*RH
y_prime = 0
r1 = sqrt(x_prime + mu2)**2.0
r2 = sqrt(x_prime - mu1)**2.0
U = (n**2.0)*(x_prime**2.0 + y_prime**2.0)/2.0 + G*m1/r1 + G*m2/r2
vx_prime = 0
vy_prime = random.uniform(sqrt(2.0*(-U_L1 + U)),sqrt(2.0*(-U_L2 + U))) + (n * x_prime)
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('U_L1 < H < U_L2 (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('U_L1 < H < U_L2 (Rotating RF)')

#v
times = np.linspace(0.0, 700.0 * 2 *pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = x_m2 + 0.1*RH
y_prime = 0
r1 = sqrt(x_prime + mu2)**2.0
r2 = sqrt(x_prime - mu1)**2.0
U = (n**2.0)*(x_prime**2.0 + y_prime**2.0)/2.0 + G*m1/r1 + G*m2/r2
vx_prime = 0
vy_prime = random.uniform(sqrt(2.0*(-U_L2 + U)),sqrt(2.0*(-U_L3 + U))) + (n * x_prime)
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('U_L2 < H < U_L3 (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('U_L2 < H < U_L3 (Rotating RF)')

#vi
times = np.linspace(0.0, 700.0 * 2 *pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = (a/2.0) + x_m1 + 0.001
y_prime = sqrt((a**2) - (a/2.0)**2)
r1 = sqrt(((x_prime + mu2)**2.0) + y_prime**2.0)
r2 = sqrt(((x_prime - mu1)**2.0) + y_prime**2.0)
U = (n**2.0)*(x_prime**2.0 + y_prime**2.0)/2.0 + G*m1/r1 + G*m2/r2
vx_prime = -n*y_prime
vy_prime = n*x_prime
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Small Tadpole (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Small Tadpole (Rotating RF)')

#vii
times = np.linspace(0.0, 700.0 * 2 *pi, 10000)
RH = a * (m2/(3*m1))**(1./3.)
x_prime = (a/2.0) + x_m1 + 0.25
y_prime = sqrt(a**2 - x_prime**2)
r1 = sqrt(((x_prime + mu2)**2.0) + y_prime**2.0)
r2 = sqrt(((x_prime - mu1)**2.0) + y_prime**2.0)
U = (n**2.0)*(x_prime**2.0 + y_prime**2.0)/2.0 + G*m1/r1 + G*m2/r2
vx_prime = -n*y_prime
vy_prime = n*x_prime
x_parray, y_parray, vx_parray, vy_parray = orb_fun(times, x_prime, y_prime, vx_prime, vy_prime)

r = np.zeros([len(times)])
x = np.zeros([len(times)])
y = np.zeros([len(times)])
theta = np.zeros([len(times)])
for i in range(len(times)):
    r[i] = (x_parray[i]**2. + y_parray[i]**2.)**(1./2.)
    theta[i] = np.arctan2(y_parray[i], x_parray[i])
    x[i] = r[i]*cos(theta[i] - n * times[i])
    y[i] = r[i]*sin(theta[i] - n * times[i])

plt.figure(figsize = [10,10])
plt.plot(x_parray, y_parray)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Small Tadpole (Inertial RF)')

plt.figure(figsize = [10,10])
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Small Tadpole (Rotating RF)')
plt.show()


