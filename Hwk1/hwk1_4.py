import matplotlib.pyplot as plt
import rebound as rb
import numpy as np
from math import *


aa = 0.051
G = 1.0
msun = 1.0
mJ = 1.0e-3

#msun = 2.0e30
#mJ = 1.0e27
#a =7.621e9

sim = rb.Simulation()
#sim.G = 6.674e-11
#sim
sim.add(m = 0.0477*mJ/sin(0.1))
sim.add(m=msun, a=aa, e=0.0, inc = 0.1)
sim.move_to_com()
times = np.linspace(0.1, 2.0*pi, 1000)
particles = sim.particles


vr = np.zeros(1000)
for i,time in enumerate(times):
	sim.integrate(time, exact_finish_time=1)
	for o in sim.calculate_orbits(heliocentric=True):
		print o.inc, o.a, o.omega, o.f, o.e
        vr[i] = mJ/(mJ+msun) * sin(o.inc) * o.a / sqrt(1-(o.e)**2) * (cos(o.omega + o.f) + o.e*cos(o.omega))

fig = plt.figure(figsize = (10,10))
plt.plot(times, vr*(1.5e11))
plt.title('Radial Velocity of 51 Peg b')
plt.xlabel('Time')
plt.ylabel('V_r (m/s)')
plt.show()


