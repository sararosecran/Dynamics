from math import *
import numpy as np

def carttoels(x, y, z, vx, vy, vz, mu):
	R = ((x**2) + (y**2) + (x**2))**(0.5)
	v2 = (vx**2) + (vy**2) + (vz**2)
	hx = (y*vz - z*vy)
	hy = (z*vx - x*vz)
	hz = (x*vy - y*vx)
	h = ((hx**2) + (hy**2) + (hz**2)) ** 0.5

	a = ((2 / R) - (v2 / mu))**(-1.0)

	e = (1 - (h**2)/(mu * a))**0.5

	i = acos(hz/h)

	if hz > 0:
		Omega = asin(hx / (h * sin(i)))
	if hz < 0:
		Omega = asin(-hx / (h * sin(i)))

	f = acos((1/e) * (((a*(1 - (e**2)))/R) - 1))
    
	omega = asin(z / (R * sin(i))) - f

	return (a, e, i, omega, Omega, f)


def elstocart(a, e, i, omega, Omega, f, mu):
	R = (a * (1 + (e**2))) / (1 + (e * cos(f)))

	z = sin(omega + f) * (R * sin(i))

	x = R * ((cos(omega+f)/sec(Omega)) - sin(Omega) * sin(omega+f) * cos(i))

	Y = ((R**2) - (x**2) + (z**2))**(0.5)

	hx = (y*vz - z*vy)
	hy = (z*vx - x*vz)
	hz = (x*vy - y*vx)
	A = [[0 , -z, y],[z, 0, -x],[y, x, 0]]
	b = [[hx] , [hy] , [hz]]
	solved = np.linalg.solve(A,b)
	vx = solved[0]
	vy = solved[1]
	vz = solved[2]

	return(x, y, z, vx, vy, vz)