import numpy as np
from math import sqrt
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


#function finds pseudo-potential of Lagrangian points
def pot(X, Y, n, m1, m2, x_m1, x_m2):
	G = 1.0
	U = (n**2)*(X**2 + Y**2)/2 + G*m1/np.sqrt((X - x_m1)**2 + (Y**2)) + G*m2/np.sqrt((X - x_m2)**2 + (Y**2))
	return U


m_sun = 1.0
m1 = 1.0*m_sun
m2 = 1.0e-3 * m_sun
G = 1.0
a = 1.0
mu = G * (m1 + m2)
n = sqrt((a**3.0) / mu)

#find x_m1 and x_m2 from barycenter
x_m1 = -(m2 / (m1+m2)) * a
x_m2 = (m1 / (m1+m2)) * a

#Hill radius
RH = a * (m2/(3.0*m1))**(1.0/3.0)

#Finding L1 position
func1 = lambda xL: (n**2)*xL - G*m1/(xL-x_m1)**2 + G*m2/(-xL+x_m2)**2
init_guess1 = x_m2 - RH
x_L1 = fsolve(func1, init_guess1)
print 'x_L1', x_L1

#Finding L2 position
func2 = lambda xL: (n**2)*xL - G*m1/(xL-x_m1)**2 - G*m2/(xL-x_m2)**2
init_guess2 = x_m2 + RH
x_L2 = fsolve(func2, init_guess2)
print 'x_L2', x_L2

#Finding L3 position
func3 = lambda xL: (n**2)*xL + G*m1/(-xL+x_m1)**2 + G*m2/(-xL+x_m2)**2
init_guess3 = -1.0
x_L3 = fsolve(func3, init_guess3)
print 'x_L3', x_L3

#L4
x_L4 = (a/2.0) - x_m1
y_L4 = sqrt((a**2.0) - ((a/2.0)**2.0))
print 'x_L4,y_L4', x_L4, y_L4

#L5
x_L5 = x_L4
y_L5 = -y_L4
print 'x_L5,y_L5', x_L5, y_L5





#Draw contours
x = np.linspace(-2.0, 2.0, 1000)
y = np.linspace(-2.0, 2.0, 1000)
X, Y = np.meshgrid(x, y)
U = pot(X, Y, n, m1, m2, x_m1, x_m2)
U_L1 = pot(x_L1, 0, n , m1, m2, x_m1, x_m2)
U_L2 = pot(x_L2, 0, n , m1, m2, x_m1, x_m2)
U_L3 = pot(x_L3, 0, n , m1, m2, x_m1, x_m2)
U_L4 = pot(x_L4, y_L4, n , m1, m2, x_m1, x_m2)
U_L5 = pot(x_L5, y_L5, n , m1, m2, x_m1, x_m2)
print U_L1, U_L2, U_L3, U_L4

plt.contour(X, Y, U, levels=np.sort([U_L5+0.002, U_L4+0.001, U_L3, U_L2, U_L1]))
plt.plot([x_m1, x_m2], [0,0], '.')
plt.plot([0], [0], '+')
plt.axis('image')
plt.axis([-1.5, 1.5, -1.5, 1.5])
plt.show()
