from orbconvert import *
from newt_ra import *
from scipy.optimize import newton
import matplotlib.pyplot as plt
from math import pi, cos

def func(x, ecc):
	return x - (ecc * sin(x))

# uses carttoels from orbconvert to convert form cartesian to orbital space
#mu = input('What is mu? ' )
#x = input('What is x? ')
#y = input('What is y? ')
#z = input('What is z? ')
#vx = input('What is x component of v? ')
#vy = input('What is y component of v? ')
#vz = input('What is z component of v? ')
#cart_convert = carttoels(x, y, z, vx, vy, vz, mu)
#print("a, e, i, omega, Omega, f", cart_convert)

# uses elstocart from orbconvert to converto orbital space to cartesian space
#a = input('What is a? ')
#e = input('What is e? ')
#i = input('What is i? ')
#omega = input('What is omega? ')
#Omega = input('What is Omega? ')
#f = input('What is f? ')
#orb_convert = elstocart(a, e, i, omega, Omega, f)
#print("x, y, z, vx, vy, vz", orb_convert)

#Solve Kepler's equation
#M = input('What is your guess for M? ')
#ecc =input('What is the eccentricity? ')
#E = newton(func, M, args = ((ecc),))
#print("Eccentric Anomoly: ",E)

#Solves Kepler's equation for Jupiter and Sun over full period
ecc =input('What is the eccentricity? ')
a = 1.5e11 #in meters
G = 6.674e-11 #graviational constant in m^3/kg*s^2
mJ = 1.9e27 #mass of Jupiter in kg
msun = 2.0e30 #mass of sun in kg
mu = G * (mJ +msun)
T = (4.0 * (pi**2) * (a**3))/mu #period
n = (2.0*pi)/T

period_space = np.linspace(0.000001, T, 1000)
M_space = n * period_space

E_array = np.zeros([1000])
for i in range(len(M_space)):
	E_array[i] = newton(func, M_space[i], args = ((ecc),))
	print(E_array[i])
	#E_array = np.append(E_array, newton(func, M_space[i], args = ((ecc),)))
r_array = np.array([])
for j in range(len(E_array)):
	r_array = np.append(r_array, a*(1 - (e * cos(E_array[j]))))

print(r_array)
plt.plot(period_space, r_array)
plt.show()

