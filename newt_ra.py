from math import *

def nr_func(M, ecc):
	sol = newton(E - (ecc * sin(E)), M)
	return sol