# coding: utf-8
"""
## Auxiliary functions for numerical calculation
"""

import numpy as np
import math
from scipy.integrate import quad
from constants import *


def Phi(x):
	r""""Almost" error function (f.24-25)"""

	if (abs(x) >= 15): # asymptotic testing for VERY LARGE x (f.25)
		f = 1 / 2 / x
		m = 1
		while True:
			s = (-1)**m * math.prod(range(2*m-1, 0, -2)) / 2**(m+1) / x**(2*m+1)

			if abs(s/f) < eps: break 
				
			f += s
			m += 1
	else: # (f.24)
		f = x / 2 * np.exp(x**2) * E_12(x**2)

	return f


def E_12(x):
	r"""Integral exponential function E_1/2 (f.18,20-22)"""

	if (abs(x) >= 5): # asymptotic testing for LARGE x (f.20-22)
		n = 1
		P = [0, 1]  # [P_0, P_1]
		Q = [1, x]  # [Q_0, Q_1]
		R1 = P[1] / Q[1]
		while True:
			n += 1
			if (n % 2 == 0):
				a = 1; b = n/2 - 1/2
			else:
				a = x; b = n // 2

			p = b * P[0] + a * P[1]
			q = b * Q[0] + a * Q[1]

			R = p / q 

			if abs((R-R1)/R1) < eps: break

			R1 = R
			P[0] = P[1]; P[1] = p
			Q[0] = Q[1]; Q[1] = q

		f = R1 * np.exp(-x)

	else: # (f.18)
		f = np.sqrt(np.pi / x)
		m = 0	
		while True:
			s = (-1)**m * x**m / math.factorial(m) / (m + 0.5)

			if abs(s/f) < eps: break 
				
			f -= s
			m += 1

	return f


def U_0(a):
	r"""Profile normalization constant (f.38)"""

	return 2 * Phi(a) / np.pi


def a_1(a):
	r"""The first moment (f.40)"""

	f = np.sqrt(2) * Phi(a*np.sqrt(2)) / np.pi / U_0(a)

	return f


def U(a, x):
	r"""Voigt function (f.37,39)"""

	if (a == 0): # Doppler profile
		f = np.exp(-x**2) / np.sqrt(np.pi)
	else:
		if (abs(x) >= 1e2): # testing for asymptoticism (f.39)
			
			f = 1
			n = 1
			while True:
				s0 = 0
				for k in range(n+1):
					s0 += (-1)**k * a**(2*n) / math.factorial(2*k+1) / math.factorial(n-k) / 2**(2*(n-k))

				s = math.factorial(2*n+1) * s0 / x**(2*n)

				if abs(s/f) < eps: break 
								
				f += s
				n += 1

			f = f * a / np.pi / x**2

		else: # (f.37)
			f, err = quad(lambda y, a, x: np.exp(-y**2) / ((x-y)**2 + a**2), -np.inf, np.inf, args=(a,x))

			if (abs(err/f) > eps):
				print(f"The error in calculating the Voigt function U({a},{x}) is too large. U={f * a / np.pi**(3/2)}, err={err}")
				
			f = f * a / np.pi**(3/2)

	return f


def a_l(a, l):
	r"""Calculation of the profile moment a_l (f.12)"""

	if l == 0:
		f = a_0
	elif l == 1:
		f = a_1(a)
	else:
		U_0_a = U_0(a)

		f, err = quad(lambda x, a, l: U(a,x)**(l+1), -np.inf, np.inf, args=(a,l))

		if (abs(err/f) > eps):
			print(f"The error in calculating the moment of the a_{l} profile is too large. a_{l}={f / U_0_a**l}, err={err}")

		f /= U_0_a**l

	return f


def a_wave_beta(a, l, b):
	r"""Calculation of the profile moment a_wave_l_beta (f.15)"""

	if (abs(b) >= 1e9): # asymptotics 
		f = a_l(a, l) * np.log(b)

	else:
		U_0_a = U_0(a)

		f, err = quad(lambda x, a, l, b: U(a,x)**(l+1) * np.log(U(a,x)/U_0_a+b), -np.inf, np.inf, args=(a,l,b))

		if (abs(err/f) > eps):
			print(f"The error in calculating the moment of the a_wave_{l} profile is too large. a_wave_{l}={f / U_0_a**l}, err={err}")

		f /= U_0_a**l

	return f


def delta_analyt(a, b):
	r"""Calculation of the profile moment \delta(\beta) (f.14')"""

	if (abs(b) >= 1e9): # asymptotics 
		f = a_1(a)
	else:
		U_0_a = U_0(a)

		f, err = quad(lambda x, a, b: U(a,x) / (U(a,x)/U_0_a+b), -np.inf, np.inf, args=(a,b))

		if (abs(err/f) > eps):
			print(f"The error in calculating the moment of the \\delta(\\beta) profile is too large. \\delta(\\beta)={f}, err={err}")

	return f


def delta(a, b):
	r"""Numerical calculation of the profile moment \delta(\beta) (f.41)"""

	t0 = 100

	if (abs(b) >= 1e9): # asymptotics 
		f = a_1(a)
	else:
		U_0_a = U_0(a)

		f1, err1 = quad(lambda x, a, b: U(a,x) / (U(a,x) + b*U_0_a), 0, 10, args=(a,b))

		if (abs(err1/f1) > eps):
			print(f"The error in calculating f1 (первое слагаемое \\delta2) is too large. f1={f1}, err={err1}")

		f2, err2 = quad(lambda t, a, b: U(a,10*np.exp(t)) * np.exp(t) / (U(a,10*np.exp(t)) + b*U_0_a), 0, t0, args=(a,b))

		if (abs(err2/f2) > eps):
			print(f"The error in calculating f2 (второе слагаемое \\delta2) is too large. f2={f2}, err={err2}")

		f = 2 * U_0_a * (f1 + 10*f2 + np.sqrt(a/b/U_0_a) * np.arctan(np.sqrt(a) / np.sqrt(b*U_0_a) / 10 / np.exp(t0)))

	return f


def c_S(a, j, b):
	r"""Coefficient c^S_j for K and L (f.34)"""

	global a_l_arr
	
	f = 0
	for m in range(j+2):
		l = f"a_{j+2-m}"
		if l in a_l_arr:
			al = a_l_arr.get(l)
		else:
			al = a_l(a, j+2-m)
			a_l_arr[l] = al
			
		f += al * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f


def c_E(a, j, b):
	r"""Coefficient c^E_j for K_0 and L_0 (f.34)"""

	global a_l_arr

	f = 0
	for m in range(j+2):
		l = f"a_{j+2-m}"
		if l in a_l_arr:
			al = a_l_arr.get(l)
		else:
			al = a_l(a, j+1-m)
			a_l_arr[l] = al

		f += al * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f

