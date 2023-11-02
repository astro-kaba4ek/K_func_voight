# coding: utf-8

import numpy as np
import math
from scipy.integrate import quad
from constants import *



# почти функция ошибок (ф.24-25)
def Phi(x):
	if (abs(x) >= 15): # проверка на асимптотичность для ОЧЕНЬ БОЛЬШИХ х (ф.25)
		f = 1 / 2 / x
		m = 1
		while True:
			s = (-1)**m * math.prod(range(2*m-1, 0, -2)) / 2**(m+1) / x**(2*m+1)

			if abs(s/f) < eps: break 
				
			f += s
			m += 1
	else: # (ф.24)
		f = x / 2 * np.exp(x**2) * E_12(x**2)

	return f


# интегрально-показательная функция E_1/2 (ф.18,20-22)
def E_12(x):
	if (abs(x) >= 5): # проверка на асимптотичность для БОЛЬШИХ х (ф.20-22)
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

	else: # (ф.18)
		f = np.sqrt(np.pi / x)
		m = 0	
		while True:
			s = (-1)**m * x**m / math.factorial(m) / (m + 0.5)

			if abs(s/f) < eps: break 
				
			f -= s
			m += 1

	return f


# постоянная нормировки профиля (ф.38)
def U_0(a):
	return 2 * Phi(a) / np.pi


# функция Фойгта (ф.37,39)
def U(a, x):
	if (a == 0): # доплервоский профиль
		f = np.exp(-x**2) / np.sqrt(np.pi)
	else:
		if (abs(x) >= 1e2): # проверка на асимптотичность (ф.39)
			
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

		else: # (ф.37)
			f, err = quad(lambda y, a, x: np.exp(-y**2) / ((x-y)**2 + a**2), -np.inf, np.inf, args=(a,x))

			if (abs(err/f) > eps):
				print(f"Ошибка при вычислении функции Фойгта U({a},{x}) слишком велика. U={f * a / np.pi**(3/2)}, err={err}")
				
			f = f * a / np.pi**(3/2)

	return f


# численный расчет моментов профиля (ф.12)
def a_l(a, l):
	U_0_a = U_0(a)

	f, err = quad(lambda x, a, l: U(a,x)**(l+1), -np.inf, np.inf, args=(a,l))

	if (abs(err/f) > eps):
		print(f"Ошибка при вычислении момента профиля a_{l} слишком велика. a_{l}={f / U_0_a**l}, err={err}")

	f /= U_0_a**l

	return f


# численный расчет моментов профиля (ф.15)
def a_wave_beta(a, l, b):
	if (abs(b) >= 1e9): # асимптотика 
		f = a_l(a, l) * np.log(b)

	else:
		U_0_a = U_0(a)

		f, err = quad(lambda x, a, l, b: U(a,x)**(l+1) * np.log(U(a,x)/U_0_a+b), -np.inf, np.inf, args=(a,l,b))

		if (abs(err/f) > eps):
			print(f"Ошибка при вычислении момента профиля a_wave_{l} слишком велика. a_wave_{l}={f / U_0_a**l}, err={err}")

		f /= U_0_a**l

	return f


def c_S(a, j, b):
	'''
	Coefficient c^S_j for K and L (f.34)
	'''
	
	f = 0
	for m in range(j+2):
		f += a_l(a, j+2-m) * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f


def c_E(a, j, b):
	'''
	Coefficient c^E_j for K_0 and L_0 (f.34)
	'''

	f = 0
	for m in range(j+2):
		f += a_l(a, j+1-m) * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f


# первый момент (ф.40)
a_1 = np.sqrt(2) * Phi(a*np.sqrt(2)) / np.pi / U_0(a)


def K_wave(a, t, b):
	'''
	Series expansion of the K = K_11 function (f.30)
	Returns: 
		K(\\tau,\\beta) + a_1 * ln(\\tau)
	'''

	f = -a_1 * C_E - a_wave_beta(a, 1, b)
	j = 0

	while True:
		s = c_S(a, j, b) * t**(j+1)

		if abs(s/f) < eps: break 
								
		f += s
		j += 1

	return f


def K_0_wave(a, t, b):
	'''
	Series expansion of the K_0 = K_10 function (f.31).
	Returns: 
		K_0(\\tau,\\beta) + ln(\\tau)
	'''

	f = -C_E - a_wave_beta(a, 0, b)
	j = 0

	while True:
		s = c_E(a, j, b) * t**(j+1)

		if abs(s/f) < eps: break 
								
		f += s
		j += 1

	return f
