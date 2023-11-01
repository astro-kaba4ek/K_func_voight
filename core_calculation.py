# coding: utf-8

import numpy as np
import math
from scipy.integrate import quad
import mpmath as mpm
from datetime import datetime 
from matplotlib import pyplot as plt

eps = 0.00001 # погрешность -- 0.001%
a_0 = 1 # нулевой момент
a = 0.1 # параметр фойгтовского профиля 

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


# коэффициенты c^S_j (ф.34)
def c_S(a, j, b):
	f = 0
	for m in range(j+2):
		f += a_l(a, j+2-m) * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f


# коэффициенты c^E_j (ф.34)
def c_E(a, j, b):
	f = 0
	for m in range(j+2):
		f += a_l(a, j+1-m) * b**m / math.factorial(m) / math.factorial(j+1-m)

	f *= (-1)**j / (j+1)

	return f


# первый момент (ф.40)
a_1 = np.sqrt(2) * Phi(a*np.sqrt(2)) / np.pi / U_0(a)

# l=2
# x =  [i for i in range(4,20,1)]
# # x = [i for i in range(1,30,1)]
# for xx in x:
# 	an = Phi_anal(xx)
# 	asy18 = Phi_asymp(xx)
# 	# asy20 = E_12_asymp_20(xx)
# 	print(xx)
# 	print(an)
# 	print(asy18)
# 	# print(asy20)
# 	print(abs(an-asy18))
# 	# print(abs(an-asy20))
# 	print("-----------------")



# x = [i/100 for i in range(0, 1001, 25)]
t = [i/10 for i in range(0, 21, 2)]

# print(t)

for xx in 10*np.exp(t):
	print(xx, round(U(a,xx),10))
	# print(round(xx,3))
	# print(U(a,xx))
	# print(U2(a,xx))
	print("-----------------")

# x = [i for i in range(1,100)]
# u = [] # асимптотика ...
# u2 = [] # точная ---
# for xx in x:
# 	u.append(U(a,xx))
# 	u2.append( U2(a,xx))

# plt.plot(x, u, ".", x, u2, "-")
# plt.show()
# f = list(range(6,0,-2))
# print(math.prod(range(-1,0,-2)))
# print(math.prod(range(5,0,-2)))
# print(E_12(10))
# # x= 5
# now0 = datetime.now()
# print(Phi(2000))
# now = datetime.now()
# print(now-now0)

# now0 = datetime.now()
# print(Phi2(2000))
# now = datetime.now()
# print(now-now0)
# # print(funk(20))
# # print(quad((ff)**3, 1, 5, 4))
# print(a_l(1,3))
# # print(a_l2(1,3))
# # print(quad(a_l_int, -5, 3, args=(2,3)))
# x = 0.1
# print(a_0)
# print(a_l(a,0))
# print("------------")
# print(a_1)
# print(a_l(a,1))
# print("------------")
# print("------------")
# now0 = datetime.now()
# print(a_wave_beta(a,0,1000005))
# now = datetime.now()
# print(now-now0)
# print("------------")
# now0 = datetime.now()
# print(a_l(a,0)*np.log(1000005))
# now = datetime.now()
# print(now-now0)