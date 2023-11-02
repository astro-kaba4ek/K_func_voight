# coding: utf-8

import numpy as np
from scipy.integrate import dblquad
from constants import *
from calculation_functions import U_0, U

def K_analyt(a, t, b):
	'''Analytical function K = K_11 (f.2)'''

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x)**2 * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	if (abs(err/f) > eps):
		print(f"Ошибка при вычислении аналитической функции K_analyt слишком велика. K_analyt={f / U_0_a}, err={err}")

	f /= U_0_a

	return f	


def K_0_analyt(a, t, b):
	'''Analytical function K_0 = K_10 (f.3)'''

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x) * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	if (abs(err/f) > eps):
		print(f"Ошибка при вычислении аналитической функции K_0_analyt слишком велика. K_0_analyt={f / U_0_a}, err={err}")

	return f	


def L_analyt(a, t, b):
	'''Analytical function L = K_21 (f.4)'''

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x)**2 / (U(a,x)/U_0_a + b) * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	if (abs(err/f) > eps):
		print(f"Ошибка при вычислении аналитической функции L_analyt слишком велика. L_analyt={f / U_0_a}, err={err}")

	f /= U_0_a
	
	return f	


def L_0_analyt(a, t, b):
	'''Analytical function L_0 = K_20 (f.5)'''

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x) / (U(a,x)/U_0_a + b) * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	if (abs(err/f) > eps):
		print(f"Ошибка при вычислении аналитической функции L_0_analyt слишком велика. L_0_analyt={f / U_0_a}, err={err}")

	return f	