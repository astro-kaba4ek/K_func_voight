# coding: utf-8
"""
## "Analytical" core functions

General formulas that are unsuitable for numerical calculations.\
Solved by built-in Python integrators.

### Contains the functions К, K_0, L, L_0
"""

import numpy as np
from scipy.integrate import dblquad
from constants import *
from calculation_functions import U_0, U


def K_analyt(a, t, b):
	r"""Analytical function K = K_11 (f.2).

	Parameters
	----------
	a : _float_
		`a` Task parameter
	t : _float_
		`\tau` The optical distance. Task parameter
	b : _float_
		`\beta` The ratio of absorption coefficients in the continuum and the center of the line. Task parameter

	Returns
	-------
	float
		`K(\tau,\beta)`
	"""

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: (U(a, x)/U_0_a)**2 * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	f *= U_0_a
	err *= U_0_a

	# f, err = dblquad(lambda y, x, a, t, b: U_0_a * (U(a, x)/U_0_a)**2 * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))


	if (abs(err/f) > eps):
		print(f"# The error in calculating the analytical function K_analyt is too large. K_analyt={f}, err={err}")


	return f	


def K_0_analyt(a, t, b):
	r"""Analytical function K_0 = K_10 (f.3).

	Parameters
	----------
	a : _float_
		`a` Task parameter
	t : _float_
		`\tau` The optical distance. Task parameter
	b : _float_
		`\beta` The ratio of absorption coefficients in the continuum and the center of the line. Task parameter

	Returns
	-------
	float
		`K_0(\tau,\beta)`
	"""

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x)/U_0_a * np.exp(-(U(a,x)/U_0_a + b) * t / y) / y, -np.inf, np.inf, 0, 1, args=(a,t,b))

	f *= U_0_a
	err *= U_0_a

	if (abs(err/f) > eps):
		print(f"# The error in calculating the analytical function K_0_analyt is too large. K_0_analyt={f}, err={err}")

	return f	


def L_analyt(a, t, b):
	r"""Analytical function L = K_21 (f.4).

	Parameters
	----------
	a : _float_
		`a` Task parameter
	t : _float_
		`\tau` The optical distance. Task parameter
	b : _float_
		`\beta` The ratio of absorption coefficients in the continuum and the center of the line. Task parameter

	Returns
	-------
	float
		`L(\tau,\beta)`
	"""

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: (U(a, x)/U_0_a)**2 / (U(a,x)/U_0_a + b) * np.exp(-(U(a,x)/U_0_a + b) * t / y), -np.inf, np.inf, 0, 1, args=(a,t,b))

	f *= U_0_a
	err *= U_0_a

	# f, err = dblquad(lambda y, x, a, t, b: U_0_a * (U(a, x)/U_0_a)**2 / (U(a,x)/U_0_a + b) * np.exp(-(U(a,x)/U_0_a + b) * t / y), -np.inf, np.inf, 0, 1, args=(a,t,b))

	if (abs(err/f) > eps):
		print(f"# The error in calculating the analytical function L_analyt is too large. L_analyt={f}, err={err}")

	return f	


def L_0_analyt(a, t, b):
	r"""Analytical function L_0 = K_20 (f.5).

	Parameters
	----------
	a : _float_
		`a` Task parameter
	t : _float_
		`\tau` The optical distance. Task parameter
	b : _float_
		`\beta` The ratio of absorption coefficients in the continuum and the center of the line. Task parameter

	Returns
	-------
	float
		`L_0(\tau,\beta)`
	"""

	U_0_a = U_0(a)

	f, err = dblquad(lambda y, x, a, t, b: U(a, x)/U_0_a / (U(a,x)/U_0_a + b) * np.exp(-(U(a,x)/U_0_a + b) * t / y), -np.inf, np.inf, 0, 1, args=(a,t,b))
	
	f *= U_0_a
	err *= U_0_a

	if (abs(err/f) > eps):
		print(f"# The error in calculating the analytical function L_0_analyt is too large. L_0_analyt={f}, err={err}")

	return f	