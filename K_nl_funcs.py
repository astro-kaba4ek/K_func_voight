# coding: utf-8
"""
## "Numerical" core functions

The functions are calculated using optimized formulas (series expansions, asymptotic approximations, etc.).

### Contains the functions К, K_0, L, L_0, К_wave, K_0_wave, L_wave, L_0_wave
"""

import numpy as np
from constants import *
from calculation_functions import *


def K_wave(a, t, b):
	r"""
	Series expansion of the K = K_11 function (f.30).\
	Used when \tau = 0.

	Returns: 
		K(\tau,\beta) + a_1 * ln(\tau)
	"""

	f = -a_1(a) * C_E - a_wave_beta(a, 1, b)
	
	return f


def K(a, t, b):
	r"""
	Series expansion of the K = K_11 function (f.30).

	Returns: 
		K(\tau,\beta)
	"""

	f = -a_1(a) * (np.log(t) + C_E) - a_wave_beta(a, 1, b)
	j = 0

	while True:
		s = c_S(a, j, b) * t**(j+1)

		if abs(s/f) < eps: break 
								
		f += s
		j += 1

	return f


def K_0_wave(a, t, b):
	r"""
	Series expansion of the K_0 = K_10 function (f.31).\
	Used when \tau = 0.

	Returns: 
		K_0(\tau,\beta) + ln(\tau)
	"""

	f = -C_E - a_wave_beta(a, 0, b)
	
	return f


def K_0(a, t, b):
	r"""
	Series expansion of the K = K_10 function (f.31).

	Returns: 
		K_0(\tau,\beta)
	"""

	f = -np.log(t) - C_E - a_wave_beta(a, 0, b)
	j = 0

	while True:
		s = c_S(a, j, b) * t**(j+1)

		if abs(s/f) < eps: break 
								
		f += s
		j += 1

	return f


def L_wave(a, t, b):
	r"""
	Series expansion of the L = K_21 function (f.32).\
	Used when \beta = 0.

	Returns: 
		L(\tau,\beta)
	"""

	if t == 0:
		f = 1
	else:
		f = 1 + (a_1(a)*(np.log(t) + C_E - 1) + a_wave_beta(a, 1, b)) * t
		j = 0

		while True:
			s = c_S(a, j, b) * t**(j+2) / (j+2)

			if abs(s/f) < eps: break 
									
			f -= s
			j += 1

	return f


def L(a, t, b):
	r"""
	Series expansion of the L = K_21 function (f.32).

	Returns: 
		L(\tau,\beta)
	"""

	if t == 0:
		f = 1 - b * delta(a, b)
	else:
		f = 1 - b * delta(a, b) + (a_1(a)*(np.log(t) + C_E - 1) + a_wave_beta(a, 1, b)) * t
		j = 0

		while True:
			s = c_S(a, j, b) * t**(j+2) / (j+2)

			if abs(s/f) < eps: break 
									
			f -= s
			j += 1

	return f


def L_0_wave(a, t, b):
	r"""
	Series expansion of the L_0 = K_20 function (f.33).\
	Used when \beta = 0.

	Returns: 
		L_0(\tau,\beta) - \delta(\beta)
	"""

	if t == 0:
		f = 0
	else:
		f = (np.log(t) + C_E - 1 + a_wave_beta(a, 0, b)) * t
		j = 0

		while True:
			s = c_E(a, j, b) * t**(j+2) / (j+2)

			if abs(s/f) < eps: break 
									
			f -= s
			j += 1

	return f


def L_0(a, t, b):
	r"""
	Series expansion of the L_0 = K_20 function (f.33).

	Returns: 
		L_0(\tau,\beta)
	"""

	if t == 0:
		f = delta(a, b)
	else:
		f = delta(a, b) + (np.log(t) + C_E - 1 + a_wave_beta(a, 0, b)) * t
		j = 0

		while True:
			s = c_E(a, j, b) * t**(j+2) / (j+2)

			if abs(s/f) < eps: break 
									
			f -= s
			j += 1

	return f
