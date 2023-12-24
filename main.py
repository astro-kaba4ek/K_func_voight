# coding: utf-8

from constants import *
from analytical_functions import *
from calculation_functions import *
from K_nl_funcs import *
from utils import *


parser = create_parser()
args = parser.parse_args()

a_array, b_array, t_array = create_a_b_t_arrays(args)

print(r"# a_array:    ", a_array)
print(r"# \beta_array:", b_array)
print(r"# \tau_array: ", t_array, end="")

res_fun_array = {}
for fun in args.core_func:
	res_a_array = {}

	for a in a_array:
		res_b_array = {}

		for b in b_array:
			res_t_array = {}
			print()

			for t in t_array:
				if fun == "K":
					if t == 0:
						print(r"# The value of the optical distance \tau = 0, therefore K_wave = K(\tau,\beta) + a_1*ln(\tau) is calculated.")
						res = K_wave(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, K_wave   = {res:.16f}", end="")
					else:
						res = K(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, K        = {res:.16f}", end="")

					if args.analytical: 
						res_an = K_analyt(a, t, b)
						print(rf", K_an   = {res_an:.16f}, dK   = {abs(res_an-res)}", end="")

				elif fun == "K_0":
					if t == 0: 
						print(r"# The value of the optical distance \tau = 0, therefore K_0_wave = K_0(\tau,\beta) + ln(\tau) is calculated.")
						res = K_0_wave(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, K_0_wave = {res:.16f}", end="")
					else:
						res = K_0(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, K_0      = {res:.16f}", end="")

					if args.analytical: 
						res_an = K_0_analyt(a, t, b)
						print(rf", K_0_an = {res_an:.16f}, dK_0 = {abs(res_an-res)}", end="")

				elif fun == "L":
					if b == 0: 
						# print(r"# The value of the coefficient \beta = 0, therefore L is calculated.")
						res = L_wave(a, t, b)
					else:
						res = L(a, t, b)
					print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, L        = {res:.16f}", end="")

					if args.analytical: 
						res_an = L_analyt(a, t, b)
						print(rf", L_an   = {res_an:.16f}, dL   = {abs(res_an-res)}", end="")

				elif fun == "L_0":
					if b == 0: 
						print(r"# The value of the coefficient \beta = 0, therefore L_0_wave = L_0(\tau,\beta) - \delta(\beta) is calculated.")
						res = L_0_wave(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, L_0_wave = {res:.16f}", end="")
					else:
						res = L_0(a, t, b)
						print(f"\n# a = {a:6.4f}, \\beta = {b:6.4f}, \\tau = {t:9.6f}, L_0      = {res:.16f}", end="")

					if args.analytical: 
						res_an = L_0_analyt(a, t, b)
						print(rf", L_0_an = {res_an:.16f}, dL_0 = {abs(res_an-res)}", end="")

				res_t_array[t] = res

			res_b_array[b] = res_t_array
			if args.plot_graphs: plot_results(res_t_array, t_array, fun, a, b)

		res_a_array[a] = res_b_array	

	res_fun_array[fun] = res_a_array	
				
print()
save_results(res_fun_array, a_array, b_array, t_array)
