# coding: utf-8

import argparse
import numpy as np
import subprocess
import matplotlib.pyplot as plt


def create_parser():
	r"""Creating an argument parser.

	Returns
	-------
	ArgumentParser
		A function for creating an object for parsing arguments of this program
	"""

	parser = argparse.ArgumentParser(
		description="""The program for calculating the core functions K, K_0, L and L_0 for the Voigt absorption profile \
					depending on the parameters a, \\beta and \\tau. \
			 		The program is not able to calculate values for large (depending on a and \\beta) \\tau, due to the peculiarities \
					of the formulas used and the overflow of variables. Use carefully when \\tau is greater than ~10. Accuracy ~1e-5.""",
    	epilog="""Program written as part of the course "Special Workshop in Theoretical Astrophysics" under the guidance of Dr. Nagirner. Author S.I. Laznevoi.""")
	
	parser.add_argument("-in_f", "--input", help="the full name of the input-file")
	parser.add_argument("-a", default=0.1, type=float, help="the value of parameter a. Default = 0.1")
	parser.add_argument("-b", "--beta", default=1, type=float, help="the value of parameter \\beta \n\
						(the ratio of absorption coefficients in the continuum and the center of the line). \n\
					 	Default = 1")
	parser.add_argument("-t", "--tau", default=1, type=float, help="the value of parameter \\tau (optical distance). Default = 1")
	parser.add_argument("-K_nl", "--core_func", nargs="+", choices=["K", "K_0", "L", "L_0"], default=["K"],
					  	help="the core functions that need to be calculated")
	parser.add_argument("-p_g", "--plot_graphs", action="store_const", const=True, 
					 	help="if this flag is enabled, graphs will be plotted")
	parser.add_argument("-a_f", "--analytical", action="store_const", const=True, 
					 	help='if this flag is enabled, the "analytical" functions will also be calculated. \
								They use general formulas for K, K_0, L and L_0, solved by the built-in numerical integration of NumPy. \
								At \\tau = 0 or \\beta = 0, strong increases in calculation time and looping are possible, use with caution')
	
	return parser


def create_a_b_t_arrays(args):
	r"""Reading the input-file, setting parameters from the console, or using default parameters.

	Parameters
	----------
	args : _Namespace_
		Arguments of main program

	Returns
	-------
	Array[float], Array[float], Array[float]
		`a_array, \beta_array, \tau_array` Arrays of task parameters 
	"""

	try:
		with open(args.input) as in_file:
			marker_01 = int(in_file.readline().split()[0])
			if not marker_01:
				a_start, a_finish, a_N = map(float, in_file.readline().split()[:3])
				a_array = np.linspace(a_start, a_finish, int(a_N))
			else:
				a_array = np.array(list(map(float, in_file.readline().split()[:marker_01])))

			marker_01 = int(in_file.readline().split()[0])
			if not marker_01:
				b_start, b_finish, b_N = map(float, in_file.readline().split()[:3])
				b_array = np.linspace(b_start, b_finish, int(b_N))
			else:
				b_array = np.array(list(map(float, in_file.readline().split()[:marker_01])))
			
			marker_01 = int(in_file.readline().split()[0])
			if not marker_01:
				t_start, t_finish, t_N = map(float, in_file.readline().split()[:3])
				t_array = np.linspace(t_start, t_finish, int(t_N))
			else:
				t_array = np.array(list(map(float, in_file.readline().split()[:marker_01])))

		print("# The input-file has been read successfully:\n")
	except:
		print("# The input-file was not found or specified.",
				"The values entered in the console or the default values are used:\n")
		a_array = np.array([args.a])
		b_array = np.array([args.beta])
		t_array = np.array([args.tau])

	return a_array, b_array, t_array


def get_stdout(command):
	"""Running third-party commands and saving the standard output.

	Parameters
	----------
	command : _str_
		The command to run in the console

	Returns
	-------
	str
		Standard output after processing the command
	"""

	res = subprocess.run(command, 
		      			stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
						encoding="utf-8", shell=True, cwd=".")
	
	return res.stdout


def save_results(res_dict, a_array, b_array, t_array):
	r"""Saving the results to files.

	Parameters
	----------
	res_dict : _dict[float]_
		Dictionary of solutions
	a_array : _Array[float]_
		`a_array` Arrays of task parameter
	b_array : _Array[float]_
		`\beta_array` Arrays of task parameter
	t_array : _Array[float]_
		`\tau_array` Arrays of task parameter
	"""

	print()
	print("# Create the results directory:")
	mkdir_res = get_stdout("mkdir results")
	print(f"# {mkdir_res}")

	print("# Saving results.")

	for k in res_dict.keys():
		with open(f"results/result_{k}.txt", "w") as out_file:
			out_file.write(f"Core function: {k} \n")
			if k == "K": 
				out_file.write(f"If the value of the optical distance \\tau = 0, then K_wave = K(\\tau,\\beta) + a_1*ln(\\tau) is calculated. \n \n")
			elif k == "K_0":
				out_file.write(f"If the value of the optical distance \\tau = 0, then K_0_wave = K_0(\\tau,\\beta) + ln(\\tau) is calculated. \n \n")
			elif k == "L_0":
				out_file.write(f"If the value of the optical distance \\beta = 0, then L_0_wave = L_0(\\tau,\\beta) - \delta(\\beta) is calculated. \n \n")				

			for a in a_array:
				out_file.write(f"Parameter_a: {a} \n")
				b_print = " ".join(str(b) for b in b_array)
				out_file.write(f"\\tau_\_\\beta {b_print} \n")
				for t in t_array:
					res_fun = res_dict.get(k)
					res_a = res_fun.get(a)

					res_t = []
					for b in b_array:
						res_b = res_a.get(b)
						res_t.append(res_b.get(t))

					t_print = " ".join(str(t) for t in res_t)

					out_file.write(f"{t} {t_print} \n")
				out_file.write("\n")


def plot_results(res_dict, t_array, fun, a, b):
	r"""Plotting and saving graphs.

	Parameters
	----------
	res_dict : _dict[float]_
		Dictionary of solutions
	t_array : _Array[float]_
		`\tau_array` Arrays of task parameter
	fun : _str_
		The type of function displayed
	a : _float_
		`a` Task parameter
	b : _float_
		`\beta` The ratio of absorption coefficients in the continuum and the center of the line. Task parameter
	"""

	# print()
	# print("# Create the results directory:")
	mkdir_gra = get_stdout("mkdir graphs")
	# print(f"# {mkdir_gra}")

	plt.close()
	plt.plot(t_array, res_dict.values(), ".")
	plt.grid()
	plt.xlabel(r"$\tau$") 
	plt.ylabel(f"${fun}$") 
	plt.title(f"$K_{{nl}} = {fun}. a = {a}. b = {b}$")
	plt.savefig(f"graphs/{fun}__a={a}__b={b}.png")