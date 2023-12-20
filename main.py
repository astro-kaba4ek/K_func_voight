# coding: utf-8



from constants import *
from analytical_functions import *
from calculation_functions import *
from K_nl_funcs import *
import time
import os



print("kek")


# t = 1
# b = 1




with open("input.txt") as in_file:
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

	marker_01 = int(in_file.readline().split()[0])
	if not marker_01:
		a_start, a_finish, a_N = map(float, in_file.readline().split()[:3])
		a_array = np.linspace(a_start, a_finish, int(a_N))
	else:
		a_array = np.array(list(map(float, in_file.readline().split()[:marker_01])))

# print(b_array)
# print(t_array)
# print(a_array)

# if marker_01: print(marker_01)

# start = time.perf_counter()
# Kk = L(0.1, t, b)
# finish = time.perf_counter()
# print(Kk, finish-start)

# start = time.perf_counter()
# Kk1 = L_analyt(0.1, t, b)
# finish = time.perf_counter()
# print(Kk1, finish-start)
# ------------------------------------------------

# # ok
# start = time.perf_counter()
# Kk = K_wave(a, t, b)
# finish = time.perf_counter()
# print("K   :", Kk-a_1*np.log(t), "time:", finish-start)

# start = time.perf_counter()
# Kk1 = K_analyt(a, t, b)
# finish = time.perf_counter()
# print("K_an:", Kk1, "time:", finish-start)
# # ------------------------------------------------

# # ok
# start = time.perf_counter()
# Kk = K_0_wave(a, t, b)
# finish = time.perf_counter()
# print("K0   :", Kk-np.log(t), "time:", finish-start)

# start = time.perf_counter()
# Kk1 = K_0_analyt(a, t, b)
# finish = time.perf_counter()
# print("K0_an:", Kk1, "time:", finish-start)
# # ------------------------------------------------


# b=0.00000001
# # почти..
# start = time.perf_counter()
# Kk = delta(a, b)
# finish = time.perf_counter()
# print("L0   :", Kk, "time:", finish-start)

# start = time.perf_counter()
# Kk1 = delta2(a, b)
# finish = time.perf_counter()
# print("L0_an:", Kk1, "time:", finish-start)
# os.system('spd-say "биип"')
# # ------------------------------------------------

# # ok
# start = time.perf_counter()
# Kk = K_wave(a, t, b)
# finish = time.perf_counter()
# print("K   :", Kk-a_1*np.log(t), "time:", finish-start)

# start = time.perf_counter()
# Kk1 = K_analyt(a, t, b)
# finish = time.perf_counter()
# print("K_an:", Kk1, "time:", finish-start)
# # ------------------------------------------------


# t=3
# b=1
# start = time.perf_counter()
# Kk = L_0_wave(a, t, b)
# finish = time.perf_counter()
# print("L0   :", Kk, "time:", finish-start)
# print(delta(a,b))
# start = time.perf_counter()
# Kk1 = L_analyt(a, t, b)
# finish = time.perf_counter()
# print("L0_an:", Kk1, "time:", finish-start)
# os.system('spd-say "биип"')
# ------------------------------------------------


# 
# start = time.perf_counter()
# Kk = L(a, t, b)
# finish = time.perf_counter()
# print("L  :", Kk, "time:", finish-start)

# start = time.perf_counter()
# Kk1 = L_analyt(a, t, b)
# finish = time.perf_counter()
# print("L_an:", Kk1, "time:", finish-start)
# os.system('spd-say "биип"')
# # ------------------------------------------------

# Kk = K_wave(a, t, b)
# print(Kk -a_1*np.log(t))

# Kk1 = K_analit(a, t, b)
# print(Kk1)
# print(delta(a,1))
# print(delta(a,0.1)*0.1)
# print(delta(a,0.01)*0.01)
# print(delta(a,0.001)*0.001)
# print(delta(a,0.0001)*0.0001)
# print(delta(a,0.00001)*0.00001)
# print(delta(a,0.000001)*0.000001)
# print(delta(a,0.0000001)*0.0000001)
# print(delta(a,0.00000001)*0.00000001)
# print(delta(a,0.000000001)*0.000000001)
# print(delta(a,0.0000000001)*0.0000000001)
# print(delta(a,0.00000000001)*0.00000000001)
# print(delta(a,0.000000000001)*0.000000000001)
# print(delta(a,0)*0)
print("lol")
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



# # x = [i/100 for i in range(0, 1001, 25)]
# t = [i/10 for i in range(0, 21, 2)]

# # print(t)

# for xx in 10*np.exp(t):
# 	print(xx, round(U(a,xx),10))
# 	# print(round(xx,3))
# 	# print(U(a,xx))
# 	# print(U2(a,xx))
# 	print("-----------------")

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