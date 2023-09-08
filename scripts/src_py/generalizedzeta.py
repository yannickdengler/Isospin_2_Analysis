"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-5
    @Last Modified by: Yannick Dengler
    
    Calculation of the function Zeta^0_00 (1,q^2)
 """

import numpy as np
import scipy.integrate as scpint
import matplotlib.pyplot as plt
from decimal import Decimal
from decimal import getcontext

getcontext().prec = 40

def y_elem_Z3(order):
    y = []
    for i in [x for x in range(-order, order+1)]:
        for j in [x for x in range(-order, order+1) if abs(x)+abs(i) < order]:
            for k in [x for x in range(-order, order+1) if abs(x)+abs(i)+abs(j) < order]:
                y.append((i,j,k))
    return y

def y_squared(order):
    res = []
    for y in y_elem_Z3(order):
        res.append(y[0]*y[0]+y[1]*y[1]+y[2]*y[2])
    return res
    
def term1(q, order=10):
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-(y_2-q))/(y_2-q)
    return res

def factorial(x):
    res = 1
    for i in range(1, x+1):
        res = res*i
    return res
def term2(q, order=23):
    res = 0
    for k in range(order):
        res += (q**k)/(factorial(k)*(k-0.5))
    return np.pi**(3/2.)*res

def exp_3(u, order):
    res = 0
    for y_2 in y_squared(order):
        res += np.exp(-np.pi**2*y_2/u)
    return res-1

def e_u_q2(u, q, order):
    return np.exp(u*q)*(np.pi/u)**(3/2.)*exp_3(u, order)

def term3(q, order=7):
    return scpint.quad(func=e_u_q2, a=0, b=1, args=(q, order,))[0]


def Zeta(q, orders=(10,23,7)):
    return (term1(q, orders[0]) + term2(q, orders[1]) + term3(q, orders[2]))/(np.sqrt(4*np.pi))

# print(Zeta(q=-0.095900719))                                     # Test for zero-points of Zeta function
# print(Zeta(q=0.472894247))
# print(Zeta(q=1.441591313))
# print(Zeta(q=2.627007612))

# xarr = []
# yarr = []

# for q in np.linspace(-1, 10.2, 400):
#     print(q)
#     func = Zeta
#     xarr.append(q)
#     yarr.append(func(q))

# with open("Zeta.dat", "w") as f:
#     for i in range(len(xarr)):
#         f.write("%e\t%e\n"%(xarr4[i], yarr[i]))

# plt.plot(xarr, yarr, label = "Zeta")
# plt.ylim((-10,10))
# plt.grid()
# plt.xlabel("$q^2$")
# plt.ylabel("$Z(1,q^2)$")
# plt.title("Zeta function")
# plt.savefig("Zeta.pdf")
# plt.show()                                                  # Function is not independant on order. Check for the error somewhere!!!

