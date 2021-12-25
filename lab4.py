import math
import random
import sympy
from sympy import *
from PIL import Image

eps = [1e-3, 1e-6, 1e-9]
root = 7.57079579576335176657004106
a, b = 0, 10


def f(x):
    return math.atan(math.tanh(x)) - 0.5 * x + 3


z = Symbol('z')
func = sympy.atan(sympy.tanh(z)) - 0.5 * z + 3
derivative_f = diff(func)


def dichotomy_method(a0=a, b0=b, epsilon=eps[0]):
    n = 0
    left, right = a0, b0
    while abs(right - left) > epsilon:
        c = (left + right) / 2
        if f(left)*f(c) < 0:
            right = c
        else:
            left = c
        n += 1
    return c, n


def newton_method(a0=a, b0=b, function=f, epsilon=eps[0]):
    n = 0
    x = (a0 + 0.456*b0) / 2
    while abs(function(x)) > epsilon:
        derivative = derivative_f.subs(z, x)
        x = x - function(x) / derivative
        n += 1
    x = float(x)
    return x, n


def chord_method(a0=a, b0=b, function=f, epsilon=eps[0]):
    n = 0
    k1 = random.random()
    k2 = random.random()
    x1 = (a0 + k1 * b0) / 2
    x2 = (a0 + k2 * b0) / 2
    while abs(function(x2)) > epsilon:
        n += 1
        x_prev = x2
        x2 -= f(x2) * (x2 - x1) / (f(x2) - f(x1))
        x1 = x_prev

    return x2, n


# Вычисления
# eps = 1e-3
print('eps = 1e-3:')
dichotomy_root0, dichotomy_n0 = dichotomy_method()
newton_root0, newton_n0 = newton_method()
print('Dichotomy method: n = ', dichotomy_n0, '; approximate root ≈ ', dichotomy_root0, sep='')
print('Newton method: n = ', newton_n0, '; approximate root ≈ ', newton_root0, sep='')
print()

# eps = 1e-6
print('eps = 1e-6:')
dichotomy_root1, dichotomy_n1 = dichotomy_method(epsilon=eps[1])
newton_root1, newton_n1 = newton_method(epsilon=eps[1])
print('Dichotomy method: n = ', dichotomy_n1, '; approximate root ≈ ', dichotomy_root1, sep='')
print('Newton method: n = ', newton_n1, '; approximate root ≈ ', newton_root1, sep='')
print()

# eps = 1e-9
print('eps = 1e-9:')
dichotomy_root2, dichotomy_n2 = dichotomy_method(epsilon=eps[2])
newton_root2, newton_n2 = newton_method(epsilon=eps[2])
print('Dichotomy method: n = ', dichotomy_n2, '; approximate root ≈ ', dichotomy_root2, sep='')
print('Newton method: n = ', newton_n2, '; approximate root ≈ ', newton_root2, sep='')
print()

# exact root
print('Exact root =', root)
print()

# chord method
print('eps = 1e-3:')
chord_root0, chord_n0 = chord_method()
print('Chord method: n = ', chord_n0, '; approximate root ≈ ', chord_root0, sep='')
print()

print('eps = 1e-6:')
chord_root1, chord_n1 = chord_method(epsilon=eps[1])
print('Chord method: n = ', chord_n1, '; approximate root ≈ ', chord_root1, sep='')
print()

print('eps = 1e-9:')
chord_root2, chord_n2 = chord_method(epsilon=eps[2])
print('Chord method: n = ', chord_n2, '; approximate root ≈ ', chord_root2, sep='')
print()


print('Chord method\'s order of convergence = 1')
img = Image.open(r'D:\Downloads\order.jpg')
img.show()


# # Order of convergence
# arg1 = 0
# arg2 = 9.54
# arg = list([arg2 - f(arg2) * (arg2 - arg1) / (f(arg2) - f(arg1))])
# arg.append(arg[0] - f(arg[0]) * (arg[0] - arg2) / (f(arg[0]) - f(arg2)))
# for i in range(2, 5):
#     arg.append(arg[i-1] - f(arg[i-1]) * (arg[i-1] - arg[i-2]) / (f(arg[i-1]) - f(arg[i-2])))
#
# arg = list(set(arg))
# j = len(arg) - 1
# q = math.log(abs((arg[j] - arg[j-1]) / (arg[j-1] - arg[j-2]))) / math.log(abs((arg[j-1] - arg[j-2]) / (arg[j-2] -
#                                                                                                        arg[j-3])))
#
# print(q)
