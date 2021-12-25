import math
from sympy import *
import matplotlib.pyplot as plt

h = 0.1  # grid step
n = int(2 / h)
knots = [(-1 + i * h) for i in range(n + 1)]
m_knots = knots[0:(len(knots)-1)]


# Function definition
def f(x):  # x ∈ [-1, 1]
    return 1 / (x ** 3 + x + 10)


# Trapeze method
def rectangles_method(x):  # for [x(k), x(k+1)]
    return f(x) * h


# Rectangles method
def trapeze_method(x):  # for [x(k), x(k+1)]
    return (f(x) + f(x+h)) * h / 2


# Simpson method
def simpson_method(x):  # for [x(k), x(k+1)]
    return (f(x) + 4*f(x + h/2) + f(x+h)) * h / 6


# Antiderivative
def antiderivative(x):
    return math.log(abs(x+2))/13 - ln(x**2-2*x+5)/26 + (3*math.atan((2*x-2)/4))/26


# Exact integral
def exact_integral(x):  # for [x(k), x(k+1)]
    return antiderivative(x+h) - antiderivative(x)


trapeze = [trapeze_method(knot) for knot in m_knots]

rectangle = [rectangles_method(knot) for knot in m_knots]

simpson = [simpson_method(knot) for knot in m_knots]

exact = [exact_integral(knot) for knot in m_knots]

standard_deviation_trapeze = [abs((trapeze[i] - exact[i])) for i in range(len(m_knots))]
standard_deviation_rectangle = [abs((rectangle[i] - exact[i])) for i in range(len(m_knots))]
standard_deviation_simpson = [abs((simpson[i] - exact[i])) for i in range(len(m_knots))]

# Graphs
plt.figure('Approximations and their deviations')

plt.subplot(231)
plt.plot(m_knots, rectangle, 'red', m_knots, exact, 'green')
plt.title('Rectangles')
plt.grid(True)

plt.subplot(232)
plt.plot(m_knots, trapeze, 'red', m_knots, exact, 'green')
plt.title('Trapezes')
plt.grid(True)

plt.subplot(233)
plt.plot(m_knots, simpson, 'red', m_knots, exact, 'green')
plt.title('Simpson')
plt.grid(True)

plt.subplot(234)
plt.plot(m_knots, standard_deviation_rectangle)
plt.title('Rectangles deviation')
plt.grid(True)

plt.subplot(235)
plt.plot(m_knots, standard_deviation_trapeze)
plt.title('Trapezes deviation')
plt.grid(True)

plt.subplot(236)
plt.plot(m_knots, standard_deviation_simpson)
plt.title('Simpson deviation')
plt.grid(True)
plt.show()


# Ошибка
list_err_rectangle = []
list_err_trapeze = []
list_err_simpson = []
list_h = []
for i in range(1, 16):
    mean_err_trapeze = 0
    mean_err_rectangle = 0
    mean_err_simpson = 0
    h1 = 0.005*i
    n1 = int(2 / h1)
    test_knots = [(-1 + i * h1) for i in range(n1)]

    test_exact = [exact_integral(knot) for knot in test_knots]
    test_trapeze = [trapeze_method(knot) for knot in test_knots]
    test_rectangle = [rectangles_method(knot) for knot in test_knots]
    test_simpson = [simpson_method(knot) for knot in test_knots]

    for j in range(len(test_rectangle)):
        mean_err_trapeze += abs((test_trapeze[j] - test_exact[j]))
        mean_err_rectangle += abs((test_rectangle[j] - test_exact[j]))
        mean_err_simpson += abs((test_simpson[j] - test_exact[j]))

    mean_err_rectangle /= len(test_rectangle)
    mean_err_trapeze /= len(test_trapeze)
    mean_err_simpson /= len(test_simpson)

    list_err_rectangle.append(mean_err_rectangle)
    list_err_trapeze.append(mean_err_trapeze)
    list_err_simpson.append(mean_err_simpson)
    list_h.append(h1)


list_err_rectangle = list(map(lambda x: abs(math.log(x)), list_err_rectangle))
list_err_trapeze = list(map(lambda x: abs(math.log(x)), list_err_trapeze))
list_err_simpson = list(map(lambda x: abs(math.log(x)), list_err_simpson))
list_h = list(map(lambda x: abs(math.log(x)), list_h))

list_err_rectangle = list(map(lambda x: x - list_err_rectangle[len(list_err_rectangle)-1], list_err_rectangle))
list_err_trapeze = list(map(lambda x: x - list_err_trapeze[len(list_err_trapeze)-1], list_err_trapeze))
list_err_simpson = list(map(lambda x: x - list_err_simpson[len(list_err_simpson)-1], list_err_simpson))
list_h = list(map(lambda x: x - list_h[len(list_h)-1], list_h))


plt.figure('Error')
plt.subplot(131)
plt.title('Rectangle error')
plt.plot(list_h, list_err_rectangle)
plt.grid(True)

plt.subplot(132)
plt.title('Trapezes error')
plt.plot(list_h, list_err_trapeze)
plt.grid(True)

plt.subplot(133)
plt.title('Simpson error')
plt.plot(list_h, list_err_simpson)
plt.grid(True)
plt.show()

# min_err_rectangle = min(list_err_rectangle)
# min_err_trapeze = min(list_err_trapeze)
# min_err_simpson = min(list_err_simpson)
#
# index_min_err_rectangle = list_err_rectangle.index(min_err_rectangle)
# index_min_err_trapeze = list_err_trapeze.index(min_err_trapeze)
# index_min_err_simpson = list_err_simpson.index(min_err_simpson)
#
# opt_h_rectangle = list_h[index_min_err_rectangle]
# opt_h_trapeze = list_h[index_min_err_trapeze]
# opt_h_simpson = list_h[index_min_err_simpson]
#
# print('Optimal h for rectangle method = ', opt_h_rectangle)
# print('Optimal h for trapeze method = ', opt_h_trapeze)
# print('Optimal h for simpson method = ', opt_h_simpson)
