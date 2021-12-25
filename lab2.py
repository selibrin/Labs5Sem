import math
from sympy import *
import sympy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

h = 0.25  # шаг сетки
n = int(10 / h)
args = [(-5 + i * h) for i in range(n + 1)]

# Точная производная
z = symbols('z')
func = sympy.atan(sympy.log(z ** 2 + 1) + 1) ** 2
der_exact = diff(func)
second_der_exact = diff(der_exact)


# Определния функции и разностей
def f(x):  # x ∈ [-5, 5]
    return math.atan(math.log(x ** 2 + 1) + 1) ** 2


def derivative_central(x):
    return (f(x+h) - f(x-h)) / (2 * h)


def derivative_right(x):
    return (f(x + h) - f(x)) / h


# Первая производная через правую и центральные разности
list_der = [derivative_central(args[i]) for i in range(1, len(args) - 1)]
list_der.append((f(args[28]) - 4 * f(args[29]) + 3 * f(args[30])) / (2 * h))
list_der.insert(0, (-3 * f(args[0]) + 4 * f(args[1]) - f(args[2])) / (2 * h))

list_der_exact = [der_exact.subs(z, x) for x in args]
standard_deviation_1 = [math.sqrt(((list_der[i] - list_der_exact[i]) ** 2)) for i in range(len(list_der))]


# Вторая производная через центральные разности
# со вторым порядком точности в узлах и точная вторая производная
def second_derivative_central_2(x):
    return (f(x + h) - 2 * f(x) + f(x - h)) / h ** 2


list_second_der_2 = [second_derivative_central_2(args[i]) for i in range(1, len(args) - 1)]
list_second_der_2.insert(0, (f(-5) - 2 * f(-5+h) + f(-5+2*h)) / (h ** 2))
list_second_der_2.append((f(5) + f(5-2*h) - 2 * f(5 - h)) / h ** 2)


list_second_der_exact = [second_der_exact.subs(z, x) for x in args]
standard_deviation_2_2 = [abs((list_second_der_2[i] - list_second_der_exact[i]))
                          for i in range(len(list_second_der_2))]


# Вторая производная через центральные разности
# с четвертым порядком точности в узлах
def second_derivative_central_4(x):
    return (-f(x-2*h) + 16*f(x-h) - 30*f(x) + 16*f(x+h) - f(x+2*h)) / (12 * h ** 2)


new_args4 = [args[0]-2*h] + [args[0]-h] + args + [args[len(args)-1]+h] + [args[len(args)-1]+2*h]
list_second_der_4 = [second_derivative_central_4(new_args4[i]) for i in range(len(new_args4))]


list_second_der_exact_4 = [second_der_exact.subs(z, x) for x in new_args4]
standard_deviation_2_4 = [abs((list_second_der_4[i] - list_second_der_exact_4[i]))
                          for i in range(len(list_second_der_4))]


# # Вывод
# print('Arguments = ', args)
# print()
# print('NEW_Arguments = ', new_args4)
# print()
# print('First derivative in terms of right differences and central differences = ', list_der)
# print()
# print('First exact derivative = ', list_der_exact)
# print()
# print('Standard deviation of the second derivative in terms of central differences with 2nd grade of accuracy = ',
#       standard_deviation_2_2)
# print()
# print('Standard deviation of the second derivative in terms of central differences with 4th grade of accuracy = ',
#       standard_deviation_2_4)
# print()
# print('Second exact derivative = ', list_second_der_exact)


# Графики
plt.figure('First derivative and Standard deviation', figsize=(7, 7))

plt.subplot(211)
plt.plot(args, list_der, 'red', args, list_der_exact, 'green')
plt.title('First derivative')
plt.grid(True)

plt.subplot(212)
plt.title('Standard deviation')
plt.plot(args, standard_deviation_1)
plt.grid(True)


plt.figure('Second derivative and Standard deviation', figsize=(15, 10))

plt.subplot(221)
plt.plot(args, list_second_der_2, 'r', args, list_second_der_exact, 'green')
plt.title('Second derivative, 2nd accuracy')
plt.grid(True)

plt.subplot(222)
plt.plot(new_args4, list_second_der_4, 'r', new_args4, list_second_der_exact_4, 'green')
plt.title('Second derivative, 4th accuracy')
plt.grid(True)

plt.subplot(223)
plt.title('Standard deviation of the 2nd acc')
plt.plot(args, standard_deviation_2_2)
plt.grid(True)

plt.subplot(224)
plt.title('Standard deviation of the 4th acc')
plt.plot(new_args4, standard_deviation_2_4)
plt.grid(True)


# Ошибка
# Первая производная
list_err = []
list_h = []
for i in range(1, 26):
    mean_err = 0
    h = 0.01*i
    n = int(10 / h)
    test_args1 = [(-5 + i * h) for i in range(n + 1)]
    test_list_der = [derivative_central(test_args1[i]) for i in range(1, len(test_args1) - 1)]
    test_list_der.append((f(test_args1[len(test_args1) - 3]) - 4 * f(test_args1[len(test_args1) - 2])
                          + 3 * f(test_args1[len(test_args1) - 1])) / (2 * h))
    test_list_der.insert(0, (-3 * f(test_args1[0]) + 4 * f(test_args1[1]) - f(test_args1[2])) / (2 * h))
    test_list_der_exact = [der_exact.subs(z, x) for x in test_args1]

    for j in range(len(test_list_der)):
        mean_err += math.sqrt(((test_list_der[j] - test_list_der_exact[j]) ** 2))

    mean_err /= len(test_list_der)
    list_err.append(mean_err)
    list_h.append(h)


list_err = list(map(lambda x: abs(math.log(x)), list_err))
list_h = list(map(lambda x: abs(math.log(x)), list_h))

list_err = list(map(lambda x: x - list_err[len(list_err)-1], list_err))
list_h = list(map(lambda x: x - list_h[len(list_h)-1], list_h))

plt.figure('Error', figsize=(5, 5))
plt.subplot(111)
plt.plot(list_h, list_err)
plt.grid(True)


# Вторая производная 2 степени точности
list_err = []
list_h = []
for i in range(1, 26):
    mean_err = 0
    h = 0.01 * i
    n = int(10 / h)
    test_args1 = [(-5 + i * h) for i in range(n + 1)]

    test_list_second_der_2 = [second_derivative_central_2(test_args1[i]) for i in range(1, len(test_args1) - 1)]
    test_list_second_der_2.insert(0, (f(test_args1[0]) - 2 * f(test_args1[1]) + f(test_args1[2])) / (h ** 2))
    test_list_second_der_2.append((f(test_args1[len(test_args1) - 1]) + f(test_args1[len(test_args1) - 3])
                                   - 2 * f(test_args1[len(test_args1) - 2])) / h ** 2)

    test_list_second_der_exact = [second_der_exact.subs(z, x) for x in test_args1]

    for j in range(len(test_list_second_der_2)):
        mean_err += math.sqrt(((test_list_second_der_2[j] - test_list_second_der_exact[j]) ** 2))

    mean_err /= len(test_list_second_der_exact)
    list_err.append(mean_err)
    list_h.append(h)


list_err = list(map(lambda x: abs(math.log(x)), list_err))
list_h = list(map(lambda x: abs(math.log(x)), list_h))

list_err = list(map(lambda x: x - list_err[len(list_err)-1], list_err))
list_h = list(map(lambda x: x - list_h[len(list_h)-1], list_h))

plt.figure('Error1', figsize=(5, 5))
plt.subplot(111)
plt.plot(list_h, list_err)
plt.grid(True)


# Вторая производная 4 степень точности
list_err = []
list_h = []
for i in range(1, 26):
    mean_err = 0
    h = 0.01 * i
    n = int(10 / h)
    test_args2 = [(-5 + i * h) for i in range(n + 1)]

    new_test_args4 = [-5 - 2 * h] + [-5 - h] + test_args2 + [5 + h] + [5 + 2 * h]
    test_list_second_der_4 = [second_derivative_central_4(new_test_args4[i]) for i in range(0, len(new_test_args4))]

    test_list_second_der_exact = [second_der_exact.subs(z, x) for x in new_test_args4]

    for j in range(len(test_list_second_der_exact)):
        mean_err += math.sqrt(((test_list_second_der_4[j] - test_list_second_der_exact[j]) ** 2))

    mean_err /= len(test_list_second_der_exact)
    list_err.append(mean_err)
    list_h.append(h)

list_err = list(map(lambda x: abs(math.log(x)), list_err))
list_h = list(map(lambda x: abs(math.log(x)), list_h))

list_err = list(map(lambda x: x - list_err[len(list_err)-1], list_err))
list_h = list(map(lambda x: x - list_h[len(list_h)-1], list_h))

plt.figure('Error2', figsize=(5, 5))
plt.subplot(111)
plt.plot(list_h, list_err, 'g')
plt.grid(True)
plt.show()
