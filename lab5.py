import math
from matplotlib import pyplot as plt


# Точное решение
def u(x):
    return x * math.atan(x) - math.exp(x)


# Края отрезка, начальные условия, шаг сетки, сетка
a, b = 0, 1
u_0 = -1
derivative_u_0 = -1
h = 0.05
n = int(1 / h)
list_x = [i*h for i in range(n + 1)]


# Уравнение может быть разрешено относительно старшей производной ( u" = f(x,u,u') )
def f(x_, u_, der_u_):
    return (-2 * x_ * der_u_ + 2 * u_ + 2 + (1 - 2 * x_ - x_ ** 2) * math.exp(x_)) / (1 + x_ ** 2)


# Метод Эйлера, порядок точности 1
def euler(list_x_, u_0_, derivative_u_0_):
    h_ = h
    list_u_ = [u_0_]
    list_der_u_ = [derivative_u_0_]

    for i in range(len(list_x_) - 1):
        der_u_ = list_der_u_[i] + h_ * f(list_x_[i], list_u_[i], list_der_u_[i])
        u_ = list_u_[i] + h_ * list_der_u_[i]

        list_der_u_.append(der_u_)
        list_u_.append(u_)

    return list_u_


# Метод Адамса, порядок точности 3, D:\Downloads\adams.png
def adams(list_x_, u_0_, derivative_u_0_):
    h_ = h
    list_u_ = [u_0_]
    list_der_u_ = [derivative_u_0_]

    # # для Адамса 3-го порядка точности нужны еще 2 доп значения u и 2 доп значения der_u. Используем метод Эйлера
    # for i in range(2):
    #     der_u_ = list_der_u_[i] + h_ * f(list_x_[i], list_u_[i], list_der_u_[i])
    #     u_ = list_u_[i] + h_ * list_der_u_[i]
    #
    #     list_der_u_.append(der_u_)
    #     list_u_.append(u_)

    # для Адамса 3-го порядка точности нужны еще 2 доп значения u и 2 доп значения der_u. Используем метод Рунге-Кутта 4
    for i in range(2):
        k1_der_u = f(list_x_[i], list_u_[i], list_der_u_[i])
        k1_u = list_der_u_[i]

        k2_der_u = f(list_x_[i] + h / 2, list_u_[i] + h / 2 * k1_u, list_der_u_[i] + h / 2 * k1_der_u)
        k2_u = list_der_u_[i] + h / 2 * k1_der_u

        k3_der_u = f(list_x_[i] + h / 2, list_u_[i] + h / 2 * k2_u, list_der_u_[i] + h / 2 * k2_der_u)
        k3_u = list_der_u_[i] + h / 2 * k2_der_u

        k4_der_u = f(list_x_[i] + h, list_u_[i] + h * k3_u, list_der_u_[i] + h * k3_der_u)
        k4_u = list_der_u_[i] + h * k3_der_u

        list_der_u_.append(list_der_u_[i] + h / 6 * (k1_der_u + 2 * (k2_der_u + k3_der_u) + k4_der_u))
        list_u_.append(list_u_[i] + h / 6 * (k1_u + 2 * (k2_u + k3_u) + k4_u))

    # сам метод Адамса
    for i in range(len(list_x_) - 3):
        der_u_ = (list_der_u_[i + 2] + h_ * (23 / 12 * f(list_x_[i + 2], list_u_[i + 2], list_der_u_[i + 2]) -
                                             4 / 3 * f(list_x_[i + 1], list_u_[i + 1], list_der_u_[i + 1]) +
                                             5 / 12 * f(list_x_[i], list_u_[i], list_der_u_[i])))
        u_ = list_u_[i + 2] + h_ * (23 / 12 * list_der_u_[i + 2] - 4 / 3 * list_der_u_[i + 1] + 5 / 12 * list_der_u_[i])

        list_der_u_.append(der_u_)
        list_u_.append(u_)

    return list_u_


# Метод Рунге-Кутта, порядок точности 4, D:\Downloads\runge4.png , D:\Downloads\runge_general.png
def runge_kutte(list_x_, u_0_, derivative_u_0_):
    h_ = h
    list_u_ = [u_0_]
    list_der_u_ = [derivative_u_0_]

    for i in range(len(list_x_) - 1):
        k1_der_u = f(list_x_[i], list_u_[i], list_der_u_[i])
        k1_u = list_der_u_[i]

        k2_der_u = f(list_x_[i] + h_ / 2, list_u_[i] + h_ / 2 * k1_u, list_der_u_[i] + h_ / 2 * k1_der_u)
        k2_u = list_der_u_[i] + h_ / 2 * k1_der_u

        k3_der_u = f(list_x_[i] + h_ / 2, list_u_[i] + h_ / 2 * k2_u, list_der_u_[i] + h_ / 2 * k2_der_u)
        k3_u = list_der_u_[i] + h_ / 2 * k2_der_u

        k4_der_u = f(list_x_[i] + h_, list_u_[i] + h_ * k3_u, list_der_u_[i] + h_ * k3_der_u)
        k4_u = list_der_u_[i] + h_ * k3_der_u

        der_u_ = list_der_u_[i] + h_ / 6 * (k1_der_u + 2 * (k2_der_u + k3_der_u) + k4_der_u)
        u_ = list_u_[i] + h_ / 6 * (k1_u + 2 * (k2_u + k3_u) + k4_u)

        list_der_u_.append(der_u_)
        list_u_.append(u_)

    return list_u_


def runge_corrections(list_u_h_, list_u_h1_, p):  # list_u_h1_ - на половине h, list_u_h_ - на h
    corrections = []

    for i in range(len(list_u_h_)):
        n = 2 * i
        correction = (list_u_h1_[n] - list_u_h_[i]) / (2 ** p - 1)
        corrections.append(correction)

    return corrections


# Вычисления приближенных решений в узлах
list_u_exact = [u(x) for x in list_x]

list_euler = euler(list_x, u_0, derivative_u_0)
list_adams = adams(list_x, u_0, derivative_u_0)
list_runge_kutte = runge_kutte(list_x, u_0, derivative_u_0)
print(list_u_exact)
print(list_euler)
print(list_adams)
print(list_runge_kutte)
print()

# Логарифмы ошибок без поправок для фиксированного шага
list_err_euler = [abs(e - u) for e, u in zip(list_euler, list_u_exact)]
# list_err_euler = list(map(math.log, list_err_euler[1:]))
list_err_adams = [abs(e - u) for e, u in zip(list_adams, list_u_exact)]
# list_err_adams = list(map(math.log, list_err_adams[1:]))
list_err_runge_kutte = [abs(e - u) for e, u in zip(list_runge_kutte, list_u_exact)]
# list_err_runge_kutte = list(map(math.log, list_err_runge_kutte[1:]))

# print(list_err_euler)
# print(list_err_adams)
# print(list_err_runge_kutte)
# print()

# Графики решений
plt.figure('Solutions')
plt.subplot(131)
plt.title('Euler')
plt.plot(list_x, list_adams, 'red', list_x, list_u_exact, 'green')
plt.grid(True)

plt.subplot(132)
plt.title('Adams')
plt.plot(list_x, list_adams, 'red', list_x, list_u_exact, 'green')
plt.grid(True)

plt.subplot(133)
plt.title('Runge-Kutta')
plt.plot(list_x, list_runge_kutte, 'red', list_x, list_u_exact, 'green')
plt.grid(True)

# Поправки для Эйлера, Адамса, Рунге-Кутта
h1 = 0.1
list_x1 = [i * h for i in range(int((b - a) / h1) + 1)]
list_x2 = []

for j in range(len(list_x1) - 1):
    list_x2.append(list_x1[j])
    list_x2.append((list_x1[j] + list_x1[j + 1]) / 2)
list_x2.append(list_x1[len(list_x1) - 1])

list_u_exact_1 = [u(x) for x in list_x1]

list_euler_1 = euler(list_x1, u_0, derivative_u_0)
list_euler_2 = euler(list_x2, u_0, derivative_u_0)

list_adams_1 = adams(list_x1, u_0, derivative_u_0)
list_adams_2 = adams(list_x2, u_0, derivative_u_0)

list_runge_kutte_1 = runge_kutte(list_x1, u_0, derivative_u_0)
list_runge_kutte_2 = runge_kutte(list_x2, u_0, derivative_u_0)

euler_corrections = runge_corrections(list_euler_1, list_euler_2, 1)
adams_corrections = runge_corrections(list_adams_1, list_adams_2, 3)
runge_kutte_corrections = runge_corrections(list_runge_kutte_1, list_runge_kutte_2, 4)

list_euler_corrected = [list_euler_2[2 * j] + euler_corrections[j] for j in range(len(euler_corrections))]
list_adams_corrected = [list_adams_2[2 * j] + adams_corrections[j] for j in range(len(adams_corrections))]
list_runge_kutte_corrected = [list_runge_kutte_2[2 * j] + runge_kutte_corrections[j]
                              for j in range(len(runge_kutte_corrections))]

# Логарифмы ошибок до и после поправок для различных шагов
list_h = []

errors_runge_kutte_1 = []
errors_runge_kutte_2 = []

errors_adams_1 = []
errors_adams_2 = []

errors_euler_1 = []
errors_euler_2 = []

for j in range(1, 11):
    h2 = 0.005 * j
    number = int((b-a)/h2)
    list_x0 = [i * h for i in range(number + 1)]

    list_x1 = []
    for k in range(len(list_x0) - 1):
        list_x1.append(list_x0[k])
        list_x1.append((list_x0[k] + list_x0[k + 1]) / 2)
    list_x1.append(list_x0[len(list_x0) - 1])

    _list_u_exact = [u(x) for x in list_x0]

    list_runge_kutte_1_ = runge_kutte(list_x0, u_0, derivative_u_0)
    list_runge_kutte_2_ = runge_kutte(list_x1, u_0, derivative_u_0)
    runge_kutte_corrections_ = runge_corrections(list_runge_kutte_1_, list_runge_kutte_2_, 4)
    list_runge_kutte_corrected = [list_runge_kutte_2_[2 * k] + runge_kutte_corrections_[k]
                                  for k in range(len(runge_kutte_corrections_))]

    list_adams_1_ = adams(list_x0, u_0, derivative_u_0)
    list_adams_2_ = adams(list_x1, u_0, derivative_u_0)
    adams_corrections_ = runge_corrections(list_adams_1_, list_adams_2_, 3)
    list_adams_corrected = [list_adams_2_[2 * k] + adams_corrections_[k]
                            for k in range(len(adams_corrections_))]

    list_euler_1_ = euler(list_x0, u_0, derivative_u_0)
    list_euler_2_ = euler(list_x1, u_0, derivative_u_0)
    euler_corrections_ = runge_corrections(list_euler_1_, list_euler_2_, 1)
    list_euler_corrected = [list_euler_2_[2 * k] + euler_corrections_[k]
                            for k in range(len(euler_corrections_))]

    max_error_runge_kutte_1_ = max([abs(e_ - u) for e_, u in zip(list_runge_kutte_1_, _list_u_exact)])
    max_error_runge_kutte_2_ = max([abs(e_ - u) for e_, u in zip(list_runge_kutte_corrected, _list_u_exact)])

    max_error_adams_1_ = max(abs(e_ - u) for e_, u in zip(list_adams_1_, _list_u_exact))
    max_error_adams_2_ = max([abs(e_ - u) for e_, u in zip(list_adams_corrected, _list_u_exact)])

    max_error_euler_1_ = max([abs(e_ - u) for e_, u in zip(list_euler_1_, _list_u_exact)])
    max_error_euler_2_ = max([abs(e_ - u) for e_, u in zip(list_euler_corrected, _list_u_exact)])

    list_h.append(h2)
    errors_runge_kutte_1.append(max_error_runge_kutte_1_)
    errors_runge_kutte_2.append(max_error_runge_kutte_2_)
    errors_adams_1.append(max_error_adams_1_)
    errors_adams_2.append(max_error_adams_2_)
    errors_euler_1.append(max_error_euler_1_)
    errors_euler_2.append(max_error_euler_2_)

ln_list_h = list(map(lambda x: abs(math.log(x)), list_h))
ln_errors_euler_1 = list(map(lambda x: math.log(x) - math.log(errors_euler_1[len(errors_euler_1)-1]), errors_euler_1))
ln_errors_euler_2 = list(map(lambda x: math.log(x) - math.log(errors_euler_2[len(errors_euler_2)-1]), errors_euler_2))
ln_errors_adams_1 = list(map(lambda x: math.log(x) - math.log(errors_adams_1[len(errors_adams_1)-1]), errors_adams_1))
ln_errors_adams_2 = list(map(lambda x: math.log(x) - math.log(errors_adams_2[len(errors_adams_2)-1]), errors_adams_2))
ln_errors_runge_kutte_1 = list(map(lambda x: math.log(x) - math.log(errors_runge_kutte_1[len(errors_runge_kutte_1)-1]),
                                   errors_runge_kutte_1))
ln_errors_runge_kutte_2 = list(map(lambda x: math.log(x) - math.log(errors_runge_kutte_2[len(errors_runge_kutte_2)-1]),
                                   errors_runge_kutte_2))

plt.figure("Corrections")
plt.title("Зависимость ошибки от шага до и после поправок Рунге")

plt.subplot(131)
plt.title('Euler')
plt.plot(ln_list_h, ln_errors_euler_1, 'red')
plt.plot(ln_list_h, ln_errors_euler_2, 'blue')
plt.grid(True)

plt.subplot(132)
plt.title('Adams')
plt.plot(ln_list_h, ln_errors_adams_1, 'red')
plt.plot(ln_list_h, ln_errors_adams_2, 'blue')
plt.grid(True)

plt.subplot(133)
plt.title('Runge-Kutta')
plt.plot(ln_list_h, ln_errors_runge_kutte_1,  'red')
plt.plot(ln_list_h, ln_errors_runge_kutte_2, 'blue')
plt.grid(True)

plt.show()

