import math
from matplotlib import pyplot as plt
import numpy as np


# Точные решение и его производная
def u(x):
    return x * math.atan(x) - math.exp(x)


# Условия
a, b = 0, 1
h = 0.05
n = int(1 / h)
list_x = np.linspace(0, 1, n + 1)
alpha1, beta1, gamma1 = 0, 1, -1
alpha2, beta2, gamma2 = 4, 1, -9.1644


def p(x):
    return (2 * x) / (1 + x * x)


def q(x):
    return -2 / (1 + x * x)


def f(x):
    return (2 + (1 - 2 * x - x * x) * math.exp(x)) / (1 + x * x)


# approximation (находим коэффициенты слау Au=f)
def find_constants(list_x_, alpha1_, alpha2_, beta1_, beta2_, gamma1_, gamma2_, p_, q_, f_, order=1):
    n_ = len(list_x_)
    h_ = list_x_[1] - list_x_[0]
    list_a_, list_b_, list_c_, list_f_ = [], [], [], []

    a0 = 0
    b0 = alpha1_ - beta1_ / h_
    c0 = beta1_ / h_
    f0 = gamma1_

    an = (-1) * beta2_ / h_
    bn = alpha2_ + beta2_ / h_
    cn = 0
    fn = gamma2_

    if order == 2:
        if beta1_ != 0 and beta2_ != 0:
            b0 = -2 + 2 * h_ * alpha1_ / beta1_ - p_(list_x_[0]) * np.power(h_, 2) * alpha1_ / beta1_ + q_(
                list_x_[0]) * np.power(
                h_, 2)
            c0 = 2
            f0 = f_(list_x_[0]) * np.power(h_, 2) + 2 * gamma1_ * h_ / beta1_ - p_(list_x_[0]) * np.power(h_,
                                                                                                          2) * gamma1_ / beta1_

            an = 2
            bn = -2 + (-2) * h_ * alpha2_ / beta2_ - p_(list_x_[n_ - 1]) * np.power(h_, 2) * alpha2_ / beta2_ + q_(
                list_x_[n_ - 1]) * np.power(h_, 2)
            fn = f_(list_x_[n_ - 1]) * np.power(h_, 2) - (
                    2 * h_ * gamma2_ / beta2_ + p_(list_x_[n_ - 1]) * np.power(h_, 2) * gamma2_ / beta2_)

    else:
        if beta1_ != 0:
            b0 = -2 + 2 * h_ * alpha1_ / beta1_ - p_(list_x_[0]) * np.power(h_, 2) * alpha1_ / beta1_ + q_(
                list_x_[0]) * np.power(h_, 2)
            c0 = 2
            f0 = f_(list_x_[0]) * np.power(h_, 2) + 2 * gamma1_ * h_ / beta1_ - p_(list_x_[0]) * np.power(h_,
                                                                                                          2) * gamma1_ / beta1_

        elif beta2_ != 0:
            an = 2
            bn = -2 + (-2) * h_ * alpha2_ / beta2_ - p_(list_x_[n_ - 1]) * np.power(h_, 2) * alpha2_ / beta2_ + q_(
                list_x_[n_ - 1]) * np.power(
                h_, 2)
            fn = f_(list_x_[n_ - 1]) * np.power(h_, 2) - (2 * h_ * gamma2_ / beta2_ + p_(list_x_[n_ - 1]) *
                                                          np.power(h_, 2) * gamma2_ / beta2_)

    list_a_.append(a0)
    list_b_.append(b0)
    list_c_.append(c0)
    list_f_.append(f0)

    for k in range(1, n_ - 1):
        ak = 1 / (h_ ** 2) - p_(list_x_[k]) / (2 * h_)
        bk = (-2) / (h_ ** 2) + q_(list_x_[k])
        ck = 1 / (h_ ** 2) + p_(list_x_[k]) / (2 * h_)
        fk = f_(list_x_[k])

        list_a_.append(ak)
        list_b_.append(bk)
        list_c_.append(ck)
        list_f_.append(fk)

    list_a_.append(an)
    list_b_.append(bn)
    list_c_.append(cn)
    list_f_.append(fn)

    return list_a_, list_b_, list_c_, list_f_


# tridiagonal matrix algorithm (решает слау)
def tma(a_plot, b_plot, c_plot, f_plot):
    n_ = len(a_plot)
    list_u_ = []
    list_A_ = []
    list_B_ = []

    list_A_.append(-c_plot[0] / b_plot[0])  # A0
    list_B_.append(f_plot[0] / b_plot[0])  # B0

    for i in range(1, n_ - 1):
        list_A_.append(-c_plot[i] / (b_plot[i] + a_plot[i] * list_A_[i - 1]))

    for i in range(1, n_):
        list_B_.append((f_plot[i] - a_plot[i] * list_B_[i - 1]) / (b_plot[i] + a_plot[i] * list_A_[i - 1]))

    list_A_.append(0)  # An
    list_u_.append(list_B_[n_ - 1])

    for i in range(n_ - 2, -1, -1):
        yi = list_B_[i] + list_A_[i] * list_u_[0]
        list_u_.insert(0, yi)

    return list_u_


# Вычисления для 1 и 2 порядков точности
list_a_order1, list_b_order1, list_c_order1, list_f_order1 = find_constants(list_x, alpha1, alpha2,
                                                                            beta1, beta2, gamma1, gamma2, p, q, f, 1)
list_a_order2, list_b_order2, list_c_order2, list_f_order2 = find_constants(list_x, alpha1, alpha2,
                                                                            beta1, beta2, gamma1, gamma2, p, q, f, 2)
u_exact = [u(node) for node in list_x]
u_order1 = tma(list_a_order1, list_b_order1, list_c_order1, list_f_order1)
u_order2 = tma(list_a_order2, list_b_order2, list_c_order2, list_f_order2)

# Графики решений
plt.figure('Solutions')
plt.subplot(121)
plt.title('Order = 1')
plt.plot(list_x, u_order1, 'red', list_x, u_exact, 'green')
plt.grid(True)

plt.subplot(122)
plt.title('Order = 2')
plt.plot(list_x, u_order2, 'red', list_x, u_exact, 'green')
plt.grid(True)
# plt.show()

# print(u_order1)
# print(u_exact)
# Ошибки
list_err1, list_err2, list_h = [], [], []
for i in range(5, 110, 10):
    h_ = 1 / i
    list_x_ = np.linspace(0, 1, i)
    list_a_order1, list_b_order1, list_c_order1, list_f_order1 = find_constants(list_x_, alpha1, alpha2,
                                                                                beta1, beta2, gamma1, gamma2, p, q, f,
                                                                                1)
    list_a_order2, list_b_order2, list_c_order2, list_f_order2 = find_constants(list_x_, alpha1, alpha2,
                                                                                beta1, beta2, gamma1, gamma2, p, q, f,
                                                                                2)
    u_exact = [u(node) for node in list_x_]
    u_order1 = tma(list_a_order1, list_b_order1, list_c_order1, list_f_order1)
    u_order2 = tma(list_a_order2, list_b_order2, list_c_order2, list_f_order2)

    err1 = max([abs(u_exact[i] - u_order1[i]) for i in range(len(u_exact))])
    err2 = max([abs(u_exact[i] - u_order2[i]) for i in range(len(u_exact))])

    list_h.append(np.log(h_))
    list_err1.append(np.log(err1))
    list_err2.append(np.log(err2))

print(list_err1)
print(list_err2)
print(list_h)


plt.figure('Errors')
plt.title('Логарифмы ошибок от логарифма шага')
plt.plot(list_h, list_err1, 'red', list_h, list_err2, 'green')
plt.grid(True)
plt.show()
