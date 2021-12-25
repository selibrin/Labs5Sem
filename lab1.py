import math
import matplotlib.pyplot as plt


def f(x):
    return math.sin(x/3 + math.exp(math.sin(x/3) ** 2))


h = 0.01  # шаг интерполяции
n = 25  # для узлов чебышева
a = 10
b = 0
upper_bound = int(10 / h) + 1

args = [i*h for i in range(upper_bound)]
values = list(map(f, args))


def polynomial(x, i, list_args):
    result = 1
    for j in range(len(list_args)):
        if i != j:
            result *= (x - list_args[j])/(list_args[i] - list_args[j])
    return result


def lagrange(x, list_args, list_values):
    result = 0
    for i in range(len(list_values)):
        result += list_values[i] * polynomial(x, i, list_args)
    return result


new_args = [(args[i] + args[i+1]) / 2 for i in range(len(args) - 1)]
new_values = [lagrange(x, args, values) for x in new_args]


list_diff = []
for i in range(len(args)):
    for j in range(len(args)):
        if i != j:
            list_diff.append(args[i] - args[j])

list_diff.sort()
print(min(map(abs, list_diff)))

print('args = ', args)
print('values = ', values)
print('new_args = ', new_args)
print('new_values = ', new_values)

for i in range(len(new_values)):
    if new_values[i] > 5:
        print(polynomial(new_args[i], i, args))


def cheb_node(k):
    return (a + b) / 2 + (b - a) / 2 * math.cos(math.pi * (2 * k + 1) / (2 * n))


cheb_args = [cheb_node(k) for k in range(n)]
cheb_values = list(map(f, cheb_args))
new_args1 = [(cheb_args[i] + cheb_args[i+1]) / 2 for i in range(len(cheb_args) - 1)]
new_values1 = [lagrange(x, cheb_args, cheb_values) for x in new_args1]


print('\n')

print('cheb_args = ', cheb_args)
print('cheb_values = ', cheb_values)
print('new_args1 = ', new_args1)
print('new_values1 = ', new_values1)

plt.subplot(2, 2, 1)
plt.plot(args, values, 'green', new_args, new_values, 'red')
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(cheb_args, cheb_values, 'green', new_args1, new_values1, 'red')
plt.grid(True)


# погрешность

error_values = list(map(lambda x, y: abs(x - y), values, new_values))
print('error_values = ', error_values)
plt.subplot(2, 2, 3)
args.pop()
plt.plot(args, error_values)
plt.grid(True)

error_cheb_values = list(map(lambda x, y: abs(x - y), cheb_values, new_values1))
print('error_cheb_values = ', error_cheb_values)
plt.subplot(2, 2, 4)
cheb_args.pop()
plt.plot(cheb_args, error_cheb_values)
plt.grid(True)

print('max_error = ', max(error_values))
print('max_cheb_error = ', max(error_cheb_values))

plt.show()



