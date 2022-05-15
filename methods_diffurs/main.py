from tabulate import tabulate

from DifferentialEquation import DifferentialEquation
from methods import euler_method, modified_euler_method, runge_kutta_4_method, adams_method, predictor_corrector_method, \
    miln_method

f = lambda x, y: y + (1 + x) * y**2
diff_eq = DifferentialEquation(f, 1, -1, 1, 1.5, 0.1, 0.001)
print("-------------------------------")
table, headers = euler_method(diff_eq)
print("Метод Эйлера")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")
table, headers = modified_euler_method(diff_eq)
print("Модифицированный метод Эйлера")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")
table, headers = runge_kutta_4_method(diff_eq)
print("Метод Рунге-Кутта четвертого порядка")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")
table, headers = adams_method(diff_eq)
print("Метод Адамса")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")
table, headers = predictor_corrector_method(diff_eq)
print("Метод прогноза и коррекции")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")
table, headers = miln_method(diff_eq)
print("Метод Милна")
print(tabulate(table, headers, tablefmt="github"))
print("-------------------------------")