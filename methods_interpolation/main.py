from IntervalMethods import linear_interpolation, square_interpolation, langrange_interpolation, newton_interpolation, \
    gauss_interpolation, stirling_interpolation, bessel_interpolation
from TableFunction import TableFunction

table_function = TableFunction(5, [0.1, 0.2, 0.3, 0.4, 0.5], [1.25, 2.38, 3.79, 5.44, 7.14])
table_function_lagrange = TableFunction(3, [100, 121, 144], [10, 11, 12])
table_function_linear_square = TableFunction(5, [0.1, 0.2, 0.3, 0.4, 0.5], [1.25, 2.38, 3.79, 5.44, 7.14])
table_bessel = TableFunction(6, [0.6, 0.7, 0.8, 0.9, 1, 1.1], [1.8221, 2.0138, 2.2255, 2.4596, 2.7183, 3.0042])
print("Линейная интерполяция: ", linear_interpolation(table_function_linear_square, 0.35))
print("Квадратичная интерполяция: ", square_interpolation(table_function_linear_square, 0.35))
print("Лагранж: ", langrange_interpolation(table_function_lagrange, 105))
print("Ньютон: ", newton_interpolation(table_function, 0.47))
print("Гаусс: ", gauss_interpolation(table_function, 0.47))
print("Стирлинг: ", stirling_interpolation(table_bessel, 0.98))
print("Бессель: ", bessel_interpolation(table_bessel, 0.95))
