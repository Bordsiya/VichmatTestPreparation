from Math import kramer_method
from TableFunction import TableFunction


def linear_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        return None
    #определяем в какой интервал попадает x
    y_left, y_right, x_left, x_right = None, None, None, None
    for i in range(0, table_function.n - 1):
        if (x > table_function.x_arr[i]) and (x < table_function.x_arr[i + 1]):
            x_left = table_function.x_arr[i]
            x_right = table_function.x_arr[i + 1]
            y_left = table_function.y_arr[i]
            y_right = table_function.y_arr[i + 1]
            break
        elif x == table_function.x_arr[i]:
            return table_function.y_arr[i]
        elif x == table_function.x_arr[i + 1]:
            return table_function.y_arr[i + 1]

    if x_left is None:
        return None
    a = (y_right - y_left) / (x_right - x_left)
    b = y_left - a * x_left
    return a * x + b


def square_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        return None
    #находим три ближайших к x узла
    #сначала находим интервал в котором лежит x
    x_left_ind, x_right_ind = None, None
    for i in range(0, table_function.n - 1):
        if (x > table_function.x_arr[i]) and (x < table_function.x_arr[i + 1]):
            x_left_ind = i
            x_right_ind = i + 1
            break
        elif x == table_function.x_arr[i]:
            return table_function.y_arr[i]
        elif x == table_function.x_arr[i + 1]:
            return table_function.y_arr[i + 1]

    if x_left_ind is None:
        return None

    x_middle = (table_function.x_arr[x_left_ind] + table_function.x_arr[x_right_ind]) / 2
    if (x <= x_middle) and x_left_ind - 1 >= 0:
        #третья точка слева от интервала
        x0_ind = x_left_ind - 1
        x1_ind = x_left_ind
        x2_ind = x_right_ind
    elif x_right_ind + 1 < table_function.n:
        #третья точка справа от интервала
        x0_ind = x_left_ind
        x1_ind = x_right_ind
        x2_ind = x_right_ind + 1
    else:
        #не смогли найти три точки
        return None

    #составляем матрицу для вычисления коэффициентов a, b, c
    matrix_coefficients = [[table_function.x_arr[x0_ind]**2, table_function.x_arr[x0_ind], 1],
                           [table_function.x_arr[x1_ind]**2, table_function.x_arr[x1_ind], 1],
                           [table_function.x_arr[x2_ind]**2, table_function.x_arr[x2_ind], 1]]
    matrix_answers = [table_function.y_arr[x0_ind],
                      table_function.y_arr[x1_ind],
                      table_function.y_arr[x2_ind]]

    coeff_arr = kramer_method(matrix_coefficients, matrix_answers)
    if coeff_arr is None:
        return None
    a = coeff_arr[0]
    b = coeff_arr[1]
    c = coeff_arr[2]
    return a * (x**2) + b * x + c


def langrange_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        return None
    L_arr = []
    for i in range(table_function.n):
        comp_1 = 1
        comp_2 = 1
        for j in range(table_function.n):
            if i != j:
                comp_1 *= (x - table_function.x_arr[j])
                comp_2 *= (table_function.x_arr[i] - table_function.x_arr[j])
        l = comp_1 / comp_2
        L_arr.append(l * table_function.y_arr[i])

    return sum(L_arr)


def get_f(ind: list[int], table_function: TableFunction):
    if len(ind) == 1:
        return table_function.y_arr[ind[0]]
    if len(ind) == 2:
        return (table_function.y_arr[ind[1]] - table_function.y_arr[ind[0]])\
               / (table_function.x_arr[ind[1]] - table_function.x_arr[ind[0]])
    return (get_f(ind[1:], table_function) - get_f(ind[:len(ind) - 1], table_function))\
           / (table_function.x_arr[ind[len(ind) - 1]] - table_function.x_arr[ind[0]])


def get_factorial(n: int):
    if n == 0 or n == 1:
        return 1
    return n * get_factorial(n - 1)


def get_lambda(k: int, ind: int, table_function: TableFunction):
    if k == 0:
        #print("//-",k,  table_function.y_arr[ind])
        return table_function.y_arr[ind]
    if k == 1:
        #print("//+",k,  table_function.y_arr[ind + 1],  table_function.y_arr[ind])
        return table_function.y_arr[ind + 1] - table_function.y_arr[ind]

    #print("/////", k, get_lambda(k - 1, ind + 1, table_function), get_lambda(k - 1, ind, table_function))
    return get_lambda(k - 1, ind + 1, table_function) - get_lambda(k - 1, ind, table_function)


def newton_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        return None
    is_equally_spaced = True
    h = round(table_function.x_arr[1] - table_function.x_arr[0], 3)
    #проверяем, равноотстоящие ли узлы
    for i in range(1, table_function.n - 1):
        if round(table_function.x_arr[i + 1] - table_function.x_arr[i], 3) != h:
            is_equally_spaced = False
            break
    if not is_equally_spaced:
        #неравноотстоящие узлы
        #ищем индекс x_i >= x
        x_right_ind = None
        for i in range(table_function.n):
            if (x > table_function.x_arr[i]) and (x < table_function.x_arr[i + 1]):
                x_right_ind = i + 1
                break
            elif x == table_function.x_arr[i]:
                return table_function.y_arr[i]
            elif x == table_function.x_arr[i + 1]:
                return table_function.y_arr[i + 1]

        if x_right_ind is None:
            return None

        n = 0
        args = []
        for i in range(x_right_ind + 1):
            if i == 0:
                n += get_f([i], table_function)
                args.append(0)
            else:
                comp = 1
                for j in range(i):
                    comp *= (x - table_function.x_arr[j])
                args.append(i)
                comp *= get_f(args, table_function)
                n += comp
        return n
    else:
        #равноотстоящие узлы
        x_middle = (table_function.x_arr[0] + table_function.x_arr[table_function.n - 1]) / 2
        if x > x_middle:
            #интерполирование назад
            t = (x - table_function.x_arr[table_function.n - 1]) / h
            n = get_lambda(0, table_function.n - 1, table_function)
            for i in range(1, table_function.n):
                comp_t = 1
                curr_t = t
                for j in range(i):
                    comp_t *= curr_t
                    curr_t += 1
                n += comp_t * get_lambda(i, table_function.n - i - 1, table_function) / get_factorial(i)
            return n
        else:
            #интерполирование вперед
            #ищем левую границу интервала для x
            x_left_ind = None
            for i in range(0, table_function.n - 1):
                if (x > table_function.x_arr[i]) and (x < table_function.x_arr[i + 1]):
                    x_left_ind = i
                    break
                elif x == table_function.x_arr[i]:
                    return table_function.y_arr[i]
                elif x == table_function.x_arr[i + 1]:
                    return table_function.y_arr[i + 1]
            if x_left_ind is None:
                return None

            t = (x - table_function.x_arr[x_left_ind]) / h
            n = get_lambda(0, x_left_ind, table_function)
            for i in range(1, table_function.n - x_left_ind):
                comp_t = 1
                curr_t = t
                for j in range(i):
                    comp_t *= curr_t
                    curr_t -= 1
                n += comp_t * get_lambda(i, x_left_ind, table_function) / get_factorial(i)
            return n


def gauss_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        return None

    is_equally_spaced = True
    h = round(table_function.x_arr[1] - table_function.x_arr[0], 3)
    # проверяем, равноотстоящие ли узлы
    for i in range(1, table_function.n - 1):
        if round(table_function.x_arr[i + 1] - table_function.x_arr[i], 3) != h:
            is_equally_spaced = False
            break
    if not is_equally_spaced:
        return None

    if table_function.n % 2 == 0:
        x_middle_ind = int(table_function.n / 2 - 1)
    else:
        x_middle_ind = int(table_function.n / 2)
    t = (x - table_function.x_arr[x_middle_ind]) / h
    if x > table_function.x_arr[x_middle_ind]:
        #первая формула гаусса
        n = get_lambda(0, x_middle_ind, table_function)
        for i in range(1, table_function.n):
            comp_t = t
            pred_sign ='+'
            pred_number = 0
            for j in range(0, i - 1, 1):
                if pred_sign == '-':
                    comp_t *= (t + pred_number)
                    pred_sign = '+'
                else:
                    pred_number += 1
                    comp_t *= (t - pred_number)
                    pred_sign = '-'
            arg = int((table_function.n - i) / 2)
            n += (comp_t / get_factorial(i)) * get_lambda(i, arg, table_function)
        return n

    else:
        #вторая формула гаусса
        n = get_lambda(0, x_middle_ind, table_function)

        for i in range(1, table_function.n):
            comp_t = t
            pred_sign = '-'
            pred_number = 0
            for j in range(0, i - 1, 1):
                if pred_sign == '-':
                    pred_number += 1
                    comp_t *= (t + pred_number)
                    pred_sign = '+'
                else:
                    comp_t *= (t - pred_number)
                    pred_sign = '-'
            if (table_function.n - i) % 2 == 0:
                arg = int((table_function.n - i) / 2) - 1
            else:
                arg = int((table_function.n - i) / 2)
            n += (comp_t / get_factorial(i)) * get_lambda(i, arg, table_function)
        return n


def stirling_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        print("количество узлов слишком мало")
        return None

    if table_function.n % 2 == 0:
        print("нужно нечетное число узлов")
        return None

    is_equally_spaced = True
    h = round(table_function.x_arr[1] - table_function.x_arr[0], 3)
    # проверяем, равноотстоящие ли узлы
    for i in range(1, table_function.n - 1):
        if round(table_function.x_arr[i + 1] - table_function.x_arr[i], 3) != h:
            is_equally_spaced = False
            break
    if not is_equally_spaced:
        print("узлы неравноотстоящие")
        return None

    # ищем узел x0, максимально близкий к x
    if table_function.n % 2 == 0:
        x_0 = int(table_function.n / 2 - 1)
    else:
        x_0 = int(table_function.n / 2)

    t = (x - table_function.x_arr[x_0]) / h
    if abs(t) > 0.25:
        print("значение t не удовлетворяет условиям применимости, будет большая погрешность")
        return None

    n = get_lambda(0, x_0, table_function)
    comp_t1 = t
    comp_t2 = t**2
    pred_number = 0
    print(n)
    for i in range(1, table_function.n):
        if i % 2 == 0:
            print("//", get_lambda(i, x_0 - (i // 2), table_function))
            print(i, (comp_t2 / get_factorial(i)) * get_lambda(i, x_0 - (i // 2), table_function))
            n += (comp_t2 / get_factorial(i)) * get_lambda(i, x_0 - (i // 2), table_function)
            comp_t2 *= (t**2 - pred_number**2)
        else:
            print("++", get_lambda(i, x_0 - ((i + 1) // 2), table_function), get_lambda(i, x_0 - (((i + 1) // 2) - 1),
                                                                               table_function))
            print(i, (comp_t1 / get_factorial(i)) * \
                 ((get_lambda(i, x_0 - ((i + 1) // 2), table_function) + get_lambda(i, x_0 - (((i + 1) // 2) - 1),
                                                                               table_function)) / 2))
            n += (comp_t1 / get_factorial(i)) * \
                 ((get_lambda(i, x_0 - ((i + 1) // 2), table_function) + get_lambda(i, x_0 - (((i + 1) // 2) - 1),
                                                                               table_function)) / 2)
            pred_number += 1
            comp_t1 *= (t**2 - pred_number**2)
        #comp_t = t
        #pred_number = 0
        #for j in range(i - 1):
        #    pred_number -= 1
        #    comp_t *= (t**2 - pred_number**2)
        #n += (comp_t / get_factorial(2 * i - 1)) * \
        #         ((get_lambda(2 * i - 1, x_0 - i, table_function) + get_lambda(2 * i - 1, x_0 - i + 1, table_function)) / 2)
        #n += (comp_t * t / get_factorial(2 * i)) * get_lambda(2 * i, x_0 - i, table_function)
    return n


def bessel_interpolation(table_function: TableFunction, x: float):
    if table_function.n < 2:
        print("количество узлов слишком мало")
        return None

    if table_function.n % 2 != 0:
        print("нужно четное число узлов")
        return None

    is_equally_spaced = True
    h = round(table_function.x_arr[1] - table_function.x_arr[0], 3)
    # проверяем, равноотстоящие ли узлы
    for i in range(1, table_function.n - 1):
        if round(table_function.x_arr[i + 1] - table_function.x_arr[i], 3) != h:
            is_equally_spaced = False
            break
    if not is_equally_spaced:
        print("неравноотстоящие узлы")
        return None

    #ищем узел x0, максимально близкий к x
    if table_function.n % 2 == 0:
        x_0 = int(table_function.n / 2 - 1)
    else:
        x_0 = int(table_function.n / 2)

    t = (x - table_function.x_arr[x_0]) / h

    if not 0.25 <= abs(t) <= 0.75:
        print("t не подходит под метод")

    n = (get_lambda(0, x_0, table_function) + get_lambda(0, x_0 + 1, table_function)) / 2
    n += (t - 0.5) * get_lambda(1, x_0, table_function)
    comp_t = t
    last_number = 0
    for i in range(2, table_function.n):
        if i % 2 == 0:
            last_number += 1
            comp_t *= (t - last_number)
            n += (comp_t / get_factorial(i)) * ((get_lambda(i, x_0 - i // 2, table_function)
                                                 + get_lambda(i, x_0 - ((i//2) - 1), table_function)) / 2)
        else:
            n += (comp_t * (t - 0.5) / get_factorial(i)) * get_lambda(i, x_0 - ((i - 1)//2), table_function)
            comp_t *= (t + last_number)
    return n
