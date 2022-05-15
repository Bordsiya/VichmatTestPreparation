from DifferentialEquation import DifferentialEquation


def euler_method(differential_equation: DifferentialEquation):
    h = differential_equation.h
    f = differential_equation.f
    table = []
    headers = ["x", "y", "f(x, y)"]
    y_pred = differential_equation.y_0
    x_pred = differential_equation.x_0
    table.append([x_pred, y_pred, f(x_pred, y_pred)])

    x_curr = differential_equation.a
    while x_curr <= differential_equation.b:
        x_curr += h
        y_curr = y_pred + h * f(x_pred, y_pred)
        table.append([x_curr, y_curr, f(x_curr, y_curr)])
        y_pred = y_curr
        x_pred = x_curr

    return table, headers


def modified_euler_method(differential_equation: DifferentialEquation):
    f = differential_equation.f
    h = differential_equation.h
    table = []
    headers = ["x", "y", "f(x, y)", "y_exemplary", "f(x, y_exemplary"]
    y_pred = differential_equation.y_0
    x_pred = differential_equation.x_0
    table.append([x_pred, y_pred, f(x_pred, y_pred)])

    x_curr = differential_equation.a
    while x_curr <= differential_equation.b:
        x_curr += h
        y_exemplary = y_pred + h * f(x_pred, y_pred)
        y_curr = y_pred + (h / 2) * (f(x_pred, y_pred) + f(x_curr, y_exemplary))
        table.append([x_curr, y_curr, f(x_curr, y_curr), y_exemplary, f(x_curr, y_exemplary)])
        x_pred = x_curr
        y_pred = y_curr

    return table, headers


def runge_kutta_4_method(differential_equation: DifferentialEquation):
    f = differential_equation.f
    h = differential_equation.h
    table = []
    headers = ["x", "y", "k1", "k2", "k3", "k4"]
    y_pred = differential_equation.y_0
    x_pred = differential_equation.x_0

    x_curr = differential_equation.a
    while x_curr <= differential_equation.b:
        x_curr += h
        k1 = h * f(x_pred, y_pred)
        k2 = h * f(x_pred + (h / 2), y_pred + (k1 / 2))
        k3 = h * f(x_pred + (h / 2), y_pred + (k2 / 2))
        k4 = h * f(x_pred + h, y_pred + k3)
        y_curr = y_pred + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        table.append([x_pred, y_pred, k1, k2, k3, k4])
        y_pred = y_curr
        x_pred = x_curr
    table.append([x_pred, y_pred])

    return table, headers


def adams_method(differential_equation: DifferentialEquation):
    if ((differential_equation.b - differential_equation.a) / differential_equation.h) + 1 <= 4:
        return runge_kutta_4_method(differential_equation)

    table =[]
    headers = ["x", "y", "fi", "lambda_fi", "lambda2_fi", "lambda3_fi"]

    h = differential_equation.h
    f = differential_equation.f
    table_runge, headers_runge = runge_kutta_4_method(differential_equation)
    x_precalculated = []
    [x_precalculated.append(l[:1][0]) for l in table_runge[:4]]
    y_precalculated = []
    [y_precalculated.append(l[1:2][0]) for l in table_runge[:4]]
    for i in range(4):
        table.append([x_precalculated[i], y_precalculated[i]])

    x_curr = x_precalculated[3]
    while x_curr <= differential_equation.b:
        x_curr += h
        lambda_fi = f(x_precalculated[3], y_precalculated[3]) - f(x_precalculated[2], y_precalculated[2])
        lambda2_fi = f(x_precalculated[3], y_precalculated[3]) - 2 * f(x_precalculated[2], y_precalculated[2]) \
                    + f(x_precalculated[1], y_precalculated[1])
        lambda3_fi = f(x_precalculated[3], y_precalculated[3]) - 3 * f(x_precalculated[2], y_precalculated[2]) \
                    + 3 * f(x_precalculated[1], y_precalculated[1]) - f(x_precalculated[0], y_precalculated[0])
        y_curr = y_precalculated[3] + h * f(x_precalculated[3], y_precalculated[3]) \
                + ((h**2) / 2) * lambda_fi + ((5 * h**3) / 12) * lambda2_fi + ((3 * h**4) / 8) * lambda3_fi
        table.append([x_curr, y_curr, f(x_precalculated[3], y_precalculated[3]), lambda_fi, lambda2_fi, lambda3_fi])
        for i in range(1, 4):
            x_precalculated[i - 1] = x_precalculated[i]
            y_precalculated[i - 1] = y_precalculated[i]
        x_precalculated[3] = x_curr
        y_precalculated[3] = y_curr

    return table, headers


def predictor_corrector_method(differential_equation: DifferentialEquation):
    MAX_N_AMOUNT = 2*20

    if ((differential_equation.b - differential_equation.a) / differential_equation.h) + 1 <= 4:
        return runge_kutta_4_method(differential_equation)

    table = []
    headers = ["x", "y"]

    h = differential_equation.h
    f = differential_equation.f
    table_runge, headers_runge = runge_kutta_4_method(differential_equation)
    x_precalculated = []
    [x_precalculated.append(l[:1][0]) for l in table_runge[:4]]
    y_precalculated = []
    [y_precalculated.append(l[1:2][0]) for l in table_runge[:4]]
    for i in range(4):
        table.append([x_precalculated[i], y_precalculated[i]])

    x_curr = x_precalculated[3]
    while x_curr <= differential_equation.b:
        x_curr += h
        fi = f(x_precalculated[3], y_precalculated[3])
        fi_1 = f(x_precalculated[2], y_precalculated[2])
        fi_2 = f(x_precalculated[1], y_precalculated[1])
        fi_3 = f(x_precalculated[0], y_precalculated[0])
        y_predicted = y_precalculated[3] + (h / 24) * (55 * fi - 59 * fi_1 + 37 * fi_2 - 9 * fi_3)

        y_pred = y_predicted
        fi_corrector = f(x_curr, y_predicted)
        y_curr = y_precalculated[3] + (h / 24) * (9 * fi_corrector - 19 * fi - 5 * fi_1 + fi_2)

        iteration_amount = 1
        while abs(y_curr - y_pred) > differential_equation.e:
            y_pred = y_curr
            fi_corrector = f(x_curr, y_pred)
            y_curr = y_precalculated[3] + (h / 24) * (9 * fi_corrector - 19 * fi - 5 * fi_1 + fi_2)
            iteration_amount += 1
            if iteration_amount > MAX_N_AMOUNT:
                return None
        table.append([x_curr, y_curr])
        for i in range(1, 4):
            x_precalculated[i - 1] = x_precalculated[i]
            y_precalculated[i - 1] = y_precalculated[i]
        x_precalculated[3] = x_curr
        y_precalculated[3] = y_curr

    return table, headers


def miln_method(differential_equation: DifferentialEquation):
    MAX_N_AMOUNT = 2 * 20

    if ((differential_equation.b - differential_equation.a) / differential_equation.h) + 1 <= 4:
        return runge_kutta_4_method(differential_equation)

    table = []
    headers = ["x", "y"]

    h = differential_equation.h
    f = differential_equation.f
    table_runge, headers_runge = runge_kutta_4_method(differential_equation)
    x_precalculated = []
    [x_precalculated.append(l[:1][0]) for l in table_runge[:4]]
    y_precalculated = []
    [y_precalculated.append(l[1:2][0]) for l in table_runge[:4]]
    for i in range(4):
        table.append([x_precalculated[i], y_precalculated[i]])

    x_curr = x_precalculated[3]
    while x_curr <= differential_equation.b:
        x_curr += h
        fi = f(x_precalculated[3], y_precalculated[3])
        fi_1 = f(x_precalculated[2], y_precalculated[2])
        fi_2 = f(x_precalculated[1], y_precalculated[1])
        y_predicted = y_precalculated[0] + (4 * h / 3) * (2 * fi_2 - fi_1 + 2 * fi)

        y_pred = y_predicted
        fi_corrector = f(x_curr, y_predicted)
        y_curr = y_precalculated[2] + (h / 3) * (fi_1 + 4 * fi + fi_corrector)

        iteration_amount = 1
        while abs(y_curr - y_pred) > differential_equation.e:
            y_pred = y_curr
            fi_corrector = f(x_curr, y_pred)
            y_curr = y_precalculated[2] + (h / 3) * (fi_1 + 4 * fi + fi_corrector)
            iteration_amount += 1
            if iteration_amount > MAX_N_AMOUNT:
                return None
        table.append([x_curr, y_curr])
        for i in range(1, 4):
            x_precalculated[i - 1] = x_precalculated[i]
            y_precalculated[i - 1] = y_precalculated[i]
        x_precalculated[3] = x_curr
        y_precalculated[3] = y_curr

    return table, headers
