import sympy as sp

# example input y' + 2y = e^x

def parse_ode(equation):
    equation = equation.replace(" ", "")
    print(equation)

    left, right = split_equation(equation)

    parsed_left = 0
    for expression in left.split("+"):
        print(expression)
        if "'" in expression:
            parsed_left += parse_derivative(expression)
        else:
            parsed_left += parse_expression(expression)

    parsed_right = 0
    for expression in right.split("+"):
        if "'" in expression:
            parsed_right += parse_derivative(expression)
        else:
            parsed_right += parse_expression(expression)

    return sp.Eq(parsed_left, parsed_right)


def split_equation(equation):
    split_equation = equation.split("=")
    print(split_equation)
    return (split_equation[0], split_equation[1])

def parse_expression(expression):
    x = sp.symbols('x')
    y = sp.Function('y')
    expression = expression.replace("^", "**")
    return sp.sympify(expression)

def parse_derivative(derivative):
    x = sp.symbols('x')
    y = sp.Function('y')

    order = derivative.count("'")
    parse_derivative = sp.diff(y(x), x, order)
    return parse_derivative


ode = "y' + 2*y = e^x"
parsed_ode = parse_ode(ode)
print(parsed_ode)