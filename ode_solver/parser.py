import re
import sympy as sp
 
# example input y' + 2y = e^x
 
def parse_ode(equation):
    equation = equation.replace(" ", "")
    print(equation)
 
    left, right = split_equation(equation)
 
    parsed_left = 0
    for expression in split_terms(left):
        print(expression)
        if "'" in expression:
            parsed_left += parse_derivative(expression)
        else:
            parsed_left += parse_expression(expression)
 
    parsed_right = 0
    for expression in split_terms(right):
        if "'" in expression:
            parsed_right += parse_derivative(expression)
        else:
            parsed_right += parse_expression(expression)
 
    return sp.Eq(parsed_left, parsed_right)
 
 
def split_equation(equation):
    split_equation = equation.split("=")
    print(split_equation)
    return (split_equation[0], split_equation[1])


def split_terms(expression):
    return re.findall(r'[+-]?[^+-]+', expression)

def parse_expression(expression):
    expression = expression.replace("^", "**")
    return parse_coefficient(expression)


def parse_coefficient(text):
    if text in ("", "+"):
        return 1
    if text == "-":
        return -1
    text = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', text)
    return sp.sympify(text)
 
 
def strip_derivative(term):
    order = term.count("'")
    coefficient_part = term.replace("'", "")
    coefficient_part = coefficient_part.replace("y", "")
    coefficient_part = coefficient_part.replace("*", "")
    return coefficient_part, order
 

def parse_derivative(derivative):
    x = sp.symbols('x')
    y = sp.Function('y')
 
    coefficient_part, order = strip_derivative(derivative)
    deriv_expr = sp.diff(y(x), x, order)
    coefficient = parse_coefficient(coefficient_part)
 
    return coefficient * deriv_expr


ode = "y' + 2*y = e^x"
parsed_ode = parse_ode(ode)
print(parsed_ode)