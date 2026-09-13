import sympy as sp
import parser

def normalize_first_order(eq):
    x = sp.symbols('x')
    y = sp.Function('y')
    dydx = sp.Derivative(y(x), x)

    solved = sp.solve(eq, dydx)
    if len(solved) == 0:
        raise ValueError("Could not isolate dy/dx")
    if len(solved) > 1:
        raise ValueError("Equation is nonlinear in dy/dx - not solvable by these methods")
    return solved[0]

def is_linear(rhs):
    x = sp.symbols('x')
    y = sp.Function('y')

    try:
        poly = sp.Poly(rhs, y(x))
    except sp.PolynomialError:
        # will fail for things like 1/y and sin(y)
        return False, None   

    degree = poly.degree()

    if degree > 1:
        return False, None

    coeffs = poly.all_coeffs()
    for c in coeffs:
        if c.has(y(x)):
            return False, None

    if degree == 1: # degree = 1 - solve with intergrating factor
        P = -coeffs[0]
        Q = coeffs[1]
    else: # degree = 0 - normal intergral
        P = 0
        Q = coeffs[0]

    return True, {"P": P, "Q": Q}



test_odes = [
    "y' + 2*y = e^x",
    "y' - 3*y = 0",
    "y' + y = x**2",
    "y' = x**2",
    "y' = y**2",
]
 
for ode_str in test_odes:
    eq = parser.parse_ode(ode_str)
    rhs = normalize_first_order(eq)
    result = is_linear(rhs)
    print(f"{ode_str}  ->  rhs = {rhs}  ->  is_linear = {result}")
    print()
