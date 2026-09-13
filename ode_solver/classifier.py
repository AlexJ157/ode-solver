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


ode = "y' + 2*y = e^x"
parsed_ode = parser.parse_ode(ode)
rhs = normalize_first_order(parsed_ode)
print(rhs)