import sympy as sp
import parser

x = sp.symbols('x')
y = sp.Function('y')

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

def is_seperable(rhs):
    combined = sp.together(rhs)
    factored = sp.factor(combined)
    factors = factored.as_ordered_factors()

    f_x_parts = []
    g_y_parts = []

    for factor in factors:
        has_y = factor.has(y(x))

        if not has_y:
            f_x_parts.append(factor)
            continue

        placeholder = sp.Dummy('Y')
        swapped = factor.subs(y(x), placeholder)

        if swapped.has(x):
            return False, None   # genuinely mixed — x remains even without y(x)

        g_y_parts.append(factor)

    if len(g_y_parts) == 0:
        return True, {"f(x)": factored, "g(x)": 1}

    f_x = sp.Mul(*f_x_parts)
    g_y = sp.Mul(*g_y_parts)

    return True, {"f(x)": f_x, "g(x)": g_y}

linear_test_odes = [
    "y' + 2*y = e^x",
    "y' - 3*y = 0",
    "y' + y = x**2",
    "y' = x**2",
    "y' = y**2",
]

seperable_test_odes = [
    "y' = x*y",
    "y' = x**2",
    "y' = y**2",
    "y' = e^x/y",
    "y' = x*y**2 + x",
    "y' = x+y",
    "y' = x**2 + y**2",
]
 
for ode_str in seperable_test_odes:
    eq = parser.parse_ode(ode_str)
    rhs = normalize_first_order(eq)
    result = is_seperable(rhs)
    print(f"{ode_str}  ->  rhs = {rhs}  ->  seperable = {result}")
    print()
