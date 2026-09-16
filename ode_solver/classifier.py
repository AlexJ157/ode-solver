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
        return True, {"f(x)": factored, "g(y)": sp.Integer(1)}

    f_x = sp.Mul(*f_x_parts)
    g_y = sp.Mul(*g_y_parts)

    return True, {"f(x)": f_x, "g(y)": g_y}

def is_homogeneous(rhs):  # TODO lots of cases will break parser as it doesnt rcognise fractions with brackets
    z = sp.Symbol('z')
    substituted = rhs.subs(y(x), z*x)
    simplified = sp.simplify(substituted)

    if simplified.has(x):
        return False, None

    return True, {"G_z": simplified}
