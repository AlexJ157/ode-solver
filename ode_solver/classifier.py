from dataclasses import dataclass
import sympy as sp

@dataclass
class Classification:
    order: int
    method: str
    data: dict


def classify_first_order(rhs):
    is_lin, lin_data = is_linear(rhs)
    if is_lin:
        return Classification(order=1, method="linear", data=lin_data)

    is_sep, sep_data = is_seperable(rhs)   # renamed from is_seperable to is_sep
    if is_sep:
        return Classification(order=1, method="separable", data=sep_data)

    is_hom, hom_data = is_homogeneous(rhs)   # check this one too — same risk
    if is_hom:
        return Classification(order=1, method="homogeneous", data=hom_data)

    return Classification(order=1, method="numerical", data={})

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
        # TODO will fail for things like 1/y and sin(y)
        return False, None   

    degree = poly.degree()

    if degree > 1:
        return False, None

    coeffs = poly.all_coeffs()
    for c in coeffs:
        if c.has(y(x)):
            return False, None

    if degree == 1:
        P = -coeffs[0]
        Q = coeffs[1]
    else:
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
            return False, None

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
