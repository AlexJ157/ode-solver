import sympy as sp
from . import classifier

x = sp.symbols('x')
y = sp.Function('y')
z = sp.symbols('z')


def intergrating_factor(rhs):
    is_linear, data = classifier.is_linear(rhs)
    print("ATTEMPTING TO SOLVE USING INTERGRATING FACTOR:\n")
 
    if not is_linear:
        print("ODE cant be solved using the intergrating factor.")
        return None, []
 
    P, Q = data["P"], data["Q"]
 
    steps = []
 
    ode = sp.Eq(sp.diff(y(x), x, 1) + P * y(x), Q)
    steps.append(("Standard form", f"y'(x) + {P}y(x) = {Q}"))
    steps.append(("Identified P, Q", f"P(x) = {P}, Q(x) = {Q}"))
 
    integral = sp.integrate(P, x)
    I_factor = sp.exp(integral)
    steps.append(("Integrating factor I(x)", I_factor))
 
    new_lhs = ode.lhs * I_factor
    new_rhs = ode.rhs * I_factor
    new_ode = sp.Eq(new_lhs, new_rhs)
    steps.append(("Multiply through by I(x)", f"{new_lhs} = {new_rhs}"))
 
    steps.append(("Collapses to", f"d/dx({I_factor * y(x)}) = {new_rhs}"))
 
    solution_lhs = I_factor * y(x)
    solution_rhs = sp.integrate(new_ode.rhs, x)
    steps.append(("Integrate both sides", f"{solution_lhs} = {solution_rhs} + C"))
 
    C = sp.symbols('C')
    solution = sp.Eq(solution_lhs, solution_rhs + C)
    return solution, steps

def seperable(rhs):
    is_seperable, data = classifier.is_seperable(rhs)
    print("ATTEMPTING TO SOLVE BY SEPERABLE METHOD:\n")

    if not is_seperable:
        print("ODE isn't seperable and therefore can't be solved this way.")
        return None, []

    f_x, g_y = data["f(x)"], data["g(y)"]
    steps = []

    steps.append(("Standard form", f"dy/dx = ({f_x})({g_y})"))
    steps.append(("Identified f(x), g(y)", f"f(x) = {f_x}, g(y) = {g_y}"))
    steps.append(("Separated variables", f"dy/({g_y}) = ({f_x})dx"))

    Y = sp.symbols('Y')
    g_y_sub = g_y.subs(y(x), Y)
    solution_lhs = sp.integrate(1/g_y_sub, Y)
    solution_lhs = solution_lhs.subs(Y, y(x))

    solution_rhs = sp.integrate(f_x, x)

    steps.append(("Integrate LHS (w.r.t. y)", solution_lhs))
    steps.append(("Integrate RHS (w.r.t. x)", solution_rhs))

    C = sp.symbols('C')
    implicit_solution = sp.Eq(solution_lhs, solution_rhs + C)
    steps.append(("Implicit general solution", f"{solution_lhs} = {solution_rhs} + C"))

    try:
        explicit = sp.solve(implicit_solution, y(x))
    except NotImplementedError:
        explicit = []

    if explicit:
        steps.append(("Solved for y", explicit))
        return explicit, steps

    return implicit_solution, steps

def homogeneous(rhs):
    is_homogeneous, data = classifier.is_homogeneous(rhs)
    print("ATTEMPTING TO SOLVE HOMOGENEOUS ODE:\n")
 
    if not is_homogeneous:
        print("ODE isn't homogeneous and therefore can't be solved this way.")
        return None, []
 
    g_z = data["G_z"]
    steps = []
 
    steps.append(("Standard form", f"dy/dx = {rhs}"))
    steps.append(("Substitution", "Let z = y/x, so y = z*x"))
    steps.append(("Differentiate y=zx", "dy/dx = x*dz/dx + z"))
    steps.append(("Substituted RHS", f"G(z) = {g_z}"))
    steps.append(("Substituted equation", f"x*dz/dx + z = {g_z}"))
 
    rearranged = (g_z - z) / x
    steps.append(("Rearranged for dz/dx", f"dz/dx = {rearranged}"))

    if sp.simplify(g_z - z) == 0:
        steps.append(("Special case: dz/dx = 0", "z is constant, so z = C"))
        C = sp.symbols('C')
        solution = sp.Eq(y(x), C * x)
        steps.append(("Solved for y", f"y(x) = {C}*x"))
        return solution, steps
    
    steps.append(("Separated variables", f"dz/({g_z - z}) = dx/x"))

    solution_lhs = sp.integrate(1/(g_z - z), z)
    solution_rhs = sp.integrate(1/x, x)
    steps.append(("Integrate LHS (w.r.t. z)", solution_lhs))
    steps.append(("Integrate RHS (w.r.t. x)", solution_rhs))

    C = sp.symbols('C')
    z_implicit_solution = sp.Eq(solution_lhs, solution_rhs + C)
    steps.append(("Implicit general solution", f"{solution_lhs} = {solution_rhs} + C"))

    implicit_solution = z_implicit_solution.subs(z, y(x)/x)
    steps.append(("Substitute back z = y/x", f"{implicit_solution.lhs} = {implicit_solution.rhs}"))

    try:
        explicit = sp.solve(implicit_solution, y(x))
    except NotImplementedError:
        explicit = []

    if explicit:
        steps.append(("Solved for y", explicit))
        return explicit, steps

    return implicit_solution, steps