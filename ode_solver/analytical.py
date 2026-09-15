import sympy as sp
import classifier
import parser

x = sp.symbols('x')
y = sp.Function('y')

def intergrating_factor(rhs):
    is_linear, data = classifier.is_linear(rhs)
 
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
 
for ode_str in linear_test_odes:
    eq = parser.parse_ode(ode_str)
    rhs = classifier.normalize_first_order(eq)
    solution, steps = intergrating_factor(rhs)
    print(f"{ode_str}  ->  solution = {solution}  ->  steps = {steps}")
    print()

for ode_str in seperable_test_odes:
    eq = parser.parse_ode(ode_str)
    rhs = classifier.normalize_first_order(eq)
    solution, steps = seperable(rhs)
    print(f"{ode_str}  ->  solution = {solution}  ->  steps = {steps}")
    print()