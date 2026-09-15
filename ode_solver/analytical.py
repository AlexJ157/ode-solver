import sympy as sp
import classifier

x = sp.symbols('x')
y = sp.Function('y')

def intergrating_factor(rhs):
    is_linear, data = classifier.is_linear(rhs)
 
    if not is_linear:
        print("ODE cant be solved using the intergrating factor.")
        return None, []
 
    P, Q = data["p"], data["Q"]
 
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
 
    solution = sp.Eq(solution_lhs, solution_rhs)
    return solution, steps

