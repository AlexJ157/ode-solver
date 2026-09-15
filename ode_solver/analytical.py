import sympy as sp
import classifier

x = sp.symbols('x')
y = sp.Function('y')

def intergrating_factor(rhs):
    is_linear, data = classifier.is_linear(rhs)

    if not is_linear:
        print("ODE cant be solved using the intergrating factor")
        return

    P = data.get("P")
    Q = data.get("Q")

    ode = sp.Eq((sp.diff(y(x), x, 1) + P*y(x)), Q)

    print(f"ODE: y'(x) + {P}y(x) = {Q}\n")
    print(f"P(x) = {P}, Q(x) = {Q}\n")

    integral = sp.integrate(P, x)
    I_factor = sp.exp(integral)
    new_lhs = ode.lhs * I_factor
    new_rhs = ode.rhs * I_factor
    new_ode = sp.Eq(new_lhs, new_rhs)
    print(f"Intergrating factor - I(x): {I_factor}\n")
    print(f"New ode:  y'(x){I_factor} + {P*I_factor}y(x) = {Q*I_factor}\n")
    print(f"=> d/dx({I_factor*y(x)}) = {Q*I_factor}\n")

    solution_lhs = I_factor*y(x)
    solution_rhs = sp.integrate(new_ode.rhs, x)
    solution = sp.Eq(solution_lhs, solution_rhs)
    print(f"Solution: {solution_lhs} = {solution_rhs} + C")

    return solution



