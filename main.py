from ode_solver import solve_ode

test_odes = [
    "y' + 2*y = e^x",
    "y' - 3*y = 0",
    "y' + y = x**2",
    "y' = x**2",
    "y' = y**2",
    "y' = x*y",
    "y' = e^x/y",
    "y' = x*y**2 + x",
]

for ode_str in test_odes:
    solution, steps = solve_ode(ode_str)
    print(f"{ode_str}  ->  solution = {solution}")