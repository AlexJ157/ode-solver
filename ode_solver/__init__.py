from . import parser, classifier, analytical

def solve_ode(ode_string):
    equation = parser.parse_ode(ode_string)
    rhs = classifier.normalize_first_order(equation)
    classification = classifier.classify_first_order(rhs)

    if classification.method == "linear":
        return analytical.intergrating_factor(rhs)
    
    elif classification.method == "separable":
        return analytical.seperable(rhs)
    
    elif classification.method == "homogeneous":
        return analytical.homogeneous(rhs)
    
    else:
        print("No analytical method found — numerical methods not yet implemented")
        return None, []
