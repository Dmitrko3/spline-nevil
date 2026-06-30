import sympy as sp

def create_mash(left_boundary, right_boundary, n):
    """
    Creates a mesh of n intervals (n + 1 points) between left and right boundaries.
    Uses eq: h = (b - a) / n
    """
    h = (right_boundary - left_boundary) / n
    # Utilizing float to ensure numerical operations
    return [float(left_boundary + i * h) for i in range(n + 1)]

def calculate_integral(mash, f):
    """
    Calculates the integral over the mesh using Simpson's 1/3 rule (Eq 1.1):
    I = (h / 3) * [f(x_0) + 4*f(x_1) + 2*f(x_2) + 4*f(x_3) + ... + f(x_n)]
    """
    n = len(mash) - 1
    h = (mash[-1] - mash[0]) / n
    
    # Initialize with the boundary terms
    integral_sum = f(mash[0]) + f(mash[-1])
    
    # Add the inner terms with alternating weights of 4 and 2
    for i in range(1, n):
        if i % 2 == 1:
            integral_sum += 4 * f(mash[i])
        else:
            integral_sum += 2 * f(mash[i])
            
    return (h / 3) * integral_sum

def sympson_method(left_boundary, right_boundary, polynomial, epsilon=1e-6):
    """
    Approximates the integral using Simpson's method. Refines the mesh 
    iteratively by adding 10 intervals each step until the error is below epsilon.
    """
    # Convert SymPy polynomial/expression into a callable function
    free_symbols = list(polynomial.free_symbols)
    if len(free_symbols) >= 1:
        var = free_symbols[0]
        f = sp.lambdify(var, polynomial, 'numpy')
    else:
        # Fallback for constant expressions
        f = lambda x: float(polynomial)

    n, num_of_iterations = 10, 0
    
    # Initial estimation
    mash = create_mash(left_boundary, right_boundary, n)
    old_integral = calculate_integral(mash, f)
    
    # Second estimation with n + 10 intervals
    n += 10
    mash = create_mash(left_boundary, right_boundary, n)
    new_integral = calculate_integral(mash, f)
    
    # Iterative refinement loop based on tolerance (epsilon)
    while abs(new_integral - old_integral) > epsilon:
        num_of_iterations += 1
        n += 10
        old_integral = new_integral
        mash = create_mash(left_boundary, right_boundary, n)
        new_integral = calculate_integral(mash, f)
        
    print(f"Convergence reached after {num_of_iterations} loop iterations with n = {n} intervals.")
    return new_integral

# Example execution to demonstrate the implementation
if __name__ == "__main__":
    # Define SymPy symbols
    x = sp.Symbol('x')
    
    # Example from slide 3: Integral of sin(x) from 0 to pi
    # The actual analytical value is 2.0
    expr = sp.sin(x)
    a = 0
    b = sp.pi
    
    print("--- Running Simpson's Method ---")
    print(f"Integrating f(x) = {expr} from {a} to {b}")
    
    # Calculating with an epsilon of 1e-5
    result = sympson_method(a, b, expr, epsilon=1e-5)
    print(f"Numerical result: {result:.6f}")