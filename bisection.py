from os import name

#Gemini chat for menu,calculation,etc.
#https://www.google.com/search?sca_esv=19ca41183bcca24c&sxsrf=ANbL-n4Sp8StCY3HPYXlFrVXM4D6iXsicw%3A1778048309859&fbs=ADc_l-aN0CWEZBOHjofHoaMMDiKpaEWjvZ2Py1XXV8d8KvlI3vWUtYx0DZdicpfE1faGYemg2KC4yuMPyQlIvlWqq2AtcdVMJmMDffRprXURy79lwWYbYYV2Zn-_3zemGi58ZoYQT-gc_TefnhK7O6tgh3cJbnRhsfD9J2bnFesALGy1NtfkyaEXUF5U34PYt9QRs96v3tWGXCPutdqk722sMCshCQjZQw&ntc=1&sa=X&ved=2ahUKEwijuIOXgqSUAxUdg_0HHZ9HCqQQoo4PegYIAQgAEAM&biw=1440&bih=765&dpr=2&q=make+function+that+does+the+next+thing%3A%0A+1.sets+polynom%28calls+a+fuction+for+this%29%0A2.sets+the+area+where+the+answer+known%0A3.sets+Number+of+segments+to+divide+the+range+defined+in+the+previous+section+or+the+size+of+the+small+segment.+%28You+can+divide+into+segments+of+0.1%29%0A4.the+program+will+ask+which+function+use+to+calc+the+roots%2Cfor+this+use+while%2C+3+options+0+to+exit+and+default+option%3Aerror+try+again%0A%0Awrite+this+code+in+python&mstk=AUtExfBG63hb7mf4-xgfLsICHKuGnlpT0S7jVY--23wuhuVaZGkXzfvmS8aAkUSepdyNZF6boXDBbDQTHsETSg3g8gepW6lH4nmGbWhDLf-zIaaLDugJwZBxEJ4Ef1-YM_hhe5tbmonFqoHPkSOfCvMk9ExO5iqbvAiRx5EqqfRNIplgHMUDKk0EOd6-YVVeCouOVUVY3oKD4P9hnTDw92JT0Rkl-qJQ9BUVIxlIowyr-jai9Caz_YSp44-STppPMWHy6MixWr_emFwMuO6Jax4rb4UWWZI_zAKiYBE&csuir=1&udm=50
#Main function
def calculate_roots():
    #Main of the program: gets data, sets data, and executes the choice
    # 1. Set polynomial
    f=lambda x:8*x**3-12*x**2+15*x-1
    # 2. Set the range/area
    print("--- Define Range ---")
    start= float(input("Range start"))
    end = float(input("Range End"))

    # 3. Set segments
    # Dividing the range into steps (e.g., 0.1)
    step = float(input("segment size (e.g., 0.1) "))
    # 4. Menu loop
    while True:
        print("\n--- Root Calculation Menu ---")
        print("1. Method A (e.g., Bisection)")
        print("2. Method B (e.g., Newton-Raphson)")
        print("3. Method C (e.g., Secant)")
        print("0. Exit")

        choice = input("Select an option: ")

        if choice == "1":
            print("1. Method A (e.g., Bisection)")
            Bisection(f,start,end)
            break

        elif choice == "2":
            print("2. Method B (e.g., Newton-Raphson)")
            Newton_Raphson(f,start,end)
            break

        elif choice == "3":
            print("3. Method C (e.g., Secant)")
            secant_method(f,start,end)
            break

        elif choice == "0":
            print("Exiting program.")
            break

        else:
            print("Error: Invalid input. Try again.")
#=========================================================================

#=========================================================================
##Root calculation methods
#https://gemini.google.com/app/89cd4cfd6cce7b75
def Bisection(poly_expression, point_start, point_end, epsilon=0.0001):
    # Helper to evaluate the polynomial at start and end
    f_start = poly_expression(point_start)
    f_end = poly_expression(point_end)

    # Check the Intermediate Value Theorem requirement
    if f_start * f_end > 0:
        print("Bisection requirement failed: f(start) and f(end) must have opposite signs.")
        return None

    iteration_count = 0
    mid = point_start

    # Loop until the interval is smaller than epsilon
    while (point_end - point_start) / 2 > epsilon:
        iteration_count += 1
        mid = (point_start + point_end) / 2
        f_mid = poly_expression(mid)

        # If we hit the root exactly or are within epsilon
        if abs(f_mid) < epsilon:
            break

        # Determine which side to keep
        if poly_expression(point_start) * f_mid < 0:
            point_end = mid
        else:
            point_start = mid

    print(f"Root found at: {mid}")
    print(f"Total iterations (convergences): {iteration_count}")
    return mid









#===========================================================================
## Newton Raphson Method
#Gemini chat for menu,calculation,etc.
#https://gemini.google.com/share/fe600a2d640f
def Newton_Raphson(poly_expression, point_start, point_end, epsilon=0.0001):

    """
    Finds a root of a polynomial within a given range using the Newton Raphson method.

    Parameters:
    poly_expression (sympy.Expr): The symbolic polynomial expression.
    point_start (float): The start of the interval.
    point_end (float): The end of the interval.
    epsilon (float): The desired precision (tolerance).

    Returns:
    float: The approximated root if found, otherwise None.
    """

    # Initial guess: midpoint of the range
    x_curr = (point_start + point_end) / 2

    iteration_count = 0
    max_iterations = 100

    # Loop until the method converges or reaches max iterations
    while iteration_count < max_iterations:
        iteration_count += 1

        # Evaluate f(x)
        f_val = poly_expression(x_curr)

        # Calculate derivative using a small difference
        h = 0.000001
        f_prime_val = (poly_expression(x_curr + h) - poly_expression(x_curr - h)) / (2 * h)

        # Avoid division by zero
        if abs(f_prime_val) < 1e-12:
            print("Newton-Raphson failed: derivative is zero.")
            return None

        # Newton Raphson formula
        x_next = x_curr - f_val / f_prime_val

        # Check if the answer is close enough
        if abs(x_next - x_curr) < epsilon:
            print(f"Root found at: {x_next}")
            print(f"Number of iterations: {iteration_count}")
            return x_next

        # Update x for the next iteration
        x_curr = x_next

    print("Method did not converge.")
    return None




#===========================================================================
# secant method
# claude link:
# https://claude.ai/share/a4e1def6-2e75-4274-abfe-01e8faff6c6d

def secant_method(polynomial, start_point, end_point, epsilon=0.0001, max_iterations=1000):
    """
    Finds a root of a polynomial using the Secant Method.

    Parameters:
        polynomial    : A callable f(x) representing the polynomial.
        start_point   : The left boundary of the search interval (x0).
        end_point     : The right boundary of the search interval (x1).
        epsilon       : Convergence tolerance (default: 0.0001).
        max_iterations: Safety cap on iterations before declaring non-convergence.

    Returns:
        The approximate root, or None if the method does not converge.
    """
    x0 = start_point
    x1 = end_point

    for iteration in range(1, max_iterations + 1):
        f0 = polynomial(x0)
        f1 = polynomial(x1)

        # Guard against division by zero (flat secant line)
        if abs(f1 - f0) < 1e-12:
            print("Method failed: denominator (f(x1) - f(x0)) is zero — the secant line is flat.")
            return None

        # Core secant formula: x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

        # Check whether the new point has left the original interval
        if x2 < min(start_point, end_point) or x2 > max(start_point, end_point):
            print(f"Warning: iterate x{iteration+1} = {x2:.6f} stepped outside "
                  f"the interval [{min(start_point, end_point)}, {max(start_point, end_point)}]. "
                  "Continuing anyway.")

        # Convergence check: root found when |f(x2)| < epsilon
        if abs(polynomial(x2)) < epsilon:
            print(f"Root found after {iteration} iteration(s).")
            print(f"Approximate root : x ≈ {x2:.6f}")
            print(f"Residual |f(x)|  : {abs(polynomial(x2)):.2e}")
            return x2

        # Slide the window forward
        x0, x1 = x1, x2

    print(f"Method did not converge within {max_iterations} iterations.")
    return None

calculate_roots()
