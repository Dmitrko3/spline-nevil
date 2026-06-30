import numpy as np
from colors import bcolors


def romberg_integration(func, a, b, n):
    """
    Romberg Integration

    Parameters:
    func (function): The function to be integrated.
    a (float): The lower limit of integration.
    b (float): The upper limit of integration.
    n (int): The number of iterations (higher value leads to better accuracy).

    Returns:
    numpy.ndarray: The full Romberg integration table of shape (n, n).
    float: The approximate definite integral (the final diagonal entry).
    """
    h = b - a
    R = np.zeros((n, n), dtype=float)

    # Initial trapezoidal approximation
    R[0, 0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, n):
        h /= 2
        sum_term = 0

        # Sum of midpoints
        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        # Trapezoidal update
        R[i, 0] = 0.5 * R[i - 1, 0] + h * sum_term

        # Richardson extrapolation columns
        for j in range(1, i + 1):
            R[i, j] = R[i, j - 1] + (R[i, j - 1] - R[i - 1, j - 1]) / ((4 ** j) - 1)

    return R, R[n - 1, n - 1]


def print_romberg_table(R):
    """
    Prints the Romberg integration table in a structured, lower-triangular format.
    """
    n = R.shape[0]
    print("\n" + "=" * 32 + " ROMBERG TABLE " + "=" * 32)
    
    # Headers correspond to the error bounds O(h^2), O(h^4), O(h^6), etc.
    headers = [f"O(h^{2*(j+1)})" for j in range(n)]
    header_row = "Row | " + " | ".join(f"{h:^14}" for h in headers)
    print(header_row)
    print("-" * len(header_row))
    
    for i in range(n):
        row_vals = []
        for j in range(i + 1):
            row_vals.append(f"{R[i, j]:.10f}")
        for j in range(i + 1, n):
            row_vals.append("")  # Leave upper triangular entries empty
        print(f"{i:3d} | " + " | ".join(f"{val:^14}" for val in row_vals))
    print("=" * len(header_row))


def f(x):
    return 1 / (2 + x ** 4)


if __name__ == '__main__':
    a = 0
    b = 3
    n = 5
    
    # Retrieve both the full matrix and the final approximation
    R_table, integral = romberg_integration(f, a, b, n)

    print(f" Division into n={n} sections ")
    print_romberg_table(R_table)
    
    print(bcolors.OKBLUE, f"\nApproximate integral in range [{a},{b}] is {integral:.12f}", bcolors.ENDC)