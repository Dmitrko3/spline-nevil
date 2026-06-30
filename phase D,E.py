"""
Numerical Analysis - Stages D & E (Standalone Complete Version - Popup Fix)
"""

# --- POPUP FIX FOR PYCHARM ---
import matplotlib
matplotlib.use('TkAgg') # Forces an external popup window
import matplotlib.pyplot as plt
# -----------------------------

# ==========================================
# 1. Core Mathematical Functions
# ==========================================

def cubic_spline_coeff(x_nodes, y_nodes):
    """ Calculates the original Spline 'c' coefficients (Stage A) """
    n = len(x_nodes)
    h = [x_nodes[i + 1] - x_nodes[i] for i in range(n - 1)]
    alpha, mu, z = [0] * n, [0] * n, [0] * n

    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (y_nodes[i + 1] - y_nodes[i]) - (3 / h[i - 1]) * (y_nodes[i] - y_nodes[i - 1])
        l_val = 2 * (x_nodes[i + 1] - x_nodes[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l_val
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l_val

    c = [0] * n
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
    return c

def f_spline(x, x_nodes, y_nodes, c_coeffs):
    """ Evaluates the Spline function at point x """
    n = len(x_nodes)
    h = [x_nodes[i + 1] - x_nodes[i] for i in range(n - 1)]

    i = 0
    while i < n - 2 and x > x_nodes[i + 1]:
        i += 1

    b = (y_nodes[i + 1] - y_nodes[i]) / h[i] - h[i] * (c_coeffs[i + 1] + 2 * c_coeffs[i]) / 3
    d = (c_coeffs[i + 1] - c_coeffs[i]) / (3 * h[i])
    dx = x - x_nodes[i]
    return y_nodes[i] + b * dx + c_coeffs[i] * (dx ** 2) + d * (dx ** 3)

def richardson_derivative(f, x, h, order=2):
    """ Calculates derivative using Richardson Extrapolation """
    A_h = (f(x + h) - f(x - h)) / (2 * h)
    half_h = h / 2
    A_h2 = (f(x + half_h) - f(x - half_h)) / (2 * half_h)
    factor = 2 ** order
    L = (factor * A_h2 - A_h) / (factor - 1)
    return L

def gaussian_elimination(A, B):
    """ Solves linear system using Gaussian Elimination """
    n = len(B)
    A = [row[:] for row in A]
    B = B[:]
    for i in range(n):
        pivot = A[i][i]
        for j in range(i+1, n):
            factor = A[j][i] / pivot
            B[j] -= factor * B[i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
    x = [0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (B[i] - sum_ax) / A[i][i]
    return x

def gauss_seidel(A, B, iterations=100, tol=1e-8):
    """ Solves linear system using Gauss-Seidel Method """
    n = len(B)
    x = [0.0] * n
    for _ in range(iterations):
        x_new = x[:]
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (B[i] - s1 - s2) / A[i][i]
        if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
            return x_new
        x = x_new
    return x


# ==========================================
# 2. Main Execution (Stages D & E)
# ==========================================

def main():
    points = [
        (0.0, 0.500), (0.3, 0.779), (0.7, -0.054), (1.0, -0.591),
        (1.5, -0.368), (2.0, -0.715), (2.5, -0.056), (3.0, 0.976)
    ]
    x_nodes = [p[0] for p in points]
    y_nodes = [p[1] for p in points]

    c_original = cubic_spline_coeff(x_nodes, y_nodes)

    def f_spline_wrapper(x):
        return f_spline(x, x_nodes, y_nodes, c_original)

    # ---------------------------------------------------------
    # STAGE D: Derivatives & Critical Points
    # ---------------------------------------------------------
    print("\n" + "="*60)
    print(" STAGE D: Derivatives & Critical Points")
    print("="*60)

    print("{:<5} | {:<8} | {:<17} | {:<17}".format("x", "f(x)", "f'(x) [h=0.01]", "f'(x) [h=0.001]"))
    print("-" * 55)
    for x in x_nodes:
        d_01 = richardson_derivative(f_spline_wrapper, x, h=0.01, order=2)
        d_001 = richardson_derivative(f_spline_wrapper, x, h=0.001, order=2)
        y_val = f_spline_wrapper(x)
        print("{:<5.2f} | {:>8.3f} | {:>17.6f} | {:>17.6f}".format(x, y_val, d_01, d_001))

    num_points = 200
    test_x = [i * (3.0 / (num_points - 1)) for i in range(num_points)]
    f_prime_vals = [richardson_derivative(f_spline_wrapper, tx, h=0.001) for tx in test_x]

    critical_x = []
    critical_y = []
    for i in range(len(test_x) - 1):
        if f_prime_vals[i] * f_prime_vals[i+1] <= 0:
            root_x = (test_x[i] + test_x[i+1]) / 2.0
            critical_x.append(root_x)
            critical_y.append(0)

    print("\nCritical Points (Min/Max approx locations):")
    for cx in critical_x:
        print(f"-> Point near x = {cx:.4f} (Derivative is 0)")

    # שרטוט גרף 3 (שלב ד')
    plt.figure(figsize=(10, 5))
    plt.plot(test_x, f_prime_vals, label="f'(x) [Richardson Extrapolation]", color="green")
    plt.axhline(0, color='black', linewidth=1, linestyle='--')
    plt.scatter(critical_x, critical_y, color="red", zorder=5, label="Critical Points (Roots of f'(x))")
    plt.title("Stage D - Graph 3: Derivative f'(x) and Critical Points")
    plt.xlabel("x")
    plt.ylabel("f'(x)")
    plt.legend()
    plt.grid(True)
    print("\n>>> OPENING GRAPH 1: Please close the popup window to continue the code... <<<")
    plt.show(block=True) # Forces code to stop until window is closed

    # ---------------------------------------------------------
    # STAGE E: Internal Spline Linear System
    # ---------------------------------------------------------
    print("\n" + "="*60)
    print(" STAGE E: Internal Spline Linear System")
    print("="*60)

    n = len(x_nodes)
    h = [x_nodes[i + 1] - x_nodes[i] for i in range(n - 1)]

    size = n - 2
    A = [[0.0] * size for _ in range(size)]
    B = [0.0] * size

    for i in range(1, n - 1):
        idx = i - 1
        if idx > 0:
            A[idx][idx - 1] = h[i - 1]
        A[idx][idx] = 2 * (h[i - 1] + h[i])
        if idx < size - 1:
            A[idx][idx + 1] = h[i]
        B[idx] = (3 / h[i]) * (y_nodes[i + 1] - y_nodes[i]) - (3 / h[i - 1]) * (y_nodes[i] - y_nodes[i - 1])

    print("1. Tridiagonal Matrix A (6x6):")
    for row in A:
        print(["{:.3f}".format(val) for val in row])

    print("\nVector B:")
    print(["{:.3f}".format(val) for val in B])

    sol_gauss = gaussian_elimination(A, B)
    sol_seidel = gauss_seidel(A, B)
    internal_c_stage_A = c_original[1:-1]

    print("\n2 & 3. Solutions Comparison (Stage E vs. Stage A):")
    print("{:<8} | {:<15} | {:<15} | {:<15}".format("Coeff", "Gauss Method", "Seidel Method", "Stage A (Original)"))
    print("-" * 65)
    labels = []
    for i in range(len(sol_gauss)):
        label = f"c[{i+1}]"
        labels.append(label)
        print("{:<8} | {:<15.6f} | {:<15.6f} | {:<15.6f}".format(
            label, sol_gauss[i], sol_seidel[i], internal_c_stage_A[i]
        ))

    # שרטוט גרף (שלב ה')
    x_indices = range(len(labels))
    width = 0.25

    plt.figure(figsize=(10, 5))
    plt.bar([x - width for x in x_indices], sol_gauss, width=width, label='Gaussian Elimination', color='blue')
    plt.bar(x_indices, sol_seidel, width=width, label='Gauss-Seidel', color='orange')
    plt.bar([x + width for x in x_indices], internal_c_stage_A, width=width, label='Stage A (Original)', color='green')

    plt.title("Stage E - Graph: Comparison of Spline 'c' Coefficients")
    plt.xlabel("Coefficient Index")
    plt.ylabel("Calculated Value")
    plt.xticks(x_indices, labels)
    plt.legend()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    print("\n>>> OPENING GRAPH 2: Please close the popup window to finish the code... <<<")
    plt.show(block=True) # Forces code to stop until window is closed

if __name__ == "__main__":
    main()