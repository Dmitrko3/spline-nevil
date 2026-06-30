# ==========================================
# 1. YOUR SIMPSON'S FUNCTIONS
# ==========================================
def create_mash(left_boundary, right_boundary, n):
    h = (right_boundary - left_boundary) / n
    return [float(left_boundary + i * h) for i in range(n + 1)]

def calculate_integral(mash, f):
    n = len(mash) - 1
    h = (mash[-1] - mash[0]) / n
    
    integral_sum = f(mash[0]) + f(mash[-1])
    for i in range(1, n):
        if i % 2 == 1:
            integral_sum += 4 * f(mash[i])
        else:
            integral_sum += 2 * f(mash[i])
            
    return (h / 3) * integral_sum


# ==========================================
# 2. ROMBERG FUNCTION (PURE PYTHON)
# ==========================================
def romberg_integration(func, a, b, n):
    h = b - a
    R = [[0.0] * n for _ in range(n)]
    R[0][0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, n):
        h /= 2
        sum_term = 0
        for k in range(1, 2 ** i, 2):
            sum_term += func(a + k * h)

        R[i][0] = 0.5 * R[i - 1][0] + h * sum_term

        for j in range(1, i + 1):
            R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / ((4 ** j) - 1)

    return R, R[n - 1][n - 1]


# ==========================================
# 3. YOUR RECONSTRUCTED FUNCTION FROM STEP A
# ==========================================
def your_reconstructed_function(x):
    """
    Replace the return statement below with your actual Spline or Neville 
    evaluation logic from Step A.
    
    If your Step A function is already named differently (e.g., cubic_spline(x)), 
    you can directly pass that function name as 'func' in the lines below.
    """
    # Placeholder calculation:
    return x ** 2 - 2 * x + 1


# ==========================================
# 4. EXECUTION AND COMPARISON
# ==========================================
if __name__ == "__main__":
    a = 0
    b = 3

    # --- Calculate Romberg ---
    # (Set 'n_romberg' to the number of rows you used for Romberg convergence)
    n_romberg = 5  
    R_table, I_Romberg = romberg_integration(your_reconstructed_function, a, b, n_romberg)

    # --- Calculate Simpson (Exactly n = 100) ---
    n_simpson = 100
    mash = create_mash(a, b, n_simpson)
    I_Simpson = calculate_integral(mash, your_reconstructed_function)

    # --- Relative Difference Check ---
    # Formula from PDF: |I_Romberg - I_Simpson| / |I_Romberg|
    relative_difference = abs(I_Romberg - I_Simpson) / abs(I_Romberg)
    
    # Check if agreement is less than 0.01 (1%) as specified in Step C
    agreement = relative_difference < 0.01

    # --- Output results for your report ---
    print("--- STEP C INTEGRATION VALUES ---")
    print(f"Romberg Integral (I_Romberg) : {I_Romberg:.10f}")
    print(f"Simpson's Integral (I_Simpson): {I_Simpson:.10f}")
    print(f"Relative Difference          : {relative_difference:.6f} ({relative_difference * 100:.4f}%)")
    print(f"Agreement (< 1% Check)       : {agreement}")