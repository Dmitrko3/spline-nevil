def cubic_spline(x_nodes, y_nodes, target_x):
    # Sort the data points by x to ensure they are in order
    x, a = zip(*sorted(zip(x_nodes, y_nodes)))
    n = len(x)

    # h is the step size between consecutive x points
    h = [x[i + 1] - x[i] for i in range(n - 1)]

    # Forward sweep for the Thomas Algorithm
    alpha, mu, z = [0] * n, [0] * n, [0] * n
    for i in range(1, n - 1):
        # Right-hand side of the tridiagonal system
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])

        # Solving the diagonals
        l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l

    # Back substitution to find 'c' (the second derivatives)
    c = [0] * n
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]

    # Find the correct interval 'i' for our target_x
    i = 0
    while i < n - 2 and target_x > x[i + 1]:
        i += 1

    # Calculate 'b' and 'd' ONLY for the interval we care about
    b = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3
    d = (c[i + 1] - c[i]) / (3 * h[i])

    # Evaluate the cubic polynomial for this specific interval
    dx = target_x - x[i]
    return a[i] + b * dx + c[i] * (dx ** 2) + d * (dx ** 3)
