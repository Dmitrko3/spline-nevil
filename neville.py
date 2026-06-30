from typing import List,Tuple

Point = Tuple[float, float]


def validate_input(points: List[Point], target_x: float) -> float | None:
    """
    Validates that the interpolation input is legal.
    If target_x is already in the table, returns the known y-value.
    Otherwise, returns None if validation passes.
    """
    # Check if target_x is already in the table
    for x, y in points:
        if x == target_x:
            return y
def neville_interpolation(points: List[Point], target_x: float) -> float:
    validate_input(points, target_x)
    n = len(points)
    p = [point[1] for point in points] # מערך חד ממדי שמאתחל עם ערכי ה-Y

    for j in range(1, n):
        for i in range(n - j):
            xi = points[i][0]
            xj = points[i + j][0]
            p[i] = ((target_x - xj) * p[i] - (target_x - xi) * p[i + 1]) / (xi - xj)

    return p[0]
