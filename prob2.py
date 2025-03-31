# Function to subtract lambda * identity from matrix A
def subtract_lambda_identity(matrix, lambda_val, size):
    result = []
    for i in range(size):
        row = []
        for j in range(size):
            if i == j:
                row.append(matrix[i][j] - lambda_val)
            else:
                row.append(matrix[i][j])
        result.append(row)
    return result

# Function to compute determinant of a 2x2 matrix
def determinant_2x2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

# Function to compute determinant of a 3x3 matrix (Rule of Sarrus)
def determinant_3x3(matrix):
    pos = (matrix[0][0] * matrix[1][1] * matrix[2][2] +
           matrix[0][1] * matrix[1][2] * matrix[2][0] +
           matrix[0][2] * matrix[1][0] * matrix[2][1])
    neg = (matrix[0][2] * matrix[1][1] * matrix[2][0] +
           matrix[0][1] * matrix[1][0] * matrix[2][2] +
           matrix[0][0] * matrix[1][2] * matrix[2][1])
    return pos - neg

# Function to evaluate polynomial at a given x
def evaluate_polynomial(coeffs, x):
    result = 0
    for i, coeff in enumerate(coeffs[::-1]):
        result += coeff * (x ** i)
    return result

# Bisection method to find a root
def find_root_bisection(coeffs, a, b, tol=1e-6, max_iter=1000):
    fa, fb = evaluate_polynomial(coeffs, a), evaluate_polynomial(coeffs, b)
    if fa * fb > 0:
        return None  # No sign change
    if fa == 0:
        return a
    if fb == 0:
        return b
    for _ in range(max_iter):
        mid = (a + b) / 2
        f_mid = evaluate_polynomial(coeffs, mid)
        if abs(f_mid) < tol:
            return mid
        if fa * f_mid < 0:
            b = mid
            fb = f_mid
        else:
            a = mid
            fa = f_mid
    return (a + b) / 2

# Solve quadratic equation ax^2 + bx + c = 0
def solve_quadratic(a, b, c):
    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return []
    elif discriminant == 0:
        return [-b / (2 * a)]
    else:
        root1 = (-b + (discriminant ** 0.5)) / (2 * a)
        root2 = (-b - (discriminant ** 0.5)) / (2 * a)
        return [root1, root2]

# Find eigenvalues for 2x2 or 3x3 matrix
def find_eigenvalues(matrix, size):
    if size == 2:
        a, b = matrix[0][0], matrix[0][1]
        c, d = matrix[1][0], matrix[1][1]
        coeff_a = 1
        coeff_b = -(a + d)
        coeff_c = a * d - b * c
        return solve_quadratic(coeff_a, coeff_b, coeff_c)
    elif size == 3:
        trace = matrix[0][0] + matrix[1][1] + matrix[2][2]
        det = determinant_3x3(matrix)
        m1 = determinant_2x2([[matrix[1][1], matrix[1][2]], [matrix[2][1], matrix[2][2]]])
        m2 = determinant_2x2([[matrix[0][0], matrix[0][2]], [matrix[2][0], matrix[2][2]]])
        m3 = determinant_2x2([[matrix[0][0], matrix[0][1]], [matrix[1][0], matrix[1][1]]])
        minors_sum = m1 + m2 + m3
        coeffs = [1, -trace, minors_sum, -det]
        
        # Adaptive root finding
        max_val = max(abs(min([min(row) for row in matrix])), abs(max([max(row) for row in matrix])))
        bounds = [-max_val - 1, max_val + 1]
        roots = []
        step = 0.5  # Check points to find sign changes
        points = [bounds[0] + i * step for i in range(int((bounds[1] - bounds[0]) / step) + 1)]
        
        for i in range(len(points) - 1):
            a, b = points[i], points[i + 1]
            root = find_root_bisection(coeffs, a, b)
            if root is not None and not any(abs(root - r) < 1e-5 for r in roots):
                roots.append(root)
        
        return sorted(roots, reverse=True) if roots else []

# Find eigenvector for a given eigenvalue
def find_eigenvector(matrix, eigenvalue, size):
    A_minus_lambda_I = subtract_lambda_identity(matrix, eigenvalue, size)
    if size == 2:
        a, b = A_minus_lambda_I[0][0], A_minus_lambda_I[0][1]
        c, d = A_minus_lambda_I[1][0], A_minus_lambda_I[1][1]
        if b != 0 and a != 0:
            y = 1
            x = -b * y / a
            return [x, y]
        elif a == 0:
            x = 1
            return [x, 0]
        else:
            y = 1
            return [0, y]
    elif size == 3:
        a11, a12, a13 = A_minus_lambda_I[0]
        a21, a22, a23 = A_minus_lambda_I[1]
        a31, a32, a33 = A_minus_lambda_I[2]
        if a11 != 0:
            y = (a13 * a21 - a11 * a23) / (a11 * a22 - a12 * a21 + 1e-10)  # Avoid division by zero
            x = (-a12 * y - a13) / (a11 + 1e-10)
            return [x, y, 1]
        elif a12 != 0:
            x = 1
            y = -a13 / (a12 + 1e-10)
            return [x, y, 0]
        else:
            return [1, 0, 0]

# Main computation function
def compute_eigen(matrix, size):
    eigenvalues = find_eigenvalues(matrix, size)
    if not eigenvalues:
        return "No real eigenvalues found."
    result = {}
    for lambda_val in eigenvalues:
        eigenvector = find_eigenvector(matrix, lambda_val, size)
        result[lambda_val] = eigenvector
    return result

# Get matrix input from user
def get_matrix_input():
    while True:
        try:
            size = int(input("Enter matrix size (2 for 2x2, 3 for 3x3): "))
            if size not in [2, 3]:
                print("Only 2x2 or 3x3 matrices are supported.")
                continue
            break
        except ValueError:
            print("Enter 2 or 3.")
    
    print(f"Enter a {size}x{size} matrix.")
    matrix = []
    for i in range(size):
        while True:
            try:
                row = input(f"Enter row {i+1} ({size} numbers separated by space): ").strip().split()
                if len(row) != size:
                    print(f"Please enter exactly {size} numbers.")
                    continue
                row = [float(x) for x in row]
                matrix.append(row)
                break
            except ValueError:
                print("Invalid input. Enter numbers only.")
    return matrix, size

# Main execution
def main():
    print("Eigenvalue and Eigenvector Calculator for 2x2 and 3x3 Matrices")
    matrix, size = get_matrix_input()
    
    print("\nYour matrix:")
    for row in matrix:
        print(row)
    
    result = compute_eigen(matrix, size)
    if isinstance(result, str):
        print(f"\nResult: {result}")
    else:
        print("\nResults:")
        for eigenvalue, eigenvector in result.items():
            print(f"Eigenvalue: {eigenvalue:.2f}")
            print(f"Eigenvector: [{', '.join(f'{x:.2f}' for x in eigenvector)}]")
            print(f"Eigenspace: Span of [{', '.join(f'{x:.2f}' for x in eigenvector)}]")
            print()

if __name__ == "__main__":
    main()