def matrix_input(n):
    """Get matrix input from user with validation."""
    print(f"Enter {n}x{n} matrix values row by row (space-separated numbers):")
    matrix = []
    for i in range(n):
        while True:
            try:
                # Get row input and convert to float list
                row = input(f"Row {i+1}: ").strip().split()
                if len(row) != n:
                    print(f"Please enter exactly {n} numbers")
                    continue
                row = [float(x) for x in row]
                matrix.append(row)
                break
            except ValueError:
                print("Please enter valid numbers")
    return matrix

def minor_matrix(matrix, row, col):
    """Create a minor matrix by removing specified row and column."""
    return [[matrix[i][j] for j in range(len(matrix)) if j != col]
            for i in range(len(matrix)) if i != row]

def determinant(matrix):
    """Calculate determinant using Laplace expansion along first row."""
    n = len(matrix)
    
    # Base cases
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    
    # Recursive case using Laplace expansion along first row
    det = 0
    for j in range(n):
        det += matrix[0][j] * (-1) ** j * determinant(minor_matrix(matrix, 0, j))
    return det

def print_matrix(matrix):
    """Print matrix in a formatted way."""
    for row in matrix:
        print(" ".join(f"{x:8.2f}" for x in row))

def main():
    # Get matrix size
    while True:
        try:
            n = int(input("Enter matrix size n (for nxn matrix): "))
            if n <= 0:
                print("Please enter a positive number")
                continue
            break
        except ValueError:
            print("Please enter a valid integer")
    
    # Get matrix input
    matrix = matrix_input(n)
    
    # Print the input matrix
    print("\nYour matrix:")
    print_matrix(matrix)
    
    # Calculate and print determinant
    result = determinant(matrix)
    print(f"\nDeterminant = {result:.2f}")

if __name__ == "__main__":
    main()