import math
import matplotlib.pyplot as plt
import time


def decompose_lu(matrix):
    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    upper = [[0.0] * size for _ in range(size)]

    for row in range(size):
        lower[row][row] = 1.0

        for col in range(row, size):
            partial_sum = sum(lower[row][k] * upper[k][col] for k in range(row))
            upper[row][col] = matrix[row][col] - partial_sum

        for col in range(row + 1, size):
            partial_sum = sum(lower[col][k] * upper[k][row] for k in range(row))
            lower[col][row] = (matrix[col][row] - partial_sum) / upper[row][row]

    return lower, upper

def solve_lower_triangular(triangular_matrix, vector):
    length = len(vector)
    solution = [0.0] * length
    for index in range(length):
        subtracted_sum = sum(triangular_matrix[index][j] * solution[j] for j in range(index))
        solution[index] = (vector[index] - subtracted_sum) / triangular_matrix[index][index]
    return solution

def solve_upper_triangular(triangular_matrix, vector):
    length = len(vector)
    solution = [0.0] * length
    for index in reversed(range(length)):
        subtracted_sum = sum(triangular_matrix[index][j] * solution[j] for j in range(index + 1, length))
        solution[index] = (vector[index] - subtracted_sum) / triangular_matrix[index][index]
    return solution

def compute_residual_norm(system_matrix, solution_vector, right_side_vector):
    size = len(solution_vector)
    residuals = [0.0] * size
    for i in range(size):
        computed_sum = sum(system_matrix[i][j] * solution_vector[j] for j in range(size))
        residuals[i] = right_side_vector[i] - computed_sum
    norm_of_residuals = math.sqrt(sum(residual ** 2 for residual in residuals))
    return norm_of_residuals



N = 916
a1 = 5 + 3  # 8
a2 = -1
a3 = -1
f = 3

A = [[0] * N for _ in range(N)]

# Wypełnianie głównej diagonali
for i in range(N):
    A[i][i] = a1

# Wypełnianie diagonali a2
for i in range(1, N):
    A[i][i - 1] = a2
    A[i - 1][i] = a2

# Wypełnianie diagonali a3
for i in range(2, N):
    A[i][i - 2] = a3
    A[i - 2][i] = a3

# Tworzenie wektora b

b = [math.sin(i * (f + 1)) for i in range(1, N+1)]

# LU
L, U = decompose_lu(A)

# Ly = b
y = solve_lower_triangular(L, b)

# Ux = y
x = solve_upper_triangular(U, y)

residual_norm = compute_residual_norm(A, x, b)
print("residuum norm:", residual_norm)
