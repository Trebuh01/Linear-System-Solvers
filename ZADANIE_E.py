import time
import math
import matplotlib.pyplot as plt
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
def calculate_row_value(matrix, current_estimate, target_vector, row_index):
    sum_of_products = sum(matrix[row_index][col_index] * current_estimate[col_index]
                          for col_index in range(len(matrix)) if col_index != row_index)
    return (target_vector[row_index] - sum_of_products) / matrix[row_index][row_index]

def compute_residuals(matrix, next_solution, target_vector):
    return [target_vector[i] - sum(matrix[i][j] * next_solution[j] for j in range(len(matrix)))
            for i in range(len(target_vector))]

def compute_residual_norm(residuals):
    return math.sqrt(sum(res ** 2 for res in residuals))

def solve_jacobi(system_matrix, vector_b, error_tolerance=1e-9, max_loop=10000):
    dimension = len(vector_b)
    current_solution = [0.0] * dimension
    error_history = []

    for current_iteration in range(max_loop):
        next_solution = [calculate_row_value(system_matrix, current_solution, vector_b, row_index)
                         for row_index in range(dimension)]

        residuals = compute_residuals(system_matrix, next_solution, vector_b)
        norm_of_residuals = compute_residual_norm(residuals)
        error_history.append(norm_of_residuals)

        if norm_of_residuals < error_tolerance:
            break

        current_solution = next_solution

    return current_solution, current_iteration + 1, error_history



def calculate_new_value(matrix, current_estimate, target_vector, row):
    summed_terms = sum(matrix[row][col] * current_estimate[col] for col in range(len(matrix)) if col != row)
    return (target_vector[row] - summed_terms) / matrix[row][row]

def calculate_residual_norm(matrix, current_estimate, target_vector):
    residuals = [target_vector[i] - sum(matrix[i][col] * current_estimate[col] for col in range(len(matrix)))
                 for i in range(len(target_vector))]
    return math.sqrt(sum(res ** 2 for res in residuals))

def perform_gauss_seidel(matrix, target_vector, precision_threshold=1e-9, max_loops=10000):
    current_estimate = [0.0] * len(target_vector)
    errors_over_time = []

    for current_step in range(max_loops):
        for row in range(len(matrix)):
            current_estimate[row] = calculate_new_value(matrix, current_estimate, target_vector, row)

        residual_norm = calculate_residual_norm(matrix, current_estimate, target_vector)
        errors_over_time.append(residual_norm)

        if residual_norm < precision_threshold:
            break

    return current_estimate, current_step + 1, errors_over_time
sizes = [100, 500, 1000, 2000, 3000]
a1 = 5 + 3  # 8
a2 = -1
a3 = -1
f = 3

results_jacobi = []
results_gauss_seidel = []
results_lu = []

for N in sizes:

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

    b = [math.sin(i * (f + 1)) for i in range(1, N + 1)]

    # dla metody Jacobiego
    start_time = time.time()
    solve_jacobi(A, b)
    results_jacobi.append(time.time() - start_time)

    # dla metody Gaussa-Seidla
    start_time = time.time()
    perform_gauss_seidel(A, b)
    results_gauss_seidel.append(time.time() - start_time)

    # dla faktoryzacji LU
    start_time = time.time()
    L, U = decompose_lu(A)
    y = solve_lower_triangular(L, b)
    x = solve_upper_triangular(U, y)
    results_lu.append(time.time() - start_time)


plt.figure(figsize=(10, 6))
plt.plot(sizes, results_jacobi, label='Jacobi Method', marker='o')
plt.plot(sizes, results_gauss_seidel, label='Gauss-Seidel Method', marker='o')
plt.plot(sizes, results_lu, label='LU Decomposition', marker='o')

plt.xlabel('Size of Matrix N')
plt.ylabel('Time in seconds')
plt.title('Comparison of Solution Methods for Linear Systems')
plt.legend()
plt.grid(True)
plt.show()