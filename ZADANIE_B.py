import math
import matplotlib.pyplot as plt
import time

def calculate_row_value(matrix, current_estimate, target_vector, row_index):
    # Sumuje produkty z wyłączeniem elementu diagonalnego
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


def measure_performance(solver, system_matrix, vector_b):
    start = time.time()
    solution, iterations, residuals = solver(system_matrix, vector_b)
    elapsed_time = time.time() - start
    return solution, iterations, residuals, elapsed_time

jacobi_solution, jacobi_iters, jacobi_residuals, time_for_jacobi = measure_performance(solve_jacobi, A, b)

seidel_solution, seidel_iters, seidel_residuals, time_for_seidel = measure_performance(perform_gauss_seidel, A, b)

plt.figure(figsize=(10, 6))
plt.semilogy(range(jacobi_iters), jacobi_residuals, label='Residuals of Jacobi')
plt.semilogy(range(seidel_iters), seidel_residuals, label='Residuals of Gauss-Seidel')
plt.xlabel('Iteration Number')
plt.ylabel('Logarithm of Residual Norm')
plt.title('Comparison of Convergence Between Jacobi and Gauss-Seidel')
plt.legend()
plt.show()

print(f"Jacobi Method: Iterations = {jacobi_iters}, Execution Time = {time_for_jacobi:.2f} seconds")
print(f"Gauss-Seidel Method: Iterations = {seidel_iters}, Execution Time = {time_for_seidel:.2f} seconds")
