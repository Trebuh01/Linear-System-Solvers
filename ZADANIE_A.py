import math

N = 916
a1 = 5 + 3  # 8
a2 = -1
a3 = -1
f = 3

# Tworzenie macierzy A
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

print("Fragment lewego górnego rogu macierzy A:")
for row in A[:5]:
    print(row[:5])
print("Pierwsze elementy wektora b:")
print(b[:5])
