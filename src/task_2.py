# *-* coding: utf-8 *-*

from Matrix import *
from utility import *


''' Метод Гаусса-Зейделя '''

def task_2():
    A = input_matrix()
    b = input_vector()
    eps = input_eps()
    x = gauss_seidel(A, b, eps)
    print(f'Ответ: \n{x}')


# Раскладывает A на верхне- и нижнетреугольную матрицы: A = L - U
def LU_decompozition(A):
    L = Matrix.zero(A.height, A.width)
    U = Matrix.zero(A.height, A.width)
    for i in range(A.height):
        for j in range(A.width):
            L[i][j] = A[i][j] * (i >= j)
            U[i][j] = -A[i][j] * (i < j)
    return L, U

# Решает систему Lx = b, где L - нижнетреугольная
def Lx_solve(L, b):
    x = Matrix.zero(L.width, 1)
    for i in range(L.width):
        prefix = 0
        for j in range(i):
            prefix += L[i][j] * x[j][0]
        x[i][0] = (b[i][0] - prefix) / L[i][i]
    return x

#  Input: квадратная матрица A, вектор b, точность eps
# Output: вектор x такой, что ||Ax - b|| < eps
def gauss_seidel(A, b, eps = global_eps):
    # Проверка
    A.square()
    
    n = A.width
    
    # На диагонали матрицы A стоят нули
    if [A[i][i] for i in range(n)] == [0] * n:
        return 0
    
    # Случайный вектор
    x = Matrix.random(n, 1)
    # Разложение на верхне- и нижнетреугольную матрицы
    L, U = LU_decompozition(A)
    
    # Количество подряд идущих неудачных итераций
    count_fail = 0
    
    while (A * x - b).norm() >= eps:
        # Следующий вектор последовательности
        x_new = Lx_solve(L, U * x + b)
        
        # Удачная ли итерация
        if abs(x_new.norm() - x.norm()) >= 1:
            count_fail += 1
        else:
            count_fail = 0
        
        x = x_new
        
        # Последние 20 итераций были неудачными
        if count_fail > 20:
            return 0
    
    return x
