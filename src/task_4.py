# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_3 import *


''' QR-разложение при помощи вращений Гивенса '''

def task_4():
    A = input_matrix()
    Q, R = qr_decomp_givens(A)
    print(f'Ответ:\nQ = \n{Q}\nR = \n{R}')


# Нахождение QR-разложения при помощи вращений Гивенса
#  Input: квадратная матрица A
# Output: ортогональная и верхнетреугольная матрицы Q, R такие, что A = QR
def qr_decomp_givens(A):
    # Проверка
    A.square()
    
    n = A.width
    
    # Q и R будут накапливать в себе ответ
    Q = Matrix.unit(n)
    R = A.copy()
    
    for j in range(n - 1):
        # Ищем первый ненулевой элемент ниже диагонали в столбце j
        k = j
        while k < n and eq(R[k][j], 0):
            k += 1
            
        # Все элементы нулевые => пропускаем этот шаг
        if k == n:
            continue
        
        # Ставим первый ненулевой элемент на диагональ
        # То есть теперь он стоит в ячейке A[j][j]
        if k != j:
            givens_mul(Q, j, k, 0, 1) # Q = G(j, k, 0, 1) * Q
            givens_mul(R, j, k, 0, 1) # R = G(j, k, 0, 1) * R
        
        for i in range(j + 1, n):
            # Ищем подходящий угол вращения
            c = 1 / (1 + R[i][j] ** 2 / R[j][j] ** 2) ** 0.5
            s = R[i][j] / R[j][j] * c
            
            givens_mul(Q, j, i, c, s) # Q = G(j, i, c, s) * Q
            givens_mul(R, j, i, c, s) # R = G(j, i, c, s) * R
    
    Q = Q.transpose() # Q = Q^T
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    R = R.map(lambda x: 0 if eq(x, 0) else x)
    
    return Q, R


# Нахождение QR-разложения при помощи вращений Гивенса
#  Input: квадратная трёхдиагональная матрица A
# Output: ортогональная и верхнетреугольная матрицы Q, R такие, что A = QR
def qr_decomp_givens_tridiag(A):
    # Проверка
    A.square()
    
    n = A.width
    
    # Q и R будут накапливать в себе ответ
    # Q - список матриц Гивенса, то есть список кортежей вида (i, j, c, s)
    Q = []
    R = A.copy()
    
    for j in range(n - 1):
        # Матрица трёхдиагональна, поэтому числа после maxN равны 0
        maxN = min(j + 2, n)
        
        # Ищем первый ненулевой элемент ниже диагонали в столбце j
        k = j
        while k < maxN and eq(R[k][j], 0):
            k += 1
        
        # Все элементы нулевые => пропускаем этот шаг
        if k == maxN:
            continue
        
        # Ставим первый ненулевой элемент на диагональ
        # То есть теперь он стоит в ячейке A[j][j]
        if k != j:
            Q.append((j, k, 0, 1))    # Q = G(j, k, 0, 1) * Q
            givens_mul(R, j, k, 0, 1) # R = G(j, k, 0, 1) * R
        
        for i in range(j + 1, maxN):
            # Ищем подходящий угол вращения
            c = 1 / (1 + R[i][j] ** 2 / R[j][j] ** 2) ** 0.5
            s = R[i][j] / R[j][j] * c
            
            Q.append((j, i, c, s))    # Q = G(j, i, c, s) * Q
            givens_mul(R, j, i, c, s) # R = G(j, i, c, s) * R
    
    # Q = Q.transpose()
    Q = [(i, j, c, -s) for (i, j, c, s) in Q]
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    R = R.map(lambda x: 0 if eq(x, 0) else x)
    
    return Q, R