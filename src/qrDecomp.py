# *-* coding: utf-8 *-*

from Matrix import *


''' Вращения Гивенса '''

# Домножение на матрицу Гивенса слева: G(i, j, c, s) * A
# NB: изменяет переданную матрицу
def givens_mul(A, i, j, c, s):
    n = A.width
    new_i, new_j = [0] * n, [0] * n
    for k in range(n):
        new_i[k] = c * A[i][k] + s * A[j][k]
        new_j[k] = c * A[j][k] - s * A[i][k]
    A[i], A[j] = new_i[:], new_j[:]
    
# Домножение на матрицу Гивенса справа: A * G(i, j, c, s)
# NB: изменяет переданную матрицу
def givens_mul_right(A, i, j, c, s):
    n = A.width
    new_i, new_j = [0] * n, [0] * n
    for k in range(n):
        new_i[k] = c * A[k][i] - s * A[k][j]
        new_j[k] = c * A[k][j] + s * A[k][i]
    for k in range(n):
        A[k][i], A[k][j] = new_i[k], new_j[k]


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


''' Отражения Хаусхолдера '''

# Домножение на матрицу Хаусхолдера слева: (I - 2vv^T) * A
# NB: не изменяет переданную матрицу
def householder_mul(A, v):
    return A - (2 * v) * (v.transpose() * A)

# Домножение на матрицу Хаусхолдера справа: A * (I - 2vv^T)
# NB: не изменяет переданную матрицу
def householder_mul_right(A, v):
    return A - (A * (2 * v)) * v.transpose()


# Нахождение QR-разложения при помощи отражений Хаусхолдера
#  Input: квадратная матрица A
# Output: ортогональная и верхнетреугольная матрицы Q, R такие, что A = QR
def qr_decomp_householder(A):
    # Проверка
    A.square()
    
    n = A.width
    
    # Q и R будут накапливать в себе ответ
    Q = Matrix.unit(n)
    R = A.copy()
    
    for j in range(n):
        # j-ый столбец матрицы R, у которого занулили первые j чисел
        u = Matrix.vector([
            0 if i < j else R[i][j] for i in range(n)
        ]).normalize() # поделили на норму
        
        # Стандартный базисный вектор e_j
        e_j = Matrix.vector([int(i == j) for i in range(n)])
        
        # Вектор, вокруг которого будем отражать
        v = (u - e_j).normalize()
        
        Q = householder_mul(Q, v) # Q = Q * H_v
        R = householder_mul(R, v) # R = R * H_v
    
    Q = Q.transpose() # Q = Q^T
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    R = R.map(lambda x: 0 if eq(x, 0) else x)
    
    return Q, R