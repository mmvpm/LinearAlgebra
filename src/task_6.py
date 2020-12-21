# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_5 import *


''' QR-разложение при помощи отражений Хаусхолдера '''

def task_6():
    A = input_matrix()
    Q, R = qr_decomp_householder(A)
    print(f'Ответ:\nQ = \n{Q}\nR = \n{R}')


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
        
        Q = householder_mul(Q, v) # Q = H_v * Q
        R = householder_mul(R, v) # R = H_v * R
    
    Q = Q.transpose() # Q = Q^T
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    R = R.map(lambda x: 0 if eq(x, 0) else x)
    
    return Q, R