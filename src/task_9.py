# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_5 import *


''' Трёхдиагональные матрицы '''

def task_9():
    A = input_matrix()
    tA, Q = get_tridiagonal(A)
    print(f'Ответ: \nA\' = \n{tA}\nQ = \n{Q}')


#  Input: квадратная симметричная матрица A
# Output: TD - тридиагонализация матрицы A, 
#         ортогональная матрица Q такая, что Q^T A Q = TD
def get_tridiagonal(A):
    # Проверка
    A.square()
    A.symmetric()
    
    n = A.width
    
    # TD и Q будут накапливать в себе ответ
    TD = A.copy()
    Q = Matrix.unit(n)
    
    # Вместо того, чтобы явно переходить к подматрицам, 
    # будем занулять несколько первых чисел у векторов
    for j in range(n - 1):
        # j-ый столбец матрицы TD, у которого занулили первые j + 1 чисел
        u = Matrix.vector([
            0 if i < j + 1 else TD[i][j] for i in range(n)
        ])
        
        if eq(u.norm(), 0):
            continue
        
        u = u.normalize() # поделили на норму
        
        # Стандартный базисный вектор e_{j + 1}
        e = Matrix.vector([int(i == j + 1) for i in range(n)])
        
        # Вектор, вокруг которого будем отражать
        v = u - e
        
        if eq(v.norm(), 0):
            continue
        
        v = v.normalize()

        TD = householder_mul(TD, v)       # TD = TD * H_v
        TD = householder_mul_right(TD, v) # TD = H_v * TD
        Q  = householder_mul(Q , v)       # Q  = Q  * H_v
    
    Q = Q.transpose() # Q = Q^T
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    TD = TD.map(lambda x: 0 if eq(x, 0) else x)
    
    return TD, Q