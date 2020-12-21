# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_4 import *


''' QR-алгоритм '''

def task_8():
    A = input_matrix()
    eps = input_eps()
    lms, q = qr_algorithm(A, eps)
    print(f'Ответ: \nСобственные числа = \n{lms}\nQ = \n{q}')


#  Input: квадратная симметричная матрица A, точность eps
# Output: список из приближённых собственных чисел evs, 
#         матрица примерно из собственных векторов Q
def qr_algorithm(A, eps = global_eps):
    # Проверка
    A.square()
    A.symmetric()
    
    # Матрица из собственных векторов
    Q = Matrix.unit(A.width)
    
    while not A.gc_check(eps):
        # Вместо создания R_k сразу изменяем A
        Q_k, A = qr_decomp_givens(A) # Q_k, R_k = qr_decomp(A)
        A *= Q_k # A = R_k * Q_k
        Q *= Q_k # Q = Q * Q_k
    
    # Собственные числа стоят на диагонали матрицы A
    evs = [A[i][i] for i in range(A.width)]
    
    return evs, Q