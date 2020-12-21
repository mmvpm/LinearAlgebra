# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_4 import *



''' QR-алгоритм для трёхдиагональных матриц '''

def task_10():
    tA = input_matrix()
    eps = input_eps()
    lms, q = qr_algorithm_tridiag(tA, eps)
    print(f'Ответ: \nСобственные числа = \n{lms}\nQ = \n{q}')


#  Input: квадратная симметричная трёхдиагональная матрица A, точность eps
# Output: список из приближённых собственных чисел evs, 
#         матрица примерно из собственных векторов Q
def qr_algorithm_tridiag(A, eps = global_eps):
    # Проверка
    A.square()
    A.symmetric()
    A.tridiagonal()
    
    # Накапливаем ответ в матрице Q
    Q = Matrix.unit(A.width)
    
    while not A.gc_check(eps):
        # Вместо создания R_k сразу изменяем A
        Q_k, A = qr_decomp_givens_tridiag(A) # Q_k, R_k = qr_decomp(A)
        # Q_k на самом деле является списком матриц Гивенса
        # Как раз отсюда и получается, что шаг алгоритма работает за O(n^2)
        for (i, j, c, s) in Q_k:
            givens_mul_right(A, i, j, c, s) 
            givens_mul_right(Q, i, j, c, s)
    
    # Собственные числа стоят на диагонали матрицы A
    evs = [A[i][i] for i in range(A.width)]
    
    return evs, Q