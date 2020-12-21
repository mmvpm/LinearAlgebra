# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_4 import *


''' Сдвиг '''

def task_11():
    tA = input_matrix()
    eps = input_eps()
    lms, q = qr_shift(tA, eps)
    print(f'Ответ: \nСобственные числа = \n{lms}\nQ = \n{q}')


#  Input: квадратная симметричная матрица A, точность eps
# Output: список из приближённых собственных чисел evs, 
#         матрица примерно из собственных векторов Q
def qr_shift(A, eps = global_eps):
    # Проверка
    A.square()
    A.symmetric()
    A.tridiagonal()
    
    # Список собственных чисел
    evs = []
    # Матрица из собственных векторов
    Q = Matrix.unit(A.width)
    
    # Переходим к подматрице, пока в A не остался один элемент
    while A.width > 1:
        # x[0], x[1] - собственные числа нижнего блока 2 на 2 матрицы A
        cf = A[-2][-2], A[-2][-1], A[-1][-2], A[-1][-1]
        D = (cf[0] - cf[3]) ** 2 + 4 * cf[1] * cf[2]
        x = (cf[0] + cf[3] + D ** 0.5) / 2, (cf[0] + cf[3] - D ** 0.5) / 2
        
        # Ближайшее к A[-1][-1] собственное число (то есть x[0] или x[1])
        v_k = x[0] if abs(x[0] - cf[3]) < abs(x[1] - cf[3]) else x[1]
        
        # Единичная матрица нужного размера
        E = Matrix.unit(A.width)
        
        # Обычный шаг QR-алгоритма, но со сдвигом на v_k * E
        Q_k, A = qr_decomp_givens_tridiag(A - v_k * E)
        for (i, j, c, s) in Q_k:
            givens_mul_right(A, i, j, c, s)
            givens_mul_right(Q, i, j, c, s)
        
        # Сдвигаем обратно
        A += v_k * E
        
        # Радиусы последнего круга Гершгорина для матриц A и A^T
        # То есть сумма внедиагональных элементов последней строки или столбца
        bottom, right = 0, 0
        for i in range(A.width - 1):
            bottom += abs(A[-1][i])
            right += abs(A[i][-1])
        
        # Если радиусы кругов < eps
        if bottom < eps and right < eps:
            # A[-1][-1] - приближение к собственному числу A
            evs.append(A[-1][-1])
            # Переходим к подматрице на единицу меньшего размера
            A = Matrix([A[i][:-1] for i in range(A.width - 1)])
    
    # Единственный оставшийся элемент A
    evs.append(A[0][0])
    
    return evs, Q