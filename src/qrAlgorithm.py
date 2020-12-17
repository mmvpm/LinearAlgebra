# *-* coding: utf-8 *-*

from Matrix import *
from qrDecomp import *


''' QR-алгоритм '''

# Верно ли, что радиусы кругов Гершгорина у A меньше eps
def __check(A, eps = global_eps):
    for cent, rad in A.gershgorin_circles():
        if rad > eps:
            return False
    return True


#  Input: квадратная симметричная матрица A, точность eps
# Output: список из приближённых собственных чисел evs, 
#         матрица примерно из собственных векторов Q
def qr_algorithm(A, eps = global_eps):
    # Проверка
    A.square()
    A.symmetric()
    
    # Матрица из собственных векторов
    Q = Matrix.unit(A.width)
    
    while not __check(A, eps):
        # Вместо создания R_k сразу изменяем A
        Q_k, A = qr_decomp_givens(A) # Q_k, R_k = qr_decomp(A)
        A *= Q_k # A = R_k * Q_k
        Q *= Q_k # Q = Q * Q_k
    
    # Собственные числа стоят на диагонали матрицы A
    evs = [A[i][i] for i in range(A.width)]
    
    return evs, Q


''' Трёхдиагональные матрицы '''

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
        ]).normalize() # поделили на норму
        
        # Стандартный базисный вектор e_{j + 1}
        e = Matrix.vector([int(i == j + 1) for i in range(n)])
        
        # Вектор, вокруг которого будем отражать
        v = (u - e).normalize()
        v = Matrix.vector([
            0 if i < j + 1 else v[i][0] for i in range(n)
        ]).normalize()
        
        TD = householder_mul_right(TD, v) # TD = H_v * TD
        TD = householder_mul(TD, v)       # TD = TD * H_v
        Q  = householder_mul(Q , v)       # Q  = Q  * H_v
    
    Q = Q.transpose() # Q = Q^T
    
    # Зануляем те числа, которые из-за ошибок округления не равны нулю
    TD = TD.map(lambda x: 0 if eq(x, 0) else x)
    
    return TD, Q

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
    
    while not __check(A, eps):
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


''' Сдвиг '''

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
        D = (cf[0] + cf[3]) ** 2 - 4 * (cf[0] * cf[3] - cf[1] * cf[2])
        x = (-cf[1] + D ** 0.5) / 2, (-cf[1] - D ** 0.5) / 2
        
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