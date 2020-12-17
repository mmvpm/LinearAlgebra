# *-* coding: utf-8 *-*   

from Matrix import *


''' Метод простой итерации '''

#  Input: квадратная матрица A, вектор b, точность eps
# Output: вектор x такой, что ||x - Ax - b|| < eps
def simple_iteration(A, b, eps = global_eps):
    # Проверка
    A.square()
    
    # Случайный вектор
    x = Matrix.random(A.width, 1)
    
    # Наиболее удалённая от 0 точка в кругах Гершгорина
    max_ev_estimation = 0
    for (center, radius) in A.gershgorin_circles():
        max_ev_estimation = max(max_ev_estimation, abs(center) + radius)
        
    # Лежат ли круги Гершгорина вне единичного круга с центром в нуле
    outside_unit = max_ev_estimation > 1
    
    # Количество подряд идущих неудачных итераций
    count_fail = 0
    
    while (x - A * x - b).norm() >= eps:
        x_new = A * x + b
        
        # Удачная ли итерация
        if abs(x_new.norm() - x.norm()) >= 1:
            count_fail += 1
        else:
            count_fail = 0
        
        x = x_new
        
        # Последние 20 итераций были неудачными
        if count_fail > 20 and outside_unit:
            return 0
    
    return x


''' Метод Гаусса-Зейделя '''

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
def Gauss_Seidel(A, b, eps = global_eps):
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


''' Метод простой итерации для нахождения максимального собственного числа '''

#  Input: квадратная матрица A, 
#         начальное приближение к собственному вектору v,
#         точность eps
# Ouput: lm – приближение к максимальному собственному числу A,
#        вектор v такой, что ||v|| = 1 и ||Av − lm v|| < eps
def ev_simple_iteration(A, v = -1, eps = global_eps):
    # Проверка
    A.square()
    
    # Если начальное приближение не задано, берём случайный вектор
    if v == -1:
        v = Matrix.random(A.width, 1).normalize()
    
    for _ in range(100):
        v = (A * v).normalize() # v = Av / ||Av||
        lm = (v.transpose() * A * v)[0][0] # v^T * A * v
        
        # Верно ли, что ||Av − lm v|| < eps
        if eq((A * v - lm * v).norm(), 0, eps):
            return lm, v
    
    # За 100 итераций не нашли ответ
    return 0