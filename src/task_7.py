# *-* coding: utf-8 *-*

from Matrix import *
from utility import *


''' Метод простой итерации для нахождения максимального собственного числа '''

def task_7():
    A = input_matrix()
    x_0 = input_vector()
    eps = input_eps()
    x = ev_simple_iteration(A, x_0, eps)
    print('Ответ: ')
    if x == 0:
        print(0)
    else:
        ev, v = x
        print(f'lambda = \n{ev}\nv = \n{v}')


#  Input: квадратная матрица A, 
#         начальное приближение к собственному вектору v,
#         точность eps
# Ouput: lm – приближение к максимальному собственному числу A,
#        вектор v такой, что ||v|| = 1 и ||Av − lm v|| < eps
def ev_simple_iteration(A, v = -1, eps = global_eps):
    # Проверка
    A.square()
    
    # Если начальное приближение не задано, берём случайный вектор
    if type(v) != Matrix and v == -1:
        v = Matrix.random(A.width, 1).normalize()
    
    for _ in range(100):
        v = (A * v).normalize() # v = Av / ||Av||
        lm = (v.transpose() * A * v)[0][0] # v^T * A * v
        
        # Верно ли, что ||Av − lm v|| < eps
        if eq((A * v - lm * v).norm(), 0, eps):
            return lm, v
    
    # За 100 итераций не нашли ответ
    return 0