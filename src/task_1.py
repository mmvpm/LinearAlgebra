# *-* coding: utf-8 *-*

from Matrix import *
from utility import *


''' Метод простой итерации '''

def task_1():
    A = input_matrix()
    b = input_vector()
    eps = input_eps()
    x = simple_iteration(A, b, eps)
    print(f'Ответ: \n{x}')


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