# *-* coding: utf-8 *-*

from Matrix import *
from utility import *


''' Отражения Хаусхолдера '''

def task_5():
    A = input_matrix()
    v = input_vector()
    A = householder_mul(A, v)
    print(f'Ответ: \n{A}')


# Домножение на матрицу Хаусхолдера слева: (I - 2vv^T) * A
# NB: не изменяет переданную матрицу
def householder_mul(A, v):
    return A - (2 * v) * (v.transpose() * A)

# Домножение на матрицу Хаусхолдера справа: A * (I - 2vv^T)
# NB: не изменяет переданную матрицу
def householder_mul_right(A, v):
    return A - (A * (2 * v)) * v.transpose()