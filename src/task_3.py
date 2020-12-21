# *-* coding: utf-8 *-*

from Matrix import *
from utility import *


''' Вращения Гивенса '''

def task_3():
    A = input_matrix()
    print('Введите два числа через пробел - i и j: ')
    i, j = map(int, input().split())
    i, j = i - 1, j - 1
    print('Введите два числа через пробел - c и s: ')
    c, s = map(float, input().split())
    x = givens_mul(A, i, j, c, s)
    print(f'Ответ: \n{A}')


# Домножение на матрицу Гивенса слева: G(i, j, c, s) * A
# NB: изменяет переданную матрицу
def givens_mul(A, i, j, c, s):
    n = A.width
    new_i, new_j = [0] * n, [0] * n
    for k in range(n):
        new_i[k] = c * A[i][k] + s * A[j][k]
        new_j[k] = c * A[j][k] - s * A[i][k]
    A[i], A[j] = new_i[:], new_j[:]
    
# Домножение на матрицу Гивенса справа: A * G(i, j, c, s)
# NB: изменяет переданную матрицу
def givens_mul_right(A, i, j, c, s):
    n = A.width
    new_i, new_j = [0] * n, [0] * n
    for k in range(n):
        new_i[k] = c * A[k][i] - s * A[k][j]
        new_j[k] = c * A[k][j] + s * A[k][i]
    for k in range(n):
        A[k][i], A[k][j] = new_i[k], new_j[k]