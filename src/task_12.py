# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_9 import *
from task_11 import *


''' Тест графов на неизоморфность '''

def task_12():
    print('Введите матрицы смежности графов')
    A = input_matrix()
    B = input_matrix()
    res = isomorphism(A, B)
    print(f'Ответ: {res}')


#  Input: матрицы смежности графов A и B
# Output: 0 - если графы не изоморфны, 1 - если это неясно
def isomorphism(A, B):
    # Проверка
    A.square()
    A.symmetric()
    
    # В точности одинаковые графы
    if A == B:
        return True
    
    # Если количество вершин разное, то графы неизоморфны
    if A.height != B.height or A.width != B.width:
        return False
    
    # Считаем степени вершин
    n = A.width
    degA, degB = [0] * n, [0] * n
    for i in range(n):
        degA[i] = sum([A[i][j] for j in range(n)])
        degB[i] = sum([B[i][j] for j in range(n)])
    
    # Если степени разные, то графы неизоморфны
    if sorted(degA) != sorted(degB):
        return False
    
    # Тридиагонализуем, что запустить qr_shift
    triA, QA = get_tridiagonal(A)
    triB, QB = get_tridiagonal(B)
    
    # Находим собственные числа
    evsA, QA = qr_shift(triA)
    evsB, QB = qr_shift(triB)
    evsA.sort()
    evsB.sort()
    
    # Если собственные числа не равны, то графы неизоморфны
    for i in range(n):
        if not eq(evsA[i], evsB[i]):
            return False
    
    # Неясно
    return True