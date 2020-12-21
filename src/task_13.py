# *-* coding: utf-8 *-*

from Matrix import *
from utility import *
from task_9 import *
from task_11 import *


''' Экспандеры '''

def task_13():
    print('Введите натуральное число n: ')
    n = int(input())
    alpha_n = get_alpha_n(n)
    print(f'Ответ: {alpha_n}')
    
    print('Введите простое число p:')
    p = int(input())
    alpha_p = get_alpha_p(p)
    print(f'Ответ: {alpha_p}')


#  Input: натуральное число n, точность eps
# Output: alpha = max(|lm_2|, |lm_n|) / deg
def get_alpha_n(n, eps = 1e-5):
    
    # Отображает пару в число
    def num(i, j):
        i, j = i % n, j % n
        return i * n + j
    
    # Граф, который будем строить
    A = Matrix.zero(n ** 2, n ** 2)
    
    for x in range(n):
        for y in range(n):
            # Соседи вершины (x, y)
            neighbours = [
                (    x + 2 * y, y),
                (    x - 2 * y, y),
                (x + 2 * y + 1, y),
                (x - 2 * y - 1, y),
                (x,     y + 2 * x),
                (x,     y - 2 * x),
                (x, y + 2 * x + 1),
                (x, y - 2 * x - 1)
            ]
            for nx, ny in neighbours:
                v, u = num(x, y), num(nx, ny)
                A[v][u] += 1
    
    # Находим собственные числа
    tA, qt = get_tridiagonal(A)
    lms, qA = qr_shift(tA, eps)
    lms.sort(reverse =  True)
    
    # Степень построенного графа
    deg = sum(A[0])
    
    return max(abs(lms[1]), abs(lms[-1])) / deg


#  Input: простое число p, точность eps
# Output: alpha = max(|lm_2|, |lm_n|) / deg
def get_alpha_p(p, eps = 1e-5):
    
    # Вершина p соответствует значению inf
    
    # (x + 1) (mod p)
    def inc(x):
        if x == p:
            return p # inf
        return (x + 1) % p
    
    # (x - 1) (mod p)
    def dec(x):
        if x == p:
            return p # inf
        return (x - 1) % p
    
    # a ** m (mod p)
    def fpow(a, m):
        res = 1
        while m > 0:
            if m & 1:
                res = res * a % p
            a = a * a % p
            m >>= 1
        return res

    # x ** (-1) == x ** (p - 2) (mod p)
    def rev(x):
        if x == 0:
            return p # inf
        if x == p:
            return 0
        return fpow(x, p - 2)

    # Граф, который будем строить
    A = Matrix.zero(p + 1, p + 1)
    
    for x in range(p + 1):
        # Соседи вершины x
        neighbours = [inc(x), dec(x), rev(x)]
        for nx in neighbours:
            A[x][nx] += 1
    
    # Находим собственные числа
    tA, qt = get_tridiagonal(A)
    lms, qA = qr_shift(tA, eps)
    lms.sort(reverse =  True)
    
    # Степень построенного графа
    deg = sum(A[0])
    
    return max(abs(lms[1]), abs(lms[-1])) / deg