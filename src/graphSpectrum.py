# *-* coding: utf-8 *-*

from Matrix import *
from qrAlgorithm import *


''' Тест графов на неизоморфность '''

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


''' Экспандеры '''

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