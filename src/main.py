# *-* coding: utf-8 *-*

from graphSpectrum import *
from iterativeMethods import *


if __name__ == '__main__':
    
    ''' Итерационные методы '''
    
    A = Matrix([
        [0.4, 0],
        [0.2, 0.01]
    ])
    E = Matrix.unit(2)
    b = Matrix.vector([0.8, 0.43])
    
    # Решаем A * x == b двумя методами
    x_it = simple_iteration(E - A, b, 1e-14)
    x_gs = Gauss_Seidel(A, b, 1e-14)
    
    # Оператор (==) реализован через сравнение с погрешностью
    assert(A * x_it == b)
    assert(A * x_gs == b)
    
    
    ''' QR-разложение '''
    
    A = Matrix([
        [0, 1, 1],
        [1, 0, 1],
        [0, 1, 0]
    ])
    QG, RG = qr_decomp_givens(A)
    QH, RH = qr_decomp_householder(A)
    
    assert(QG * RG == A)
    assert(QH * RH == A)
    
    
    ''' Простая итерация для нахождения максимального собственного числа '''
    
    # Матрица A из предыдущего пункта
    lm, v = ev_simple_iteration(A, eps = 1e-14)
    
    assert(A * v == lm * v)
    
    
    ''' QR-алгоритм '''
    
    A = Matrix([
        [4,    1,    1,    0,    0],
        [1, -1.1,    1,    2,   -1],
        [1,    1,   -3, 0.02, 0.02],
        [0,    2, 0.02,    0,    5],
        [0,   -1, 0.02,    5,   -1]
    ])
    
    # Список собственных чисел evs, матрица из собственных векторов qA
    evs, qA = qr_algorithm(A)
    
    # Собственные вектора A
    evec = []
    for j in range(qA.width):
        evec.append(Matrix.vector(
            [qA[i][j] for i in range(qA.width)]
        ))
    
    # A * v == lm * v
    for i in range(len(evec)):
        assert(A * evec[i] == evs[i] * evec[i])
    
    
    ''' Трёхдиагональные матрицы '''
    
    # Получаем трёхдиагональную матрицу, запускаем улучшенный алгоритм
    tA, qq = get_tridiagonal(A)
    tevs, tqA = qr_algorithm_tridiag(tA)
    
    # Собственные вектора A
    tevec = []
    for j in range(tqA.width):
        tevec.append(Matrix.vector(
            [tqA[i][j] for i in range(tqA.width)]
        ))
    
    # A * v == lm * v
    for i in range(len(tevec)):
        assert(tA * tevec[i] == tevs[i] * tevec[i])
    
    
    ''' Изоморфизм '''
    
    C = to_matrix(
        '{{1,1,0,0,0,0},\
          {1,0,0,0,0,0},\
          {0,0,0,1,1,0},\
          {0,0,1,0,0,1},\
          {0,0,1,0,0,1},\
          {0,0,0,1,1,0}}'
    )
    D = to_matrix(
        '{{0,1,1,0,0,0},\
          {1,0,0,1,0,0},\
          {1,0,0,0,1,0},\
          {0,1,0,0,0,1},\
          {0,0,1,0,1,0},\
          {0,0,0,1,0,0}}'
    )
    
    assert(isomorphism(C, C))
    assert(not isomorphism(C, D))
    
    
    ''' Экспандеры '''
    
    assert(eq(get_alpha_n(2), 0.5))
    assert(eq(get_alpha_p(2), 0.57735027))
