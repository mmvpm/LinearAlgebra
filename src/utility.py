# *-* coding: utf-8 *-*

from Matrix import *

def input_matrix():
    print('Введите два числа через пробел - размер матрицы: ')
    n, m = map(int, input().split())
    
    print('Введите матрицу построчно, разделяя числа пробелами: ')
    result = [0] * n
    for i in range(n):
        result[i] = list(map(float, input().split()))
    
    return Matrix(result)

def input_vector():
    print('Введите одно число - размер вектора: ')
    n = int(input())
    
    print('Введите вектор (по одному числу в строке): ')
    result = [0] * n
    for i in range(n):
        result[i] = float(input())
    
    return Matrix.vector(result)

def input_eps():
    print('Введите точность: ')
    eps = float(input())
    return eps