# *-* coding: utf-8 *-*

from random import randint, random


# Точность вычислений по умолчанию
global_eps = 1e-9

# Сравнение двух double на равенство
def eq(a, b, eps = global_eps):
    return abs(a - b) < eps


class Matrix:

    # Нулевая матрица размером n на m
    @staticmethod
    def zero(n, m):
        return Matrix([[0] * m for _ in range(n)])
    
    # Диагональная матрица с элементами списка diag на диагонали
    @staticmethod
    def diagonal(diag):
        
        def get_row(i):
            return [0] * i + [diag[i]] + [0] * (n - i - 1)
        
        n = len(diag)
        return Matrix([get_row(i) for i in range(n)])

    # Единичная матрица размером n на n
    @staticmethod
    def unit(n):
        return Matrix.diagonal([1] * n)
    
    # Вектор из списка a (матрица размером len(a) на 1)
    @staticmethod
    def vector(a):
        return Matrix([[i] for i in a])
    
    # Случаная матрица размером n на m
    # Числа из полуинтервала [0, 1)
    # precision - количество знаков после запятой
    @staticmethod
    def random(n, m, precision = 16):
        result = Matrix.zero(n, m)
        for i in range(n):
            for j in range(m):
                result[i][j] = random()
        return result.prec(precision)
    
    # Случайная целочисленная матрица размером n на m
    # Числа из отрезка [first, last]
    @staticmethod
    def randint(n, m, first = 0, last = 1):
        result = Matrix.zero(n, m)
        for i in range(n):
            for j in range(m):
                result[i][j] = randint(first, last)
        return result
    
    
    def __init__(self, data):
        if len(data) == 0:
            raise ValueError('Matrix must be non-empty')
        if len(data[0]) == 0:
            raise ValueError('Matrix rows must be non-empty')
        
        # Высота и ширина матрицы
        self.height = len(data)
        self.width = len(data[0])
        
        if [len(line) for line in data].count(self.width) != self.height:
            raise ValueError('Matrix must be rectangular')
        
        # Двумерный массив - сама матрица
        self.data = data
    
    
    # Применяет функцию fun ко всем элементам матрицы
    def map(self, fun):
        result = Matrix.zero(self.height, self.width)
        for i in range(self.height):
            result[i] = list(map(fun, self[i]))
        return result
    
    # Копия текущей матрицы
    def copy(self):
        return self.map(lambda x: x)
    
    # Копия текущей матрицы, с округлёнными до p знаков после запятой числами
    def prec(self, p = 3):
        return self.map(lambda x: int(x * 10 ** p) / 10 ** p)
    
    
    # Строковое представление
    def __str__(self):
        
        def shift(x):
            return ' ' * (offset - len(x)) + x
        
        str_data = self.map(str)
        # Выравнивание элементов
        offset = max([max(map(len, line)) for line in str_data])
        str_data = [' '.join(map(shift, line)) for line in str_data]
        return '\n'.join(str_data)
    
    def __repr__(self):
        return str(self)
    
    
    # self[key]
    def __getitem__(self, key):
        return self.data[key]
    
    # self[key] = value
    def __setitem__(self, key, value):
        if len(value) != self.width:
            raise ValueError('Matrix must be rectangular')
        
        self.data[key] = value
    
    # self == other (с точностью до global_eps)
    def __eq__(self, other):
        if self.height != other.height or self.width != other.width:
            return False
        for i in range(self.height):
            for j in range(self.width):
                if not eq(self[i][j], other[i][j]):
                    return False
        return True
    
    
    # self * n, где n - число
    def __rmul__(self, n):
        result = Matrix.zero(self.height, self.width)
        for i in range(self.height):
            for j in range(self.width):
                result[i][j] = self[i][j] * n
        return result    
    
    # -self
    def __neg__(self):
        return -1 * self
    
    # self + other
    def __add__(self, other):
        if self.height != other.height or self.width != other.width:
            raise ValueError('Dimensions of the matrices must be equal')
        
        result = Matrix.zero(self.height, self.width)
        for i in range(self.height):
            for j in range(self.width):
                result[i][j] = self[i][j] + other[i][j]
        return result
    
    # self - other
    def __sub__(self, other):
        return self + -other    
    
    # self += other
    def __iadd__(self, other):
        self = self + other
        return self
    
    # self -= other
    def __isub__(self, other):
        self = self - other
        return self
    
    # self * other
    def __mul__(self, other):
        if self.width != other.height:
            raise ValueError('Dimensions of the matrices must be suitable')
        
        result = Matrix.zero(self.height, other.width)
        for i in range(result.height):
            for j in range(result.width):
                for k in range(self.width):
                    result[i][j] += self[i][k] * other[k][j]
        return result
    
    # self *= other
    def __imul__(self, other):
        self = self * other
        return self
    
    # self ** n
    def __pow__(self, n):
        if type(n) != int:
            raise TypeError('Argument must be of the type int')
        if n < 0:
            raise ValueError('Negative degree has not implemented yet')
        
        result = Matrix.unit(self.height)
        while n > 0:
            if n % 2 == 1:
                result *= self
            self *= self
            n //= 2
        return result
    
    # Евклидова норма
    def norm(self):
        result = 0
        for i in range(self.height):
            for j in range(self.width):
                result += abs(self[i][j]) ** 2
        return result ** 0.5
    
    # Нормализация, приведение к единичной норме
    def normalize(self):
        norm = self.norm()
        if eq(norm, 0):
            return self
        return (1 / norm) * self
    
    # Транспонирование
    def transpose(self):
        result = Matrix.zero(self.width, self.height)
        for i in range(self.height):
            for j in range(self.width):
                result[j][i] = self[i][j]
        return result
    
    # Круги Гершгорина
    def gershgorin_circles(self):
        if self.height != self.width:
            raise ValueError('Gershgorin circles are defined only for square matrices')
        result = [(0, 0)] * self.height
        for i in range(self.height):
            result[i] = (self[i][i], sum(map(abs, self[i])) - abs(self[i][i]))
        return result
    
        
    # Является ли матрица квадратной
    def square(self):
        if self.height != self.width:
            raise ValueError('Matrix must be square')
    
    # Является ли матрица симметричной
    def symmetric(self):
        if self != self.transpose():
            raise ValueError('Matrix must be symmetric')
    
    # Является ли матрица трёхдиагональной
    def tridiagonal(self):
        fail = False
        for i in range(self.height):
            for j in range(self.width):
                if abs(i - j) > 1 and not eq(self[i][j], 0):
                    fail = True
                    break
        if fail:
            raise ValueError('Matrix must be tridiagonal')


# Матрица из строки вида s = '{{..}, .., {..}}'
def to_matrix(s):
    
    def remove_all(s, trash):
        return ''.join([c for c in s if c not in trash])
    
    result = []
    s = ''.join(s.split()).split(',')
    for num in s:
        if '{' in num or '[' in num:
            result.append([])
        num = remove_all(num, '{}[]')
        result[-1].append(float(num))
    
    return Matrix(result)