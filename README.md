# Класс Matrix: Реализация матричных операций на C++

## Описание проекта
Этот проект представляет собой реализацию класса `Matrix` на C++, который поддерживает основные матричные операции, включая решение систем линейных уравнений методами Гаусса и LU-разложения. Проект разработан для учебных целей и демонстрирует работу с матрицами без использования сторонних библиотек.

## Основные возможности
- Создание матриц различных размеров
- Базовые операции: сложение, вычитание, умножение
- Транспонирование и нахождение обратной матрицы
- Вычисление определителя и миноров (рекурсией)
- Решение СЛАУ методами:
  - Гаусса
  - LU-разложения
- Поддержка пользовательского ввода и вывода матриц

## Структура класса
### Публичные методы
```cpp
// Конструкторы
Matrix();              
Matrix(int);                
Matrix(int, int);
Matrix(int, int, double);

// Геттеры
int getRow();
int getColumn();
double getElement(int, int) const;
double getMinor(int, int);
double det(); 
void print();

// Сеттеры
void setElement(int, int, double);
void set();

// Операции с матрицами
Matrix transpose();
Matrix inverse();

// Решение СЛАУ
Matrix solveGauss(const Matrix& b) const; 
bool LU(Matrix& L, Matrix& U) const;
Matrix solveLU(const Matrix& b) const;

// Перегрузка операторов
const Matrix operator+ (const Matrix&);
const Matrix operator- (const Matrix&);
const Matrix operator* (const Matrix&);
const Matrix operator* (const double&);
