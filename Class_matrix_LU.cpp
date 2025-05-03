#include <iostream>
#include <chrono>
#include <cmath>    // понадобилось для операции (-1)^(i+j)
#define N 10

class Matrix
{
public:
    Matrix();              
    Matrix(int);                
    Matrix(int, int);
    Matrix(int, int, double);

    int getRow();
    int getColumn();
    double getElement(int, int) const;
    double getMinor(int, int);
    double det(); 
    void print();

    void setElement(int, int, double);
    void set();                                   // запускает программу ввода всех элементов 

    Matrix transpose();                           // возвращает транспонированную матрицу
    Matrix inverse();                             // возвращает обратную матрицу

    Matrix solveGauss(const Matrix& b) const; 

    bool LU(Matrix& L, Matrix& U) const;
    Matrix solveLU(const Matrix& b) const;

    const Matrix operator+ (const Matrix&) ;
    const Matrix operator- (const Matrix&) ;
    const Matrix operator* (const Matrix&) ;
    const Matrix operator* (const double&) ; 

private:
    int row;
    int column;
    double data[N][N];

    void initMatrix(int, int, double);
    Matrix ret();                                 // матрица возвращает сама себя
    Matrix sum(double, double, Matrix, Matrix);   // k1*M1 + k2*M2
    Matrix product(Matrix, Matrix);               // M1*M2
    Matrix small(int, int);                       // возвращает матрицу без i-ой строки и j-того столбца
    Matrix dopV();                                // возвращает присоединенную матрицу

    int maxRow(int, int);                         // получает на вход столбец А и строку B.
    //Возвращает номер строки, которая >= B, значение в столбце A у которой по модулю больше остальных
    
    Matrix swapRows(int, int);                    // меняет строки местами
    Matrix splitRow(int, double);                 // делит строку на число
    Matrix minusRow(int, int, double);            // из одной строки вычитает другую, умноженную на k
    Matrix inverseMove();                         // обратный ход в методе Гаусса для упрощенной матрицы

    Matrix forwardSubstitution(const Matrix& b) const;   //Ly=b
    Matrix backwardSubstitution(const Matrix& y) const;  //Ux=y
};

///////////////////////////////////////////////////////////////////////////////////////////////


void Matrix::initMatrix(int row, int column, double value)
{
    this->row = row;
    this->column = column;


    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            data[i][j] = value;
        }
    }
}
Matrix Matrix::ret()
{
    Matrix result(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            result.setElement(i, j, data[i][j]);
        }
    }
    return result;
}
Matrix Matrix::sum(double k1, double k2, Matrix m1, Matrix m2)
{
    double value;
    Matrix result(m1.getRow(), m1.getColumn());
    for (int i = 0; i < m1.getRow(); i++) {
        for (int j = 0; j < m1.getRow(); j++) {
            value = k1 * m1.getElement(i, j) + k2 * m2.getElement(i, j);
            result.setElement(i, j, value);
        }
    }
    return result;
}
Matrix Matrix::product(Matrix m1, Matrix m2)
{
    Matrix result(m1.getRow(), m2.getColumn());
    double count;
    for (int i = 0; i < m1.getRow(); i++) {
        for (int j = 0; j < m2.getColumn(); j++) {
            count = 0;
            for (int i1 = 0; i1 < m1.getColumn(); i1++)
            {
                count += m1.getElement(i, i1) * m2.getElement(i1, j);
            }
            result.setElement(i, j, count);

        }
    }
    return result;
}
Matrix Matrix::small(int a, int b)
{

    Matrix result(getRow() - 1, getColumn() - 1);
    int i1 = 0, j1 = 0;
    for (int i = 0; i < result.getRow(); i++) {
        for (int j = 0; j < result.getColumn(); j++) {
            while (i1 == a || j1 == b) {      
                j1++;
                if (j1 == getRow()) {                         
                    j1 = 0; i1++;
                }
            }
            result.setElement(i, j, getElement(i1, j1));
            j1++;
            if (j1 == getRow()) {
                j1 = 0; i1++;
            }
            

           
        }
    }
    return result;
}
Matrix Matrix::dopV()
{
    double n;
    Matrix result(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            n = getMinor(i, j) * pow(-1, (i + j));
            result.setElement(i, j, n);
        }
    }
    return result;
}

int Matrix::maxRow(int a, int i0)
{
    Matrix result(row, column);
    double maxvalue = getElement(i0, a) * getElement(i0, a);
    int maxnum = i0;
    for (int i = i0; i < row; i++) {
        if (getElement(i, a) * getElement(i, a) > maxvalue) {
            maxvalue = getElement(i, a) * getElement(i, a);
            maxnum = i;
        }
    }
    return maxnum;
}
Matrix Matrix::swapRows(int a , int b) 
{
    Matrix result(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            if (i == a) result.setElement(i, j, getElement(b, j));
            if (i == b) result.setElement(i, j, getElement(a, j));
            if (i != a && i != b) result.setElement(i, j, getElement(i, j));
        }
    }
    return result;

}
Matrix Matrix::splitRow(int a, double c)
{
    Matrix result(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            if (i == a) result.setElement(i, j, (getElement(i, j)/c) );
            if (i != a) result.setElement(i, j, getElement(i, j));
        }
    }
    return result;
}
Matrix Matrix::minusRow(int a, int b, double c)
{
    Matrix result(row, column);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            if (i == a) result.setElement(i, j, (getElement(i, j) - getElement(b, j) * c) );
            if (i != a) result.setElement(i, j, getElement(i, j));
        }
    }
    return result;
}
Matrix Matrix::inverseMove()
{
    Matrix result(row, 1, 0);
    result.setElement(row - 1, 0, getElement(row - 1, column - 1));

    for (int i = row - 2; i >= 0; i--) {
        double x = getElement(i, column - 1); //сначала делал int, потом долго не мог найти почему на 1 меньше результаты некоторые
        double count = 0;
        for (int j = 0; j < row; j++) {
            count += (result.getElement(j, 0)) * (getElement(i, j));
        }
        x = x - count;
        result.setElement(i, 0, x);

    }
    return result;

}

Matrix Matrix::forwardSubstitution(const Matrix& b) const {
    Matrix y(row, 1);
    for (int i = 0; i < row; i++) {
        y.data[i][0] = b.data[i][0];
        for (int j = 0; j < i; j++) {
            y.data[i][0] -= data[i][j] * y.data[j][0];
        }
        y.data[i][0] /= data[i][i];
    }
    return y;
}
Matrix Matrix::backwardSubstitution(const Matrix& y) const {
    Matrix x(row, 1);
    for (int i = row - 1; i >= 0; i--) {
        x.data[i][0] = y.data[i][0];
        for (int j = i + 1; j < row; j++) {
            x.data[i][0] -= data[i][j] * x.data[j][0];
        }
        x.data[i][0] /= data[i][i];
    }
    return x;
}

///////////////////////////////////////////////////////////////////////////////////////////////

Matrix::Matrix()
{
    initMatrix(3, 3, 0);
}
Matrix::Matrix(int n)
{
    initMatrix(n, n, 0);
}
Matrix::Matrix(int n, int m)
{
    initMatrix(n, m, 0);
}
Matrix::Matrix(int row, int column, double c)
{
    initMatrix(row, column, c);
}

int Matrix::getRow()
{
    return row;
}
int Matrix::getColumn()
{
    return column;
}
double Matrix::getElement(int a, int b) const
{
    return data[a][b];
}
double Matrix::getMinor(int a, int b) 
{
    return small(a, b).det();
}
double Matrix::det()  //реализовано разложением по первой строке с помощью рекурсии
{
    if (getRow() == 1 && getColumn() == 1) return getElement(0, 0);
    if (getRow()==2 && getColumn()==2) {
        double c = getElement(0, 0) * getElement(1, 1) - getElement(0, 1) * getElement(1, 0);
        return c;
    }
    else {
        double count = 0;
        for (int j = 0; j < column; j++) {
            count += getElement(0, j) * pow(-1, 0 + j) * getMinor(0, j);
        }
        return count;
    }
}      
void Matrix::print()
{
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            std::cout << data[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void Matrix::setElement(int a, int b, double val)
{
    data[a][b] = val;
}
void Matrix::set()
{
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            std::cin >> data[i][j];
        }
    }
}

Matrix Matrix::transpose()
{
    Matrix result(column, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            result.setElement(j, i, data[i][j]);
        }
    }
    return result;
}
Matrix Matrix::inverse() {

    Matrix result(column, row);
    result = dopV();
    result = result.transpose();
    return sum((1 / det()), 0, result, ret());
}

const Matrix Matrix::operator+ (const Matrix& matrix)
{
    return sum(1, 1, ret(), matrix);
}
const Matrix Matrix::operator- (const Matrix& matrix) 
{
    return sum(1, -1, ret(), matrix);
}
const Matrix Matrix::operator* (const Matrix& matrix) 
{
    return product(ret(), matrix);
}
const Matrix Matrix::operator* (const double& a) 
{
    return sum(a, 0, ret(), ret());
}

Matrix Matrix::solveGauss(const Matrix& b) const {

    Matrix Ab(row, column + 1);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            Ab.setElement(i, j, data[i][j]);
        }
        Ab.setElement(i, column, b.data[i][0]);
    }

    for (int i = 0; i < row; i++) {
        
        int maxRow = Ab.maxRow(i, i);
        if (maxRow != i) {
            Ab = Ab.swapRows(i, maxRow);
        }

        double pivot = Ab.getElement(i, i);

        Ab = Ab.splitRow(i, pivot);

        for (int k = i + 1; k < row; k++) {
            double factor = Ab.getElement(k, i);
            Ab = Ab.minusRow(k, i, factor);
        }
    }

    Matrix x(row, 1);
    for (int i = row - 1; i >= 0; i--) {
        x.setElement(i, 0, Ab.getElement(i, column));
        for (int j = i + 1; j < row; j++) {
            x.setElement(i, 0, x.getElement(i, 0) - Ab.getElement(i, j) * x.getElement(j, 0));
        }
    }

    return x;
}
bool Matrix::LU(Matrix& L, Matrix& U) const {
    if (row != column) return false;

    L = Matrix(row, row, 0);
    U = Matrix(row, row, 0);

    for (int i = 0; i < row; i++) {

        for (int k = i; k < row; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += L.data[i][j] * U.data[j][k];
                U.data[i][k] = data[i][k] - sum;
        }

        for (int k = i; k < row; k++) {
            if (i == k)
                L.data[i][i] = 1;
            else {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L.data[k][j] * U.data[j][i];

                if (U.data[i][i] == 0) return false;
                L.data[k][i] = (data[k][i] - sum) / U.data[i][i];
            }
        }
    }
    return true;
}
Matrix Matrix::solveLU(const Matrix& b) const {
    Matrix L, U;
    if (!this->LU(L, U)) {
        std::cout << "LU decomposition failed!" << std::endl;
        return Matrix();
    }

    Matrix y = L.forwardSubstitution(b);
    return U.backwardSubstitution(y);
}

///////////////////////////////////////////////////////////////////////////////////////////

void testMethods() {

    Matrix A(3);
    A.setElement(0, 0, 2); A.setElement(0, 1, 1); A.setElement(0, 2, -1);
    A.setElement(1, 0, -3); A.setElement(1, 1, -1); A.setElement(1, 2, 2);
    A.setElement(2, 0, -2); A.setElement(2, 1, 1); A.setElement(2, 2, 2);

    Matrix b1(3, 1);
    b1.setElement(0, 0, 8);
    b1.setElement(1, 0, -11);
    b1.setElement(2, 0, -3);

    Matrix b2(3, 1);
    b2.setElement(0, 0, 5);
    b2.setElement(1, 0, -7);
    b2.setElement(2, 0, -1);

    A.print(); b1.print(); b2.print();

    auto start = std::chrono::high_resolution_clock::now();
    Matrix x1_gauss = A.solveGauss(b1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_gauss1 = end - start;
    
    start = std::chrono::high_resolution_clock::now();
    Matrix x1_lu = A.solveLU(b1);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_lu1 = end - start;

    start = std::chrono::high_resolution_clock::now();
    Matrix x2_gauss = A.solveGauss(b2);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_gauss2 = end - start;
    
    start = std::chrono::high_resolution_clock::now();
    Matrix x2_lu = A.solveLU(b2);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_lu2 = end - start;

    
    std::cout << "b1 (Gauss)\n";
    x1_gauss.print();
    std::cout << "Time: " << elapsed_gauss1.count() << " s\n\n";

    std::cout << "b1 (LU):\n";
    x1_lu.print();
    std::cout << "Time: " << elapsed_lu1.count() << " s\n\n";

    std::cout << "b2 (Gauss):\n";
    x2_gauss.print();
    std::cout << "Time: " << elapsed_gauss2.count() << " s\n\n";

    std::cout << "b2 (LU):\n";
    x2_lu.print();
    std::cout << "Time: " << elapsed_lu2.count() << " s\n\n";

}

int main() {
    testMethods();
    return 0;
}
