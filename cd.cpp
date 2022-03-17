#include<iostream>
#include<cmath>
using namespace std;

class Vector;

class Matrix
{
private:
   static unsigned int instances;
   static const unsigned int dimension = 3;
   double component[dimension][dimension];
public:
   Matrix(); // standard constructor
   Matrix(const Matrix&); // copy constructor
   ~Matrix();
   static unsigned int get_instances();
   unsigned int get_dimension();
   friend ostream& operator<<(ostream& out, Matrix&); // binary operator << as friend function for a short and elegant output of a matrix
   friend istream& operator>>(istream& in, Matrix&); // binary operator >> as friend function for a short and elegant input of a matrix
   Vector operator*(Vector&); // binary operator * as member operator for multiplication of a matrix with a vector
   friend Matrix operator+(Matrix&, Matrix&);
   friend Matrix operator-(Matrix&, Matrix&);
   friend Matrix operator*(Matrix&, Matrix&);
//   double& operator[](unsigned int i, unsigned int j); // not possible since operator[] is unary
   double& operator()(unsigned int i, unsigned int j);
   // ... possible further matrix operations
};

class Vector
{
   friend class Matrix; // allows class Matrix to access the private members of class Vector
private:
   static unsigned int instances;
   static const unsigned int dimension = 3;
   double component[dimension];
public:
   Vector(); // standard constructor
   Vector(const Vector&); // copy constructor
   ~Vector();
   static unsigned int get_instances();
   unsigned int get_dimension();
   friend istream& operator>>(istream&, Vector&);
   friend ostream& operator<<(ostream&, Vector&);
   double& get_component(unsigned int); // access to component of vector
   double& operator[](unsigned int);
   double& operator()(unsigned int);
   friend Vector operator+(Vector&, Vector&);
   friend Vector operator-(Vector&, Vector&);
   friend double operator*(Vector&, Vector&);
   // ... possible further vector operations
};

unsigned int Matrix::instances = 0;
unsigned int Vector::instances = 0;

Vector::Vector()
{
   instances++;
   for(unsigned int k = 0; k < dimension; k++)
      component[k] = 0.0;
}

Vector::Vector(const Vector &v)
{
   instances++;
   for(unsigned int k = 0; k < dimension; k++)
      component[k] = v.component[k];
}

Vector::~Vector()
{
   instances--;
}

inline unsigned int Vector::get_instances()
{
   return instances;
}

inline unsigned int Vector::get_dimension()
{
   return dimension;
}

inline double& Vector::get_component(unsigned int k)
{
   return component[k];
}

inline double& Vector::operator[](unsigned int k)
{
   return component[k];
}

double& Vector::operator()(unsigned int k)
{
//   return component[k-1]; // without index transformation
   return component[k-1]; // with index transformation
}

Vector operator+(Vector& x, Vector& y)
{
   Vector z; // z = x + y;
   for(unsigned int i = 0; i < x.dimension; i++)
      z.component[i] = x.component[i] + y.component[i];
//      z[i] = x[i] + y[i]; // alternatively uses double& operator[](unsigned int)
   return z;
}

Vector operator-(Vector& x, Vector& y)
{
   Vector z; // z = x - y;
   for(unsigned int i = 0; i < x.dimension; i++)
      z.component[i] = x.component[i] - y.component[i];
//      z[i] = x[i] - y[i]; // alternatively uses double& operator[](unsigned int)
   return z;
}

double operator*(Vector& x, Vector& y)
{
   double d = 0.0; // d = x^T * y
   for(unsigned int i = 0; i < x.dimension; i++)
      d += x.component[i] * y.component[i];
   return d;
}

istream& operator>>(istream& in, Vector& v)
{
   for(unsigned int i = 0; i < v.dimension; i++)
   {
      cout << "vector component[" << i+1 << "]: "; // only useful with keyboard input, generally to omit
      in >> v.component[i];
   }
   return in;
}

ostream& operator<<(ostream& out, Vector& v)
{
   for(unsigned int i = 0; i < v.dimension; i++)
      out << v.component[i] << endl;
   return out;
}

Matrix::Matrix()
{
   instances++;
   for(unsigned int i = 0; i < dimension; i++)
      for(unsigned int j = 0; j < dimension; j++)
         component[i][j] = 0.0;
}

Matrix::Matrix(const Matrix &M)
{
   instances++;
   for(unsigned int i = 0; i < dimension; i++)
      for(unsigned int j = 0; j < dimension; j++)
         component[i][j] = M.component[i][j];
}

Matrix::~Matrix()
{
   instances--;
}

inline unsigned int Matrix::get_instances()
{
   return instances;
}

unsigned int Matrix::get_dimension()
{
   return dimension;
}

ostream& operator<<(ostream& out, Matrix& A)
{
   for(unsigned int i = 0; i < A.dimension; i++)
   {
      for(unsigned int j = 0; j < A.dimension; j++)
         out << '\t' << A.component[i][j];
      out << endl;
   }
   return out;
}

istream& operator>>(istream& in, Matrix& A)
{
   for(unsigned int i = 0; i < A.dimension; i++)
      for(unsigned int j = 0; j < A.dimension; j++)
      {
         cout << "matrix component[" << i+1 << ',' << j+1 << "]: ";
         in >> A.component[i][j];
      }
   return in;
}

Vector Matrix::operator*(Vector& x)
{
   Vector y; // y = Matrix * x
   for(unsigned int i = 0; i < dimension; i++)
      for(unsigned int j = 0; j < dimension; j++)
//         y.component[i] += component[i][j] * x.component[j];
//         y.get_component(i) += component[i][j] * x.get_component(j); // x.get_component(j) is a function call
         y[i] += component[i][j] * x[j]; // alternatively x[j] is an operator[] call
   return y;
}

Matrix operator+(Matrix &X, Matrix &Y)
{
   Matrix Z; // Z = X + Y
   for(unsigned int i = 0; i < X.dimension; i++)
      for(unsigned int j = 0; j < X.dimension; j++)
         Z.component[i][j] = X.component[i][j] + Y.component[i][j];
   return Z;
}

Matrix operator-(Matrix &X, Matrix &Y)
{
   Matrix Z; // Z = X - Y
   for(unsigned int i = 0; i < X.dimension; i++)
      for(unsigned int j = 0; j < X.dimension; j++)
         Z.component[i][j] = X.component[i][j] - Y.component[i][j];
   return Z;
}

Matrix operator*(Matrix &X, Matrix &Y)
{
   Matrix Z; // Z = X * Y
   for (unsigned int i = 0; i < X.dimension; i++)
      for (unsigned int j = 0; j < X.dimension; j++)
         for (unsigned int k = 0; k < X.dimension; k++)
            Z.component[i][j] += X.component[i][k] * Y.component[k][j];
   return Z;
}

double& Matrix::operator()(unsigned int i, unsigned int j)
{
//   return component[i][j]; // without index transformation
   return component[i-1][j-1]; // with index transformation
}

int main(void)
{
   Vector x, y, z;
   Matrix A, B, C1, C2, C3;
   cout << "please input vector x:" << endl;
   cin >> x;
   cout << "please input matrix A:" << endl;
   cin >> A;
   cout << "matrix A is:" << endl << A << endl;
   cout << "vector x is:" << endl << x << endl;
   y = A * x;
   cout << "vector y = A * x is:" << endl << y << endl;
   z = x + y;
   cout << "vector z = x + y is:" << endl << z << endl;
   z = x - y;
   cout << "vector z = x - y is:" << endl << z << endl;
   double s = x * y;
   cout << "scalar product s = x * y = " << s << endl;
   cout << "please input matrix B:" << endl;
   cin >> B;
   cout << "matrix A is:" << endl << A << endl;
   cout << "matrix B is:" << endl << B << endl;
   C1 = A + B;
   cout << "matrix C1 = A + B is:" << endl << C1 << endl;
   C2 = A - B;
   cout << "matrix C2 = A - B is:" << endl << C2 << endl;
   C3 = A * B;
   cout << "matrix C3 = A * B is:" << endl << C3 << endl;

   Vector v = x;
   v(1) = -1.5; // decide yourself whether indices from 1 to n or 0 to n-1 shall be used
   cout << "please input v(2): ";
   cin >> v(2); // decide yourself whether indices from 1 to n or 0 to n-1 shall be used
   v(3) = 2 * v(1) * v(2); // decide yourself whether indices from 1 to n or 0 to n-1 shall be used
   cout << "vector v is:" << endl << v << endl;

   Matrix D = A;
   cout << "please input D(2,3): "; // as example
   cin >> D(2,3); // decide yourself whether indices from 1 to n or 0 to n-1 shall be used
   D(1,1) = -0.5 * D(3, 3); // as example // decide yourself whether indices from 1 to n or 0 to n-1 shall be used
   cout << "matrix D is:" << endl << D << endl;

   cout << "number of matrix objects: " << Matrix::get_instances() << endl;
   cout << "number of vector objects: " << Vector::get_instances() << endl << endl;

   return 0;
}
