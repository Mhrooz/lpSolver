#include <iostream>
#include <exception>
#include <utility>
#include <vector>
using namespace std;

typedef long double Element;
const bool debug = false;

class Matrix{
public:
    // conductors
    // default matrix:
    // 1,0
    // 0,1
    Matrix():row(2),column(2){
        Element e[4]= {1,0,0,1};
        element = new Element[4];
        element = e;
    }
    // generate zero matrix with dimension r*c;
    Matrix(int r,int c):row(r),column(c){
        element = new Element[row*column];
    }
    // generate square matrix
    Matrix(int r):row(r),column(r){
        element = new Element[r*r];
    }
    // generate matrix with existing element
    Matrix(int r,int c, Element * e):row(r),column(c){
        element = new Element[row*column];
        for(int i = 0 ; i < r*c ;i ++){
            element[i] = e[i];
            if(debug)
                cout<<element[i]<<' ';
        }
        if(debug)
            cout<<endl;
    };
    // generate matrix
    Matrix(int r,int c,Element x):row(r),column(c){
        element = new Element[row*column];
        for(int i = 0 ; i < r*c;i++){
            element[i]=x;
        }
    }
    //merge A and B to [A,B]
    Matrix(Matrix A,Matrix B){
        try {
            if (A.row != B.row) {
                throw (sizeMismatch(A, B, "Merge Matrix"));
            }
        }catch(sizeMismatch &e){
            e.what();
        }
        row = A.getRow();
        column = A.column+B.column;
        element = new Element (row*column);
        for(int i = 0;i<row;i++)
            for(int j=0;j<column;j++){
                if(j<A.column){
                    (*this)(i,j)=A(i,j);
                }
                else{
                    (*this)(i,j)=B(i,j-A.column);
                }
            }
    }
    /*
     * The full matrix after insert function:
     * [M1]
     * [M2]
     */
    Matrix insert(Matrix &M1,Matrix &M2){
        try {
            if (M1.column!= M2.column) {
                throw (sizeMismatch(M1, M2, "Insert Matrix"));
            }
        }catch(sizeMismatch &e){
            e.what();
        }
        auto rlt = new Matrix(M1.getRow()+M2.getRow(),M1.getColumn());
        for(int i = 0 ; i < M1.getRow();i++){
            for(int j = 0 ; j < M1.getColumn();j++){
                (*rlt)(i,j) = M1(i,j);
            }
        }
        for(int i = 0 ; i < M2.getRow();i++){
            for(int j = 0 ; j < M1.getColumn();j++){
                (*rlt)(i+M1.getRow(),j) = M2(i,j);
            }
        }
        return *rlt;
    }
    //overloaded operator + for Matrix;
    friend Matrix operator +(const Matrix & M1,const Matrix & M){
        try{
            if(M1.column!=M.column||M1.row!=M.row)
                throw(sizeMismatch( M1, M));
        }
        catch(sizeMismatch &e){
            cout<<e.what()<<endl;
        }
        Matrix rlt = Matrix(M.getRow(),M.getColumn());

        for(int i = 0 ; i < M.row;i++)
            for(int j = 0 ; j < M.column;j++){
                rlt(i,j) = M1(i,j)+M(i,j);
            }
        return rlt;
    }

    Matrix operator + (int x) const{
        auto * rlt = new Matrix(this->column,this->row,x);
        return *rlt+*this;
    }

    friend Matrix operator * (const Matrix & M1,const Matrix & M2){
        try{
            if(M1.column != M2.row)
                throw(sizeMismatch(M1,M2,"This is Matrix multiplication"));
        }catch(sizeMismatch &e){
            e.what();
        }
        auto rlt = new Matrix(M1.row,M2.column);
        for(int i = 0 ; i < M1.row ; i++)
            for(int j = 0 ; j < M2.row ; j++){
                for(int k = 0 ; k < M1.column ;k++){
                    (*rlt)(i,j) += M1(i,k)*M2(k,j);
                }
            }
        return *rlt;
    }

    friend Matrix operator * (const Matrix& M1, const int &x) {
        auto rlt = new Matrix(M1);
        for(int i = 0 ; i < M1.row;i++){
            for(int j = 0 ;j < M1.column ;j ++){
                (*rlt)(i,j)=-1*M1(i,j);
            }
        }
        return *rlt;
    }

    friend Matrix operator * (const int &x , const Matrix& M1){
        return M1*x;
    }

    friend Matrix operator -(const Matrix & M1,const Matrix & M2) {
        auto rlt = -1 * M2;
        return M1+rlt;
    }
    //overloaded operator [] for get the number of i(i/column = x,i%column = y)
    Element & operator [](int i)const{
        try{
            if(i > column * row){
                throw(overSize(*this,column*row));
            }
        }catch(sizeMismatch &e){
            e.what();
        }
        return element[i];
    }

    //overloaded operator () for get the number of x,y
    Element & operator () (int x,int y)const {
        try{
            if(x>row||y>column)
                throw(overSize(*this,row*column));
        }
        catch(overSize &e){
            e.what();
        }
        return (*this)[x*column+y];
    }

    int getRow() const{return this->row;}

    int getColumn() const{return this->column;}

    Element * getElement() const{return this->element;}

    class overSize :public exception{ ;

    public:
        overSize(Matrix a,int id):row(a.getRow()),col(a.getColumn()),i(id){}

        const char * what() const noexcept override{
            return "the size is not correct";
        }

    private:
        int row{};int col{};int i{};
    };

    class sizeMismatch:public exception{
    public:
        sizeMismatch(Matrix a, Matrix b):aRow(a.getRow()),bRow(b.getRow()),aCol(a.getColumn()),bCol(b.getColumn()){};
        sizeMismatch(Matrix a, Matrix b,string ss):aRow(a.getRow()),bRow(b.getRow()),aCol(a.getColumn()),bCol(b.getColumn()),s(std::move(ss)){};
        const char * what()const noexcept override{
            if(!s.empty())cout<<s<<endl;
            cout<<"The scale of left Matrix"<<aRow<<'x'<<aCol;
            cout<<"The scale of right Matrix"<<bRow<<'x'<<bCol;
            return "The scale of two matrix does not match";
        }
    private:
        string s; int aRow{}; int bRow{}; int aCol{}; int bCol{};
    };

    void matrixPrint(){
        for(int i = 0 ; i <row;i++) {
            for (int j = 0; j < column; j++) {
                cout << (*this)(i,j) << ' ';
            }
            cout << endl;
        }
    }

private:
    int row,column;// the scale of Matrix;
    Element * element;
};

//Transfer the matrix into loose matrix:
Matrix getLooseMatrix(Matrix M){
    int row = M.getRow();
    int column = M.getColumn();
    auto looseMatrix = new Matrix(row,row+column);
    for(int i = 0 ; i < row ;i ++){
        for(int j = 0 ; j < column;j++)
            (*looseMatrix)(i,j) = M(i,j);
        (*looseMatrix)(i,i+column) = 1;
    }
    return *looseMatrix;
}

//Join three Matrix into one
Matrix joinMatrix(Matrix coefficientMatrix,Matrix constrainsMatrix,Matrix minMatrix){
    int row = coefficientMatrix.getRow();
    int column = coefficientMatrix.getColumn();
    Matrix *rlt = new Matrix(row+1,column+1);
    for(int i = 0 ; i < row ; i++)
        for(int j = 0 ; j < column ; j++){
            (*rlt)(i+1,j+1) = coefficientMatrix(i,j);//coefficientMatrix at right downside
        }
    for(int i = 0 ; i < row ; i++)
        (*rlt)(i+1,0) = constrainsMatrix(0,i);//constrainsMatrix at left side
    for(int i = 1 ; i <= minMatrix.getColumn();i++){
        (*rlt)(0,i) = minMatrix(0,i-1);//min Matrix at upper side
    }
    return *rlt;
}
//pivot the matrix
Matrix pivot_matrix(Matrix &matrix, int k, int j){
//    cout<<"----------"<<endl;
//    matrix.matrixPrint();
//    cout<<"----------"<<endl;
    double tmp = matrix(k,j);
    for(int i = 0; i < matrix.getColumn();i++) {
        matrix(k, i) = matrix(k, i) / tmp;
    }
    for(int i = 0 ; i < matrix.getRow();i++) {
        if(i!=k) {
            Element tmp = matrix(i,j);
            for (int m = 0; m < matrix.getColumn(); m++) {
                matrix(i, m) = matrix(i, m) - matrix(k, m) * tmp;
            }
        }
    }
    return matrix;
}
bool judge(Matrix rlt){
    for(int i = 1 ; i <rlt.getColumn();i++)
        if(rlt(0,i)<0)
            return true;
   return false;
}
Matrix simplex(Matrix matrix,vector<int>&base_ids){

    Matrix rlt = Matrix(matrix);
    while(judge(rlt)){
        int j;
        for(j=1;j<rlt.getRow();j++) {
            if (rlt(0, j) < 0)
                break;
        }
        Element minn = 0x7fff;
        int k ;//k is the constriction
        for(int i = 1 ; i<rlt.getRow();i++){
            double x;
            if(rlt(i,j)>0){
                x = rlt(i,0)/rlt(i,j);
            }
            else
                x = 0x7fff;
            if(x<minn) {
                k = i; minn = x;
            }
        }
        if(rlt(k,j)<0){
            cout<<"the problem has no limited"<<endl;
            Element *tmp = nullptr;
            return Matrix(0,0,tmp);
        }
        rlt = pivot_matrix(rlt,k,j);
        base_ids[k-1]=j-1;
    }
    return rlt;
}
Matrix Laux(Matrix matrix, vector<int>&base_ids){
    auto rlt = Matrix(matrix);
    //Matrix(Matrix A,Matrix B) will merge two matrices as [A,B]
    rlt = Matrix(rlt,Matrix(rlt.getRow(),1,-1));// Matrix(rlt.getRow(),1,-1) will generate a matrix of assistant element x which use -1
    for(int i = 0 ; i<rlt.getColumn()-1;i++)rlt(0,i)=0;
    rlt(0,rlt.getColumn()-1)=1;
    int k ; //k is substitution out variable
    int j = rlt.getColumn()-1; //j is substitution in variable
    int minn = 0x3f3f3f3f;
    for(int i = 1 ; i < rlt.getRow();i++){
        if(minn>rlt(i,0)){minn = rlt(i,0);k=i;}
    }
    rlt = pivot_matrix(rlt,k,j); // pivot the matrix so that all constrainsMatrix bigger than 0
    base_ids[k-1]=j;
    rlt = simplex(rlt,base_ids);
    auto it = find(base_ids.begin(),base_ids.end(),rlt.getColumn()-1);
    if(it!=base_ids.end()){//if substitution element in vector
        int firstNotZero;
        int k;
        for(int i = 1;i<rlt.getRow();i++)
            if(rlt(0,i)!=0){
                firstNotZero = i;
                break;
            }
        it = find(base_ids.begin(),base_ids.end(),rlt.getColumn()-1);
        rlt = pivot_matrix(rlt, distance(base_ids.begin(),it),j);
        *it = j;
    }
    return rlt;

}
Matrix restorFromLaux(Matrix l_matrix,Matrix z, vector<int>&base_ids){
    vector<int>z_ids;
    for(int i=0;i<z.getColumn();i++){
        if(z(0,i)!=0)
            z_ids.push_back(i-1);
    }

    Matrix *restore_matrix = new Matrix(l_matrix.getRow(),l_matrix.getColumn()-1);
    for(int i = 0 ; i < l_matrix.getRow();i++) {
        for(int j = 0 ; j < l_matrix.getColumn()-1;j++){
            (*restore_matrix)(i,j) = l_matrix(i,j);
        }
    }
    for(int i = 0 ; i <restore_matrix->getColumn();i++){
        (*restore_matrix)(0,i)=z(0,i);
    }
    int cnt = 0;
    for(auto i : base_ids){
        auto f = find(z_ids.begin(),z_ids.end(),i );
        if(f!=z_ids.end()){
            Element t = (*restore_matrix)(0,i+1);
            for(int j = 0 ; j <restore_matrix->getColumn();j++){
                (*restore_matrix)(0,j) -= t*(*restore_matrix)(cnt+1,j);
            }
        }
        cnt++;
    }
    return *restore_matrix;
}
vector<Element> get_base_solution(Matrix matrix,vector<int>&base_ids){
    vector<Element>X(matrix.getColumn());
    for(int i = 0 ; i < base_ids.size();i++){
        if(base_ids[i]==-1)break;
        X[base_ids[i]]=matrix(i+1,0);
    }
    return X;
}

Matrix solve(Matrix coefficientMatrix, Matrix constrainsMatrix, Matrix minMatrix,Matrix equal) {
    auto looseMatrix = getLooseMatrix(coefficientMatrix); //get Loose Matrix
    looseMatrix = looseMatrix.insert(equal,looseMatrix);
    Matrix matrix = joinMatrix(looseMatrix,constrainsMatrix, minMatrix);//join three matrix into one matrix
    matrix.matrixPrint();
    vector<int> base_ids(matrix.getColumn());
    for(int i = 0 ; i < matrix.getColumn();i++)
        base_ids[i]=-1;
    for(int i = minMatrix.getColumn();i<constrainsMatrix.getColumn()+minMatrix.getColumn();i++)
        base_ids[i-minMatrix.getColumn()]=i;
    bool flag = false; // to judge if constrains Matrix has negative constrains;
    for(int i = 0 ; i < matrix.getRow() ;i++)
        if(matrix(i,0)<0){
            flag = true;
            break;
        }
    if(flag){
        cout<<"Constructing solution-assisted linear programming functions"<<endl;
        Matrix l_matrix = Laux(matrix,base_ids);
        if(l_matrix.getColumn()!=0&&l_matrix.getRow()!=0){
            matrix = restorFromLaux(l_matrix,matrix,base_ids);
        }
        else{
            cout<<"No solution"<<endl;
            Element a[3] = {5};
            return Matrix(0,0,a);
        }

    }
    auto ret_matrix = simplex(matrix,base_ids);
    auto X = get_base_solution(ret_matrix,base_ids);
    cout<<"the solution is"<<endl;
    for(int i = 0 ; i < minMatrix.getColumn();i++)
        cout<<X[i]<<' ';
    cout<<endl;
    cout<<"the value is "<<endl;
    cout<<-1*ret_matrix(0,0)<<endl;
}
Matrix solve(Matrix coefficientMatrix, Matrix constrainsMatrix, Matrix minMatrix){
    auto looseMatrix = getLooseMatrix(coefficientMatrix); //get Loose Matrix
    Matrix matrix = joinMatrix(looseMatrix,constrainsMatrix, minMatrix);//join three matrix into one matrix
    vector<int> base_ids(matrix.getColumn());
    for(int i = 0 ; i < matrix.getColumn();i++)
        base_ids[i]=-1;
    for(int i = minMatrix.getColumn();i<constrainsMatrix.getColumn()+minMatrix.getColumn();i++)
        base_ids[i-minMatrix.getColumn()]=i;
    bool flag = false; // to judge if constrains Matrix has negative constrains;
    for(int i = 0 ; i < matrix.getRow() ;i++)

        if(matrix(i,0)<0){
            flag = true;
            break;
        }
    if(flag){
        cout<<"Constructing solution-assisted linear programming functions"<<endl;
        Matrix l_matrix = Laux(matrix,base_ids);
        if(l_matrix.getColumn()!=0&&l_matrix.getRow()!=0){
            matrix = restorFromLaux(l_matrix,matrix,base_ids);
        }
        else{
            cout<<"No solution"<<endl;
            Element a[3] = {5};
            return Matrix(0,0,a);
        }

    }
    auto ret_matrix = simplex(matrix,base_ids);
    auto X = get_base_solution(ret_matrix,base_ids);
    cout<<"the solution is"<<endl;
    for(int i = 0 ; i < minMatrix.getColumn();i++)
        cout<<X[i]<<' ';
    cout<<endl;
    cout<<"the value is "<<endl;
    cout<<-1*ret_matrix(0,0)<<endl;
}
int main() {
//    Element l[6] = {1,-1,-1.5,1,50,20};
//    auto *a = new Matrix(3,2,l);//coefficientMatrix
//    Element s[3]={0,0,2000};
//    auto *b = new Matrix(1,3,s);//constrainMatrix
//    Element q[3]={-1,-1};
//    auto *c = new Matrix(1,2,q);//minMatrix
//    Element l[12] = {-21,0,0,0,0,-28,-4,0,-2,-1,-10,-11};
//    auto *a = new Matrix(3,4,l);//coefficientMatrix
//    Element s[3]={-500,-600,-250};
//    auto *b = new Matrix(1,3,s);//constrainMatrix
//    Element    q[4]={1,1,1,1};
//    auto *c = new Matrix(1,4,q);//minMatrix
//    double l[12] = {1,1,-1,-1};
//    auto *a = new Matrix(2,2,l);//coefficientMatrix
//    double s[3]={2,-1};
//    auto *b = new Matrix(1,2,s);//constrainMatrix
//    double q[4]={1,2};
//    auto *c = new Matrix(1,2,q);//minMatrix
    Element l[12] = {1,-2,1,4,-1,-2};
    auto *a = new Matrix(2,3,l);//coefficientMatrix
    Element e[5] = {-2,0,1,0,0};
    auto *equal = new Matrix(1,a->getRow()+3,e);
    Element s[3]={1,11,-3};
    auto *b = new Matrix(1,3,s);//constrainMatrix
    Element q[4]={-3,1,1};
    auto *c = new Matrix(1,3,q);//minMatrix

    solve(*a,*b,*c,*equal);
//    *a = getLooseMatrix(*a);
//    cout<<a->getRow()<<' '<<a->getColumn();
//    a->matrixPrint();
//    cout<<endl<<"============================";
//    *a = joinMatrix(*a,*b,*c);
//    cout<<a->getRow()<<' '<<a->getColumn()<<endl;
//    a->matrixPrint();



    return 0;
}

