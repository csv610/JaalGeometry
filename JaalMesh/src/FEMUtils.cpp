#include "FEMUtils.hpp"
//
int  getTexCoords1(const JFacePtr &face, const Point2D &xy, Point2D &uv)
{
    ////////////////////////////////////////////////////////////////////////////////////
    //     x = a0 + a1*xi + a2*eta + a3*xi*eta;
    //     y = b0 + b1*xi + b2*eta + b3*xi*eta;
    //
    //     This following is the inteverse of
    //              1    -1   -1   1
    //              1     1   -1  -1
    //              1     1    1   1
    //              1    -1    1  -1
    //
    ////////////////////////////////////////////////////////////////////////////////////
    cout << "Use iterative method " << endl;
    exit(0);
    Matrix4d  mat;
    mat(0,0) = 0.25;
    mat(0,1) = 0.25;
    mat(0,2) = 0.25;
    mat(0,3) = 0.25;
    mat(1,0) = -0.25;
    mat(1,1) =  0.25;
    mat(1,2) = 0.25;
    mat(1,3) = -0.25;
    mat(2,0) = -0.25;
    mat(2,1) = -0.25;
    mat(2,2) = 0.25;
    mat(2,3) =  0.25;
    mat(3,0) =  0.25;
    mat(3,1) = -0.25;
    mat(3,2) = 0.25;
    mat(3,3) = -0.25;

    Vector4d  x, y;
    for( int i = 0; i < 4; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        x(i) = xyz[0];
        y(i) = xyz[1];
    }
    Vector4d a = mat*x;
    Vector4d b = mat*y;

    // Calculate the coefficient of    A*y^2 + By + C = 0;
    //
    double x0 = xy[0] - a[0];
    double y0 = xy[1] - b[0];
    double A  = a[3]*b[2]  - a[2]*b[3];
    double B  = (x0*b[3] + a[1]*b[2])  - (y0*a[3] + a[2]*b[1]);
    double C  = x0*b[1]  - y0*a[1];

    double K = B*B - 4*A*C;
    assert( K >= 0);

    double u = 2.0, v = 2.0;

    if(fabs(A) < 1.0E-15) {
        assert(fabs(B) > 1.0E-10);
        v = -C/B;
    } else {
        double v1 = (-B + sqrt(K))/(2.0*A);
        double v2 = (-B - sqrt(K))/(2.0*A);
        if( v1 >= -1.0 && v1 <= 1.0) v = v1;
        if( v2 >= -1.0 && v2 <= 1.0) v = v2;
    }

    assert( v >= -1.0 && v <= 1.0);

    if( fabs(a[1] + v*a[3]) > 1.0E-10) {
        u = (x0 - a[2]*v)/(a[1] + a[3]*v);
        uv[0] = u;
        uv[1] = v;
        return 0;
    }

    if( fabs(a[3]) > 1.0E-10) {
        v = -a[1]/a[3];
        u = (y0*a[3] + a[1]*b[2])/(a[3]*b[1]-a[1]*b[3]);
        uv[0] =  u;
        uv[1] =  v;
        return 0;
    }

    if( fabs(a[2]) > 1.0E-10) {
        v = x0/a[2];
        u = (y0*a[2] + x0*b[2])/(a[2]*b[1]+x0*b[3]);
        uv[0] =  u;
        uv[1] =  v;
        return 0;
    }

    cout << "Error: Invalid reverse mapping: " << endl;
    exit(0);
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int  getTexCoords2(const JFacePtr &face, const Point2D &xy, Point2D &uv)
{
    Matrix4d  mat;
    mat(0,0) = 0.25;
    mat(0,1) = 0.25;
    mat(0,2) = 0.25;
    mat(0,3) = 0.25;
    mat(1,0) = -0.25;
    mat(1,1) =  0.25;
    mat(1,2) = 0.25;
    mat(1,3) = -0.25;
    mat(2,0) = -0.25;
    mat(2,1) = -0.25;
    mat(2,2) = 0.25;
    mat(2,3) =  0.25;
    mat(3,0) =  0.25;
    mat(3,1) = -0.25;
    mat(3,2) = 0.25;
    mat(3,3) = -0.25;

    Vector4d  x, y;
    for( int i = 0; i < 4; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        x(i) = xyz[0];
        y(i) = xyz[1];
    }
    Vector4d a = mat*x;
    Vector4d b = mat*y;

    Matrix2d  A1;
    A1(0,0) = a[1];
    A1(0,1) = a[2];
    A1(1,0) = b[1];
    A1(1,1) = b[2];
    Matrix2d A = A1.inverse();

    Vector2d param, rhs;
    param[0] = 0.0;
    param[1] = 0.0;
    double xg, yg, dx, dy, dl;
    for( int i = 0; i < 100; i++) {
        rhs[0] =  xy[0] - a[0] - a[3]*param[0]*param[1];
        rhs[1] =  xy[1] - b[0] - b[3]*param[0]*param[1];
        param  =  A*rhs;
        xg = a[0] + a[1]*param[0] + a[2]*param[1] + a[3]*param[0]*param[1];
        yg = b[0] + b[1]*param[0] + b[2]*param[1] + b[3]*param[0]*param[1];
        dx = xy[0] - xg;
        dy = xy[1] - yg;
        dl = sqrt(dx*dx + dy*dy);
        if( dl < 1.0E-12) break;
    }
    if( dl > 1.E0-12)
        cout << "Warning: Iteration did not converge " << dl << endl;

    uv[0] = param[0];
    uv[1] = param[1];
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool isSymmetric( const Eigen::SparseMatrix<double> &A)
{
    double eps = 1.0E-10;
    int nrows = A.rows();
    int ncols = A.cols();
    for( int i = 0;   i < nrows; i++) {
        for( int j = i+1; j < ncols; j++) {
            if( fabs(A.coeff(i,j) - A.coeff(j,i)) > eps ) {
                return 0;
            }
        }
    }
    return 1;
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//
bool isSymmetric( const MatrixXd &A)
{
    double eps = 1.0E-10;
    int nrows = A.rows();
    int ncols = A.cols();
    for( int i = 0;   i < nrows; i++) {
        for( int j = i+1; j < ncols; j++) {
            if( fabs(A.coeff(i,j) - A.coeff(j,i)) > eps ) {
                return 0;
            }
        }
    }
    return 1;
}
//
//
/////////////////////////////////////////////////////////////////////////////////////////
//

void latexMatrix( const MatrixXd &A)
{
    // Print a dense matrix in Latex (useful for the documentation);
    cout << fixed;
    int numRows = A.rows();
    int numCols = A.cols();
    cout << setprecision(5);
    cout << "\\begin{equation}" << endl;
    cout << "\\left[" << endl;
    cout << "\\begin{array}{";
    for( int j = 0; j < numCols; j++)
        cout << "r";
    cout << "}" << endl;
    for( int i = 0; i < numRows; i++) {
        for( int j = 0; j < numCols; j++) {
            cout << A(i,j);
            if( j == numCols-1 )
                cout << "  \\\\ ";
            else
                cout << " & ";
        }
        cout << endl;
    }
    cout << "\\end{array}" << endl;
    cout << "\\right]" << endl;
    cout << "\\end{equation}" << endl;
}
//
///////////////////////////////////////////////////////////////////////////////
//
void latexVector( const vector<double> &v)
{
    //
    // Write a vector in the Latex format (useful for the documentation)
    //
    cout << fixed;
    cout << "\\begin{array}{r}" << endl;
    for( size_t i = 0; i < v.size(); i++)
        cout << v[i] <<  " \\\\" << endl;
    cout << "\\end{array}" << endl;
}
//
///////////////////////////////////////////////////////////////////////////////
