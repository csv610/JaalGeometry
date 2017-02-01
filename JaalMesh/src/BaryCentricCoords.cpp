#include "BarycentricCoords.hpp"


////////////////////////////////////////////////////////////////////////////////

int JBarycentricCoordinates :: setFieldValues(const vector<Point2D> &cageCoords,
        const vector<double> &cageValues,
        const JMeshPtr &mesh, const  string &attribname)
{
    if( mesh == nullptr ) return 1;
    size_t numnodes = mesh->getSize(0);

    Point2D qPoint;
    double val;
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        qPoint  = vtx->getXYCoords();
        val = this->getFieldValueAt( cageCoords, cageValues, qPoint);
        vtx->setAttribute(attribname, val);
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JBarycentricCoordinates :: getXYCoords( const vector<Point2D> &cageCoords,
        const vector<double> &baryCoords, Point2D &xyCoord)
{
    if( cageCoords.size() !=  baryCoords.size() ) return 1;

    int nSize = cageCoords.size();

    double xsum = 0.0;
    double ysum = 0.0;
    for( int i = 0; i < nSize; i++) {
        xsum += baryCoords[i]*cageCoords[i][0];
        ysum += baryCoords[i]*cageCoords[i][1];
    }
    xyCoord[0] = xsum;
    xyCoord[1] = ysum;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JBarycentricCoordinates :: setXYCoords( const JMeshPtr &mesh, const vector<Point2D> &cageCoords)
{
    vector<double> baryCoords;
    size_t numnodes = mesh->getSize(0);
    Point2D xyCoord;

    for( size_t i  = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("BaryCoords", baryCoords);
            if( !err)  {
                getXYCoords( cageCoords, baryCoords, xyCoord);
                Point3D p3d = vtx->getXYZCoords();
                p3d[0] = xyCoord[0];
                p3d[1] = xyCoord[1];
                vtx->setXYZCoords(p3d);
            }
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeanValueCoordinates :: getW( const vector<Point2D> &cageCoords,
                                   const Point2D &queryCoord,
                                   vector<double> &w)
{
    //////////////////////////////////////////////////////////////////////////////////
    // Input :
    //      1. polyCoords :  Coordinates of closed polygon in the Counter
    //                       clockwise direction. The input is not tested inside.
    //
    //       2. queryCoord:   the xyCoords of the query Point
    // Output:
    //       1:  baryCoords: baryCentric Coords of the query Point.
    //
    // Reference: Mean Value Coordinates for Arbitrary Planar Polygons:
    //            Kai Hormann and Michael Floater;
    // Written by:
    //            Chaman Singh Verma
    //            University of Wisconsin at Madison.
    //            18th March, 2011.
    /////////////////////////////////////////////////////////////////////////////////

    int nSize = cageCoords.size();
    if( nSize == 0) return 1;
    w.resize(nSize);

    double dx, dy;

    vector<Point2D>  s(nSize);
    for( int i = 0; i < nSize; i++)
    {
        dx  =   cageCoords[i][0] - queryCoord[0];
        dy  =   cageCoords[i][1] - queryCoord[1];
        s[i][0]  = dx;
        s[i][1]  = dy;
        w[i]     = 0.0;
    }

    int ip, im;      // (i+1) and (i-1)
    double ri, rp, Ai, Di, dl, mu;  // Distance
    double eps = 10.0*std::numeric_limits<double>::min();

    // First check if any coordinates close to the cage point or
    // lie on the cage boundary. These are special cases.
    for( int i = 0; i < nSize; i++)
    {
        ip = (i+1)%nSize;
        ri = sqrt( s[i][0]*s[i][0] + s[i][1]*s[i][1] );
        Ai = 0.5*(s[i][0]*s[ip][1] - s[ip][0]*s[i][1]);
        Di = s[ip][0]*s[i][0] + s[ip][1]*s[i][1];
        if( ri <= eps)
        {
            w[i] = 1.0;
            return 0;
        }
        else if( fabs(Ai) <= 1.0E-10 && Di < 0.0)
        {
            dx = cageCoords[ip][0] - cageCoords[i][0];
            dy = cageCoords[ip][1] - cageCoords[i][1];
            dl = sqrt(dx*dx + dy*dy);
            assert( dl > 1.0E-06);
            assert(dl > eps);
            dx = queryCoord[0] - cageCoords[i][0];
            dy = queryCoord[1] - cageCoords[i][1];
            mu = sqrt(dx*dx + dy*dy)/dl;
            assert( mu >= 0.0 && mu <= 1.0);
            w[i]  = 1.0-mu;
            w[ip] = mu;
            return 0;
        }
    }

    // Page #12, from the paper
    vector<double> tanalpha(nSize); // tan(alpha/2)
    for( int i = 0; i < nSize; i++)
    {
        ip = (i+1)%nSize;
        im = (nSize-1+i)%nSize;
        ri = sqrt( s[i][0]*s[i][0] + s[i][1]*s[i][1] );
        rp = sqrt( s[ip][0]*s[ip][0] + s[ip][1]*s[ip][1] );
        Ai = 0.5*(s[i][0]*s[ip][1] - s[ip][0]*s[i][1]);
        if( fabs(Ai) < 1.0E-15) {
            return 1;
        }
        Di = s[ip][0]*s[i][0] + s[ip][1]*s[i][1];
        tanalpha[i] = (ri*rp - Di)/(2.0*Ai);
    }

    // Equation #11, from the paper
    for( int i = 0; i < nSize; i++)
    {
        im = (nSize-1+i)%nSize;
        ri = sqrt( s[i][0]*s[i][0] + s[i][1]*s[i][1] );
        if( ri < 1.0E-15) return 1;
        w[i] = (tanalpha[i] + tanalpha[im] )/ri;
    }

    double sum = 0.0;
    for( int i = 0; i < nSize; i++) sum += w[i];
    for( int i = 0; i < nSize; i++) w[i] /= sum;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeanValueCoordinates :: getGradW( const vector<Point2D> &cageCoords,
                                       const Point2D &qPoint, vector<Vec2D> &gradW)
{
    int n = cageCoords.size();

    vector<Vec2D>  e(n), s(n);
    vector<Vec2D>  c(n);
    vector<double>  t(n);
    vector<double>  sina(n);
    vector<double>  w(n);
    vector<double>  r(n);

    for( int i = 0; i < n; i++) {
        double dx = cageCoords[i][0] - qPoint[0];
        double dy = cageCoords[i][1] - qPoint[1];
        r[i]    = sqrt(dx*dx + dy*dy);
        if( r[i]  < 1.0E-15) return 1;
        s[i][0]  =  dx;
        s[i][1]  =  dy;
        e[i][0]  = dx/r[i];
        e[i][1]  = dy/r[i];
    }

    double ri, rp, Ai, Di;
    for( int i = 0; i < n; i++)  {
        int ip = (i+1)%n;
        sina[i] =  e[i][0]*e[ip][1] -  e[i][1]*e[ip][0];
        ri = sqrt( s[i][0]*s[i][0] + s[i][1]*s[i][1] );
        rp = sqrt( s[ip][0]*s[ip][0] + s[ip][1]*s[ip][1] );
        Ai = 0.5*(s[i][0]*s[ip][1] - s[ip][0]*s[i][1]);
        if( fabs(Ai) < 1.0E-15) return 2;
        Di = s[ip][0]*s[i][0] + s[ip][1]*s[i][1];
        t[i] = (ri*rp - Di)/(2.0*Ai);
    }

    for( int i = 0; i < n; i++) {
        int ip = (i+1)%n;
        c[i][0] = e[i][0]/r[i] - e[ip][0]/r[ip];
        c[i][1] = e[i][1]/r[i] - e[ip][1]/r[ip];
    }

    for( int i = 0; i < n; i++) {
        double a1 = c[i][0];
        double a2 = c[i][1];
        c[i][0]   = -a2;
        c[i][1]   =  a1;
    }

    for( int i = 0; i < n; i++) {
        int im = (n-1+i)%n;
        w[i] =  (t[i] + t[im])/r[i];
    }

    gradW.resize(n);
    for( int i = 0; i < n; i++) {
        int im = (n-1+i)%n;
        gradW[i][0] =  (t[im]*c[im][0]/sina[im] + t[i]*c[i][0]/sina[i] + w[i]*e[i][0])/r[i];
        gradW[i][1] =  (t[im]*c[im][1]/sina[im] + t[i]*c[i][1]/sina[i] + w[i]*e[i][1])/r[i];
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////

int JMeanValueCoordinates :: getCoords( const vector<Point2D> &cageCoords,
                                        const Point2D &queryCoord,
                                        vector<double> &baryCoords)
{
    vector<double>  w;
    getW( cageCoords, queryCoord, w);

    int nSize = cageCoords.size();

    double wsum = 0.0;
    for( int i = 0; i < nSize; i++) {
        wsum += w[i];
    }
    if( fabs(wsum)  < 1.0E-6) {
        cout << w.size() << endl;
        cout << "Error in MVC " << endl;
        cout << "Cage Coords " <<  endl;
        cout << setprecision(15) << fixed;
        for( int i = 0;  i < cageCoords.size(); i++)
            cout << cageCoords[i][0] << " " << cageCoords[i][1] << endl;
        cout << "Point " << queryCoord[0] << "  " << queryCoord[1] << endl;
        cout << "Weight " << endl;
        for( int i = 0;  i < cageCoords.size(); i++)
            cout << w[i] << endl;
        exit(0);
    }

    baryCoords.resize(nSize);
    for( int i = 0; i < nSize; i++) baryCoords[i] = w[i]/wsum;

    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////

int JMeanValueCoordinates :: getShapeGradients( const vector<Point2D> &cageCoords,
        const Point2D &qPoint, vector<Vec2D> &gradB)
{
    size_t nSize = cageCoords.size();

    vector<double>  w;
    getW( cageCoords, qPoint, w );
    if( w.size() != nSize ) return 1;

    vector<Vec2D>  gradW;
    getGradW( cageCoords, qPoint, gradW);
    if( gradW.size() != nSize ) return 2;

    double wsum = 0.0;
    double gradSum[] = {0.0, 0.0};

    for( int i = 0; i < nSize; i++)  {
        wsum += w[i];
        gradSum[0] +=  gradW[i][0];
        gradSum[1] +=  gradW[i][1];
    }

    // Simple quotion rule ...
    gradB.resize(nSize);
    for( int i = 0; i < nSize; i++)  {
        gradB[i][0] =  (gradW[i][0]*wsum - w[i]*gradSum[0])/( wsum*wsum);
        gradB[i][1] =  (gradW[i][1]*wsum - w[i]*gradSum[1])/( wsum*wsum);
    }

    /*
         double h = 1.0E-05;
         vector<double>  f0, f1, f2;
         Point2D qp;
         vector<Vec2D> gradBB;
         gradBB.resize(nSize);

         getCoords( cageCoords, qPoint, f0);

         qp[0] = qPoint[0] + h;
         qp[1] = qPoint[1];
         getCoords( cageCoords, qp, f1);
         for( int i = 0; i < nSize; i++)  gradBB[i][0] = (f1[i] - f0[i])/h;

         qp[0] = qPoint[0];
         qp[1] = qPoint[1] + h;
         getCoords( cageCoords, qp, f1);
         for( int i = 0; i < nSize; i++)  gradBB[i][1] = (f1[i] - f0[i])/h;

         double maxerror = 0.0;
         for( int i = 0; i < nSize; i++) {
              maxerror = max( maxerror, fabs( gradB[i][0] - gradBB[i][0] ));
              maxerror = max( maxerror, fabs( gradB[i][1] - gradBB[i][1] ));
         }
         cout << "Max error in grad calculattion " << maxerror << endl;
    */


    return 0;
}

////////////////////////////////////////////////////////////////////////////////

double JMeanValueCoordinates :: getFieldValueAt( const vector<Point2D> &cageCoords,
        const vector<double> &scalarField, const Point2D &queryCoord)
{
    int nSize = cageCoords.size();
    assert( (int) scalarField.size() == nSize);

    vector<double> baryCoords;
    getCoords( cageCoords, queryCoord, baryCoords);

    double sum = 0.0;
    for( int i = 0; i < nSize; i++) sum += baryCoords[i]*scalarField[i];
    return sum;
}

////////////////////////////////////////////////////////////////////////////////

void JMeanValueCoordinates :: setBaryCoords( const JMeshPtr &mesh, const vector<Point2D> &cageCoords)
{
    vector<double> baryCoords;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point2D &queryCoord = vtx->getXYCoords();
            getCoords( cageCoords,queryCoord, baryCoords);
            vtx->setAttribute("BaryCoords", baryCoords);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

#ifdef CSV

int  JGreenCoordinates :: getCoords( const vector<Point2D> &cageCoords,
                                     const Point2D  &queryCoord,
                                     vector<double>  &baryCoords)
{
    ////////////////////////////////////////////////////////////////////////////
    // Input :
    //      1. cageCoords : Coordinates of closed polygon in the
    //                      Clockwise direction. The input is not tested inside.
    //      2. queryCoord:  the xyCoords of the query Point
    // Output:
    //       1:  baryCoords: baryCentric Coords of the query Point with respect
    //           to cageCoords.
    //
    // Reference: Green Coordinates
    //            Yaron Lipman, David Levin and Daniel Cohen-or
    // Written by:
    //            Chaman Singh Verma
    //            University of Wisconsin at Madison.
    //            18th March, 2011.
    ///////////////////////////////////////////////////////////////////////////
    //
    int nSize = cageCoords.size();

    Point2D  a, b, nj;
    double   V, Q, S, R, BA, SRT, L0, L1, A0, A1, A10, L10;

    vector<double> phi( nSize );
    vector<double> psi( nSize );

    for( int i = 0; i < nSize; i++) phi[i] = 0.0;
    for( int i = 0; i < nSize; i++) psi[i] = 0.0;

    for( int j = 0; j < nSize; j++)
    {
        const Point2D &vj1 =  cageCoords[j];
        b[0] = vj1[0] - queryCoord[0];
        b[1] = vj1[1] - queryCoord[1];
        if( b[0]*b[0] + b[1]*b[1]  <= 1.0E-20)
        {
            cout << " TOO ClOSED " << endl;
            return 1;
        }
    }

    for( int j = 0; j < nSize; j++)
    {
        int jp = (j+1)%nSize;
        const Point2D &vj1 =  cageCoords[j];
        const Point2D &vj2 =  cageCoords[jp];

        a[0] = vj2[0] - vj1[0];
        a[1] = vj2[1] - vj1[1];

        nj[0] =  a[1];
        nj[1] = -a[0];

        b[0] = vj1[0] - queryCoord[0];
        b[1] = vj1[1] - queryCoord[1];

        Q     = a[0]*a[0] + a[1]*a[1];
        S     = b[0]*b[0] + b[1]*b[1];
        R     = 2.0*( a[0]*b[0] + a[1]*b[1]);
        BA    = b[0]*nj[0] + b[1]*nj[1];   // Remember do not normalize Normal vector here.
        V     = 4.0*S*Q - R*R;
        assert( V >  0.0);
        SRT   = sqrt( V );
        assert(SRT > 0.0 );
        L0    = log(S);
        L1    = log( S + Q + R );
        A0    = atan( R/SRT)/SRT;
        A1    = atan((2*Q+R)/SRT)/SRT;
        A10   = A1-A0;
        L10   = L1-L0;
        psi[j]    = sqrt(Q)/(4.0*M_PI)*(( 4.0*S-(R*R/Q))*A10 + (R/(2.0*Q))*L10 + L1 - 2);
        phi[jp]  -=  (BA/(2.0*M_PI))*( (L10/(2.0*Q)) - A10*R/Q );
        phi[j]   +=  (BA/(2.0*M_PI))*( (L10/(2.0*Q)) - A10*(2.0 + R/Q) );
    }

    double sum = 0.0;
    for( int i = 0; i < nSize; i++)
        sum += phi[i];

    for( int i = 0; i < nSize; i++)
        phi[i] /= sum;

    sum = 0.0;
    for( int i = 0; i < nSize; i++)
        sum += phi[i];
    assert(fabs(sum -1.0) < 1.0E-06);

    baryCoords.resize(2*nSize);
    for( int i = 0; i < nSize; i++)
    {
        baryCoords[i] = phi[i];
        baryCoords[i+nSize] = psi[i];
    }

    return 0;
}

double JGreenCoordinates :: getFieldValueAt( const vector<Point2D> &cageCoords,
        const vector<double> &scalarField, const Point2D &queryCoord)
{
    int nSize = cageCoords.size();
    assert( (int) scalarField.size() == nSize);

    vector<double> baryCoords;
    getCoords( cageCoords, queryCoord, baryCoords);

    double sum = 0.0;
    for( int i = 0; i < nSize; i++) sum += baryCoords[i]*scalarField[i];
    return sum;
}
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <algorithm>
#include <math.h>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include "StopWatch.hpp"

typedef mtl::compressed2D<double> SMatrix;
typedef std::vector<Vertex*> Segment;

double rellength = 0.001;

struct SparseMatrix
{
public:
    typedef std::map<int, double>  SparseRow;

    int setSize( int n, int m)
    {
        nRows = n;
        nCols = m;
        data.resize(nRows);
    }

    int  getNumOfRows() const
    {
        return nRows;
    }

    int  getNumOfCols() const
    {
        return nCols;
    }

    bool isSymmetric() const
    {
        SparseRow row;
        SparseRow::const_iterator it;
        for( int irow = 0; irow < nRows; irow++)
        {
            row = data[irow];
            for( it = row.begin(); it != row.end(); ++it)
            {
                int jcol = it->first;
                if( getVal(irow,jcol) != getVal( jcol, irow) ) return 0;
            }
        }
        return 1;
    }

    size_t getNumOfNonZeros() const
    {
        size_t ncount = 0;
        for( int i = 0; i < nRows; i++)
            for( int j = 0; j < nCols; j++)
                if( getVal(i,j) != 0.0) ncount++;
        return ncount;
    }

    const SparseRow  &getRow(int i) const
    {
        assert( i < nRows);
        return data[i];
    }

    void setVal(int i, int j, const double &val)
    {
        assert( i >= 0 && i < nRows );
        assert( j >= 0 && j < nCols );

        if( val == 0.0)
        {
            data[i].erase(j);
            return;
        }

        data[i][j] = val;
    }

    double getVal( int i, int j ) const
    {
        assert( i >=0 && i < nRows );
        assert( j >=0 && j < nCols );

        SparseRow::const_iterator it = data[i].find(j);
        if( it == data[i].end() ) return 0.0;

        return it->second;
    }

    void clear()
    {
        data.clear();
    }

    void print() const;

private:
    typedef std::vector<SparseRow> Matrix;

    int nRows, nCols;
    Matrix data;
};

using namespace std;

////////////////////////////////////////////////////////////////////////////////

void MatMat_Mult ( const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &C)
{
    int nRows = A.getNumOfRows();
    int nCols = B.getNumOfCols();
    int nK    = A.getNumOfCols();

    assert( nK == B.getNumOfRows() );

    C.setSize( nRows, nCols);
    SparseMatrix::SparseRow::const_iterator it;

    for( int i = 0; i < nRows; i++)
    {
        for( int j = 0; j < nCols; j++)
        {
            double sum = 0.0;
            SparseMatrix::SparseRow row = A.getRow(i);
            for( it = row.begin(); it != row.end(); ++it)
            {
                int col = it->first;
                sum += (it->second)*B.getVal(col, j);
            }
            C.setVal(i, j, sum );
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////

void MatVec_Mult(const SparseMatrix &A, const vector<double> &x, vector<double> &y)
{
    size_t nRows = A.getNumOfRows();
    size_t nCols = A.getNumOfCols();

    assert( x.size() == nCols );
    y.resize( nRows );

    SparseMatrix::SparseRow::const_iterator it;

    for( int i = 0; i < nRows; i++)
    {
        double sum = 0.0;
        const SparseMatrix::SparseRow &row = A.getRow( i );
        for( it = row.begin(); it != row.end(); ++it)
        {
            int j = it->first;
            sum += A.getVal(i,j)*x[j];
        }
        y[i] = sum;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

int Jacobi( const SparseMatrix &A, const vector<double> &b, vector<double> &x,
            double &maxerror )
{
    vector<double> tmp;

    size_t nRows = A.getNumOfRows();
    size_t nCols = A.getNumOfCols();

    assert( b.size() == nRows );

    x.resize(nRows);
    for( size_t i = 0; i < nRows; i++)
        x[i] = 0.0;

    SparseMatrix::SparseRow::const_iterator it;
    size_t iter = 0;
    while(1)
    {
        for( size_t i = 0; i < nRows; i++)
        {
            double sum = b[i];
            const SparseMatrix::SparseRow &row = A.getRow(i);
            for ( it = row.begin(); it != row.end(); ++it)
            {
                int j = it->first;
                if( i != j) sum += - (it->second)*x[j];
            }
            x[i] = sum/A.getVal(i, i);
        }

        MatVec_Mult(A, x, tmp);
        maxerror = 0.0;
        for( size_t i = 0; i < nRows; i++)
            maxerror = max( maxerror, fabs(tmp[i] - b[i] ));
        if( maxerror < 1.0E-08) break;
        if( iter++ == 1000000) break;
        cout << "iteration " << iter << ": resid " << maxerror << endl;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

double getBoundaryLength(const vector<Vertex*> &nodes)
{
    int nSize = nodes.size();
    double sum = 0.0;
    for( int i = 0; i < nSize; i++)
    {
        Point2D &p0 = nodes[i]->getCoords();
        Point2D &p1 = nodes[(i+1)%nSize]->getCoords();
        double  dx  = p1[0] - p0[0];
        double  dy  = p1[1] - p0[1];
        sum  += sqrt( dx*dx + dy*dy );
    }
    return sum;
}
///////////////////////////////////////////////////////////////////////////////

void refinedBoundary( const vector<Vertex*> &nodes, double minlen,
                      vector<Segment> &refinedSegments)
{
    size_t numNodes = nodes.size();

    if( numNodes < 2 ) return;

    refinedSegments.resize(numNodes);

    Point2D p;
    for( int i = 0; i < numNodes; i++)
    {
        Vertex  *v0 = nodes[i];
        Vertex  *v1 = nodes[(i+1)%numNodes];
        Point2D &p0 = v0->getCoords();
        Point2D &p1 = v1->getCoords();
        double  dx  = p1[0] - p0[0];
        double  dy  = p1[1] - p0[1];
        double  dl  = sqrt( dx*dx + dy*dy );
        int n = (ceil)(dl/minlen);  // Higher value better.
        refinedSegments[i].reserve(n+1);
        refinedSegments[i].push_back( v0 );
        if( n )
        {
            double dt = 1.0/(double)n;
            for( int j = 1; j < n; j++)
            {
                double t = j*dt;
                p[0]  = (1.0-t)*p0[0] + t*p1[0];
                p[1]  = (1.0-t)*p0[1] + t*p1[1];
                Vertex *vm = new Vertex;
                vm->setCoords( p );
                refinedSegments[i].push_back( vm );
            }
        }
        refinedSegments[i].push_back( v1 );
    }
}

///////////////////////////////////////////////////////////////////////////////

FaceMesh *delaunay_triangulation( const vector<Segment> &cageSegments,
                                  const vector<Segment> &modelSegments)
{
    int numCageSegments   = cageSegments.size();
    int numModelSegments  = modelSegments.size();

    // Mark all Nodes Visit = 0;
    for( int i = 0; i < numCageSegments; i++)
    {
        int nSub = cageSegments[i].size();
        for( int j = 0; j < nSub; j++)
        {
            Vertex *v = cageSegments[i][j];
            v->setVisitMark(0);
        }
    }

    for( int i = 0; i < numModelSegments; i++)
    {
        int nSub = modelSegments[i].size();
        for( int j = 0; j < nSub; j++)
        {
            Vertex *v = modelSegments[i][j];
            v->setVisitMark(0);
        }
    }

    vector<Vertex*>  trinodes;

    int numEdges = 0;
    int nodeid   = 0;

    for( int i = 0; i < numModelSegments; i++)
    {
        int nSub = modelSegments[i].size();
        numEdges += nSub-1;
        for( int j = 0; j < nSub; j++)
        {
            Vertex *v = modelSegments[i][j];
            if( !v->isVisited() )
            {
                v->setBoundaryMark(2);
                v->setID( nodeid++);
                v->setVisitMark(1);
                trinodes.push_back(v);
            }
        }
    }

    for( int i = 0; i < numCageSegments; i++)
    {
        int nSub = cageSegments[i].size();
        numEdges += nSub-1;
        for( int j = 0; j < nSub; j++)
        {
            Vertex *v = cageSegments[i][j];
            if( !v->isVisited() )
            {
                v->setBoundaryMark(1);
                v->setID( nodeid++);
                v->setVisitMark(1);
                trinodes.push_back(v);
            }
        }
    }

    ofstream ofile( "data.poly", ios::out);
    ofile << trinodes.size() << " 2  0  1 " << endl;
    int numNodes = trinodes.size();
    for( int i = 0; i < numNodes; i++)
    {
        Vertex *v = trinodes[i];
        Point2D p = v->getCoords();
        ofile << v->getID() << setprecision(25) << " "
              << p[0] << " " << p[1] << " " << v->isBoundary() << endl;
    }
    ofile << numEdges << " 0 " << endl;

    int edgeid  = 1;
    for( int i = 0; i < numModelSegments; i++)
    {
        int nSub = modelSegments[i].size();
        for( int j = 0; j < nSub-1; j++)
        {
            Vertex *v0 = modelSegments[i][j];
            Vertex *v1 = modelSegments[i][j+1];
            ofile <<  edgeid++  << " " << v0->getID() << " " << v1->getID() << endl;
        }
    }

    for( int i = 0; i < numCageSegments; i++)
    {
        int nSub = cageSegments[i].size();
        for( int j = 0; j < nSub-1; j++)
        {
            Vertex *v0 = cageSegments[i][j];
            Vertex *v1 = cageSegments[i][j+1];
            ofile <<  edgeid++  << " " << v0->getID() << " " << v1->getID() << endl;
        }
    }

    ofile << " 0 " << endl;
    ofile.close();

    int err = system( "triangle -pga0.001q30.0Y data.poly");

    ifstream infile( "data.1.off", ios::in);
    assert( !infile.fail() ) ;

    string str;
    infile >> str;
    assert( str == "OFF");

    int numBoundNodes = trinodes.size();
    assert( numBoundNodes );

    int numFaces;
    infile >> numNodes >> numFaces >> numEdges;

    assert( numNodes >= numBoundNodes);

    assert( numNodes >=3 );

    trinodes.resize( numNodes );

    FaceMesh *mesh = new FaceMesh;
    assert( mesh );

    mesh->nodes.resize( numNodes );

    double  z;
    Point2D p;
    for( int i = 0; i < numNodes; i++)
    {
        infile >> p[0] >> p[1] >> z;
        if( i >= numBoundNodes )
        {
            Vertex *vm = new Vertex;
            vm->setCoords( p );
            trinodes[i] = vm;
            trinodes[i]->setID(i);
        }
        mesh->nodes[i] = trinodes[i];
    }

    mesh->faces.resize( numFaces );

    int nnodes, vid[3];
    vector<Vertex*> triConnect(3);
    for( int i = 0; i < numFaces; i++)
    {
        infile >> nnodes >> vid[0] >> vid[1] >> vid[2];
        assert( nnodes == 3 );
        triConnect[0] = trinodes[ vid[0] ];
        triConnect[1] = trinodes[ vid[1] ];
        triConnect[2] = trinodes[ vid[2] ];
        Face *face = new Face;
        face->setNodes( triConnect );
        mesh->faces[i] = face;
    }

    for( int i = 0; i < numModelSegments; i++)
    {
        int nSub = modelSegments[i].size();
        for( int j = 0; j < nSub; j++)
        {
            Vertex *v = modelSegments[i][j];
            v->setBoundaryMark(0);
        }
    }

    return mesh;

}

///////////////////////////////////////////////////////////////////////////////

int build_laplacian_matrix( FaceMesh *mesh,  SparseMatrix &matrix)
{
    mesh->buildRelations(0,0);

    int  numNodes = mesh->getSize(0);
    assert( numNodes );
    matrix.setSize( numNodes, numNodes);

    for( int i = 0; i < numNodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        v->setID(i);
    }

    vector<Vertex*> vneighs;
    for( int i = 0; i < numNodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        int irow  = v->getID();
        if( v->isBoundary() )
        {
            matrix.setVal(irow,irow,1.0);
        }
        else
        {
            vneighs = v->getRelations0();
            for( size_t j = 0; j < vneighs.size(); j++)
            {
                int jcol = vneighs[j]->getID();
                matrix.setVal(irow, jcol, 1.0);
            }
            matrix.setVal(irow, irow, -(double)vneighs.size() );
        }
    }
    return 0;
}

SMatrix* mtl_matrix(const SparseMatrix &matrix, bool deleteOrg = 1)
{
    int numRows = matrix.getNumOfRows();

    SMatrix *smat = new SMatrix(numRows, numRows );
    assert( smat );
    mtl::matrix::inserter<SMatrix>  matrix_inserter(*smat);

    SparseMatrix::SparseRow::const_iterator it;
    for( int irow = 0; irow < numRows; irow++)
    {
        SparseMatrix::SparseRow row = matrix.getRow(irow);
        for( it = row.begin(); it != row.end(); ++it)
        {
            int jcol   = it->first;
            double aij = it->second;
            matrix_inserter(irow, jcol) <<  (double)aij;
        }
    }
    return smat;
}


//////////////////////////////////////////////////////////////////////////////////

int apply_boundary_conditions( Vertex *cageNode,
                               const vector<Segment> &cageSegments,
                               vector<double> &b)
{
    int nSize = cageSegments.size();
    int next_seg = -1;
    for( int i = 0; i < nSize; i++)
    {
        if( cageSegments[i].front() == cageNode )
        {
            next_seg = i;
            break;
        }
    }

    assert( next_seg >= 0);
    int prev_seg = (nSize-1 + next_seg)%nSize;
    assert( cageSegments[prev_seg].back() == cageNode );

    for( int i = 0; i < b.size(); i++)
        b[i] = 0.0;

    int id, nSub;
    double dt;

    nSub = cageSegments[next_seg].size();
    assert( nSub >= 2);
    dt = 1.0/(nSub-1);
    for( int i = 0; i < nSub; i++ )
    {
        id = cageSegments[next_seg][i]->getID();
        b[id]  = 1.0-i*dt;
    }

    nSub = cageSegments[prev_seg].size();
    assert( nSub >= 2);
    dt = 1.0/(nSub-1);
    for( int i = 0; i < nSub; i++ )
    {
        id = cageSegments[prev_seg][i]->getID();
        b[id]  = i*dt;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void modify_linear_system( SparseMatrix &A, vector<double> &b)
{
    //
    // Modify the System so that it becomes symmstric. This is important for many
    // CG type methods.
    //
    size_t nRows = A.getNumOfRows();
    vector<double>  bb( nRows);

    for( int i = 0; i < nRows; i++)
        bb[i] = 0.0;

    std::set<int> bound;
    for( int i = 0; i < nRows; i++)
    {
        SparseMatrix::SparseRow row = A.getRow(i);
        if( row.size() == 1) bound.insert(i);
    }
    assert( bound.size() > 0);

    SparseMatrix::SparseRow::const_iterator it;
    for( int irow = 0; irow < nRows; irow++)
    {
        SparseMatrix::SparseRow row = A.getRow(irow);
        if( row.size() > 1)
        {
            for( it = row.begin(); it != row.end(); ++it)
            {
                int jcol = it->first;
                if( bound.find(jcol) != bound.end() )
                {
                    double aij = it->second;
                    bb[irow]  -= aij*b[jcol];
                    A.setVal( irow, jcol, 0.0);
                }
            }
        }
    }

    for( int i = 0; i < nRows; i++)
        b[i] += bb[i];

// Now the matrix must be symmetric...
    if( !A.isSymmetric()  )
    {
        cout << "Fatal Error: The matix must be symmetric now " << endl;
        exit(0);
    }

}
//////////////////////////////////////////////////////////////////////////////////

int solve_laplace( const SparseMatrix &A, const vector<double> &b,
                   vector<double> &x, double &maxerror)
{
    int nrows = A.getNumOfRows();
    assert( b.size() == nrows);
    x.resize( nrows );

    int err = Jacobi( A, b, x, maxerror);
    for( int i = 0; i < x.size(); i++)
    {
        if( x[i] < 0.0 || x[i] > 1.0)
        {
            cout << "Warning: Wrong Laplace Solution " << endl;
            return 1;
        }
    }
    cout << "Solver Converged : MaxError  " << maxerror << endl;
    return 0;

    SMatrix *smat = mtl_matrix(A);

    int N = nrows;

    // Create an ILU(0) preconditioner
    itl::pc::ilu_0<SMatrix>        P(*smat);
//  itl::pc::diagonal<SMatrix>     P(*smat);
//  itl::pc::identity<SMatrix>     P(*smat);

    // Set b such that x == 1 is solution; start with x == 0
    mtl::dense_vector<double>   xx(N), bb(N);
    for( int i = 0; i < N; i++)
    {
        bb[i] = b[i];
        xx[i] = 0.0;
    }

    // Termination criterion: r < 1e-6 * b or N iterations
    itl::noisy_iteration<double>  iter(bb, 1000, 1.0E-10);

    // Solve Ax == b with left preconditioner P
    itl::cg( *smat, xx, bb, P, iter);

    for( int i = 0; i < N; i++)
        x[i]  = xx[i];

    delete smat;
}

//////////////////////////////////////////////////////////////////////////////////

int writeVTK( const string &filename, FaceMesh *mesh, const vector<double> &attrib)
{
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    size_t numNodes = mesh->getSize(0);

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << " Jaal Mesh Generator " << endl;

    ofile << "ASCII " << endl;

    ofile << "DATASET UNSTRUCTURED_GRID " << endl;
    ofile << "POINTS " << numNodes << " float " << endl;
    ofile << std::setiosflags(ios::fixed);

    for( size_t i = 0; i < numNodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        const Point2D &xyz = vertex->getCoords();
        ofile << xyz[0] << " " << xyz[1] << "  0.0  " << endl;
    }

    size_t numFaces = mesh->getSize(2);
    size_t numConn = 0;

    for( size_t i = 0; i < numFaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        numConn += face->getSize(0);
    }

    ofile << "CELLS " << numFaces << "  " << numFaces + numConn << endl;

    for( size_t i = 0; i < numFaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        ofile << nsize << " ";
        for(int i = 0; i < nsize; i++)
        {
            Vertex *v = face->getNodeAt(i);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }

    ofile << "CELL_TYPES " << numFaces << endl;
    for( size_t i = 0; i < numFaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        int type  = 7;
        if( nsize == 3 ) type = 5;
        if( nsize == 4 ) type = 9;
        ofile << type << " ";
    }

    ofile << "POINT_DATA " << numNodes << endl;
    ofile << "SCALARS BaryCoord float  1" << endl;
    ofile << "LOOKUP_TABLE default" << endl;
    for( int i = 0; i < numNodes; i++)
        ofile << attrib[i] << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////

int loadPoly ( const string &s, vector<Vertex*> &nodes )
{
    ifstream infile( s.c_str(), ios::in);
    if( infile.fail() ) return 1;

    int numNodes, dim , nattribs, boundmark;
    infile >> numNodes >> dim >> nattribs >> boundmark;

    assert( numNodes > 0);
    assert( dim  == 2);
    assert( nattribs == 0);

    nodes.resize( numNodes );
    int vid, bid;
    Point2D p;
    for( int i = 0; i < numNodes; i++)
    {
        infile >> vid >> p[0] >> p[1];
        Vertex *v = new Vertex;
        v->setCoords(p);
        if( boundmark )
        {
            infile >> bid;
            v->setBoundaryMark(bid);
        }
        nodes[i] = v;
    }
    return 0;
}


///////////////////////////////////////////////////////////////////////////////////////

void HarmonicCoordinates( const vector<Vertex*> & modelNodes,
                          const vector<Vertex*> & cageNodes, bool genfield)
{
    cout << " Harmonic Coordinates" << endl;
    StopWatch  stopwatch;
    stopwatch.start();

    int numCageNodes   = cageNodes.size();
    int numModelNodes  = modelNodes.size();

    cout << " #of Cage  Points      " << numCageNodes  << endl;
    cout << " #of Model Points      " << numModelNodes << endl;

    double len0 = getBoundaryLength( modelNodes );
    vector<Segment> modelSegments;
    refinedBoundary( modelNodes,  rellength*len0, modelSegments);

    // Find the length of the cage ...
    double len1 = getBoundaryLength( cageNodes );
    vector<Segment> cageSegments;
    refinedBoundary( cageNodes,  rellength*len1, cageSegments);

    FaceMesh *trimesh = delaunay_triangulation( cageSegments, modelSegments);

    assert( trimesh );

    size_t numNodes = trimesh->getSize(0);

    cout << "#Nodes for HC mesh " << numNodes << endl;

    vector<double>  b, x;
    double maxerror = 1.0E-06;

    b.resize(numNodes);

    SparseMatrix A, spmat;
    build_laplacian_matrix(trimesh, spmat);

    vector< vector<double> > modelBaryCoords;
    modelBaryCoords.resize( numModelNodes );

    b.resize(numNodes);
    for( int i = 0; i < numModelNodes; i++)
        modelBaryCoords[i].resize(numCageNodes);

    for( int j = 0; j < cageNodes.size(); j++)
    {
        apply_boundary_conditions(cageNodes[j], cageSegments, b);
        A = spmat;
        modify_linear_system(A, b);
        int err = solve_laplace( A, b, x, maxerror);
        if( genfield ) {
            ostringstream oss;
            oss << "hc" << j << ".vtk";
            writeVTK( oss.str(), trimesh, x );
        }
        for( size_t i = 0; i < numModelNodes; i++)
        {
            int id = modelNodes[i]->getID();
            modelBaryCoords[i][j] =  x[id];
        }
    }
    stopwatch.stop();

    cout << "# Samples " << numNodes << endl;
    cout << "HC Execution Time : " << stopwatch.getSeconds() << endl;

    ofstream ofile("hc.dat", ios::out);
    ofile << "H " << numModelNodes << "  " << numCageNodes << endl;
    for( int i = 0; i < numModelNodes; i++)
    {
        Vertex *v = modelNodes[i];
        v->setBaryCoords( modelBaryCoords[i] );
        for( int j = 0; j < numCageNodes; j++)
            ofile <<  modelBaryCoords[i][j] << " ";
        ofile << endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////

void MeanValueCoordinates( const vector<Vertex*> & modelNodes,
                           const vector<Vertex*> & cageNodes, bool genfield )
{
    cout << " Mean Value Coordinates" << endl;
    int numCageNodes   = cageNodes.size();
    int numModelNodes  = modelNodes.size();

    cout << " #of Cage  Points      " << numCageNodes  << endl;
    cout << " #of Model Points      " << numModelNodes << endl;

    vector<Point2D> cageCoords( numCageNodes );
    vector<double> x(numCageNodes);
    vector<double> y(numCageNodes);
    for( int j = 0; j < numCageNodes; j++)
    {
        const Point2D &p = cageNodes[j]->getCoords();
        x[j] = p[0];
        y[j] = p[1];
        cageCoords[j] = p;
    }

    int dir = 1;
    double area = polygon_area( &x[0], &y[0], numCageNodes);
    if( area < 0)
    {
        cout << " Reverse the direction " << endl;
        std::reverse( cageCoords.begin(), cageCoords.end() );
        dir = -1;
    }

    vector<double> baryCoords;
    for( int i = 0; i < numModelNodes; i++)
    {
        Point2D &queryCoord = modelNodes[i]->getCoords();
        mean_value_coordinates( cageCoords, queryCoord, baryCoords);
        if( dir == -1)
            std::reverse( baryCoords.begin(), baryCoords.end() );
        modelNodes[i]->setBaryCoords( baryCoords );
    }

    if( genfield )
    {
        StopWatch stopwatch;

        stopwatch.start();
        double len0 = getBoundaryLength( modelNodes );
        vector<Segment> modelSegments;
        refinedBoundary( modelNodes,  rellength*len0, modelSegments);

        // Find the length of the cage ...
        double len1 = getBoundaryLength( cageNodes );
        vector<Segment> cageSegments;
        refinedBoundary( cageNodes,  rellength*len1, cageSegments);

        FaceMesh *trimesh = delaunay_triangulation( cageSegments, modelSegments);

        assert( trimesh );

        size_t numNodes = trimesh->getSize(0);

        vector<double>  b, x;
        double maxerror = 1.0E-06;

        b.resize(numNodes);

        vector<vector<double> > triBaryCoords;

        triBaryCoords.resize( numNodes );

        for( int i = 0; i < numNodes; i++)
        {
            Point2D &queryCoord = trimesh->getNodeAt(i)->getCoords();
            int id = trimesh->getNodeAt(i)->getID();
            mean_value_coordinates( cageCoords, queryCoord, baryCoords);
            if( dir == -1)
                std::reverse( baryCoords.begin(), baryCoords.end() );
            triBaryCoords[id] =  baryCoords;
        }
        stopwatch.stop();
        cout << "# Samples " << numNodes << endl;
        cout << " Mean Value Calculation Time : " << stopwatch.getSeconds() << endl;

        x.resize(numNodes);
        for( int j = 0; j < numCageNodes; j++)
        {
            for( int i = 0; i < numNodes; i++)
            {
                int id = trimesh->getNodeAt(i)->getID();
                x[id]   = triBaryCoords[id][j];
            }
            ostringstream oss;
            oss << "mvc" << j << ".vtk";
            writeVTK( oss.str(), trimesh, x );
        }
    }

}

//////////////////////////////////////////////////////////////////////////////////

void GreenCoordinates( const vector<Vertex*>  & modelNodes,
                       const vector<Vertex*>    & cageNodes, bool genfield)
{
    cout << "Green Coordinates " << endl;
    int numCageNodes   = cageNodes.size();
    int numModelNodes  = modelNodes.size();

    cout << " #of Cage  Points      " << numCageNodes  << endl;
    cout << " #of Model Points      " << numModelNodes << endl;

    // For Green's Coordinates the cage must be clockwise oriented:
    vector<double> x(numCageNodes);
    vector<double> y(numCageNodes);
    vector<Point2D> cageCoords( numCageNodes );
    for( int j = 0; j < numCageNodes; j++)
    {
        const Point2D &p = cageNodes[j]->getCoords();
        x[j] = p[0];
        y[j] = p[1];
        cageCoords[j] = p;
    }

    int dir = 1;
    double area = polygon_area( &x[0], &y[0], numCageNodes);
    /*
        if( area > 0) {
            cout << " Reverse the direction " << endl;
            std::reverse( cageCoords.begin(), cageCoords.end() );
            dir = -1;
        }
    */

    vector<double> baryCoords;
    for( int i = 0; i < numModelNodes; i++)
    {
        Point2D &queryCoord = modelNodes[i]->getCoords();
        green_coordinates( cageCoords, queryCoord, baryCoords);
        if( dir == -1)
            std::reverse( baryCoords.begin(), baryCoords.end() );
        modelNodes[i]->setBaryCoords( baryCoords );
    }

    if( genfield )
    {
        StopWatch  stopwatch;
        stopwatch.start();

        double len0 = getBoundaryLength( modelNodes );
        vector<Segment> modelSegments;
        refinedBoundary( modelNodes,  rellength*len0, modelSegments);

        // Find the length of the cage ...
        double len1 = getBoundaryLength( cageNodes );
        vector<Segment> cageSegments;
        refinedBoundary( cageNodes,  rellength*len1, cageSegments);

        FaceMesh *trimesh = delaunay_triangulation( cageSegments, modelSegments);

        assert( trimesh );

        size_t numNodes = trimesh->getSize(0);

        vector<vector<double> > triBaryCoords;

        triBaryCoords.resize( numNodes );
        for( int i = 0; i < numNodes; i++)
        {
            triBaryCoords[i].resize( 2*numCageNodes );
            for( int j = 0; j < 2*numCageNodes; j++ )
                triBaryCoords[i][j] = 0.0;
        }

        for( int i = 0; i < numNodes; i++)
        {
            Point2D &queryCoord = trimesh->getNodeAt(i)->getCoords();
            int id = trimesh->getNodeAt(i)->getID();
            if( !trimesh->getNodeAt(i)->isBoundary() )
            {
                green_coordinates( cageCoords, queryCoord, baryCoords);
                if( dir == -1)
                    std::reverse( baryCoords.begin(), baryCoords.end() );
                triBaryCoords[id] =  baryCoords;
            }
        }
        stopwatch.stop();

        cout << "# Samples " << numNodes << endl;
        cout << "Green Coordinates Calculation Time : " << stopwatch.getSeconds()<< endl;

        vector<double>  x, y;
        x.resize(numNodes);
        y.resize(numNodes);
        for( int j = 0; j < numCageNodes; j++)
        {
            for( int i = 0; i < numNodes; i++)
            {
                int id = trimesh->getNodeAt(i)->getID();
                x[id]   = triBaryCoords[id][j];
                y[id]   = triBaryCoords[id][j+ numCageNodes];
            }
            ostringstream oss1;
            oss1 << "gcphi" << j << ".vtk";
            writeVTK( oss1.str(), trimesh, x );

            ostringstream oss2;
            oss2 << "gcpsi" << j << ".vtk";
            writeVTK( oss2.str(), trimesh, y );
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////


void setCageNormals( const vector<Vertex*> &cageNodes, vector<Point2D> &normals)
{
    int nSize = cageNodes.size();

    normals.resize(nSize);

    Point2D  a, nj;
    double    len;
    for( int j = 0; j < nSize; j++)
    {
        const Point2D &vj1 =  cageNodes[j]->getCoords();
        const Point2D &vj2 =  cageNodes[(j+1)%nSize]->getCoords();
        a[0]  = vj2[0] - vj1[0];
        a[1]  = vj2[1] - vj1[1];

        len   =  sqrt( a[0]*a[0] + a[1]*a[1] );
        nj[0] =  a[1]/len;
        nj[1] = -a[0]/len;
        normals[j] = nj;
    }
}
//////////////////////////////////////////////////////////////////////////////////

void setGreenBary2GlobalCoords( const vector<Vertex*> &cageNodes,
                                const vector<double>  &newSegLen,
                                const vector<double>  &orgSegLen,
                                const vector<Point2D> &cageNormals,
                                Vertex *anyNode)
{
    int nSize = cageNodes.size();

    double  xsum = 0;
    double  ysum = 0;

    const vector<double> &b =  anyNode->getBaryCoords();
    assert( b.size() == 2*nSize);

    for( int i = 0; i < nSize; i++)
    {
        Point2D &pi = cageNodes[i]->getCoords();
        double si  = newSegLen[i]/orgSegLen[i];
        double phi =  b[i];
        double psi =  b[i+ nSize];
        xsum += phi*pi[0] + psi*si*cageNormals[i][0];
        ysum += phi*pi[1] + psi*si*cageNormals[i][1];
    }
    anyNode->setCoords( xsum, ysum );
}
//////////////////////////////////////////////////////////////////////////////////

void setBary2GlobalCoords( const vector<Vertex*> &cageNodes, Vertex *anyNode )
{
    double xsum = 0.0;
    double ysum = 0.0;
    int nSize = cageNodes.size();

    vector<double> b = anyNode->getBaryCoords();
    double wsum = 0.0;
    for( int i = 0; i < nSize; i++)
    {
        Point2D &p = cageNodes[i]->getCoords();
        xsum +=  b[i]*p[0];
        ysum +=  b[i]*p[1];
        wsum +=  b[i];
    }
    xsum /= wsum;
    ysum /= wsum;
    anyNode->setCoords( xsum, ysum );
}

//////////////////////////////////////////////////////////////////////////////////
int JBarycentricCoords :: getQuadShapeFunc( const vector<Point2D> &quadPoints, const Point2D &uvquad, vector<double> &shapeFunc)
{
    vector<Point2D> uvPoints(4);
    uvPoints[0][0] = -1.0;
    uvPoints[0][1] = -1.0;

    uvPoints[1][0] = -1.0;
    uvPoints[1][1] =  1.0;

    uvPoints[2][0] =  1.0;
    uvPoints[2][1] =  1.0;

    uvPoints[3][0] = -1.0;
    uvPoints[3][1] =  1.0;

    // There must be four points defining the quadrilateral shape ( in 2D xy-plane).
    assert( quadPoints.size() == 4);

    // The parametric domain of the quad is (-1,-1) to (1,1).
    assert( uvquad[0] >= -1.0 && uvquad[0] <= 1.0);
    assert( uvquad[1] >= -1.0 && uvquad[1] <= 1.0);

    // in  which triangle the UV lies ?
    Point tripoints1 =  { Point(-1.0,-1.0),  Point(1.0,-1.0), Point( 1.0,1.0)};
    Point tripoints2 =  { Point(-1.0,-1.0),  Point(1.0, 1.0), Point(-1.0,1.0)};

    Point  queryPoint = Point( uvquad[0], uvquad[1] );

    int stat;
    Point2D uvtri, xytri;

    stat = CGAL::bounded_side_2( qPoint, tripoints1, tripoints1+3, K());
    if( stat == CGAL::ON_UNBOUNDED_SIDE)  {
        stat = CGAL::bounded_side_2( qPoint, tripoints2, tripoints2+3, K());
        if( stat == CGAL::ON_UNBOUNDED_SIDE)   return 1;
        tripnts[0] = 0;
        tripnts[1] = 1;
        tripnts[2] = 2;
        TriGeometry::getUVCoords( &uvPoints[0][0], &uvPoints[1][0], &uvPoints[2][0], &uvquad[0], &uvtri[0]);
    } else {
        tripnts[0] = 0;
        tripnts[1] = 2;
        tripnts[2] = 3;
        TriGeometry::getUVCoords( &uvPoints[0][0], &uvPoints[2][0], &uvPoints[3][0], &uvquad[0], &uvtri[0]);
    }

    double trishapefunc[3];
    trishapfunc[0] = 1.0-uvtri[0] - uvtri[1];
    trishapfunc[1] = uvtri[0];
    trishapfunc[2] = uvtri[1];
    xytri[0] = 0.0;
    xytri[1] = 0.0;
    for( int i = 0; i < 3; i++)  {
        int tp = tripnts[i];
        xytri[0] += trishapefunc[i]*quadPoints[tp][0];
        xytri[1] += trishapefunc[i]*quadPoints[tp][1];
    }

    getCoords( quadPoints, xytri, shapefunc);
}
#endif


