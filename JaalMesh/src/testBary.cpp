#include "Mesh.hpp"
#include <sstream>

using namespace Jaal;

int getTriShapeFunc( const Point2D &uv, double *N)
{
    N[0] = 1 - uv[0] - uv[1];
    N[1] = uv[0];
    N[2] = uv[1];
}

//////////////////////////////////////////////////////////////////////////////////////////////////
int getQuadShapeFunc( const Point2D &uv, double *N)
{
    double u = uv[0];
    double v = uv[1];
    N[0] = 0.25*(1 - u)*(1-v);
    N[1] = 0.25*(1 + u)*(1-v);
    N[2] = 0.25*(1 + u)*(1+v);
    N[3] = 0.25*(1 - u)*(1+v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void  linearInterp(const Point2D &p0, const Point2D &p1, const Point2D &p2,  const Point2D &uv, Point2D &xy)
{
    double N[3];
    getTriShapeFunc(uv, N);

    double x[3], y[3];
    x[0] = p0[0];
    x[1] = p1[0];
    x[2] = p2[0];

    y[0] = p0[1];
    y[1] = p1[1];
    y[2] = p2[1];

    double xsum = 0.0;
    double ysum = 0.0;
    for( int i = 0; i < 3; i++) {
        xsum += N[i]*x[i];
        ysum += N[i]*y[i];
    }
    xy[0] = xsum;
    xy[1] = ysum;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
double bilinearInterp(const vector<double> &field, const Point2D &uv)
{
    double N[4];
    getQuadShapeFunc(uv, N);

    double sum = 0.0;
    for( int i = 0; i < 4; i++)
        sum += N[i]*field[i];
    return sum;
}
//////////////////////////////////////////////////////////////////////////////////////////////////

int getTriangle( const vector<Point2D> &tpoints, const Point2D &pquery)
{
    int side;
    side = FaceGeometry::getBoundedSide( &tpoints[0][0], &tpoints[1][0], &tpoints[2][0], &pquery[0] );
    if( side != GeomOrient::OUTSIDE) return 0;

    side = FaceGeometry::getBoundedSide( &tpoints[0][0], &tpoints[2][0], &tpoints[3][0], &pquery[0] );
    if( side != GeomOrient::OUTSIDE) return 1;

    cout << "Fatal error : Point location " << endl;
    cout << "QuadPoints " << endl;
    cout << tpoints[0][0] << " " << tpoints[0][1] << endl;
    cout << tpoints[1][0] << " " << tpoints[1][1] << endl;
    cout << tpoints[2][0] << " " << tpoints[2][1] << endl;
    cout << tpoints[3][0] << " " << tpoints[3][1] << endl;
    cout << "Query point " << endl;
    cout << pquery[0] << " " << pquery[1] << endl;
    exit(0);

    return 2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

double getInterpValue( const vector<Point2D> &xypoints, vector<double> &scalarfield, const Point2D &xy)
{
    int triid =  getTriangle(xypoints, xy);

    vector<Point2D>  qpoints(4);
    qpoints[0][0] = -1.0;
    qpoints[0][1] = -1.0;

    qpoints[1][0] =  1.0;
    qpoints[1][1] = -1.0;

    qpoints[2][0] =  1.00;
    qpoints[2][1] =  1.00;

    qpoints[3][0] = -1.0;
    qpoints[3][1] =  1.0;

    Point2D tuv, quv;
    if( triid == 0) {
        TriGeometry::getUVCoords( &xypoints[0][0], &xypoints[1][0], &xypoints[2][0], &xy[0], &tuv[0]);
        linearInterp(qpoints[0], qpoints[1], qpoints[2], tuv, quv);
    } else {
        TriGeometry::getUVCoords( &xypoints[0][0], &xypoints[2][0], &xypoints[3][0], &xy[0], &tuv[0]);
        linearInterp( qpoints[0], qpoints[2], qpoints[3], tuv, quv);
    }

    double val = bilinearInterp(scalarfield, quv);
    return val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
double getInterpValue( const vector<Point2D> &xypoints, vector<double> &scalarfield, const Point2D &xy, int triid)
{
    vector<Point2D>  qpoints(4);
    qpoints[0][0] = -1.0;
    qpoints[0][1] = -1.0;

    qpoints[1][0] =  1.0;
    qpoints[1][1] = -1.0;

    qpoints[2][0] =  1.00;
    qpoints[2][1] =  1.00;

    qpoints[3][0] = -1.0;
    qpoints[3][1] =  1.0;

    Point2D tuv, quv;
    if( triid == 0) {
        TriGeometry::getUVCoords( &xypoints[0][0], &xypoints[1][0], &xypoints[2][0], &xy[0], &tuv[0]);
        linearInterp(qpoints[0], qpoints[1], qpoints[2], tuv, quv);
    } else {
        TriGeometry::getUVCoords( &xypoints[0][0], &xypoints[2][0], &xypoints[3][0], &xy[0], &tuv[0]);
        linearInterp( qpoints[0], qpoints[2], qpoints[3], tuv, quv);
    }

    double val = bilinearInterp(scalarfield, quv);
    return val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

JMeshPtr getField( vector<Point2D> &quad, int trino, vector<double> &qnodevalues)
{
    ostringstream oss;
    ofstream ofile("sq.poly", ios::out);
    ofile << "3 2 0 0" << endl;

    if( trino == 0) {
        ofile << "0  " << p[0][0] << " " << p[0][1] << endl;
        ofile << "1  " << p[1][0] << " " << p[1][1] << endl;
        ofile << "2  " << p[2][0] << " " << p[2][1] << endl;
        ofile << "3  " << p[3][0] << " " << p[3][1] << endl;
        ofile << "4 1 " << endl;
        ofile << "1 0 1  1" << endl;
        ofile << "2 1 2  2" << endl;
        ofile << "3 2 3  3" << endl;
        ofile << "4 3 0  4" << endl;
        ofile << "0" << endl;
    } else {
        ofile << "0  " << p[0][0] << " " << p[0][1] << endl;
        ofile << "1  " << p[1][0] << " " << p[1][1] << endl;
        ofile << "2  " << p[2][0] << " " << p[2][1] << endl;
        ofile << "3  " << p[3][0] << " " << p[3][1] << endl;
        ofile << "4 1 " << endl;
        ofile << "1 0 1  1" << endl;
        ofile << "2 1 2  2" << endl;
        ofile << "3 2 3  3" << endl;
        ofile << "4 3 0  4" << endl;
        ofile << "0" << endl;
    }
    ofile.close();

    system( "triangle -pzq30.0a0.005 sq.poly");

    JMeshPtr mesh = JMeshIO::readFile("sq.1.ele");
    cout << "File read " << endl;

    int numNodes = mesh->getSize(0);
    for( int i = 0; i < numNodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        Point2D  p2d = vtx->getXYCoords();
        xy[0] = p3d[0];
        xy[1] = p3d[1];
        double val = getInterpValue(p, qnodevalues, xy, trino);
        vtx->setAttribute("ScalarField", val);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
    exactinit();
    Point2D xy;
    vector<double> scalarfield(4);
    cout << "Scalar field at the four corners" << endl;
    cin >> scalarfield[0] >> scalarfield[1] >> scalarfield[2] >> scalarfield[3];

    JMeshVTKExporter mexp;
    mexp.addNodeAttribute("ScalarField");

    vector<Point2D> p(4);
    p[0][0] = 0.0;
    p[0][1] = 0.0;

    p[1][0] = 1.0;
    p[1][1] = 0.0;

    p[2][0] = 0.25;
    p[2][1] = 0.25;
//  cout << "Give second Coordinate" << endl;
//   cin >> p[2][0] >> p[2][1];

    p[3][0] = 0.0;
    p[3][1] = 1.0;

    double du = 0.01;
    for( int j = 0; j < 100; j++) {
        cout << j << endl;
        p[2][0] = 0.15 + j*du;
        p[2][1] = 0.15 + j*du;
        oss << "bary" << j << ".vtk";
        mexp.writeFile(mesh, oss.str() );
    }

}
