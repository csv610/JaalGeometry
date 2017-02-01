#include <iostream>
#include <iomanip>

#include <list>
#include "TriElasticity.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                     Point_2;
typedef CGAL::Polygon_2<Kernel>             Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>  Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>     Polylist;

typedef CGAL::Box_intersection_d::Box_d<double,2> Box;
typedef CGAL::Bbox_2  Bbox;

std::multimap<int,int> boxPairs;

/////////////////////////////////////////////////////////////////////////////////////////
void BoxBoxIntersectCallback( const Box& a, const Box& b )
{
    int minid = min(a.id(), b.id());
    int maxid = max(a.id(), b.id());

    if( minid < maxid )
        boxPairs.insert(std::pair<int,int>(minid, maxid));
}

//
/////////////////////////////////////////////////////////////////////////////////////////
//
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
template<class Kernel, class Container>
void getPolyPoints(const CGAL::Polygon_with_holes_2<Kernel, Container> &pwh,
                   vector<Point2D> &points)
{
    //
    // Collect the points on the polygon from the datastructure of CGAL...
    //
    typename CGAL::Polygon_2<Kernel, Container> poly;
    typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;

    Point2D xy;
    points.clear();
    if (!pwh.is_unbounded())  {
        poly = pwh.outer_boundary();
        if( poly.size() ) {
            points.reserve(poly.size() );
            for (vit = poly.vertices_begin(); vit != poly.vertices_end(); ++vit) {
                Point_2 pnt = *vit;
                xy[0] = CGAL::to_double(pnt.x());
                xy[1] = CGAL::to_double(pnt.y());
                points.push_back( xy);
            }
        }
    }
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//

bool hasIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints)
{
    //
    // Check whether tiangle A and triangle B overlap, if they do, then return the
    // points of the polygon formed from the intersection ...
    // The points must be counter-clockwisely oriented. The outpoint points are in the
    // counter clockwise direction ....
    //
    int nsize;

    Polygon_2 P;
    nsize = aPoints.size();
    for( int i = 0; i < nsize; i++)
        P.push_back (Point_2 (aPoints[i][0], aPoints[i][1]));

    Polygon_2 Q;
    nsize = bPoints.size();
    for( int i = 0; i < nsize; i++)
        Q.push_back (Point_2 (bPoints[i][0], bPoints[i][1]));

    Polylist polylist;
    return CGAL::do_intersect(P, Q);
}
//
/////////////////////////////////////////////////////////////////////////////////////////

void getTriTriIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints,
                            vector<Point2D> &cPoints)
{
    //
    // Check whether tiangle A and triangle B overlap, if they do, then return the
    // points of the polygon formed from the intersection ...
    // The points must be counter-clockwisely oriented. The outpoint points are in the
    // counter clockwise direction ....
    //

    int nsize;
    cPoints.clear();

    Polygon_2 P;
    nsize = aPoints.size();
    for( int i = 0; i < nsize; i++)
        P.push_back (Point_2 (aPoints[i][0], aPoints[i][1]));

    Polygon_2 Q;
    nsize = bPoints.size();
    for( int i = 0; i < nsize; i++)
        Q.push_back (Point_2 (bPoints[i][0], bPoints[i][1]));

    Polylist polylist;
    CGAL::intersection (Q, P, std::back_inserter(polylist));

    Polylist::const_iterator  it;
    for (it = polylist.begin(); it != polylist.end(); ++it)
        getPolyPoints(*it, cPoints);
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//
int isInside( const vector<Point2D> &polyPoints, const Point2D &qPoint)
{
    vector<Point_2> points;
    int np = polyPoints.size();

    points.resize(np);
    for( int i = 0; i < np; i++)
        points[i] = Point_2(polyPoints[i][0], polyPoints[i][1]);
    Point_2 testPoint( qPoint[0], qPoint[1] );

    int stat;
    stat = CGAL::bounded_side_2( points.begin(), points.end(), testPoint, Kernel() );

    if( stat == CGAL::ON_BOUNDED_SIDE) return 1;

    return 0;
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//
int isInside( const double *p0, const double *p1, const double *p2, const double *qPoint)
{
    vector<Point_2> points(3);

    points[0] = Point_2( p0[0], p0[1] );
    points[1] = Point_2( p1[0], p1[1] );
    points[2] = Point_2( p2[0], p2[1] );

    Point_2 testPoint( qPoint[0], qPoint[1] );

    int stat;
    stat = CGAL::bounded_side_2( points.begin(), points.end(), testPoint, Kernel() );

    if( stat == CGAL::ON_BOUNDED_SIDE) return 1;

    return 0;
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//
void latexMatrix( const MatrixXd &A)
{
    // Print a dense matrix in Latex (useful for the documentation);
    cout << fixed;
    int numRows = A.rows();
    int numCols = A.cols();
    cout << setprecision(2);
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
//
TriElasticity::TriElasticity()
{
    mesh   = nullptr;
    problem = PLAIN_STRESS;

    E = 2.0E11;
    nu = 0.33;

    F[0] = 0.0;
    F[1] = 0.0;
    solverMethod = 2;
    dofPerNode   = 2;
    numDOF       = 0;
    debugStep    = 0;
    verbose      = 1;
    tangleDetect = 1;

    setOrder(1);
    computeDMatrix();
}

///////////////////////////////////////////////////////////////////////////////

void TriElasticity :: computeDMatrix()
{
    double C;

    if( problem == PLAIN_STRAIN ) {
        C =   E/((1+nu)*(1.0-2*nu));
        D(0,0) =   C*(1-nu);
        D(0,1) =   C*nu;
        D(0,2) =   0.0;

        D(1,0) =   C*nu;
        D(1,1) =   C*(1.0-nu);
        D(1,2) =   0.0;

        D(2,0) =   0.0;
        D(2,1) =   0.0;
        D(2,2) =   0.5*C*(1.0-2.0*nu);
    }

    if( problem == PLAIN_STRESS) {
        C =   E/(1-nu*nu);
        D(0,0) =   C;
        D(0,1) =   C*nu;
        D(0,2) =   0.0;

        D(1,0) =   C*nu;
        D(1,1) =   C;
        D(1,2) =   0.0;

        D(2,0) =   0.0;
        D(2,1) =   0.0;
        D(2,2) =   0.5*C*(1.0-nu);
    }

    if( debugStep == 1) {
        cout << "D Matrix " << endl;
        latexMatrix(D);
    }
}

/////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: setOrder(int order)
{
    shapeOrder = order;

    femshape.setOrder(order);
    nodesPerEdge = femshape.getNumNodesPerEdge();
    nodesPerFace = femshape.getNumNodesPerFace();
    dofPerEdge       = dofPerNode*nodesPerEdge;
    dofPerFace       = dofPerNode*nodesPerFace;
    numFaceGaussPnts = femshape.getNumFaceGaussPnts();
    numEdgeGaussPnts = femshape.getNumEdgeGaussPnts();

    femshape.getGaussPoints(exi, eweight);
    femshape.getGaussPoints(fxi, feta, fweight);

}

////////////////////////////////////////////////////////////////////////////////

bool TriElasticity :: isTangled()
{
    if( mesh == nullptr ) return 0;

    meshtangle.setMesh(mesh);
    if( meshtangle.getNumInvertedElements() == 0) return 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: setMesh( Mesh *m)
{
    mesh = m;
    tangledK = 0;
    initMesh();
}

////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: getDisplacements( vector<double> &unode, vector<double> &vnode)
{

    unode.clear();
    vnode.clear();
    if( mesh == nullptr || u.empty() || v.empty() ) return;

    int numnodes = mesh->getSize(0);
    unode.resize(numnodes);
    vnode.resize(numnodes);
    for( int i = 0; i < numnodes; i++) {
        unode[i] = u[i];
        vnode[i] = v[i];
    }
}

////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: getPoints( const Face *face, vector<Point2D> &points)
{
    //
    // Get the points of the face in counter clockwise direction. If the mesh
    // is tangled then some of the elements will have negative orientation,
    // and having positive orientation is prerequisite for overlapping detection.
    //
    short int fsign;
    int err = face->getAttribute("Orient", fsign);
    if( err ) {
        fsign  = FaceGeometry::getOrientation2D(face);
        (const_cast<Face*>(face))->setAttribute("Orient", fsign);
    }

    int np = face->getSize(0);
    points.resize(np);

    if( fsign > 0) {
        for( int i = 0; i < np; i++) {
            const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
            points[i][0] = xyz[0];
            points[i][1] = xyz[1];
        }
        return;
    }

    for( int i = 0; i < np; i++) {
        const Point3D &xyz = face->getNodeAt(np-i-1)->getXYZCoords();
        points[i][0] = xyz[0];
        points[i][1] = xyz[1];
    }

    return;

}
////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: getNodes( const Face *face, vector<int> &nodes)
{
    // Give the nodes of elements including "Steiner nodes" on edges ...
    nodes.clear();

    for( int i = 0; i < face->getSize(0); i++)
        nodes.push_back(face->getNodeAt(i)->getID());

    if( shapeOrder == FEMShape::QUADRATIC) {
        int vid;
        for( int i = 0; i < face->getSize(1); i++) {
            Edge *edge = face->getEdgeAt(i);
            edge->getAttribute("QuadNode", vid );
            nodes.push_back(vid);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: setXForce( JEdgeSequence &boundaryEdges, double force)
{
    // Specify the force on "geometric edge". The mesh edges on a geometric edge
    // gets uniform force.

    size_t nSize = boundaryEdges.size();

    if( nSize < 1) {
        cout << "Info: X-force boundary edges are empty " << endl;
        return;
    }

    double len = 0.0;
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        len += edge->getLength();
    }

    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        edge->setAttribute("Xforce", force/len);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: setYForce( JEdgeSequence &boundaryEdges, double force)
{
    // Specify the force on "geometric edge". The mesh edges on a geometric edge
    // gets uniform force.

    size_t nSize = boundaryEdges.size();
    if( nSize < 1) {
        cout << "Info: X-force boundary edges are empty " << endl;
        return;
    }

    double len = 0.0;
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        len += edge->getLength();
    }

    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        edge->setAttribute("Yforce", force/len);
    }
}

/////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: fixX( const Vertex *vtx, double u)
{
    int id = vtx->getID();
    fixMap[2*id] = u;
}
/////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: fixY( const Vertex *vtx, double v)
{
    int id = vtx->getID();
    fixMap[2*id+1] = v;
}
/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: fixX( JEdgeSequence &boundaryEdges, double u)
{
    // Fix the "U" component of a edge.
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        edge->setAttribute("Xfixed", u);
    }
}
/////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: fixY( JEdgeSequence &boundaryEdges, double v)
{
    // Fix the "V" component of a edge.
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        edge->setAttribute("Yfixed", v);
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: fixEdges( JEdgeSequence &boundaryEdges, double uv)
{
    // Specify which edges are completely restricted in both the directions...
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        Edge *edge = boundaryEdges[i];
        edge->setAttribute("Xfixed", uv);
        edge->setAttribute("Yfixed", uv);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: BMatrix(const MatrixXd &gradUV, const Matrix2d &invJ, MatrixXd &Bg)
{
    //------------------------------------------------------------------------------
    // Calculate :
    //   dNx_1   0     dNx_2   0         .....   dNx_p  0
    //   0     dNy_1   0       dNy_2     ....    0      dNy_p
    //   dNy_1 dNx_1   dNy_2   dNx_2     .....   dNy_p  dNx_p
    //------------------------------------------------------------------------------

    assert( gradUV.rows() == 2);
    assert( gradUV.cols() == nodesPerFace);

    // Calculate the "B" Matrix at the Guass point on a face ...
    double j00 = invJ.coeffRef(0,0);
    double j01 = invJ.coeffRef(0,1);
    double j10 = invJ.coeffRef(1,0);
    double j11 = invJ.coeffRef(1,1);

    Bg.setZero();
    for( int i= 0; i < nodesPerFace; i++) {
        double dNdU = gradUV.coeffRef(0,i);   // Partial N over Partial xi
        double dNdV = gradUV.coeffRef(1,i);   // Partial N over Partial eta
        double dNdX = j00*dNdU + j01*dNdV;    // partial N over Partial X
        double dNdY = j10*dNdU + j11*dNdV;    // partial N over Partial Y
        Bg.coeffRef(0,2*i)   = dNdX;
        Bg.coeffRef(1,2*i+1) = dNdY;
        Bg.coeffRef(2,2*i)   = dNdY;
        Bg.coeffRef(2,2*i+1) = dNdX;
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: getBMatrix(const double *uv, const Matrix2d &invJ, MatrixXd &B)
{
    // Calculate the "B" matrix at any point within the trinagle. We pass the inverse
    // of the jacobian of the face, where we need to find the matrix "B". Notice that
    // we are pssing Jacobian of the face and not the face object, because the Jacobian
    // on the linear triangle is always constant. If we had passed, then we would have
    // required calculating inverse jacobian at every "UV", Thereofore, passing the
    // inverse is big optimization step...

    MatrixXd gradN;
    femshape.getShapeDeriv( uv[0], uv[1], gradN);
    BMatrix(gradN, invJ, B);
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: getInvJacobian(const Face *face, Matrix2d &invJ, double &dJ)
{
    // Give the inverse Jacobain( 2x2) matrix and its determinant...

    if( face->getSize(0) == 3)  {
        const Point3D &p1 = face->getNodeAt(0)->getXYZCoords();
        const Point3D &p2 = face->getNodeAt(1)->getXYZCoords();
        const Point3D &p3 = face->getNodeAt(2)->getXYZCoords();

        invJ(0,0) =   p3[1] - p1[1];
        invJ(0,1) =   p1[1] - p2[1];
        invJ(1,0) =   p1[0] - p3[0];
        invJ(1,1) =   p2[0] - p1[0];
        dJ   = invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0);
        invJ  = invJ/dJ;
        return;
    }
    cout << "Error: element not supported yet " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: integrate(const Face *face, MatrixXd &Ke, vector<double> &fe)
{
    //
    // Integrate over the face and return Element stiffness matrix "Ke" and "f" vector.
    // The size of matrix and f depends on the number of nodes.

    double dJ;
    Matrix2d invJ;

    getInvJacobian(face, invJ, dJ);

    Ke.setZero();
    MatrixXd B;
    B.resize(3, dofPerFace);

    // Calculate the integral of (B'D*B*dA) on each face. The function
    // value is evaluated at the "Gauss Points". In our case, the
    // Inverse Jacobian is constant.
    for( int ig = 0; ig < numFaceGaussPnts; ig++) {
        BMatrix(gradUV[ig], invJ, B);
        double C = fweight[ig]*fabs(dJ);
        Ke = Ke + C*B.transpose()*D*B;
    }

    // Calculate the "f" vector, which depends on the body force ...
    fe.resize(dofPerFace);
    for( int i = 0; i < dofPerFace; i++)
        fe[i] = 0.0;

    // In classical method, we need to take +ve Jacobian. But when we have
    // tangled mesh, Jacobian will indicate sign of the element.
    if( !tangledK ) dJ = fabs(dJ);

    double bx = bodyForce[0];
    double by = bodyForce[1];
    for( int ip = 0; ip < nodesPerFace; ip++) {
        for( int ig = 0; ig < numFaceGaussPnts; ig++) {
            fe[2*ip]   += fweight[ig]*NFace[ig][ip]*bx*dJ;
            fe[2*ip+1] += fweight[ig]*NFace[ig][ip]*by*dJ;
        }
    }

    if( debugStep == 2 ) {
        cout << "Calculating local stiffness matrix " << endl;
        latexMatrix(Ke);
        latexVector(fe);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: getAssemblyIndices(const Face *face, vector<int> &globalPos)
{
    int vid;
    //
    // Find the global pos of the nodes on a face. We follow the convention that
    // dof are specified as (u_1, v_1, u_2, v_2 ..... u_n, v_n)
    //
    // Also, first we enumerate the nodes of the elements and then the nodes on
    // the edges ( which are higher order nodes)..
    //
    globalPos.clear();
    globalPos.reserve(dofPerFace);

    for( int i = 0; i < face->getSize(0); i++) {
        vid = face->getNodeAt(i)->getID();
        for( int j = 0; j < dofPerNode; j++)
            globalPos.push_back(2*vid+j);
    }

    if( shapeOrder == FEMShape::QUADRATIC) {
        for( int i = 0; i < nodesPerEdge; i++) {
            Edge *edge = face->getEdgeAt(i);
            edge->getAttribute("QuadNode", vid);
            for( int j = 0; j < dofPerNode; j++)
                globalPos.push_back(2*vid+j);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: getAssemblyIndices(const Edge *edge, vector<int> &globalPos)
{
    int vid;
    // Find the global pos of the nodes on a edge. We follow the convention that
    // dof are specified as (u_1, v_1, u_2, v_2 ..... u_n, v_n)
    //
    // Also, first we enumerate the end nodes of a edge  and then any higher
    // order nodes on the edge..

    globalPos.clear();
    globalPos.reserve(dofPerEdge);

    vid = edge->getNodeAt(0)->getID();
    for( int j = 0; j < dofPerNode; j++)
        globalPos.push_back(2*vid+j);

    vid = edge->getNodeAt(1)->getID();
    for( int j = 0; j < dofPerNode; j++)
        globalPos.push_back(2*vid+j);

    if( shapeOrder == FEMShape::QUADRATIC ) {
        edge->getAttribute("QuadNode", vid);
        for( int j = 0; j < dofPerNode; j++)
            globalPos.push_back(2*vid+j);
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: assembleK()
{

    if( verbose)
        cout << "Assemble classical stiffness matrix starts ... " << endl;

    if( mesh == nullptr) return;

    // Calculate the stiffness matrix of each face and assemble into a global matrix..
    int numFaces = mesh->getSize(2);

    // Since all the elements are of same shape the Guass point location and weights
    // are constant, so calcualte them once..
    femshape.getGaussPoints(fxi, feta, fweight);

    // Calculate the shape function (N) and its derivative in UV space at each Gauss
    // point...
    NFace.resize(numFaceGaussPnts);
    gradUV.resize(numFaceGaussPnts);
    for( int i = 0; i < numFaceGaussPnts; i++) {
        femshape.getShapeFunc(  fxi[i], feta[i], NFace[i]  );
        femshape.getShapeDeriv( fxi[i], feta[i], gradUV[i] );
    }

    // Global matrix ...
    K.resize(numDOF, numDOF);

    // Element matrix ( local );
    MatrixXd KElem(dofPerFace, dofPerFace);
    vector<int> globalPos;
    vector<double> fElem(dofPerFace);

    //
    // for each face, calculate the Ke and f and assemble them into global
    // matrix and vector respectively....
    //
    for( int i = 0; i < numFaces; i++) {
        Face *face = mesh->getFaceAt(i);
        integrate(face, KElem, fElem);
        getAssemblyIndices(face, globalPos);
        for( int j = 0; j < dofPerFace; j++) {
            int irow = globalPos[j];
            for( int k = 0; k < dofPerFace; k++) {
                int icol = globalPos[k];
                K.coeffRef(irow,icol) += KElem(j,k);
            }
            f[irow] += fElem[j];
        }
    }

    if( debugStep == 3 ) {
        cout << "Global stiffness matrix " << endl;
        latexMatrix(K);
        cout << "Force vector " << endl;
        latexVector(f);
    }

    if( verbose)
        cout << "Stiffness matrix assembled, " << endl;

}

////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: integrate( Edge *edge, double force, vector<double> &fb)
{
    // Integrate "force" vector which is on the right hand side of the system...
    // Remember that Guass points are assigned from the left point to the right
    // point in the interval (-1, 1). But as per convention, the end nodes are
    // assembled first and then the nodes on the edge, therefore, we need to
    // make adjustment. For this reason, we first use "ftmp" and then assemble
    // it in proper order in the vector "fb"...

    vector<double> ftmp(nodesPerEdge);
    for( int i = 0; i < nodesPerEdge; i++) ftmp[i] = 0.0;

    // For linear edges, the dJ can be calculated simply....
    double len = edge->getLength();
    double dJ = 0.5*len;

    // For each gauss point, evaluate the function f ( which is constant in
    // ou modelling).
    for( int ig = 0; ig < numEdgeGaussPnts; ig++) {
        for( int ip = 0; ip < nodesPerEdge; ip++)
            ftmp[ip] +=  eweight[ig]*NEdge[ig][ip]*force*dJ;
    }

    // Rearrange now. First, two nodes of the edge and then higher order
    // nodes...
    fb.resize( nodesPerEdge);

    fb[0] = ftmp[0];
    fb[1] = ftmp[nodesPerEdge-1];
    for( int i = 1; i < nodesPerEdge-1; i++)
        fb[i+1] = ftmp[i];
}

////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: applyNeumannBC()
{
    if( verbose )
        cout << "Info: Assplying force boundary condition " << endl;

    femshape.getGaussPoints(exi, eweight);

    NEdge.resize(numEdgeGaussPnts);
    gradU.resize(numEdgeGaussPnts);
    for( int i = 0; i < numEdgeGaussPnts; i++)
        femshape.getShapeFunc(  exi[i], NEdge[i]  );

    vector<double>  fb;
    vector<int>     globalPos;

    JEdgeSequence boundedges;

    int nCount = 0;
    double val = 0.0;
    mesh->getEntities("Xforce", boundedges);
    nCount += boundedges.size();
    for( Edge *edge : boundedges) {
        edge->getAttribute("Xforce", val);
        integrate(edge, val, fb);
        getAssemblyIndices(edge, globalPos);
        for( size_t i = 0; i < fb.size(); i++) {
            int gid  = globalPos[2*i];
            f[gid]  += fb[i];
        }
    }

    mesh->getEntities("Yforce", boundedges);
    nCount += boundedges.size();
    for( Edge *edge : boundedges) {
        edge->getAttribute("Yforce", val);
        integrate( edge, val, fb);
        getAssemblyIndices(edge, globalPos);
        for( size_t i = 0; i < fb.size(); i++) {
            int gid = globalPos[2*i+1];
            f[gid] += fb[i];
        }
    }

    if( debugStep == 4) {
        cout << " F-Vector: After Neumann BC" << endl;
        latexVector(f);
    }

    if( nCount == 0)
        cout << "Warning:  There should be atleast one force specified on the boundary" << endl;

    if( verbose )
        cout << "Info: Force boundary condition applied" << endl;
}
////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: applyDirichletBC()
{
    int    nodeID;
//  double val = 0.0;
    size_t numedges = mesh->getSize(1);

    for( size_t i = 0; i < numedges; i++) {
        Edge *edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("Xfixed") ) {
            nodeID = edge->getNodeAt(0)->getID();
            fixMap[2*nodeID] = 0.0;
            nodeID = edge->getNodeAt(1)->getID();
            fixMap[2*nodeID] = 0.0;
            if( shapeOrder == FEMShape::QUADRATIC) {
                int err = edge->getAttribute("QuadNode", nodeID);
                if( !err) fixMap[2*nodeID] = 0.0;
            }
        }
        if( edge->hasAttribute("Yfixed") ) {
            nodeID = edge->getNodeAt(0)->getID();
            fixMap[2*nodeID+1] = 0.0;
            nodeID = edge->getNodeAt(1)->getID();
            fixMap[2*nodeID+1] = 0.0;
            if( shapeOrder == FEMShape::QUADRATIC) {
                int err = edge->getAttribute("QuadNode", nodeID);
                if( !err) fixMap[2*nodeID+1] = 0.0;
            }
        }
    }

    if( fixMap.empty() )
        cout << "Warning: There are no restricted nodes in the mesh " << endl;
}

////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: assembleBC()
{
    // Apply boundary force conditions ...
    applyNeumannBC();

    // Specify the fixed boundary condition....
    applyDirichletBC();

}

////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: genQuadraticNodes()
{
    if( verbose )
        cout << "Info: Generating Quadratic mesh ... " << endl;
    // Acutally, I would not this as below, but this is how matlab finds uniques
    // edges and assign new nodes..
    int numNodes = mesh->getSize(0);
    int numEdges = mesh->getSize(1);
    int numFaces = mesh->getSize(2);

    int index = numNodes;
    mesh->delete_edge_attribute("QuadNode");

    /*
        JEdgeSequence edges;
        mesh->getTopology()->getLexicographic(edges);
        for( int i = 0; i < numEdges; i++)
            edges[i]->setAttribute("QuadNode", index++);
    */

    for( int iface =  0; iface < numFaces; iface++) {
        Face *face = mesh->getFaceAt(iface);
        for( int iedge = 0; iedge < face->getSize(1); iedge++) {
            Edge *edge = face->getEdgeAt(iedge);
            if( !edge->hasAttribute("QuadNode") ) {
                edge->setAttribute("QuadNode", index++);
            }
        }
    }

    if( verbose )
        cout << "Info: Quadratic nodes insertted" << endl;
}

////////////////////////////////////////////////////////////////////////////////////

int TriElasticity :: getUVCoords( const Face *face, const Point2D &xy, Point2D &uv)
{
    const Point3D &p0 = face->getNodeAt(0)->getXYZCoords();
    const Point3D &p1 = face->getNodeAt(1)->getXYZCoords();
    const Point3D &p2 = face->getNodeAt(2)->getXYZCoords();

    if( !isInside( &p0[0], &p1[0], &p2[0], &xy[0])  ) {
        cout << "Fatal Error: Vertex not inside face " << endl;
        exit(0);
    }

    int err = TriGeometry::getUVCoords( &p0[0], &p1[0], &p2[0], &xy[0], &uv[0] );
    assert( uv[0] >=0.0 && uv[0] <= 1.0);
    assert( uv[1] >=0.0 && uv[1] <= 1.0);

    return err;
}

////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: integrate( const Face *face1, const Face *face2,
                                 const vector<Point2D> &triPnts, MatrixXd &Ke)
{
// CSV
    assert( face1 != face2);
    assert( triPnts.size() == 3);

    double  dJ;
    Matrix2d  invJ1, invJ2, invJ;

    getInvJacobian(face1, invJ1, dJ);
    getInvJacobian(face2, invJ2, dJ);

    // Find the area of subtriangle of the region...
    invJ(0,0)  = triPnts[2][1] - triPnts[0][1];
    invJ(0,1)  = triPnts[0][1] - triPnts[1][1];
    invJ(1,0)  = triPnts[0][0] - triPnts[2][0];
    invJ(1,1)  = triPnts[1][0] - triPnts[0][0];
    dJ   = invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0);
    invJ = invJ/dJ;

    Ke.resize( dofPerFace, dofPerFace);
    Ke.setZero();

    if( fabs(dJ) < 1.0E-10) return;

    // Determine the orientation of the two faces ...
    double sign1 = FaceGeometry::getOrientation2D( face1 );
    double sign2 = FaceGeometry::getOrientation2D( face2 );

    assert(dJ > 0.0);

    // Determine the XY Coordinates of the gauss points in the triangle
    vector<Point2D> xy;
    femshape.getFaceGaussXY(triPnts, xy);

    Point2D uv;
    MatrixXd gradN, Bi, Bj;

    Bi.resize(3, dofPerFace);
    Bj.resize(3, dofPerFace);

    // Calculate the (u,v) of each gauss points in the two triangles
    // and calculate the Bi and Bj  Matrix..
    double coeff;
    for ( int ig = 0; ig < numFaceGaussPnts; ig++) {
        getUVCoords( face1, xy[ig], uv );
        femshape.getShapeDeriv( uv[0], uv[1], gradN);
        BMatrix( gradN, invJ1, Bi);

        getUVCoords( face2, xy[ig], uv );
        femshape.getShapeDeriv( uv[0], uv[1], gradN);
        BMatrix( gradN, invJ2, Bj);

        coeff = sign1*sign2*fweight[ig]*fabs(dJ);
        Ke += coeff*Bi.transpose()*D*Bj;
    }
}

////////////////////////////////////////////////////////////////////////////////////

int TriElasticity :: integrate( const Face *face1, const Face *face2, MatrixXd &Kp)
{
    assert( face1 != face2);

    Kp.setZero();

    // Check for intersection of the two faces ..

    vector<Point2D> faPoints, fbPoints, polyPoints;
    getPoints( face1, faPoints);
    getPoints( face2, fbPoints);
    assert( faPoints.size() == 3);
    assert( fbPoints.size() == 3);
    getTriTriIntersection( faPoints, fbPoints, polyPoints);

    int np = polyPoints.size();

    // If there was no intersection, no polygon will be formed...
    if( np < 3 ) return 1;

    // If the polygonal region is a triangle, no need to subdivide
    if( np == 3 ) {
        integrate(face1, face2, polyPoints, Kp);
        return 0;
    }

    // Since the intersection of two convex faces ( here triangles) is a also
    // a convex, we can simply triangulate the region by inserting a new vertex
    // at the center.

    Point2D center;
    center[0] = 0.0;
    center[1] = 0.0;
    for( int i = 0; i < np ; i++) {
        center[0] += polyPoints[i][0];
        center[1] += polyPoints[i][1];
    }
    center[0] /= (double)np;
    center[1] /= (double)np;

    vector<Point2D> nodeCoords(3);

    // Integrate in the region..

    MatrixXd Ke(dofPerFace, dofPerFace);
    for( int i = 0; i < np; i++) {
        nodeCoords[0] = center;
        nodeCoords[1] = polyPoints[i];
        nodeCoords[2] = polyPoints[(i+1)%np];
        integrate(face1, face2, nodeCoords, Ke);
        Kp += Ke;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void TriElasticity :: assembleTangledK( const Face *face1, const Face *face2)
{
    assert( face1 != face2);

    MatrixXd Kp(dofPerFace, dofPerFace);

    Kp.setZero();

    assert(Kp.rows() == dofPerFace);
    assert(Kp.cols() == dofPerFace);

    int err = integrate(face1, face2, Kp);

    // If overlapping indeed happened, we should assemble it in a global matrix..

    if( !err) {
        vector<int> adof;
        getAssemblyIndices(face1, adof);
        assert( adof.size() == dofPerFace );

        vector<int> bdof;
        getAssemblyIndices(face2, bdof);
        assert( bdof.size() == dofPerFace );

        for( int i = 0; i < dofPerFace; i++) {
            for( int j = 0; j < dofPerFace; j++) {
                Ktangle.coeffRef(adof[i],bdof[j]) += Kp(i,j);
                Ktangle.coeffRef(bdof[i],adof[j]) += Kp(j,i);
            }
        }
        tangledK = 1;
        numTangledPairs++;
    }
}
////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: assembleTangledK()
{
    tangledK = 0;
    numTangledPairs = 0;

    Ktangle.resize(numDOF, numDOF);
    Ktangle.setZero();

    if( !tangleDetect ) return;

    if( verbose)
        cout << "Detecting negative elements ... " << endl;

    // Check if there are negative faces in the mesh. If there are none, then
    // overlapping can not happen and Ktangle will be zero.

    int numfaces = mesh->getSize(2);
    int nCount = 0;
    for( int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        short int fsign  = FaceGeometry::getOrientation2D(face);
        face->setAttribute("Orient", fsign);
        if( fsign < 0.0) nCount++;
    }

    // If there is no tangle, no need to keep the attribute ....
    if( nCount == 0) {
        mesh->delete_face_attribute("Orient");
        return;
    }

    if( verbose)
        cout << "Detecting mesh tangle and assembling tangled K ... " << endl;

    boxPairs.clear();
    vector<Box> triBoxes;
    triBoxes.reserve(numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        assert( face->getID() == i);
        const BoundingBox  &bx = face->getBoundingBox();
        const Point3D  &pmin   = bx.getLower();
        const Point3D  &pmax   = bx.getUpper();
        triBoxes.push_back(Bbox(pmin[0], pmin[1], pmax[0], pmax[1]));
        assert( triBoxes[i].id() == i);
    }

    CGAL::box_self_intersection_d( triBoxes.begin(), triBoxes.end(), BoxBoxIntersectCallback);

    for( const std::pair<int,int> &boxPair: boxPairs) {
        Face *face1 = mesh->getFaceAt(boxPair.first);
        Face *face2 = mesh->getFaceAt(boxPair.second);
        assembleTangledK(face1, face2);
    }


    // Now search for overlapping faces. This is a brute-force method,
    // Checking every element with all the subsequence faces. Not an
    // efficient method, but ok for small size problems.

    /*
        for( size_t i = 0; i < numfaces; i++) {
            cout << i << " " << numfaces << endl;
            Face *face1 = mesh->getFaceAt(i);
            for( size_t j = i+1; j < numfaces; j++) {
                Face *face2 = mesh->getFaceAt(j);
                assembleTangledK(face1, face2);
            }
        }
    */

    if( verbose)
        cout << "Info: #of Tangled Pairs " << numTangledPairs << endl;
}
////////////////////////////////////////////////////////////////////////////////////
void TriElasticity :: solve()
{
    if( mesh == nullptr ) return;
    // The D matrix will remain constant for the problem, so calculate it once...
    computeDMatrix();

    int numNodes = mesh->getSize(0);
    int numEdges = mesh->getSize(1);
    assert( numNodes );
    assert( numEdges );

    numDOF   = dofPerNode*(numNodes + numEdges*(nodesPerEdge-2));

    f.resize(numDOF);

    for( int i = 0; i < numDOF; i++) f[i] = 0.0;

    if( shapeOrder == FEMShape::QUADRATIC) genQuadraticNodes();

//  First calculate the Tangled K, so that we know whether we need to modify the
//  right hand side vector "f" in the classical formulation.
    assembleTangledK();

//  Classical assembly for tangle free mesh...
    assembleK();

    //
    // If the mesh is tangled, we will modify the K Matrix, Only the
    // internal edges are tangled and no edge must cross the boundary
    // edges ..
    //
    K = K + Ktangle;

    assembleBC();

    // If there are large number of fixed value of "U" then we should reduce the
    // System for two reasons:
    //
    // (1) Run time performanec will increase with the reduced matrix size.
    // (2) Condition number will improve ...
    //
    freduced = f;
    reduceSystem(K, freduced, fixMap, freeDof);

    // Solve the Lineay sytem either with a direct or indirect solver...
    solveLinearSystem();
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity :: initMesh()
{
    if( mesh == nullptr ) return;

    if( mesh->getTopology()->getDimension() != 2) {
        cout << "Error: Tri-Elasticity is only for triangle mesh " << endl;
        return;
    }

    if( mesh->getTopology()->getElementsType(2) != 3) {
        cout << "Warning: The mesh must consists of only triangle elements: Mixed element not supported" << endl;
        return;
    }

    mesh->pruneAll();

    mesh->getTopology()->search_boundary();
    mesh->getTopology()->collect_edges();

    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity::reduceSystem( Eigen::SparseMatrix<double> &A, vector<double> &b,
                                  map<int,double> &fixMap, vector<int> &freeDof)
{
    if( verbose )
        cout << "Info: Reducing the linear system ... " << endl;

    int nrows = A.rows();
    int ncols = A.cols();
    assert( b.size() == size_t(nrows));

    freeDof.resize(nrows);
    for(int i = 0; i < nrows; i++) freeDof[i] = i;

    if( fixMap.empty() ) return;

    for( int i = 0; i < nrows; i++) {
        for( const pair<int,double> &keyVal: fixMap) {
            int j = keyVal.first;
            double val = keyVal.second;
            b[i] = b[i] - A.coeff(i,j)*val;
            A.coeffRef(i,j) = 0.0;
        }
    }

    for( const pair<int,double> &keyVal: fixMap) {
        int i = keyVal.first;
        for( int j = 0; j < ncols; j++) {
            if( A.coeff(i,j) != 0.0) A.coeffRef(i,j) = 0.0;
        }
        A.coeffRef(i,i) = 1.0;
        b[i] = keyVal.second;
    }
    A.prune(0.0);

    if( debugStep == 5) {
        cout << "Linear System after Dirichlet boundary conditions" << endl;
        latexMatrix(A);
        latexVector(b);
    }

    vector<int> fullDof(nrows);
    for(int i = 0; i < nrows; i++) fullDof[i] = i;

    vector<int> fixDof;
    fixDof.reserve(nrows);
    for( const pair<int,double> &keyVal: fixMap)
        fixDof.push_back(keyVal.first);

    freeDof.clear();
    boost::set_difference(fullDof, fixDof, back_inserter(freeDof));

    vector<int> permrow;
    permrow.reserve(nrows);
    permrow.insert( permrow.end(), freeDof.begin(), freeDof.end() );
    permrow.insert( permrow.end(), fixDof.begin(),  fixDof.end()  );

    Eigen::SparseMatrix<double> permuteMatrix(nrows,nrows);
    for( int i = 0; i < nrows; i++)
        permuteMatrix.coeffRef(i, permrow[i] ) = 1.0;

    A = permuteMatrix*A*permuteMatrix.transpose();

    int nfree = freeDof.size();
    vector<double> bb(nfree);
    for( int i = 0; i < nfree; i++)
        bb[i] = b[permrow[i]];
    b = bb;

    Eigen::SparseMatrix<double> Am(nfree,nfree);
    for( int i = 0; i < nfree; i++) {
        for( int j = 0; j < nfree; j++) {
            double val = A.coeff(i,j);
            if( val != 0.0) Am.insert(i,j) = val;
        }
    }
    A  = Am;

    if( debugStep == 6) {
        cout << "Reduced linear system after Dirichlet boundary conditions" << endl;
        latexMatrix(A);
        latexVector(b);
    }

    if( verbose )
        cout << "Info: Reducing system ready" << endl;

}
/////////////////////////////////////////////////////////////////////////////////////
int TriElasticity::solveLinearSystem()
{
    if( verbose )
        cout << "Info: Solving linear system ... " << endl;

    cout << "Info: #Nonzeros in Sparse Matrix " << K.nonZeros() << endl;

    // For small size problems, we can use Direct solver ....
    vector<double> ureduced;
    int err =  LinearSystem::solve( K, freduced, ureduced);

    cout << "Solved Vector" << endl;
    latexVector(ureduced);

    // Reassemble the vector "U" in full order (i.e. reduced + fixed value)..
    vector<double> ufull(numDOF);
    for( size_t i = 0; i < freeDof.size(); i++)
        ufull[freeDof[i]] = ureduced[i];

    for( const pair<int,double> &keyVal: fixMap)
        ufull[keyVal.first] = keyVal.second;

    // Sperate "u" and "v" into indepedent vectors, This is what matters
    // in the final stage...
    int numNodes = mesh->getSize(0) + (nodesPerEdge-2)*mesh->getSize(1);
    u.resize( numNodes);
    v.resize( numNodes);
    for( int i = 0; i < numNodes; i++) {
        u[i] = ufull[2*i];
        v[i] = ufull[2*i+1];
    }

    if( debugStep == 7) {
        cout << "Displacement vector " << endl;
        latexVector(ufull);
    }

    if( verbose )
        cout << "Info: Linear system solved ... " << endl;

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////

void TriElasticity:: computeStress( const Face *face)
{
    // Once the displacement are know, we can calculate the strain and stress
    // in each face. For Von-Mises stress is the maximum stress at any gauss point..
    int fid = face->getID();

    double   dJ;
    Matrix2d invJ;
    getInvJacobian(face, invJ, dJ);

    VectorXd un(nodesPerFace), vn(nodesPerFace), gradU(2), gradV(2);

    vector<int> nodes;
    getNodes(face, nodes);
    for( int i = 0; i < nodesPerFace; i++) {
        un[i] = u[nodes[i]];
        vn[i] = v[nodes[i]];
    }

    MatrixXd B;
    Matrix2d strain, stress;

    for( int ig = 0; ig < numFaceGaussPnts; ig++) {
        B = invJ*gradUV[ig];
        gradU = B*un;
        gradV = B*vn;
        double ux = gradU[0];
        double uy = gradU[1];
        double vx = gradV[0];
        double vy = gradV[1];
        strain(0,0) =  ux;
        strain(0,1) = 0.5*(uy+vx);
        strain(1,0) = 0.5*(uy+vx);
        strain(1,1) = vy;
        double sxx = D(0,0)*strain(0,0) + D(0,1)*strain(1,1);
        double syy = D(1,0)*strain(0,0) + D(1,1)*strain(1,1);
        double sxy = 2.0*D(2,2)*strain(0,1);
        double maxVal = sqrt(sxx*sxx + syy*syy - sxx*syy + 3.0*sxy*sxy);
        if( maxVal > faceVonMises[fid] ) {
            faceStrain[fid] = strain;
            stress(0,0) = sxx;
            stress(0,1) = sxy;
            stress(1,0) = sxy;
            stress(1,1) = syy;
            faceStress[fid] = stress;
            faceVonMises[fid] = maxVal;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void TriElasticity:: computeStress()
{
    if(mesh == nullptr) return;

    if( verbose )
        cout << "Info: Computing face stress " << endl;

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    nodeVonMises.resize(numnodes);
    for( size_t i = 0; i < numnodes; i++)
        nodeVonMises[i] = 0.0;

    faceVonMises.resize(numfaces);
    for( size_t i = 0; i < numfaces; i++)
        faceVonMises[i] = 0.0;

    faceStress.resize(numfaces);
    faceStrain.resize(numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        faceStress[i].setZero();
        faceStrain[i].setZero();
    }

    // Calculate the strain and stress at every face ....
    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        computeStress(face);
    }

    if( verbose )
        cout << "Info: Computing nodes stress " << endl;


    // Calculate the Von-Mises stress on each nodes of the
    // mesh. ( Presently only on the origial nodes and not
    // on any higher order nodes ...)

    mesh->buildRelations(0,2);
    JFaceSequence faceneighs;
    for( size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->getRelations(faceneighs);
        double sum = 0.0;
        int    numneighs = faceneighs.size();
        for( int j = 0; j < numneighs; j++) {
            int fid = faceneighs[j]->getID();
            sum += faceVonMises[fid];
        }
        nodeVonMises[i] = sum/(double)numneighs;
    }

    if( debugStep == 8) {
        cout << "VOnMises faces " << endl;
        latexVector( faceVonMises );
        cout << "VOnMises nodes " << endl;
        latexVector( nodeVonMises );
    }

    if( verbose)
        cout <<"Info: Stress calculations completed: All done " << endl;
}

/////////////////////////////////////////////////////////////////////////////////////
void TriElasticity:: saveAs( const string &filename, const string &var)
{
    if( mesh == nullptr ) return;

    if( u.empty() || v.empty() ) return;

    size_t numnodes = mesh->getSize(0);

    JMeshExporter mexp;

    if( var == "U" || var == "V") {
        MatrixXd uv;
        uv.resize(numnodes,2);

        if( !v.empty())  {
            for( size_t i = 0; i < numnodes; i++) {
                uv.coeffRef(i,0)  = u[i];
                uv.coeffRef(i,1)  = v[i];
                Vertex *vtx = mesh->getNodeAt(i);
                if( var == "U")
                    vtx->setAttribute(var, u[i]);
                else
                    vtx->setAttribute(var, v[i]);
            }
        }
        mexp.addNodeAttribute(var);
    }

    if( var == "S") {
        computeStress();
        for( size_t i = 0; i < numnodes; i++) {
            Vertex *vtx = mesh->getNodeAt(i);
            vtx->setAttribute(var, nodeVonMises[i]);
        }
        mexp.addNodeAttribute(var);
    }


    mexp.saveAs(mesh, filename);
}

/////////////////////////////////////////////////////////////////////////////////////
Mesh*  readMatlabMesh( const string filename)
{
    ifstream infile(  filename.c_str(), ios::in);
    if( infile.fail() ) return nullptr;

    int numNodes, numFaces, numEdges;
    int nOrder;

    infile >> numNodes >> numFaces >> numEdges;
    infile >> nOrder;

    Mesh *mesh = Mesh::newObject();

    Point3D xyz;
    xyz[2] = 0.0;
    for( int i = 0; i < numNodes; i++) {
        infile >> xyz[0] >> xyz[1];
        Vertex *vtx = Vertex::newObject();
        vtx->setXYZCoords(xyz);
        vtx->setID(i);
        mesh->addObject(vtx);
    }

    string line;
    int maxid = 0;
    int v[7];
    JNodeSequence nodes(3);
    for(int i = 0; i < numFaces; i++) {
        for( int j = 0; j < 3; j++) {
            infile >> v[j];
        }
        nodes[0] = mesh->getNodeAt( v[0]-1 );
        nodes[1] = mesh->getNodeAt( v[1]-1 );
        nodes[2] = mesh->getNodeAt( v[2]-1 );
        maxid    = max( maxid, v[0]-1);
        maxid    = max( maxid, v[1]-1);
        maxid    = max( maxid, v[2]-1);
        Triangle *tri = Triangle::newObject( nodes);
        mesh->addObject( tri );
        if( nOrder == 2) {
            for( int j = 0; j < 3; j++) {
                infile >> v[j];
                Edge *edge = tri->getEdgeAt(j);
                assert( edge);
                edge->setAttribute("QuadNode", v[j]-1);
            }
        }
        getline(infile, line);
    }

    if( nOrder == 1) mesh->getTopology()->collect_edges();

    for(int i = 0; i < numEdges; i++) {
        for( int j = 0; j < 4; j++)
            infile >> v[j];
        nodes[0] = mesh->getNodeAt( v[0]-1 );
        nodes[1] = mesh->getNodeAt( v[1]-1 );
        Edge *edge = Simplex::getEdgeOf( nodes[0], nodes[1] );
        assert( edge);
        edge->setAttribute("Boundary", v[3] );
    }

    for( int i = maxid+1; i < numNodes; i++) {
        Vertex *vtx = mesh->getNodeAt(i);
        vtx->setStatus(MeshEntity::REMOVE);
    }
    mesh->pruneAll();
    return mesh;
}

/////////////////////////////////////////////////////////////////////////////////////

