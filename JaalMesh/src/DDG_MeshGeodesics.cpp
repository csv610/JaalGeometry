#include "DDG_MeshGeodesics.hpp"

using namespace DDG;

void DDG::JMeshGeodesics :: setMesh( const JMeshPtr &m)
{
    jmesh = m;
    if( jmesh == nullptr) return;

    JMeshOBJExporter mexp;
    mexp.writeFile( jmesh, "tmp.obj");

    mesh.read("tmp.obj");

    HodgeStar0Form<Real>::build( mesh, star0 );

    SparseMatrix<Real> star1;
    HodgeStar1Form<Real>::build( mesh, star1 );

    SparseMatrix<Real> d0;
    ExteriorDerivative0Form<Real>::build( mesh, d0 );

    // zero Neumann boundary condition
    L  = d0.transpose() * star1 * d0;

    // make L positive-definite
    L += Real(1.0e-8)*star0;
}
///////////////////////////////////////////////////////////////////////////////
int DDG::JMeshGeodesics :: execute()
{
    // initial condiiton
    DenseMatrix<Real> u0;
    int nb = builImpulseSignal(u0);
    if( nb == 0 ) return 1.0;

    // heat flow for short interval
    dt *= sqr(mesh.meanEdgeLength());
    SparseMatrix<Real> A = star0 + Real(dt) * L;

    DenseMatrix<Real> u;
    solvePositiveDefinite(A, u, u0);

    // extract geodesic
    computeVectorField(u);

    DenseMatrix<Real> div;
    computeDivergence(div);

    DenseMatrix<Real> phi;
    solvePositiveDefinite(L, phi, div);

    setMinToZero(phi);
    assignDistance(phi);
//    return phi.norm();
}
///////////////////////////////////////////////////////////////////////////////

int DDG::JMeshGeodesics :: builImpulseSignal(DenseMatrix<Real>& x)
{
    for( VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end(); v ++ ) v->tag = 0;

    for(  const JNodePtr &vtx : srcNodes)  {
        mesh.vertices[vtx->getID()].tag = 1;
    }

    int nb = 0;
    x = DenseMatrix<Real>(mesh.vertices.size());
    for( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++ )
    {
        x(v->index) = 0.0;
        if( v->tag )
        {
            x(v->index) = 1.0;
            nb++;
        }
    }
    return nb;
}
///////////////////////////////////////////////////////////////////////////////

void DDG::JMeshGeodesics :: computeVectorField(const DenseMatrix<Real>& u)
{
    for( FaceIter f = mesh.faces.begin();
            f != mesh.faces.end();
            f++ )
    {
        if( f->isBoundary() ) continue;

        HalfEdgeIter hij = f->he;
        HalfEdgeIter hjk = hij->next;
        HalfEdgeIter hki = hjk->next;

        VertexIter vi = hij->vertex;
        VertexIter vj = hjk->vertex;
        VertexIter vk = hki->vertex;

        double ui = u(vi->index);
        double uj = u(vj->index);
        double uk = u(vk->index);

        Vector eij90 = hij->rotatedEdge();
        Vector ejk90 = hjk->rotatedEdge();
        Vector eki90 = hki->rotatedEdge();

        Vector X = 0.5 * ( ui*ejk90 + uj*eki90 + uk*eij90 ) / f->area();
        f->vector = - X.unit();
    }
}
///////////////////////////////////////////////////////////////////////////////

void DDG::JMeshGeodesics::computeDivergence(DenseMatrix<Real>& div)
{
    div = DenseMatrix<Real>(mesh.vertices.size());
    for( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        double sum = 0.0;
        HalfEdgeIter he = v->he;
        do
        {
            if( not he->onBoundary )
            {
                Vector n = he->next->rotatedEdge();
                Vector v = he->face->vector;
                sum += dot( n, v );
            }
            he = he->flip->next;
        }
        while( he != v->he );
        div(v->index) = sum;
    }
}

///////////////////////////////////////////////////////////////////////////////

void DDG::JMeshGeodesics :: assignDistance(const DenseMatrix<Real>& phi)
{
    /*
        for( VertexIter v = mesh.vertices.begin();
                v != mesh.vertices.end();
                v ++)
        {
            v->distance = phi(v->index);
        }
    */
    size_t index = 0;
    size_t numnodes = jmesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = jmesh->getNodeAt(i);
        if( vtx->isActive() ) {
            double d = phi(index++);
            vtx->setAttribute("Distance", d);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void DDG::JMeshGeodesics :: setMinToZero(DenseMatrix<Real>& phi)
{
    double minValue = 1.0e100;
    for( int i = 0; i < phi.nRows(); ++i )
        for( int j = 0; j < phi.nColumns(); ++j )
            minValue = std::min( minValue, (double) phi(i,j) );

    for( int i = 0; i < phi.nRows(); ++i )
        for( int j = 0; j < phi.nColumns(); ++j )
            phi(i,j) -= minValue;
}
///////////////////////////////////////////////////////////////////////////////
