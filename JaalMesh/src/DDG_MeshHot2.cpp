#include "DDG_MeshHot2.hpp"

using namespace DDG;
void DDG::JMeshHot2 :: setMesh( const JMeshPtr &m)
{
    jmesh = m;
    if( jmesh == nullptr) return;

    JMeshOBJExporter mexp;
    mexp.writeFile( jmesh, "tmp.obj");

    mesh.read("tmp.obj");
}

///////////////////////////////////////////////////////////////////////////
void DDG::JMeshHot2 :: optimizeWeights()
{
    SparseMatrix<Real> d0, star0, star1, Delta;
    HodgeStar1Form<Real>::build( mesh, star1 );
    HodgeStar0Form<Real>::build( mesh, star0 );
    ExteriorDerivative0Form<Real>::build( mesh, d0 );
    Delta = d0.transpose() * star1 * d0;
    Delta += Real(1e-8)*star0;

    DenseMatrix<Real> rhs;
    buildRhs(rhs);

    DenseMatrix<Real> x;
    solvePositiveDefinite(Delta, x, rhs);

    assignSolution(x);
}
///////////////////////////////////////////////////////////////////////////

void DDG::JMeshHot2 :: buildRhs(DenseMatrix<Real>& rhs)
{
    rhs = DenseMatrix<Real>(mesh.vertices.size(),1);
    for( VertexCIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        double sum = 0.0;
        HalfEdgeIter he = v->he;
        do
        {
            if( !he->onBoundary )
            {
                Vector n = he->next->rotatedEdge();
                Vector c = he->face->circumcenter();
                Vector b = he->face->barycenter();
                sum += dot( n, c-b );
            }
            he = he->flip->next;
        }
        while( he != v->he );
        rhs(v->index,0) = sum;
    }
}
///////////////////////////////////////////////////////////////////////////

void DDG::JMeshHot2 :: assignSolution(const DenseMatrix<Real>& x)
{
    for( VertexIter v = mesh.vertices.begin();
            v != mesh.vertices.end();
            v ++)
    {
        v->weight = x(v->index,0);
    }
}
///////////////////////////////////////////////////////////////////////////
JNodePtr DDG::JMeshHot2 :: getNewNode( const Vector &p)
{
    JNodePtr v = JNode::newObject();
    Point3D xyz;
    xyz[0] = p[0];
    xyz[1] = p[1];
    xyz[2] = p[2];
    v->setXYZCoords(xyz);
    return v;
}
///////////////////////////////////////////////////////////////////////////
JMeshPtr DDG::JMeshHot2 :: getDualMesh()
{
    if( jmesh == nullptr) return nullptr;
    optimizeWeights();

    JMeshPtr dmesh = JMesh::newObject();
    std::map<DDG::EdgeIter, JNodePtr>  vmap;

    JNodePtr v0, v1;


    for( FaceCIter f  = mesh.faces.begin(); f != mesh.faces.end(); f ++ )
    {
        if( f->isBoundary() ) continue;

        Vector cf = f->dualPoint();
        v0 = getNewNode(cf);
        dmesh->addObject(v0);
        HalfEdgeCIter he = f->he;
        do
        {
            if( vmap.find(he->edge) == vmap.end() ) {
                Vector ce = he->edge->dualPoint();
                v1 = getNewNode(ce);
                dmesh->addObject(v1);
            }
            v1 = vmap[he->edge];
            JEdgePtr newedge = JEdge::newObject(v0,v1);
            dmesh->addObject(newedge);
            he = he->next;
        }
        while( he != f->he );
    }
    dmesh->enumerate(0);
    dmesh->enumerate(1);
    dmesh->setName("DualGraph");
    return dmesh;
}
///////////////////////////////////////////////////////////////////////////
