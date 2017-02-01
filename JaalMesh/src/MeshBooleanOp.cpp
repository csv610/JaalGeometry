#include "MeshBoolean.hpp"


#ifdef USE_CGAL
using namespace igl::copyleft::cgal;
#endif

/////////////////////////////////////////////////////////////////////////
int JMeshBoolean :: getOp( const string &str)
{
    if( str == "Union") return MESH_UNION;
    if( str == "Difference") return MESH_DIFFERENCE;
    if( str == "Intersection") return MESH_INTERSECTION;
    if( str == "Symmetric Difference") return MESH_SYMMETRIC_DIFFERENCE;
    return -1;
}

/////////////////////////////////////////////////////////////////////////

void JMeshBoolean :: setMesh( const JMeshPtr &meshA, const JMeshPtr &meshB)
{
    setMesh1(meshA);
    setMesh2(meshB);
    if( meshA == nullptr || meshB == nullptr)  return;
}

/////////////////////////////////////////////////////////////////////////
void JMeshBoolean :: setMesh1( const JMeshPtr &meshA)
{
    if( meshA == nullptr) return;
    JMeshEigenMatrix mat;
    mat.setMesh(meshA);
    VA = mat.getNodeMatrix();
    FA = mat.getFaceMatrix();
}
/////////////////////////////////////////////////////////////////////////

void JMeshBoolean :: setMesh2( const JMeshPtr &meshB)
{
    if( meshB == nullptr) return;

    JMeshEigenMatrix mat;
    mat.setMesh(meshB);
    VB = mat.getNodeMatrix();
    FB = mat.getFaceMatrix();
}
/////////////////////////////////////////////////////////////////////////

JMeshPtr  JMeshBoolean :: getMesh( int op )
{
#ifdef USE_CGAL
    Eigen::VectorXi J;

    if( op == MESH_UNION) {
        mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_UNION,VC,FC);
    }

    if( op == MESH_INTERSECTION) {
        mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_INTERSECT,VC,FC);
    }

    if( op == MESH_DIFFERENCE)
        mesh_boolean(VA,FA,VB,FB,igl::MESH_BOOLEAN_TYPE_MINUS,VC,FC);

    JMeshEigenMatrix mat;
    JMeshPtr meshC = mat.getMesh(VC, FC);

    int Aid = 0;
    int Bid = 1;
    for(size_t i = 0; i< FC.rows(); i++)
    {
        const JFacePtr &face = meshC->getFaceAt(i);
        if(J(i)< FA.rows())
            face->setAttribute("Partition", Aid);
        else
            face->setAttribute("Partition", Bid);
    }

    return meshC;
#endif
    return nullptr;
}
/////////////////////////////////////////////////////////////////////////
void JMeshBoolean :: splitAlongIntersection()
{
//   JMeshPtr meshC = getOp(MESH_INTERSECTION);

}
/////////////////////////////////////////////////////////////////////////

#ifdef CSV
/////////////////////////////////////////////////////////////////////////
bool JMeshBoolean :: anyIntersection()
{
    if( aMesh == nullptr || bMesh == nullptr) {
        cout << "Warning:  Two mesh must be provided for the query" << endl;
        return 0;
    }

    JMeshPtr am, bm;
    if( aMesh->getSize(1) > bMesh->getSize(1) )  {
        am = aMesh;
        bm = bMesh;
    } else {
        am = bMesh;
        bm = aMesh;
    }
    JMeshIO mio;
    mio.saveAs(am, "tmp.off");

    std::ifstream input("tmp.off");
    CGALMesh mesh;
    input >> mesh;
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    bool v = CGAL::Polygon_mesh_processing::is_outward_oriented(mesh);

    if( !v ) {
        cout << "Error: Theb input mesh has normal in the negative diretion" << endl;
        return 0;
    }

    size_t numEdges = bm->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = bm->getEdgeAt(i);
        const Point3D  &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D  &p1 = edge->getNodeAt(1)->getXYZCoords();
        // constructs segment query
        Point a(p0[0], p0[1], p0[2] );
        Point b(p1[0], p1[1], p1[2] );
        Segment segment_query(a,b);
        // tests intersections with segment query
        if(tree.do_intersect(segment_query)) return 1;
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////
int JMeshBoolean :: getSkinMesh( const JMeshPtr &inmesh)
{
/*
    JMeshComponents jc;
    jc.setMesh(jc);
    jc.searchComponents();

    int nc = jc.getNumComponents();
    if( nc == 0) return 1;

    vector<JMeshPtr> submeshes(nc);
    for( int i = 0; i < nc; i++) {
        submeshes[i] = jc.getComponent(i);
        submeshes[i]->setAttribute("Component", i);
    }

    for( int i = 0; i < nc; i++) {
        if( !submeshes[i]->isEmpty() ) {
            for( int j = i+1; j < nc; j++) {
                if( !submeshes[j]->isEmpty() ) {
                    if(  doIntersect( submeshes[i], submeshes[j] ) ) {
                        setMesh(submeshs[i], submeshs[j] );
                        JMeshPtr AminusB = getOp( MESH_DIFFERENCE);
                        JMeshPtr BminusA = getOp( MESH_DIFFERENCE);
                        setMesh(AminusB, BminusA);
                        submeshes[i] = getOp(MESH_UNION);
                        submeshes[i]->setAttribute("Component", i);
                        submeshes[j]->clearAll();
                    }
                }
            }
          }
      }
*/
}

/////////////////////////////////////////////////////////////////////////
#endif

