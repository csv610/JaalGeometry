#ifdef CSV

#include "Mesh.hpp"
#include "basic_math.hpp"
#include "circumcenter.hpp"
#include "MeshTopology.hpp"

sdfsdsdfsfs

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
int MeshDualGrapher :: update_primal( JNodePtr vertex)
{
    if( vertex == NULL || vertex->isBoundary() ) return 1;

    JNode::getRelations(vertex, faceneighs);

    int numneighs = faceneighs.size();

    JNodePtr dualnode;
    Point3D newpos, xyz;
    newpos[0] = 0.0;
    newpos[1] = 0.0;
    newpos[2] = 0.0;
    for( int i = 0; i < numneighs; i++) {
        faceneighs[i]->getAttribute("DualNode", dualnode);
        xyz = dualnode->getXYZCoords();
        newpos[0] += xyz[0];
        newpos[1] += xyz[1];
        newpos[2] += xyz[2];
    }

    newpos[0] /= (double)numneighs;
    newpos[1] /= (double)numneighs;
    newpos[2] /= (double)numneighs;

    vertex->setXYZCoords(newpos);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int MeshDualGrapher :: update_dual( JNodePtr dvertex)
{
    if( dvertex == NULL || dvertex->isBoundary() ) return 1;

    JFacePtr face;
    dvertex->getAttribute("PrimalFace", face);

    Point3D xyz;
    JFaceGeometry::getDualPosition(face,xyz);
    dvertex->setXYZCoords(xyz);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MeshDualGrapher ::update_primal_from_dual()
{
    if( mesh == NULL ) return 1;

    size_t numnodes = mesh->getSize(0);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    for( size_t i = 0; i < numnodes; i++)
        update_primal( mesh->getNodeAt(i));

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int MeshDualGrapher ::update_dual_from_primal()
{
    /*
        if( mesh == NULL ) return 1;

        JMeshPtr dgraph = mesh->getDualGraph();

        if( dgraph == NULL ) return 1;
        size_t numnodes = dgraph->getSize(0);

        for( size_t i = 0; i < numnodes; i++)
            update_dual( dgraph->getNodeAt(i));
    */
    cout << "Exit " << endl;
    exit(0);

    return 0;

}
////////////////////////////////////////////////////////////////////////////////


JNodePtr MeshDualGrapher :: newDualNode( JFacePtr face)
{
    JNodePtr dvertex = JNode::newObject();
    face->setAttribute("DualNode", dvertex);
    dvertex->setAttribute("PrimalFace", face);
    update_dual(dvertex);
    return dvertex;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr MeshDualGrapher :: getGraph(const JMeshPtr &m)
{
    mesh = m;

    if( mesh == NULL ) return NULL;

    mesh->getLogger()->setInfo("Building Dual graph");


    JMeshPtr dualGraph = JMesh::newObject();
    dualGraph->setName("DualGraph");
    int bmark;
    Point3D xyz;
    JNodePtr dv0, dv1, vmid, dvertex;
    JFaceSequence faceneighs;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim == 2 ) {
        mesh->buildRelations(1,2);
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            JFacePtr face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                if( !face->hasAttribute("DualNode") ) {
                    dvertex = newDualNode(face);
                    dualGraph->addObject(dvertex);
                } else {
                    face->getAttribute("DualNode", dvertex);
                    update_dual(dvertex);
                }
            }
        }

        size_t numedges = mesh->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            JEdgePtr edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                JEdge::getRelations(edge, faceneighs);
                if( faceneighs.size() == 2 ) {
                    faceneighs[0]->getAttribute("DualNode", dv0);
                    faceneighs[1]->getAttribute("DualNode", dv1);
                    JEdgePtr dedge = JEdge::newObject(dv0,dv1);
                    if( midnodes ) {
                        vmid = JNodeGeometry::getMidNode(edge->getNodeAt(0), edge->getNodeAt(1));
                        dedge->setAttribute("Steiner", vmid);
                    }
                    dualGraph->addObject(dedge);
                } else {
                    if( !edge->hasAttribute("DualNode") ) {
                        dv0 = JNodeGeometry::getMidNode(edge->getNodeAt(0), edge->getNodeAt(1));
                        edge->setAttribute("DualNode", dv0);
                    }
                    edge->getAttribute("DualNode", dv0);
                    dualGraph->addObject(dv0);
//                    bmark = max( 1, edge->getBoundaryMark() );
//                    dv0->setBoundaryMark(bmark);
                    faceneighs[0]->getAttribute("DualNode", dv1);
                    JEdgePtr dedge = JEdge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                }
            }
        }
    }

    JCellSequence cellneighs;
    if( topDim == 3 ) {
        mesh->buildRelations(2,3);
        size_t numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                if(!cell->hasAttribute("DualNode") ) {
                    cell->getAvgXYZ(xyz);
                    JNodePtr dvertex = JNode::newObject();
                    dvertex->setXYZCoords(xyz);
                    dualGraph->addObject(dvertex);
                    cell->setAttribute("DualNode", dvertex);
                    dvertex->setAttribute("DualCell", cell);
                }
            }
        }

        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            JFacePtr face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                JFace::getRelations(face, cellneighs);
                if( cellneighs.size() == 2 ) {
                    cellneighs[0]->getAttribute("DualNode", dv0);
                    cellneighs[1]->getAttribute("DualNode", dv1);
                    JEdgePtr dedge = JEdge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                } else {
                    if(!face->hasAttribute("DualNode") ) {
                        dv0 = JNode::newObject();
                        face->setAttribute("DualNode", dv0);
                    }
                    face->getAttribute("DualNode", dv0);
                    face->getAvgXYZ(xyz);
                    dv0->setXYZCoords(xyz);
                    dualGraph->addObject(dv0);
//                    bmark = max( 1, face->getBoundaryMark() );
//                    dv0->setBoundaryMark(bmark);
                    cellneighs[0]->getAttribute("DualNode", dv1);
                    JEdgePtr dedge = JEdge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                }
            }
        }
    }

    return dualGraph;
}
#endif

