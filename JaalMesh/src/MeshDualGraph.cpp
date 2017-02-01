#include "MeshDualGraph.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshDualGraph:: getGraph()
{
    if( mesh == nullptr ) return nullptr;

    if( dualPos == DUAL_AT_HODGE_CENTER) {
        DDG::JMeshHot2  appl;
        appl.setMesh(mesh);
        graph = appl.getDualMesh();
        mesh->setAttribute("DualGraph", graph);
        return graph;
    }

    int err = mesh->getAttribute("DualGraph", graph);
    if( err) {
        graph = JMesh::newObject();
        mesh->setAttribute("DualGraph", graph);
    }

    graph->clearAll();
    int dim = mesh->getTopology()->getDimension();

    if( dim == 2 ) build2d();
    if( dim == 3 ) build3d();
    return graph;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualGraph:: build2d()
{
    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);

    JNodePtr dualvtx;
    JNodeSequence nodes;
    nodes.reserve(numfaces);

    int index = 0;
    for (size_t iface = 0; iface < numfaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            face->setID(index);
            dualvtx = JMeshDualGraph::newObject(face);
            dualvtx->setID(index);
            nodes.push_back(dualvtx);
            index++;
        }
    }
    graph->addObjects(nodes);

    JNodePtr dv0, dv1;

    size_t numedges = mesh->getSize(1);
    assert( numedges ) ;

    JFaceSequence neighs;
    for (size_t iedge = 0; iedge < numedges; iedge++) {
        const JEdgePtr &edge = mesh->getEdgeAt(iedge);
        if( edge->isActive() ) {
            JEdge::getRelations(edge, neighs);
            assert( !neighs.empty() ) ;
            if (neighs.size() == 2) {
                neighs[0]->getAttribute("DualNode", dv0);
                neighs[1]->getAttribute("DualNode", dv1);
                JEdgePtr dedge = JSimplex::getEdgeOf(dv0, dv1, 1);
                graph->addObject( dedge );
            } else {
                if( node_on_boundary) {
                    neighs[0]->getAttribute("DualNode", dv0);
                    dv1 = JMeshDualGraph::newObject(edge);
                    dv1->setAttribute("Boundary", 1);
                    graph->addObject(dv1);
                    JEdgePtr dedge = JSimplex::getEdgeOf(dv0, dv1,1);
                    graph->addObject( dedge );
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualGraph:: build3d()
{
    size_t numcells = mesh->getSize(3);

    int relexist = mesh->buildRelations(2,3);

    JNodePtr dualvtx;
    JCellPtr cell;
    int index = 0;
    for (size_t icell = 0; icell < numcells; icell++) {
        cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            cell->setID(index++);
            dualvtx = JMeshDualGraph::newObject(cell);
            dualvtx->setID(icell);
            graph->addObject(dualvtx);
        }
    }

    JNodePtr dv0, dv1;

    size_t numfaces = mesh->getSize(2);

    JCellSequence neighs;
    for (size_t iface = 0; iface < numfaces; iface++) {
        JFacePtr face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            JFace::getRelations(face, neighs);
            if (neighs.size() == 2) {
                neighs[0]->getAttribute("DualNode", dv0);
                neighs[1]->getAttribute("DualNode", dv1);
                JEdgePtr dedge = JEdge::newObject(dv0, dv1);
                graph->addObject( dedge );
            } else {
                neighs[0]->getAttribute("DualNode", dv0);
                dv0->setAttribute("Boundary", 1);
            }
        }
    }

    if (!relexist)
        mesh->clearRelations(2, 3);
}

///////////////////////////////////////////////////////////////////////////////
