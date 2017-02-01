#include "Singlet.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
bool
JSinglet::isSinglet (const JNodePtr &vertex)
{
    assert( vertex );

    if( !vertex->isActive() ) return 0;
    // A Singlet is always on boundary...
    if( !vertex->isBoundary() ) return 0;

    // A singlet has only  neigbour face ..
    int numfaces = vertex->getNumRelations(2);

    if( numfaces == 0) {
        cout << "Warning: Singlet not determined, may be vertex-face relations absent" << endl;
    }

    if (numfaces == 1) return 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int  JSinglet :: getSize() {
    searchSinglets();
    return singlets.size();
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence JSinglet::getSinglets()
{
    searchSinglets();
    return singlets;
}

///////////////////////////////////////////////////////////////////////////////

void JSinglet::searchSinglets()
{
    singlets.clear();

    if( mesh == nullptr)  return;
    //
    // A boundary singlet is a vertex which is shared by only one face.
    // They are undesirables in the quad mesh as that would mean large
    // angle on some of the edges..
    //
    // For the flat singlet ( angle closer to 180 degree ). it is easy
    // to remove the neighbouring quad from the mesh.
    //
    if( mesh->getTopology()->getDimension() != 2 ) {
        cout << "Warning: Singlet presently only for 2D mesh " << endl;
        return;
    }
    mesh->getTopology()->searchBoundary();

    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() && JSinglet::isSinglet(v))
            singlets.push_back(v);
    }
}

///////////////////////////////////////////////////////////////////////////////

int JSinglet::remove( const JNodePtr &vertex )
{
    JFaceSequence vfaces;
    JNode::getRelations( vertex, vfaces );
    if (vfaces.size() > 1) return 1;
    quadRefiner.setMesh( mesh );

    return quadRefiner.refine5(vfaces[0]);
}

///////////////////////////////////////////////////////////////////////////////
int JSinglet::removeAll()
{
    searchSinglets();
    for( const JNodePtr &vtx : singlets)
        remove(vtx);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JSinglet :: mergeOpposite(const JNodePtr &vtx)
{
    if( !vtx->isActive() ) return 1;

    JFaceSequence vfaces;
    JNode::getRelations( vtx, vfaces );
    if (vfaces.size() > 1) return 1;

    JFacePtr thisface = vfaces[0];

    int pos    = thisface->getPosOf(vtx);
    int nnodes = thisface->getSize(0);

    JNodePtr  v1, v2;
    JEdgePtr  nextEdge;

    v1 = thisface->getNodeAt(pos+1);
    if( v1->getNumRelations(2) > 2 ) {
        nextEdge = thisface->getEdgeAt(pos+1);
        v2 = thisface->getNodeAt(pos+2);
        if( v2->getNumRelations(2) > 2 )  {
            int err = edgeSwap.applyAt(nextEdge, vtx);
            if( !err) return 0;
        }
    }

    v1 = thisface->getNodeAt(pos+nnodes-1);
    if( v1->getNumRelations(2) > 2 ) {
        nextEdge = thisface->getEdgeAt(pos + nnodes-2);
        v2 = thisface->getNodeAt(pos+nnodes-2);
        if( v2->getNumRelations(2) > 2 )  {
            int err = edgeSwap.applyAt(nextEdge, vtx);
            if( !err) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void JSinglet :: mergeAll()
{
    edgeSwap.setMesh(mesh);
    searchSinglets();
    for( const JNodePtr &vtx: singlets ) mergeOpposite(vtx);
}

