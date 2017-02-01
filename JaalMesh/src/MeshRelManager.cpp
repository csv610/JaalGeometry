#include "Mesh.hpp"

using namespace Jaal;

JLogger* JMeshRelationManager::logger = JLogger::getInstance();

void JMeshRelationManager  :: printAdjTable() const
{
    cout <<  "  V E F C" << endl;
    cout << "V " << (int) adjTable[0][0] << " "
         << (int) adjTable[0][1] << " "
         << (int) adjTable[0][2] << " "
         << (int) adjTable[0][3] << endl;

    cout << "E " << (int) adjTable[1][0] << " "
         << (int) adjTable[1][1] << " "
         << (int) adjTable[1][2] << " "
         << (int) adjTable[1][3] << endl;

    cout << "F " << (int) adjTable[2][0] << " "
         << (int) adjTable[2][1] << " "
         << (int) adjTable[2][2] << " "
         << (int) adjTable[2][3] << endl;

    cout << "C " << (int) adjTable[3][0] << " "
         << (int) adjTable[3][1] << " "
         << (int) adjTable[3][2] << " "
         << (int) adjTable[3][3] << endl;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshRelationManager :: updateRelations()
{
    for( int i = 0; i < 4; i++) {
        for( int j = 0; j < 4; j++) {
            if( adjTable[i][j] ) buildRelations(i,j);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager :: removeRelations(const JNodePtr &)
{
//   LOG4CXX_ERROR( Mesh::logger, "vertex relations not removed ");
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager :: removeRelations( const JEdgePtr &)
{
//   LOG4CXX_WARN( Mesh::logger,  "Edge relations not removed" );
    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager :: removeRelations(const JFacePtr &face )
{
    int nSize;

    if( adjTable[0][2] ) {
        nSize = face->getSize(0);
        for( int i = 0; i < nSize; i++) {
            const JNodePtr &vtx = face->getNodeAt(i);
            vtx->removeRelation(face);
        }
    }

    if( adjTable[1][2] ) {
        nSize = face->getSize(0);
        for( int i = 0; i < nSize; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            edge->removeRelation(face);
        }
    }

    if( adjTable[2][1] ) {
        nSize = face->getSize(0);
        for( int i = 0; i < nSize; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            face->removeRelation(edge);
        }
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager :: removeRelations( const JCellPtr &cell )
{
    int nSize;
    if( adjTable[0][3] ) {
        nSize = cell->getSize(0);
        for( int  i = 0; i < nSize; i++) {
            const JNodePtr &vtx = cell->getNodeAt(i);
            vtx->removeRelation( cell );
        }
    }

    JEdgeSequence celledges;
    if( adjTable[1][3] ) {
        celledges = cell->getEdges();
        nSize = celledges.size();
        for( int i = 0; i < nSize; i++)
            celledges[i]->removeRelation(cell);
    }

    if( adjTable[3][1] ) {
        celledges = cell->getEdges();
        nSize = celledges.size();
        for( int i = 0; i < nSize; i++)
            cell->removeRelation(celledges[i] );
    }

    JFaceSequence cellfaces;
    if( adjTable[2][3] ) {
        cellfaces = cell->getFaces();
        nSize = cellfaces.size();
        for( int  i = 0; i < nSize; i++)
            cellfaces[i]->removeRelation(cell);
    }

    if( adjTable[3][2] ) {
        cellfaces = cell->getFaces();
        nSize = cellfaces.size();
        for( int i = 0; i < nSize; i++)
            cell->removeRelation(cellfaces[i] );
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshRelationManager:: buildRelations( const JEdgePtr &edge)
{
    if( edge == nullptr ) return;

    if( !edge->isActive() ) return;

    if( adjTable[0][0] ) {
        const JNodePtr &v0 = edge->getNodeAt(0);
        const JNodePtr &v1 = edge->getNodeAt(1);
        v0->addRelation(v1);
        v1->addRelation(v0);
    }

    if( adjTable[0][1] ) {
        const JNodePtr &v0 = edge->getNodeAt(0);
        const JNodePtr &v1 = edge->getNodeAt(1);
        v0->addRelation(edge);
        v1->addRelation(edge);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRelationManager :: buildRelations( const JFacePtr &newface)
{
    if( !newface->isActive() ) return;

    JEdgeSequence faceedges;
    int numedges;

    if (adjTable[0][0]) {
        for (int i = 0; i < newface->getSize(0); i++) {
            const JNodePtr &v0 = newface->getNodeAt(i);
            const JNodePtr &v1 = newface->getNodeAt(i+1);
            v0->addRelation(v1);
            v1->addRelation(v0);
        }
    }

    if (adjTable[0][1]) {
        faceedges = newface->getEdges();
        numedges = faceedges.size();
        for (int i = 0; i < numedges; i++) {
            const JNodePtr &v0 = faceedges[i]->getNodeAt(0);
            const JNodePtr &v1 = faceedges[i]->getNodeAt(1);
            v0->addRelation( faceedges[i] );
            v1->addRelation( faceedges[i] );
        }
    }

    if (adjTable[0][2]) {
        for (int i = 0; i < newface->getSize(0); i++) {
            const JNodePtr &vtx = newface->getNodeAt(i);
            vtx->addRelation(newface);
        }
    }

    if(adjTable[1][2] ) {
        faceedges = newface->getEdges();
        numedges = faceedges.size();
        for (int i = 0; i < numedges; i++)
            faceedges[i]->addRelation(newface);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRelationManager :: buildRelations( const JCellPtr &newcell)
{
    if( !newcell->isActive() )  return;

    int numedges, numfaces;

    JNodePtr vhash;

    JFaceSequence cellfaces;
    cellfaces = newcell->getFaces();
    numfaces = cellfaces.size();
    for( int i = 0; i < numfaces; i++) {
        const JFacePtr &cface = cellfaces[i];
        vhash = cface->getHashNode();
        if( vhash ) vhash->attach(cface );
    }

    JEdgeSequence celledges;
    celledges = newcell->getEdges();
    numedges = celledges.size();
    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &cedge = celledges[i];
        assert( cedge );
        vhash = cedge->getHashNode();
        if( vhash ) vhash->attach(cedge);
    }


    if( adjTable[0][0] ) {
        for (int i = 0; i < numedges; i++) {
            const JNodePtr &vi = celledges[i]->getNodeAt((i + 0));
            const JNodePtr &vj = celledges[i]->getNodeAt((i + 1));
            vi->addRelation(vj);
            vj->addRelation(vi);
        }
    }

    if( adjTable[0][1] ) {
        for (int i = 0; i < numedges; i++) {
            const JNodePtr &vi = celledges[i]->getNodeAt((i + 0));
            const JNodePtr &vj = celledges[i]->getNodeAt((i + 1));
            vi->addRelation( celledges[i] );
            vj->addRelation( celledges[i] );
        }
    }

    if( adjTable[0][2] ) {
        for (int i = 0; i < numfaces; i++) {
            const JFacePtr &f  = cellfaces[i];
            for( int j = 0; j < f->getSize(0); j++) {
                const JNodePtr &vtx = f->getNodeAt(j);
                vtx->addRelation(f);
            }
        }
    }

    if (adjTable[0][3]) {
        int numnodes = newcell->getSize(0);
        for (int i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = newcell->getNodeAt(i);
            vtx->addRelation(newcell);
        }
    }

    if( adjTable[2][3] ) {
        for (int i = 0; i < numfaces; i++) {
            const JFacePtr &f  = cellfaces[i];
            f->addRelation(newcell);
        }
    }

    if( adjTable[3][1] ) {
        for (int i = 0; i < numedges; i++) {
            const JEdgePtr &edge  = celledges[i];
            newcell->addRelation(edge);
        }
    }

    if( adjTable[3][2] ) {
        for (int i = 0; i < numfaces; i++) {
            const JFacePtr &face  = cellfaces[i];
            newcell->addRelation(face);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations00()
{
    if( mesh == nullptr ) return 1;

    clearRelations(0, 0);

    JLogger *logger = JLogger::getInstance();
    if( logger ) logger->setInfo("Building vertex-vertex relations ");

    JEdgeSequence celledges;

    size_t numcells = mesh->cells.size();
    size_t numfaces = mesh->faces.size();
    size_t numedges = mesh->edges.size();

    for (size_t icell = 0; icell < numcells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            celledges = cell->getEdges();
            for (size_t j = 0; j < celledges.size(); j++) {
                const JNodePtr &v0 = celledges[j]->getNodeAt(0);
                const JNodePtr &v1 = celledges[j]->getNodeAt(1);
                v0->addRelation(v1);
                v1->addRelation(v0);
            }
        }
    }

    if( numcells == 0) {
        for (size_t iface = 0; iface < numfaces; iface++) {
            const JFacePtr &face = mesh->getFaceAt(iface);
            if( face->isActive() ) {
                int  nnodes = face->getSize(0);
                for (int j = 0; j < nnodes; j++) {
                    const JNodePtr &v0 = face->getNodeAt(j);
                    const JNodePtr &v1 = face->getNodeAt(j + 1);
                    v0->addRelation(v1);
                    v1->addRelation(v0);
                }
            }
        }
    }


    if( numfaces == 0) {
        for (size_t iedge = 0; iedge < numedges; iedge++) {
            const JEdgePtr &edge = mesh->getEdgeAt(iedge);
            if( edge->isActive() ) {
                const JNodePtr &v0 = edge->getNodeAt(0);
                const JNodePtr &v1 = edge->getNodeAt(1);
                v0->addRelation(v1);
                v1->addRelation(v0);
            }
        }
    }

    adjTable[0][0] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations01()
{
    clearRelations(0, 1);

    size_t numedges = mesh->getSize(1);
    if( numedges == 0)  {
        logger->setWarn("No mesh edges present ");
        return 1;
    }

    logger->setInfo("Building new vertex-edge relations ");

    for (size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            const JNodePtr &v0 = edge->getNodeAt(0);
            const JNodePtr &v1 = edge->getNodeAt(1);
            v0->addRelation(edge);
            v1->addRelation(edge);
        }
    }

    adjTable[0][1] = 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations02()
{
    clearRelations(0,2);

    size_t numfaces = mesh->faces.size();
    if( numfaces == 0) return 1;

    logger->setInfo("Building new vertex->faces relations ");

    for (size_t iface = 0; iface < numfaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            int nf = face->getSize(0);
            for (int j = 0; j < nf; j++) {
                const JNodePtr &vtx = face->getNodeAt(j);
                vtx->addRelation(face);
            }
        }
    }

    adjTable[0][2] = 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager   ::buildRelations03()
{
    clearRelations(0,3);

    size_t numCells = mesh->getSize(3);

    if( numCells == 0) return 1;

    logger->setInfo("Building new vertex->cells relations ");

    for (size_t icell = 0; icell < numCells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            int nf = cell->getSize(0);
            for (int j = 0; j < nf; j++) {
                const JNodePtr &vtx = cell->getNodeAt(j);
                vtx->addRelation(cell);
            }
        }
    }

    adjTable[0][3] = 1;
    return 0;
}


///////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations12()
{
    clearRelations(1,2);

    size_t numfaces = mesh->getSize(2);

    if( numfaces == 0) return 1;

    logger->setInfo("Building new edge->faces relations ");

    JEdgeSequence faceedges;

    for (size_t iface = 0; iface < numfaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            faceedges = face->getEdges();
            int nedges = faceedges.size();
            for( int  i = 0; i < nedges; i++)
                faceedges[i]->addRelation(face);
        }
    }

    adjTable[1][2] = 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations13()
{
    clearRelations(1,3);

    size_t numCells = mesh->cells.size();
    if( numCells == 0) return 1;

    logger->setInfo("Building edge->cells relations ");

    JEdgeSequence celledges;
    for (size_t icell = 0; icell < numCells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            celledges = cell->getEdges();
            int nedges = celledges.size();
            for( int i = 0; i < nedges; i++)
                celledges[i]->addRelation(cell);
        }
    }
    adjTable[1][3] = 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::buildRelations21()
{
    clearRelations(2,1);

    size_t numfaces = mesh->getSize(2);
    if( numfaces == 0) return 1;

    logger->setInfo("Building new face->edges relations ");

    JEdgeSequence faceedges;
    for (size_t iface = 0; iface < numfaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            faceedges = face->getEdges();
            for( size_t i = 0; i < faceedges.size(); i++)
                face->addRelation( faceedges[i] );
        }
    }

    adjTable[2][1] = 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager  ::buildRelations23()
{
    clearRelations(2,3);

    size_t numCells = mesh->getSize(3);
    if( numCells == 0) return 1;

    logger->setInfo("Building new face->cells relations ");

    JFaceSequence cellfaces;
    for (size_t icell = 0; icell < numCells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            cellfaces = cell->getFaces();
            int numfaces = cellfaces.size();
            for( int j = 0;  j < numfaces; j++)
                cellfaces[j]->addRelation(cell);
        }
    }

    adjTable[2][3] = 1;
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager  ::buildRelations31()
{
    mesh->clearRelations(3,1);

    size_t numCells = mesh->getSize(3);
    if( numCells == 0) return 1;

    logger->setInfo("Building cell->edges relations ");

    JEdgeSequence celledges;
    for (size_t icell = 0; icell < numCells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            celledges = cell->getEdges();
            int numedges = celledges.size();
            for( int j = 0;  j < numedges; j++)
                cell->addRelation(celledges[j], JRelationManager::ORDERED_RELATIONS);
        }
    }

    adjTable[3][1] = 1;
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager  ::buildRelations32()
{
    mesh->clearRelations(3,2);

    size_t numCells = mesh->getSize(3);
    if( numCells == 0) return 1;

    logger->setInfo("Building cell->faces relations ");

    JFaceSequence cellfaces;
    for (size_t icell = 0; icell < numCells; icell++) {
        JCellPtr cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            cellfaces = cell->getFaces();
            int numfaces = cellfaces.size();
            for( int j = 0;  j < numfaces; j++)
                cell->addRelation(cellfaces[j], JRelationManager::ORDERED_RELATIONS);
        }
    }
    adjTable[3][2] = 1;
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int
JMeshRelationManager ::clearRelations(int src, int dst)
{
    assert( src >= 0 && src <= 3);
    assert( dst >= 0 && dst <= 3);

    if (src == 0 ) {
        switch(dst) {
        case 0:
            logger->setInfo("Clearning vertex->vertex relations");
            break;
        case 1:
            logger->setInfo("Clearning vertex->edge relations");
            break;
        case 2:
            logger->setInfo("Clearning vertex->face relations");
            break;
        case 3:
            logger->setInfo("Clearning vertex->cell relations");
            break;
        }
        size_t numnodes = mesh->getSize(0);
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            vtx->clearRelations(dst);
        }
        adjTable[0][dst] = 0;
    }

    if (src == 1 ) {
        switch(dst) {
        case 0:
            logger->setInfo("Clearning edge->vertex relations");
            break;
        case 1:
            logger->setInfo("Clearning edge->edge relations");
            break;
        case 2:
            logger->setInfo("Clearning edge->face relations");
            break;
        case 3:
            logger->setInfo("Clearning edge->cell relations");
            break;
        }
        size_t numedges = mesh->getSize(1);
        for (size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            edge->clearRelations(dst);
        }
        adjTable[1][dst] = 0;
    }

    if (src == 2 ) {
        switch(dst) {
        case 0:
            logger->setInfo("Clearning face->vertex relations");
            break;
        case 1:
            logger->setInfo("Clearning face->edge relations");
            break;
        case 2:
            logger->setInfo("Clearning face->face relations");
            break;
        case 3:
            logger->setInfo("Clearning face->cell relations");
            break;
        }
        size_t numfaces = mesh->getSize(2);
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            face->clearRelations(dst);
        }
        adjTable[2][dst] = 0;
    }

    if (src == 3 ) {
        switch(dst) {
        case 0:
            logger->setInfo("Clearning cell->vertex relations");
            break;
        case 1:
            logger->setInfo("Clearning cell->edge relations");
            break;
        case 2:
            logger->setInfo("Clearning cell->face relations");
            break;
        case 3:
            logger->setInfo("Clearning cell->cell relations");
            break;
        }
        size_t numcells = mesh->getSize(3);
        for (size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            cell->clearRelations(dst);
        }
        adjTable[3][dst] = 0;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshRelationManager :: buildRelations(int src, int dst)
{
    if (src == 0 && dst == 0) return buildRelations00( );
    if (src == 0 && dst == 1) return buildRelations01( );
    if (src == 0 && dst == 2) return buildRelations02( );
    if (src == 0 && dst == 3) return buildRelations03( );

    if (src == 1 && dst == 2) return buildRelations12( );
    if (src == 1 && dst == 3) return buildRelations13( );

    if (src == 2 && dst == 1) return buildRelations21( );
    if (src == 2 && dst == 3) return buildRelations23( );

    if (src == 3 && dst == 1) return buildRelations31( );
    if (src == 3 && dst == 2) return buildRelations32( );

    return 1;
}

////////////////////////////////////////////////////////////////////////////////
