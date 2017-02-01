#include "MeshTopology.hpp"

void JMeshComponents :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    dim  = mesh->getTopology()->getDimension();
}

///////////////////////////////////////////////////////////////////////////////
int  JMeshComponents :: getNumComponents() const
{
    return components.size();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshComponents :: searchFaceComponents()
{
    JFacePtr face;
    deque<JFacePtr> faceQ;
    JFaceSequence fneighs;

    size_t numfaces = mesh->getSize(2);

    for (size_t iface = 0; iface < numfaces; iface++) {
        face = mesh->getFaceAt(iface);
        face->setID(iface);
        face->setVisitBit(0);
    }

    int numComponents = 0;

    while (1) {
        faceQ.clear();
        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt(iface);
            if (face->isActive() && !face->getVisitBit()) {
                face->setAttribute("Component", numComponents);
                faceQ.push_back(face);
                break;
            }
        }

        if (faceQ.empty())
            break;

        while (!faceQ.empty()) {
            JFacePtr face = faceQ.front();
            faceQ.pop_front();
            if (!face->getVisitBit()) {
                face->setAttribute( "Component", numComponents);
                face->setVisitBit(1);
                int nedges = face->getSize(1);
                for (int j = 0; j < nedges; j++) {
                    const JEdgePtr &edge = face->getEdgeAt(j);
                    int proceed = 1;
                    if (proceed) {
                        JEdge::getRelations( edge, fneighs );
                        int numneighs = fneighs.size();
                        for( int k = 0; k < numneighs; k++) {
                            if( !fneighs[k]->getVisitBit() )
                                faceQ.push_back(fneighs[k]);
                        }
                    }
                }
            }
        } // Complete one Component
        numComponents++;
    }

    for (size_t iface = 0; iface < numfaces; iface++) {
        face = mesh->getFaceAt(iface);
        if (face->isActive() && !face->getVisitBit())
            mesh->getLogger()->setError("Component search incomplete" );
    }

    components.resize(numComponents);
    for( int i = 0; i < numComponents; i++)
        components[i] = getFaceComponent(i);

    auto comp = [] (const JMeshPtr &a, const JMeshPtr &b)
    {
        return a->getSize(2) < b->getSize(2);
    };

    sort( components.begin(), components.end(), comp);
}
///////////////////////////////////////////////////////////////////////////////////////////

void
JMeshComponents :: searchCellComponents()
{
    JCellPtr cell;
    JCellSequence cneighs;
    cell.reset();
    deque<JCellPtr> cellQ;

    size_t numcells = mesh->getSize(3);

    for (size_t icell = 0; icell < numcells; icell++) {
        cell = mesh->getCellAt(icell);
        cell->setID(icell);
        cell->setVisitBit(0);
    }

    int numComponents = 0;
    while (1) {
        cell.reset();
        cellQ.clear();
        for (size_t icell = 0; icell < numcells; icell++) {
            cell = mesh->getCellAt(icell);
            if( cell->isActive() &&!cell->getVisitBit()) {
                cell->setAttribute("Component", numComponents);
                cellQ.push_back(cell);
                break;
            }
        }

        if (cellQ.empty())
            break;

        while (!cellQ.empty()) {
            JCellPtr cell = cellQ.front();
            cellQ.pop_front();
            if (!cell->getVisitBit()) {
                cell->setAttribute( "Component", numComponents);
                cell->setVisitBit(1);
                int nfaces = cell->getSize(2);
                for (int j = 0; j < nfaces; j++) {
                    const JFacePtr &face = cell->getFaceAt(j);
                    int proceed = 1;
                    if (proceed) {
                        JFace::getRelations( face, cneighs );
                        int numneighs = cneighs.size();
                        for( int k = 0; k < numneighs; k++) {
                            if( !cneighs[k]->getVisitBit() )
                                cellQ.push_back(cneighs[k]);
                        }
                    }

                }
            }
        } // Complete one Component
        numComponents++;
    }

    for (size_t icell = 0; icell < numcells; icell++) {
        cell = mesh->getCellAt(icell);
        if( cell->isActive() && !cell->getVisitBit())
            mesh->getLogger()->setError("Component search incomplete " );
    }

    components.resize(numComponents);
    for( int i = 0; i < numComponents; i++)
        components[i] = getCellComponent(i);

    auto comp = [] ( const JMeshPtr &a, const JMeshPtr &b)
    {
        return a->getSize(3) < b->getSize(3);
    };
    sort( components.begin(),components.end(), comp);
}

////////////////////////////////////////////////////////////////////////////////////

void JMeshComponents :: searchComponents()
{
    if( mesh == nullptr) return;
    if( dim == 2 ) searchFaceComponents();
    if( dim  == 3) searchCellComponents();
}

////////////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshComponents::getFaceComponent(int id) const
{
    JFaceSequence faces;
    size_t numfaces = mesh->getSize(2);

    int compId = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err =  face->getAttribute( "Component", compId);
            if( (err == 0) && (compId == id)) faces.push_back(face);
        }
    }

    JNodeSequence nodes;
    JMeshTopology::getEntitySet(faces, nodes);
    JMeshPtr submesh = JMesh::newObject();

    submesh->addObjects(nodes);
    submesh->addObjects(faces);
    return submesh;
}
////////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshComponents::getCellComponent(int id) const
{
    JCellSequence cells;
    size_t numcells = mesh->getSize(3);

    int compId = 0;
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            int err =  cell->getAttribute( "Component", compId);
            if( err == 0 && compId == id) cells.push_back(cell);
        }
    }
    JNodeSequence nodes;
    JMeshTopology::getEntitySet(cells, nodes);
    JMeshPtr submesh = JMesh::newObject();
    submesh->addObjects(nodes);
    submesh->addObjects(cells);
    return submesh;

}

JMeshPtr JMeshComponents::getComponent(int id) const
{
    if( id >= (int)components.size() ) return nullptr;
    return components[id];
}
///////////////////////////////////////////////////////////////////////////////
