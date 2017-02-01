#include "MeshPaver.hpp"

///////////////////////////////////////////////////////////////////////////////

void JMeshPaver :: init()
{
    boundnodes.clear();
    boundedges.clear();

    if( mesh == nullptr ) return;

    mesh->getTopology()->getBoundary( boundnodes);
    mesh->getTopology()->getBoundary( boundedges);
    int topDim = mesh->getTopology()->getDimension();
    elemType = mesh->getTopology()->getElementsType(topDim);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPaver :: execute()
{
    if( mesh == nullptr ) return;

    mesh->buildRelations(0,0);

    size_t numfaces = mesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nedges = face->getSize(1);
            for( int j = 0; j < nedges; j++) {
                JEdgePtr edge = face->getEdgeAt(j);
                if( edge->isBoundary() ) {
                    if( face->getOrientation(edge) == -1) {
                        cout << "Edge has reverse orientation " << endl;
                        edge->reverse();
                    }
                }
            }
        }
    }

    init();
    for( int i = 0; i < numLayers; i++)
        addNewLayer();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPaver :: addNewLayer()
{
    int nSize = boundnodes.size();

    JNodeSequence newnodes = JNode::newObjects( nSize );
    map<JNodePtr,JNodePtr> vertexMap;
    JNodeSequence vneighs;

    double nx, ny;
    Point3D xc, xyz;

    // Just put the old vertex slightly shifted inside so that the area of the
    // new face is non-zero...
    for( int i = 0; i < nSize; i++) {
        vertexMap[ boundnodes[i] ] = newnodes[i];
        JNode::getRelations(boundnodes[i], vneighs);
        xc = boundnodes[i]->getXYZCoords();
        newnodes[i]->setXYZCoords( xc );
        /*
        if( vneighs.empty() )  return;
        for( size_t j = 0; j < vneighs.size(); j++) {
             if( vneighs[j]->isBoundary() ) {
             xyz = vneighs[j]->getXYZCoords();
             nx += xyz[0] - xc[0];
             ny += xyz[1] - xc[1];
             }
        }
        double mag = sqrt(nx*nx + ny*ny);
        nx = nx/mag;
        ny = ny/mag;
        xyz = boundnodes[i]->getXYZCoords();
        xyz[0] += 0.000001*nx;
        xyz[1] += 0.000001*ny;
        newnodes[i]->setXYZCoords( xyz );
        */
        /*
        newnodes[i]->setBoundaryMark(boundnodes[i]->getBoundaryMark());
        boundnodes[i]->setBoundaryMark(0);
        if( boundnodes[i]->isBoundary() ) {
            cout << "Warning: Old node still considered as boundary" << endl;
        }
        */
        cout << "Debug exit " << endl;
        exit(0);
    }
    mesh->addObjects( newnodes);

    int numedges = boundedges.size();

    JFaceSequence newfaces;
    JEdgeSequence newedges;
    newedges.resize(numedges);

    JFacePtr newface;
    for( int i = 0; i < numedges; i++) {
        JNodePtr v0 = boundedges[i]->getNodeAt(0);
        JNodePtr v1 = boundedges[i]->getNodeAt(1);
        JNodePtr v2 = vertexMap[v1];
        JNodePtr v3 = vertexMap[v0];
        if( elemType == JFace::QUADRILATERAL ) {
            newface = JQuadrilateral::newObject( v3,v2,v1,v0);
            newfaces.push_back(newface);
        }

        if( elemType == JFace::TRIANGLE ) {
            if( hybrid ) {
                newface = JQuadrilateral::newObject( v3,v2,v1,v0);
                newfaces.push_back(newface);
            } else {
                newface = JTriangle::newObject( v3,v2,v1);
                newfaces.push_back(newface);
                newface = JTriangle::newObject( v3,v1,v0);
                newfaces.push_back(newface);
            }
        }
        newedges[i] = JEdge::newObject(v3,v2);
        cout << "Debug exit " << endl;
        exit(0);
        /*
                 newedges[i]->setBoundaryMark(boundedges[i]->getBoundaryMark());
                 boundedges[i]->setBoundaryMark(0);
                 if( boundedges[i]->isBoundary() ) {
                     cout << "Warning: Old edge still considered as boundary" << endl;
                 }
        */
    }
    mesh->addObjects( newedges );
    mesh->addObjects( newfaces );

    int bmark;
    for( size_t i = 0; i < newnodes.size(); i++) {
        int err = boundnodes[i]->getAttribute("Constraint", bmark);
        if( !err) {
            boundnodes[i]->deleteAttribute("Constraint");
            newnodes[i]->setAttribute("Constraint", bmark);
        }
        if( !newnodes[i]->isBoundary() )
            cout << "Warning: new node is not marked as boundary " << endl;
    }

    for( size_t i = 0; i < newedges.size(); i++) {
        int err = boundedges[i]->getAttribute("Constraint", bmark);
        if( !err) {
            boundedges[i]->deleteAttribute("Constraint");
            newedges[i]->setAttribute("Constraint", bmark);
        }
        if( !newedges[i]->isBoundary() )
            cout << "Warning: new edge is not marked as boundary " << endl;
    }


    boundnodes = newnodes;
    boundedges = newedges;
}

///////////////////////////////////////////////////////////////////////////////
