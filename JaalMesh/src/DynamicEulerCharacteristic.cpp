#include "DynamicEulerCharacteristic.hpp"


///////////////////////////////////////////////////////////////////////////////

void JDynamicEulerCharacteristic :: reset()
{
    nodeSet.clear();
    edgeSet.clear();
    faceSet.clear();
    cellSet.clear();
    birth_time = 0;
}
///////////////////////////////////////////////////////////////////////////////
int JDynamicEulerCharacteristic :: addObject(const JNodePtr &v)
{
    if( nodeSet.find(v) != nodeSet.end() ) return 1;
    nodeSet.insert(v);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: addObject(const JEdgePtr &e)
{
    if( edgeSet.find(e) != edgeSet.end() ) return 1;
    edgeSet.insert(e);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: addObject(const JFacePtr &f)
{
    if( faceSet.find(f) != faceSet.end() ) return 1;
    faceSet.insert(f);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: addObject(const JCellPtr &c)
{
    if( cellSet.find(c) != cellSet.end() ) return 1;
    cellSet.insert(c);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: removeObject(const JNodePtr &vtx)
{
    nodeSet.erase(vtx);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: removeObject(const JEdgePtr &edge)
{
    JFaceSequence faceneighs;
    JEdge::getRelations(edge, faceneighs);
    if( !faceneighs.empty() ) return 1;
    edgeSet.erase(edge);

    if( edge->isActive() ) edge->setStatus(JMeshEntity::INACTIVE);
    JEdgeSequence edgeneighs;
    for( int i = 0; i < 2; i++)
        removeObject( edge->getNodeAt(i) );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: removeObject(const JFacePtr &face)
{
    JCellSequence cellneighs;
    JFace::getRelations(face, cellneighs);
    if( !cellneighs.empty() ) return 1;
    faceSet.erase(face);

    if( face->isActive() )   face->setStatus(JMeshEntity::INACTIVE);
    JEdgeSequence edgeneighs;
    size_t numedges =  face->getSize(1);
    for( int i = 0; i < numedges; i++)
        removeObject( face->getEdgeAt(i) );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JDynamicEulerCharacteristic :: removeObject(const JCellPtr &cell)
{
    cellSet.erase(cell);
    if( cell->isActive() ) cell->setStatus(JMeshEntity::INACTIVE);

    JFaceSequence faceneighs;
    int numfaces = cell->getSize(2);
    for( int i = 0; i < numfaces; i++) removeObject( cell->getFaceAt(i) );
}
///////////////////////////////////////////////////////////////////////////////////////////////
