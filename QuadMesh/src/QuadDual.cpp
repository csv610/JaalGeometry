#include "QuadDual.hpp"

JLogger*  JQuadDual::logger = JLogger::getInstance();

///////////////////////////////////////////////////////////////////////////////
void JQuadDual :: initialize()
{
    if( mesh == nullptr) return;
    mesh->buildRelations(0,1);
    mesh->buildRelations(0,2);
    mesh->buildRelations(1,2);
    mesh->getTopology()->searchBoundary();
}

///////////////////////////////////////////////////////////////////////////////

JQuadChordPtr JQuadDual :: getChord(const JEdgePtr &edge)
{
    JQuadChordPtr chord(new JQuadChord);
    chord->setMesh(mesh);
    chord->setSeed(edge);
    return chord;
}

//////////////////////////////////////////////////////////////////////////////

vector<JQuadChordPtr> JQuadDual :: getChords(JEdgeSequence &edges)
{
    vector<JQuadChordPtr> chords;

    size_t nSize = edges.size();
    if( nSize == 0) return  chords;

    chords.reserve( nSize );

    JEdgePtr seedEdge, edge = nullptr;
    int val = 0;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) edge->setAttribute("ChordID", val );
    }
    JEdgeSequence chedges;

    int  currID = 1;
    while(1) {
        numedges = edges.size();
        seedEdge = nullptr;
        for( const JEdgePtr &edge : edges) {
            if( edge->isActive() ) {
                edge->getAttribute("ChordID", val );
                if( val == 0) {
                    seedEdge = edge;
                    break;
                }
            }
        }

        if( seedEdge == nullptr ) break;

        JQuadChordPtr ch = dynamic_pointer_cast<JQuadChord>(getChord(seedEdge));
        chedges.clear();
        if( ch ) {
            chedges = ch->getEdges();
            chords.push_back(ch);
        }
        numedges = chedges.size();
        for( size_t i = 0; i < numedges; i++) {
            chedges[i]->setAttribute("ChordID", currID);
        }
        currID++;
        if( size_t(currID) == mesh->getSize(1) ) break;
    }

    mesh->deleteEdgeAttribute("ChordID");
    return chords;
}

//////////////////////////////////////////////////////////////////////////////

vector<JQuadChordPtr> JQuadDual :: getColumn(const JFacePtr &face)
{
    vector<JQuadChordPtr> chords;
    JQuadChordPtr ch1 = getChord(face->getEdgeAt(0));
    /*
        if( ch1 ) {
            chords.push_back(ch1);
             if( !ch1->hasEdge(face->getEdgeAt(1) ) ) {
                   JQuadChordPtr ch2 = getChord(face->getEdgeAt(1));
                   if( ch2 ) chords.push_back(ch2);
              }
         }
    */
    return chords;
}
//////////////////////////////////////////////////////////////////////////////

vector<JQuadChordPtr> JQuadDual :: getColumns( JFaceSequence &seeds)
{
    /*
        vector<JDualChordPtr>  chords, twoChords;

        for( size_t i = 0; i < seeds.size(); i++) {
            twoChords = getColumn( seeds[i]);
            chords.push_back( twoChords[0] );
            if( twoChords.size() > 1)
                chords.push_back( twoChords[1] );
        }
        return chords;
    */
}

///////////////////////////////////////////////////////////////////////////////

vector<JQuadChordPtr> JQuadDual :: getAllChords()
{
    vector<JQuadChordPtr> nochords;
    if( mesh == nullptr) return  nochords;
    JEdgeSequence alledges = mesh->getEdges();
    return getChords( alledges);
}

///////////////////////////////////////////////////////////////////////////////
vector<JQuadChordPtr> JQuadDual :: getAllCyclicChords()
{
    vector<JQuadChordPtr>  allchords;
    allchords = getAllChords();

    vector<JQuadChordPtr>  cyclicChords;
    for( const JQuadChordPtr &c : allchords)
        if( c->isCyclic() ) cyclicChords.push_back(c);
    return cyclicChords;
}
///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
int JQuadDual :: NoIntersection_Refinement( const JFacePtr &face, JNodeSequence &newnodes,
        JFaceSequence &newfaces)
{

    ///////////////////////////////////////////////////////////////////////////
    // Reference :
    // Quadrilateral surface mesh without self intersecting dual cycles for
    // Hexahedral mesh ...
    // Matthias Muller Hannemann.. 2002
    ///////////////////////////////////////////////////////////////////////////


    newnodes.clear();
    newfaces.clear();

    JNodeSequence v(17);

    Point3D xyz;
    face->getAvgXYZ(xyz);
    v[0] = JVertex ::newObject();
    v[0]->setXYZCoords(xyz);
    newnodes.push_back( v[0] );

    v[1] = face->getNodeAt(0);
    v[2] = face->getNodeAt(1);
    v[3] = face->getNodeAt(2);
    v[4] = face->getNodeAt(3);

    JEdgePtr edge;

    edge = face->getEdgeAt(0);
    if( !edge->hasAttribute("Steiner") ) {
        v[5] = JVertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[5]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[5]);
        newnodes.push_back( v[5] );
    }
    edge->getAttribute("Steiner", v[5] );

    edge = face->getEdgeAt(1);
    if( !edge->hasAttribute("Steiner") ) {
        v[6] = JVertex::newObject();
        edge->getAvgXYZ(xyz);
        v[6]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[6]);
        newnodes.push_back( v[6] );
    }
    edge->getAttribute("Steiner", v[6] );

    edge = face->getEdgeAt(2);
    if( !edge->hasAttribute("Steiner") ) {
        v[7] = JVertex::newObject();
        edge->getAvgXYZ(xyz);
        v[7]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[7]);
        newnodes.push_back( v[7] );
    }
    edge->getAttribute("Steiner", v[7] );

    edge = face->getEdgeAt(3);
    if( !edge->hasAttribute("Steiner") ) {
        v[8] = JVertex ::newObject();
        edge->getAvgXYZ(xyz);
        v[8]->setXYZCoords(xyz);
        edge->setAttribute("Steiner", v[8]);
        newnodes.push_back( v[8] );

    }
    edge->getAttribute("Steiner", v[8] );

    v[9]  = JNodeGeometry::getMidNode(v[0], v[1]);
    v[10] = JNodeGeometry::getMidNode(v[0], v[2]);
    v[11] = JNodeGeometry::getMidNode(v[0], v[3]);
    v[12] = JNodeGeometry::getMidNode(v[0], v[4]);

    v[13]  = JFaceGeometry::getCentroid(v[0], v[1], v[2] );
    v[14]  = JFaceGeometry::getCentroid(v[0], v[2], v[3] );
    v[15]  = JFaceGeometry::getCentroid(v[0], v[3], v[4] );
    v[16]  = JFaceGeometry::getCentroid(v[0], v[4], v[1] );

    newnodes.push_back( v[9] );
    newnodes.push_back( v[10] );
    newnodes.push_back( v[11] );
    newnodes.push_back( v[12] );
    newnodes.push_back( v[13] );
    newnodes.push_back( v[14] );
    newnodes.push_back( v[15] );
    newnodes.push_back( v[16] );

    newfaces.resize(12);
    newfaces[0] = Quadrilateral::newObject( v[0], v[9], v[13], v[10] );
    newfaces[1] = Quadrilateral::newObject( v[1], v[5], v[13], v[9] );
    newfaces[2] = Quadrilateral::newObject( v[2], v[10], v[13], v[5] );

    newfaces[3] = Quadrilateral::newObject( v[0], v[10], v[14], v[11] );
    newfaces[4] = Quadrilateral::newObject( v[2], v[6], v[14], v[10] );
    newfaces[5] = Quadrilateral::newObject( v[3], v[11], v[14], v[6] );

    newfaces[6] = Quadrilateral::newObject( v[0], v[11], v[15], v[12] );
    newfaces[7] = Quadrilateral::newObject( v[3], v[7], v[15], v[11] );
    newfaces[8] = Quadrilateral::newObject( v[4], v[12], v[15], v[7] );

    newfaces[9]  = Quadrilateral::newObject( v[0], v[12], v[16], v[9] );
    newfaces[10] = Quadrilateral::newObject( v[4], v[8], v[16], v[12] );
    newfaces[11] = Quadrilateral::newObject( v[1], v[9], v[16], v[8] );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
#endif

