#include "QuadDominant2PureQuadsMesher.hpp"

void JQuadDominant2PureQuadsMesher :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

}

/////////////////////////////////////////////////////////////////////////////
JMeshPtr JQuadDominant2PureQuadsMesher :: getCatmullClarkMesh()
{
    AllQuadMeshGenerator allquads;
    allquads.setMesh(mesh);
    JMeshPtr quadMesh = allquads.getCatmullClarkMesh();
    return quadMesh;
}
/////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsMesher :: setOnlyTrianglesAsNonQuads()
{
    if( mesh == nullptr ) return;

    size_t numfaces = mesh->getSize(2);

    int err;
    JNodeSequence  connect;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive()) {
            int nn = face->getSize(0);
            if( nn == 5) {
                connect = face->getNodes();
                JFacePtr quad = JQuadrilateral::newObject( connect[0], connect[1], connect[2], connect[3] );
                JFacePtr tri  = JTriangle::newObject( connect[0], connect[3], connect[4]);
                err = mesh->addObject(quad);
                assert(!err);
                err = mesh->addObject(tri);
                assert(!err);
                face->setStatus(JMeshEntity::REMOVE);
            }
        }
    }
    mesh->pruneFaces();
}
/////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsMesher :: searchNonQuads()
{
    if( mesh == nullptr ) return;

    nonQuadsSet.clear();

    size_t numfaces = mesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) != 4)
            nonQuadsSet.insert(face);
    }
}

/////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsMesher :: genDualGraph()
{
    dualGraph = JMesh::newObject();

    JNodePtr dualnode;
    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("DualNode", dualnode);
            dualnode = JNode::newObject();
            dualnode->setID(i);
            face->setAttribute("DualNode", dualnode);
            dualnode->setAttribute("PrimalFace", face);
            dualGraph->addObject(dualnode);
        }
    }

    JNodePtr v0, v1;
    JFaceSequence faceneighs;
    JEdgePtr dualedge;
    size_t numEdges = mesh->getSize(1);

    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            int err = edge->getAttribute("DualEdge", dualedge);
            if( err ) {
                JEdge::getRelations(edge, faceneighs);
                if( faceneighs.size() == 2 ) {
                    faceneighs[0]->getAttribute("DualNode", v0);
                    faceneighs[1]->getAttribute("DualNode", v1);
                    dualedge = JEdge::newObject(v0,v1);
                    dualGraph->addObject(dualedge);
                }
            }
        }
    }

    djk.reset( new JDijkstraShortestPath);
    djk->setMesh(dualGraph);
    dualGraph->buildRelations(0,0);
}

/////////////////////////////////////////////////////////////////////////////

JFaceSequence JQuadDominant2PureQuadsMesher :: getNewStrip()
{
    newStrip.clear();
    searchNewStrip();
    return newStrip;
}

/////////////////////////////////////////////////////////////////////////////
JFaceSequence JQuadDominant2PureQuadsMesher :: getNewStrip(const JFacePtr &src)
{
    newStrip.clear();
    searchNewStrip(src);
    return newStrip;
}
/////////////////////////////////////////////////////////////////////////////

JFaceSequence JQuadDominant2PureQuadsMesher :: getNewStrip(const JFacePtr &src, const JFacePtr &dst)
{
    newStrip.clear();
    searchNewStrip(src,dst);
    return newStrip;
}

/////////////////////////////////////////////////////////////////////////////
bool JQuadDominant2PureQuadsMesher :: checkStrip()
{
    // A Quad Strip contains minimum of two faces ..
    if( newStrip.size() < 2) return 1;

    // A Quad Strip is always open  ...
    if( newStrip.front() == newStrip.back() ) return 1;

    // Both the end nodes are irregular ,,,,
    if( newStrip.front()->getSize(0) == 4) return 1;
    if( newStrip.back()->getSize(0)  == 4) return 1;


    // All the interior faces are regular ...
    size_t numfaces = newStrip.size();
    for( size_t i = 1; i < numfaces-1; i++)
        if( newStrip[i]->getSize(0)  != 4) return 1;

    return 0;

}
/////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsMesher :: searchNewStrip()
{
    if( nonQuadsSet.empty() ) searchNonQuads();

    newStrip.clear();

    if( nonQuadsSet.size() < 2) return;

    JFacePtr nonQuad1 = *nonQuadsSet.begin();

    searchNewStrip( nonQuad1);
}
/////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsMesher :: searchNewStrip( const JFacePtr &nonQuad1)
{
    newStrip.clear();
    if( nonQuad1 == nullptr) return;
    if( nonQuad1->getSize(0) == 4) return;

    JFacePtr nonQuad2;

    // Search for the next non-quad using Breadth-first search.
    JFaceSet facesVisited;

    deque<JFacePtr> faceQ;
    faceQ.push_back(nonQuad1);
    JFaceSequence faceneighs;

    while(!faceQ.empty() ) {
        JFacePtr currface = faceQ.front();
        faceQ.pop_front();
        if( currface->getSize(0) != 4) {
            if( currface != nonQuad1) {
                if( nonQuadsSet.find(currface) != nonQuadsSet.end() ) {
                    nonQuad2 = currface;
                    break;
                }
            }
        }
        if( facesVisited.find(currface) == facesVisited.end() ) {
            facesVisited.insert(currface);
            JFace::getRelations12(currface, faceneighs);
            for( const JFacePtr &f : faceneighs) {
                if( facesVisited.find(f) == facesVisited.end() )
                    faceQ.push_back(f);
            }
        }
    }
    if( nonQuad2 == nullptr) {
        nonQuadsSet.erase(nonQuad1);
        return;
    }

    searchNewStrip( nonQuad1, nonQuad2);
}

///////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsMesher :: searchNewStrip( const JFacePtr &nonQuad1, const JFacePtr &nonQuad2)
{
    newStrip.clear();
    if( (nonQuad1 == nullptr) || (nonQuad2 == nullptr)) return;
    if( nonQuad1->getSize(0) == 4) return;
    if( nonQuad2->getSize(0) == 4) return;

    if( dualGraph == nullptr) genDualGraph();

    JNodePtr vsrc;
    nonQuad1->getAttribute("DualNode", vsrc);

    JNodePtr vdst;
    nonQuad2->getAttribute("DualNode", vdst);

    JNodeSequence path = djk->getPath(vsrc, vdst);

    JFacePtr pface;
    for(const JNodePtr &dnode: path) {
        dnode->getAttribute("PrimalFace", pface);
        newStrip.push_back( pface );
    }
    nonQuadsSet.erase(nonQuad1);
    nonQuadsSet.erase(nonQuad2);
}
/////////////////////////////////////////////////////////////////////////////

vector<JFaceSequence> JQuadDominant2PureQuadsMesher :: getAllStrips()
{
    vector<JFaceSequence> quadStrips;
    if( mesh == nullptr ) return quadStrips;

    if( nonQuadsSet.empty() ) searchNonQuads();
    if( nonQuadsSet.size() < 2) return quadStrips;

    int nCount = 0;

    int nPairs = nonQuadsSet.size()/2;

    for( int i = 0; i < nPairs; i++) {
        JFacePtr nonQuad1 = *nonQuadsSet.begin();
        searchNewStrip( nonQuad1);
        if( !newStrip.empty()) quadStrips.push_back(newStrip);
    }
    return quadStrips;
}

/////////////////////////////////////////////////////////////////////////////
/*
int JQuadDominant2PureQuadsMesher :: remeshStrip()
{
    if( quadStrips.empty() ) return 1;

    JFaceSequence currSeq = quadStrips[0];
    quadStrips.pop_front();

    int nsize = currSeq.size();
    if( nsize < 2) {
        cout << "Warning: At least two faces are required in the sequence " << endl;
        return 1;
    }

    if( currSeq.front()->getSize(0) != 3 ) {
        cout << "Warning: Begining of the sequence is not a triangle " << endl;
        return 1;
    }

    if( currSeq.back()->getSize(0) != 3 ) {
        cout << "Warning: End of the sequence is not a triangle " << endl;
        return 1;
    }

    for( int i = 1; i < nsize-1; i++) {
        if( currSeq[i]->getSize(0) != 4) {
            cout << "Warning: Intermediate elements must be quad " << endl;
            return 1;
        }
    }

    JEdgePtr commonEdge, nextEdge, newEdge;
    JEdgeSequence commonEdges;
    JNodeSequence nodes(5);

    deque<JFacePtr> faceQ;
    for( const JFacePtr &f : currSeq) faceQ.push_back(f);

    JFacePtr newQuad, newTri;

    for( int i = 0; i < nsize-2; i++) {
        const JFacePtr &tri  = faceQ[0];        // Must be triangle ...
        const JFacePtr &quad = faceQ[1];       // Must be Quad ..
        const JFacePtr &nextFace = faceQ[2];   // Can be any face ...
        faceQ.pop_front();
        faceQ.pop_front();

        // First face must be triangle and second facemust be quad.
        if( (tri->getSize(0) != 3) && (quad->getSize(0) != 4) ) break;
        JFace::getSharedEntities( tri, quad, commonEdges);
        if( commonEdges.size() != 1 ) return 1;
        commonEdge = commonEdges[0];

        // With respect to the quad, what is the position of "commonEdge"
        int pos0 = quad->getPosOf(commonEdge);
        assert( pos0 >= 0 && pos0 < 4);

        // Which edge the Quad and Next face sharing ...
        JFace::getSharedEntities( quad, nextFace, commonEdges);
        if( commonEdges.size() != 1 ) return 1;
        nextEdge = commonEdges[0];

        // What is the position of "NextEdge" relative to pos0 ? It can be 1,2, or 3
        int pos1 = -1;
        if( quad->getEdgeAt(pos0+1) == nextEdge ) pos1 = 1;
        if( quad->getEdgeAt(pos0+2) == nextEdge ) pos1 = 2;
        if( quad->getEdgeAt(pos0+3) == nextEdge ) pos1 = 3;

        assert( pos1 >= 0);

        nodes[0] = commonEdge->getNodeAt(0);
        nodes[2] = commonEdge->getNodeAt(1);
        nodes[1] = Triangle::getOppositeNode( tri, nodes[0], nodes[2] );
        int pos3 = tri->getPosOf( nodes[1] );
        nodes[2] = tri->getNodeAt(pos3+1);
        nodes[0] = tri->getNodeAt(pos3+2);
        nodes[3] = Quadrilateral::getDiagonalNode( quad, nodes[0] );
        nodes[4] = Quadrilateral::getDiagonalNode( quad, nodes[2] );

        switch(pos1)
        {
        case 1:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[4]);
            newTri  = Triangle::newObject( nodes[2], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[2], nodes[4],1);
            break;
        case 2:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[4]);
            newTri  = Triangle::newObject( nodes[2], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[2], nodes[4],1);
            break;
        case 3:
            newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[3]);
            newTri  = Triangle::newObject( nodes[0], nodes[3], nodes[4] );
            newEdge = JSimplex::getEdgeOf( nodes[0], nodes[3],1);
            break;
        default:
            return 1;
        }

        faceQ.push_front(newTri);

        mesh->addObject(newEdge);
        mesh->addObject(newTri);
        mesh->addObject(newQuad);

        tri->setStatus( JMeshEntity::REMOVE);
        quad->setStatus( JMeshEntity::REMOVE);
        commonEdge->setStatus( JMeshEntity::REMOVE);
    }

    assert( faceQ.size() == 2);
    assert( faceQ[0]->getSize(0) == 3);
    assert( faceQ[1]->getSize(0) == 3);

    JFace::getSharedEntities( faceQ[0], faceQ[1], commonEdges);
    assert( commonEdges.size() == 1);

    nodes[0] = commonEdges[0]->getNodeAt(0);
    nodes[2] = commonEdges[0]->getNodeAt(1);
    nodes[1] = Triangle::getOppositeNode( faceQ[0], nodes[0], nodes[2] );
    nodes[3] = Triangle::getOppositeNode( faceQ[1], nodes[0], nodes[2] );

    newQuad = Quadrilateral::newObject( nodes[0], nodes[1], nodes[2], nodes[3]);
    mesh->addObject(newQuad);
    faceQ[0]->setStatus( JMeshEntity::REMOVE);
    faceQ[1]->setStatus( JMeshEntity::REMOVE);
    commonEdges[0]->setStatus( JMeshEntity::REMOVE);

    mesh->pruneFaces();
    mesh->pruneEdges();
    return 0;
}
//////////////////////////////////////////////////////////////////////////////
*/
