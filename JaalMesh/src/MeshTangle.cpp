#include "MeshTangle.hpp"

///////////////////////////////////////////////////////////////////////////////
size_t JMeshTangle :: getNumInvertedElements()
{
    searchNegativeFaces();
    return negativeFaces.size();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshTangle :: searchOverlap()
{
    searchNegativeFaces();
    searchOverlapFaces();
    genIntersectPoints();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTangle:: searchNegativeFaces()
{
    negativeFaces.clear();

    if( mesh == nullptr ) return;

    double x[4], y[4];
    vector<Point3D> p;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        face->getXYZCoords(p);
        int n = face->getSize(0);
        for( int j = 0; j < n; j++) {
            x[j] = p[j][0];
            y[j] = p[j][1];
        }
        double area = JGeometry::getSignedArea(x, y, n);
        if( area < 0.0) negativeFaces.push_back( face );
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTangle ::  searchOverlapFaces()
{
    switch( method) {
    case BRUTE_SEARCH:
        brute_overlap_search();
        break;
    case BREADTH_SEARCH:
        breadth_overlap_search();
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////
vector< pair<JFacePtr,JFacePtr> > JMeshTangle :: getTangledFaces()
{
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTangle :: brute_overlap_search()
{
    searchNegativeFaces();
    edgePairs.clear();
    facePairs.clear();
    intersectPoints.clear();

    Point2D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;

    size_t numedges =  mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr ei = mesh->getEdgeAt(i);
        for( size_t j = i+1; j < numedges; j++) {
            JEdgePtr ej = mesh->getEdgeAt(j);
            int stat  = JEdgeGeometry::intersectPredicate2d(ei, ej);
            if( stat > 0) {
                JEdgeGeometry::intersectPos2d(ei, ej, &xyz[0] );
                intersectPoints.push_back(xyz);
            }
            if( stat >= 0) {
                JEdgePtr minedge = min( ei,ej);
                JEdgePtr maxedge = max( ei,ej);
                edgePairs.insert( pair<JEdgePtr,JEdgePtr>(minedge,maxedge));
            }
        }
    }

    size_t numfaces =  mesh->getSize(2);
    for( size_t jface = 0; jface < numfaces; jface++)  {
        JFacePtr face1 = mesh->getFaceAt(jface);
        for( size_t iface = jface+1; iface < numfaces; iface++)  {
            JFacePtr face2 = mesh->getFaceAt(iface);
            if(JFaceGeometry::intersectPredicate2d(face1, face2)) {
                JFacePtr minface = min(face1,face2);
                JFacePtr maxface = max(face1,face2);
                facePairs.insert( pair<JFacePtr,JFacePtr>(minface,maxface));
            }
        }
    }

    tangledFaces.clear();
    JFaceSet faceSet;
    for( std::pair<JFacePtr,JFacePtr> key : facePairs) {
        faceSet.insert( key.first  );
        faceSet.insert( key.second );
    }
    boost::copy( faceSet, back_inserter(tangledFaces));

    intersectEdges.clear();
    JEdgeSet edgeSet;
    for( std::pair<JEdgePtr,JEdgePtr> key : edgePairs) {
        edgeSet.insert( key.first  );
        edgeSet.insert( key.second );
    }

    boost::copy( edgeSet, back_inserter(intersectEdges));
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTangle :: breadth_overlap_search()
{
    facePairs.clear();
    tangledFaces.clear();

    if( mesh == nullptr ) return;
    if( mesh->getAdjTable(0,2) == 0) mesh->buildRelations(0,2);

    int numNegative = negativeFaces.size();
    if( numNegative == 0) return;

    for( JFacePtr face1: negativeFaces) {
        rangeSearch.setQueryRegion( face1 );
        JFaceSequence faces = rangeSearch.getOverlappedFaces();
        for(JFacePtr face2: faces) {
            facePairs.insert(pair<JFacePtr,JFacePtr>(face1,face2));
        }
    }

    JFaceSet faceSet;
    for( std::pair<JFacePtr,JFacePtr> key : facePairs) {
        faceSet.insert( key.first  );
        faceSet.insert( key.second );
    }
    boost::copy( faceSet, back_inserter(tangledFaces));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTangle :: searchOverlapEdges()
{

#ifdef BUGGY
    edgePairs.clear();
    JEdgeSequence edges;

    if( facePairs.size() == 0) return edges;

    for( std::pair<Face*,Face*> key : facePairs) {
        Face *face1 = key.first;
        Face *face2 = key.second;
        int numEdges1 = face1->getSize(1);
        int numEdges2 = face2->getSize(1);
        for( int i = 0; i < numEdges1; i++) {
            Edge *edge1 = face1->getEdgeAt(i);
            for( int j = 0; j < numEdges2; j++) {
                Edge *edge2 = face2->getEdgeAt(j);
                if( edge1 != edge2 ) {
                    bool val = JEdgeGeometry::intersect2d( edge1, edge2);
                    if( val ) edgePairs.insert(pair<Edge*,Edge*>(edge1,edge2));
                }
            }
        }
    }

    for( std::pair<Edge*,Edge*> key : edgePairs) {
        edges.push_back( key.first );
        edges.push_back( key.second );
    }
    return edges;
#endif


}

///////////////////////////////////////////////////////////////////////////////

void JMeshTangle :: genIntersectPoints()
{
    /*
        intersectPoints.clear();
        if( edgePairs.size() == 0) return;

        Point3D xyz;
        xyz[0] = 0.0;
        xyz[1] = 0.0;
        xyz[2] = 0.0;

        for( std::pair<Edge*,Edge*> key : edgePairs) {
            Edge *edge1 = key.first;
            Edge *edge2 = key.second;
            double u = EdgeGeometry::intersect2d( edge1, edge2, &xyz[0]);
            if( u >= 0.0 && u <= 1.0) intersectPoints.push_back(xyz);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTangle:: getIntersectionOf( const JFacePtr face1, const JFacePtr face2,
                                     vector<Point2D> &polyPnts)
{
    vector<Point2D> polyA, polyB, x, y;

    int np = face1->getSize(0);
    polyA.resize(np);
    for( int i = 0; i < np; i++) {
        const Point3D &p = face1->getNodeAt(i)->getXYZCoords();
        polyA[i][0] = p[0];
        polyA[i][1] = p[1];
    }

    np = face2->getSize(0);
    polyB.resize(np);
    for( int i = 0; i < np; i++) {
        const Point3D &p = face2->getNodeAt(i)->getXYZCoords();
        polyB[i][0] = p[0];
        polyB[i][1] = p[1];
    }

    int err = polyBool.getIntersectionOf(polyA, polyB, polyPnts);
    return err;

}
///////////////////////////////////////////////////////////////////////////////
int JMeshTangle :: random_tangle( JMeshPtr mesh, size_t ntangle, vector<size_t> &permute)
{
    int numnodes = mesh->getSize(0);

    permute.resize(numnodes);
    for( int i = 0; i < numnodes; i++)
        permute[i] = i;

    JNodeSequence internalnodes;
    for( int i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( !vertex->isBoundary() )
            internalnodes.push_back(vertex);
    }

    boost::random_shuffle( internalnodes );

    int numpairs = min(ntangle, internalnodes.size()/2);

    for( int i = 0; i < numpairs; i++) {
        int v1 = internalnodes[2*i]->getID();
        int v2 = internalnodes[2*i+1]->getID();
        const Point3D &p1 = internalnodes[2*i]->getXYZCoords();
        const Point3D &p2 = internalnodes[2*i+1]->getXYZCoords();
        internalnodes[2*i]->setXYZCoords(p2);
        internalnodes[2*i+1]->setXYZCoords(p1);
        permute[v1] = v2;
        permute[v2] = v1;
    }
    return 0;
}

