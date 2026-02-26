#include "PolyPartitioner.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JPolyPartitioner :: getPartitions(int algo)
{
    if( mesh == nullptr) return nullptr;

    partmesh  = JMesh::newObject();

    TPPLPartition polypart;
    list<TPPLPoly>  inpoly, result;

    createPolyList(inpoly);

    switch( algo )
    {
    case EAR_CLIPPED_POLYGONS:
        polypart.Triangulate_EC( &inpoly, &result);
        break;
    case EDGE_LENGTH_OPT_POLYGONS:
        polypart.Triangulate_EC( &inpoly, &result);
        break;
    case MONOTONE_POLYGONS:
        polypart.Triangulate_MONO( &inpoly, &result);
        break;
    case CONVEX_HM_POLYGONS:
        polypart.ConvexPartition_HM( &inpoly, &result);
        break;
    case CONVEX_OPT_POLYGONS:
        TPPLPoly apoly = inpoly.front();
        polypart.ConvexPartition_OPT( &apoly, &result);
        break;
    }

    JNodeSequence nodes = mesh->getNodes();
    partmesh->addObjects(nodes);

    list<TPPLPoly>::iterator it;
    for( it = result.begin(); it  != result.end(); ++it) createFace( *it );

    partmesh->getTopology()->collectEdges();

    return partmesh;
}

/////////////////////////////////////////////////////////////////////////
int JPolyPartitioner :: createPoly( const JEdgeSequence &edges, TPPLPoly &poly, bool hole)
{
    JNodeSequence nodes;
    JEdgeTopology::getChainNodes( edges, nodes);
    size_t numnodes = nodes.size();

    int ori = JEdgeGeometry::getOrientation(edges);
    if( ori < 0) boost::reverse(nodes);

    if( numnodes ) {
        poly.Init(numnodes);
        if( hole) poly.SetHole(true);
        for( int i = 0; i < numnodes; i++) {
            boundnodes.push_back(nodes[i] );
            const Point2D  xy = nodes[i]->getXYCoords();
            poly[i].x  = xy[0];
            poly[i].y  = xy[1];
        }
    }
}

/////////////////////////////////////////////////////////////////////////
const JNodePtr &JPolyPartitioner :: searchNode( const Point2D &pquery)
{
    for( const JNodePtr &vtx: boundnodes) {
        const Point2D &xy = vtx->getXYCoords();
        if( xy[0] == pquery[0] && xy[1] == pquery[1] ) return vtx;
    }
    cout << "Fatal error: node not found " << endl;
    exit(0);
}

/////////////////////////////////////////////////////////////////////////

int JPolyPartitioner :: createFace( TPPLPoly &poly)
{
    int numpoints = poly.GetNumPoints();

    JNodeSequence nodes(numpoints);

    Point2D xy;
    for( int i = 0; i < numpoints; i++) {
        xy[0] =  poly[i].x;
        xy[1] =  poly[i].y;
        nodes[i] = searchNode( xy );
    }

    JFacePtr newface;
    switch(numpoints)
    {
    case 3:
        newface = JTriangle::newObject( nodes );
        break;
    case 4:
        newface = JQuadrilateral::newObject( nodes );
        break;
    default:
        newface = JPolygon::newObject( nodes );
        break;
    }

    partmesh->addObject(newface);
}

/////////////////////////////////////////////////////////////////////////

int JPolyPartitioner :: createPolyList( list<TPPLPoly> &polylist)
{
    JEdgeSequence boundedges;

    if( mesh->getDimension() == 2 )  {
        mesh->getTopology()->getBoundary(boundedges);
    } else
        boundedges = mesh->getEdges();

    vector<JEdgeSequence> edgeloops;
    JEdgeTopology::getLoops(boundedges, edgeloops);

    boundnodes.clear();

    bool hole = 0;
    TPPLPoly poly;
    int nsize = edgeloops.size();
    for( int i = 0; i < nsize; i++)  {
        createPoly( edgeloops[i], poly, hole);
        polylist.push_back(poly);
    }

    return 0;
}
/////////////////////////////////////////////////////////////////////////

