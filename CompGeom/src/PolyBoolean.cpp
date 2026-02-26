#include "PolyBoolean.hpp"

int JPolyBoolean :: getIntersectionOf( PolyType &pA, PolyType &pB, PolyType &pC)
{
    points.clear();

    int index = 0;
    PolyPoint p;

    double signarea;

    signarea = JGeometry::getSignedArea( pB );
    if( signarea < 0.0) reverse(pB.begin(), pB.end() );

    // Orientation of node of polygon "A" with respect to Polygon "B".
    for( size_t i = 0; i < pA.size(); i++) {
        p.xy   = pA[i];
        p.sign = JGeometry::getBoundedSide(&pB[0][0], &pB[1][0], &pB[2][0], &pA[i][0] );
        p.id   = index++;
        points.push_back(p);
    }

    signarea = JGeometry::getSignedArea( pA );
    if( signarea < 0.0) reverse(pA.begin(), pA.end() );

    // Orientation of node of polygon "B" with respect to Polygon "A".
    for( size_t i = 0; i < pB.size(); i++) {
        p.xy   = pB[i];
        p.sign = JGeometry::getBoundedSide(&pA[0][0], &pA[1][0], &pA[2][0], &pB[i][0] );
        p.id   = index++;
        points.push_back(p);
    }
    assert( points.size() == 6);

    vector< vector<int> > aPoints(3), bPoints(3);
    aPoints.resize(3);
    aPoints.resize(3);

    int numEdges1 = pA.size();
    int numEdges2 = pB.size();
    for( int i = 0; i < numEdges1; i++) {
        const Point2D &p0 = pA[i];
        const Point2D &p1 = pA[(i+1)%numEdges1];
        for( int j = 0; j < numEdges2; j++) {
            const Point2D &p2 = pB[j];
            const Point2D &p3 = pB[(j+1)%numEdges2];
            int sign = JEdgeGeometry::intersectPredicate2d( &p0[0], &p1[0], &p2[0], &p3[0]);
            if( sign > 0) {
                int err = JEdgeGeometry::intersectPos2d( &p0[0], &p1[0], &p2[0], &p3[0], &p.xy[0]);
                if( !err ) {
                    p.sign = 0;
                    p.id   = index++;
                    aPoints[i].push_back( p.id );
                    bPoints[j].push_back( p.id );
                    points.push_back(p);
                }

            }
        }
    }

    /*
        int numpoints = points.size();
        for( int i = 0; i < numpoints; i++) {
            for( int j = i+1; j < numpoints; j++) {
                double dx = points[i].xy[0] - points[j].xy[0];
                double dy = points[i].xy[1] - points[j].xy[1];
                if( dx*dx + dy*dy < 1.E-10) points[j].id = points[i].id;
            }
        }
    */
    exit(0);

    vector<int> poly(3);
    poly[0] = 0;
    poly[1] = 1;
    poly[2] = 2;
    extract_segments( poly, aPoints);




    poly[0] = 3;
    poly[1] = 4;
    poly[2] = 5;
    extract_segments( poly, bPoints);
    for( int i = 0; i < edges.size(); i++)
        cout << edges[i].first << " " << edges[i].second << endl;

    create_chain();

    for( int i = 0; i < edges.size(); i++)
        cout << edges[i].first << " " << edges[i].second << endl;
    exit(0);
}

////////////////////////////////////////////////////////////////////////////////////////////////

void JPolyBoolean :: extract_segments( const vector<int> &poly, const vector< vector<int> > &edgePoints)
{
    size_t nedges = poly.size();
    assert( nedges == edgePoints.size() );

    PolyPoint p;
    vector<PolyPoint> segPoints;
    vector<int> result;

    for( int i = 0; i < nedges; i++) {
        int npoints = edgePoints[i].size();
        int id0 = poly[i];
        int id1 = poly[(i+1)%nedges];

        if( npoints > 1) {
            const Point2D &p0 = points[id0].xy;
            const Point2D &p1 = points[id1].xy;
            segPoints.resize(npoints);
            for( int j = 0; j < npoints; j++) {
                int id = edgePoints[i][j];
                p   = points[id];
                p.u = JEdgeGeometry::getU( &p0[0], &p1[0], &p.xy[0]);
                segPoints[j] = p;
            }
            boost::sort( segPoints );
        }

        result.clear();
        if( points[id0].sign <= 0)
            result.push_back(points[id0].id);

        for( int j = 0; j < npoints; j++)
            result.push_back( segPoints[j].id );

        if( points[id1].sign <= 0)
            result.push_back(points[id1].id);

        npoints = result.size();
        cout << "Ede " << i << " " << npoints << endl;
        if( npoints > 1) {
            for( int j = 0; j < npoints-1; j++) {
                edges.push_back( make_pair(result[j], result[j+1] ));
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////
void JPolyBoolean :: create_chain()
{
    int numedges = edges.size();

    for( int i = 0; i < numedges; i++) {
        int lastnode = edges[i].second;
        for( int j = i+1; j < numedges; j++) {
            int id1 = edges[j].first;
            int id2 = edges[j].second;
            if( id1  == lastnode || id2 == lastnode) {
                swap( edges[i+1], edges[j] );
                if( id2 == lastnode ) {
                    edges[i+1].first  = id2;
                    edges[i+1].second = id1;
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////

int JPolyBoolean :: unitTests()
{
    vector<Point2D>  pA(3), pB(3), pC;

    // Case-I: The oout contains six points...
    pA[0][0] = 0.0;
    pA[0][1] = 0.0;
    pA[1][0] = 1.0;
    pA[1][1] = 0.0;
    pA[2][0] = 0.0;
    pA[2][1] = 1.0;

    pB[0][0] =  0.5;
    pB[0][1] = -0.00;
    pB[1][0] = 0.50000;
    pB[1][1] = 0.5;
    pB[2][0] = 0.00000;
    pB[2][1] = 0.5;
    getIntersectionOf( pA, pB, pC);
}
////////////////////////////////////////////////////////////////////////////////////////////////
