#pragma once

#include "Mesh.hpp"
//#include "PolygonSimplify.hpp"

class JMeshContour
{
    typedef std::pair<JNodePtr,JNodePtr> NodePair;
public:
    static const int CLOSED_CONTOUR = 0;
    static const int OPEN_CONTOUR   = 1;

    static const int UNIFORM     = 0;
    static const int NON_UNIFORM = 1;

    JMeshContour()
    {
        cornerAngle = 90.0;
        dilationFactor = 2.0;
    }

    int    setMesh( const JMeshPtr &m);

    int    setSegments( const JEdgeSequence &e);
    int    setNodes( const JNodeSequence &e, int type);

    void    setDilationFactor( double d) {
        dilationFactor = d;
    }

    int     getNumComponents() const;
    int     splitAt( const JNodePtr &v, JEdgeSequence &seq1, JEdgeSequence &seq2);
    int     splitAtAngle( double t, vector<JEdgeSequence> &);
    int     getOrientation() const;
    bool    isChain() const;
    bool    startAt( const JEdgePtr &edge);
    bool    startAt( const JNodePtr &edge);
    bool    isTopologicalSimple();
    int     isClosed();
    int     isCloseable();
    void    reverse();
    double  getLength() const;
    double  getArea();
    int     makeCircular();
    bool    isConvex();
    Point2D getCenter();
    bool    isPointInside( const Point2D &p);
    void    smooth( int n = 1);
    double  getDilation( const JNodePtr &v0, const JNodePtr &v1);
    JNodePtr  getNearestNode( double t);

    vector<JEdgeSequence> getLoops();

    JEdgeSequence getChain();
    JEdgeSequence getChain(const JNodePtr &start);
    JEdgeSequence getChain(const JEdgePtr &start);
    JNodeSequence getChainNodes();
    JNodeSequence getCorners( double angle);
    JEdgeSequence getDiscretized(int N, int method = UNIFORM);
    JEdgeSequence getDiscretized( JEdgeSequence &e, int N);
    JEdgeSequence getSimplify( int algo, double tol);

    JBoundingBox   getAxisAlignedBox();
    JHexahedronPtr getMinimumBox();
private:
    double  cornerAngle;
    double  dilationFactor;

    JEdgeSequence edges;
    vector<JEdgeSequence>  components;

    int     getEdgeNumber(const JNodePtr &v);
    void    reduceDilation();
    double  getCurvedDistance( const JNodePtr &p0, const JNodePtr &p1);
    JNodePtr getMidNode(const JNodePtr &v0, const JNodePtr &v1);
};

