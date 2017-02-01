#pragma once

#include "Mesh.hpp"

using namespace Jaal;

class JCurve {
public:
    static const int OPEN_CURVE  = 0;
    static const int CLOSE_CURVE = 1;
    static JCurve *readFromFile( const string &s);

    JCurve() {
        type = 0;
    }

    void setType( int t ) {
        type = t;
    }

    void copyOf(JCurve *jc) {

        type = jc->type;
        ctrlnodes = jc->ctrlnodes;
        nodes    = jc->nodes;
        edges    = jc->edges;
    }

    int getSize(int e ) const {
        int nsize = nodes.size();
        if( e == 0)
            return nsize;
        if( e == 1) {
            if( type == 0) return nsize-1;
            if( type == 1) return nsize;
        }
        return 0;

    }

    void  mergeNodes(double d);

    void addControlNode( JNodePtr v ) {
        if( !ctrlnodes.empty() ) {
            JEdgePtr newedge = JEdge::newObject( ctrlnodes.back(), v );
            edges.push_back(newedge);
        }
        ctrlnodes.push_back(v);
        nodes.push_back(v);
    }

    void setControlNodes( JNodeSequence &s );
    void getControlNodes( JNodeSequence &s) {
        s = ctrlnodes;
    }


    void finalize() {
        if( ctrlnodes.size() > 2) {
            if( type ) {
                JEdgePtr newedge = JEdge::newObject( ctrlnodes.back(), ctrlnodes.front() );
                edges.push_back(newedge);
            }
        }
    }

    void getNodes( JNodeSequence &s) {
        s = nodes;
    }

    void getEdges( JEdgeSequence &s) {
        s = edges;
    }

    // Between two successive nodes, give the minimum Distance ...
    double getMinimumLocalDistance() const;

    // Between any two nodes, gives the minimum distance ....
    double getMinimumGlobalDistance() const;

    // Redistribue nodes so that in the parametric domain, they have are equidistant..
    void uniformDiscretize();

    // Just subdivide all the segments into same number of segments ...
    void refineSegments( int n );

    // Adaptive discretization. Set the maximum length...
    void setMaxSegmentLength( double l);

    double getLength() const;
    int getTangentAt( const JNodePtr n, Vec3D &v) const;

    void deleteAll();

private:
    int type;
    JNodeSequence ctrlnodes;
    JNodeSequence nodes;
    JEdgeSequence edges;
};

