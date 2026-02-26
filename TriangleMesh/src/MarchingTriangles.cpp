#include "MarchingTriangles.hpp"

void JMarchingTriangles :: setInsideOutside(double scalarVal)
{
    if( mesh == NULL ) return;

    double val = 0.0;
    int  scalsign = 0;
    size_t numnodes =  mesh->getSize(0);
    for( size_t i = 0; i <  numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute( attribname, val);
            if( scalarVal == val ) scalsign = 0;
            if( scalarVal <  val ) scalsign = -1;
            if( scalarVal >  val ) scalsign =  1;
            vtx->setAttribute("InOutSign", scalsign);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
void JMarchingTriangles :: getContours(  int nCount, vector<JEdgeSequence> &contours)
{
    if( mesh == NULL ) return;

    contours.clear();

    size_t numnodes = mesh->getSize(0);
    vector<double> darray;

    double scalar;

    darray.reserve(numnodes);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            vertex->getAttribute(attribname, scalar);
            darray.push_back(scalar);
        }
    }
    double minval = *min_element( darray );
    double maxval = *max_element( darray );

    double dval   = (maxval-minval)/(double)nCount;

    JEdgeSequence eseq;

    size_t  numfaces =  mesh->getSize(2);
    for(int j = 0; j < nCount; j++) {
        double scalar =  minval + j*dval;
        for( size_t i = 0; i < numfaces; i++) {
            checkTriangle( mesh->getFaceAt(i), scalar, eseq );
            contours.push_back(eseq);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//
void JMarchingTriangles :: getContour(  double scalarVal, JEdgeSequence &eseq)
{
    if( mesh == NULL ) return;

    eseq.clear();

    size_t  numfaces =  mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++)
        checkTriangle( mesh->getFaceAt(i), scalarVal, eseq );
}

////////////////////////////////////////////////////////////////////////////////

void JMarchingTriangles :: checkTriangle(  const JFacePtr face, double scalar, JEdgeSequence &eseq)
{
    JNodePtr v[3];
    v[0] = face->getNodeAt(0);
    v[1] = face->getNodeAt(1);
    v[2] = face->getNodeAt(2);

    double f[3];
    v[0]->getAttribute(attribname, f[0] );
    v[1]->getAttribute(attribname, f[1] );
    v[2]->getAttribute(attribname, f[2] );

    JNodePtr newnodes[3];
    newnodes[0] = NULL;
    newnodes[1] = NULL;
    newnodes[2] = NULL;

    Point3D  p1, p2, p3;
    for( int i = 0; i <  3; i++) {
        int i1 =   (i+1)%3;
        int i2 =   (i+2)%3;
        double minval =   std::min( f[i1], f[i2] );
        double maxval =   std::max( f[i1], f[i2] );
        if( (minval < scalar) && (maxval > scalar)) {
            double t =  (scalar - f[i1] )/( f[i2] - f[i1]);
            p1      = v[i1]->getXYZCoords();
            p2      = v[i2]->getXYZCoords();
            p3[0]   = t*p2[0] + (1-t)*p1[0];
            p3[1]   = t*p2[1] + (1-t)*p1[1];
            p3[2]   = t*p2[2] + (1-t)*p1[2];
            newnodes[i1] = JNode::newObject();
            newnodes[i1]->setXYZCoords(p3);
        }
    }

    JEdgePtr edge;
    for( int i = 0; i < 3; i++) {
        int i1 =   (i+1)%3;
        int i2 =   (i+2)%3;
        if( newnodes[i1] && newnodes[i2] ) {
            edge = JEdge::newObject(newnodes[i1], newnodes[i2] );
            eseq.push_back(edge);
        }
    }
}


