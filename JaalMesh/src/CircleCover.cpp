#include "CircleCover.hpp"

/////////////////////////////////////////////////////////////////////////////////

void JCircleCover :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    int dim = mesh->getTopology()->getDimension();

    JNodeSequence nodes;
    JDelaunayMesh2D delmesh;
    delmesh.setMesh(mesh);
    JMeshPtr medial = delmesh.getMedialAxis();
    nodes = medial->getNodes();

    boost::copy( nodes, std::back_inserter(nodesQ));
}

/////////////////////////////////////////////////////////////////////////////////

int JCircleCover :: getNewCircle( JCircle &circle)
{
    if( nodesQ.empty() ) return 1;

    JNearestNeighbours neighs;
    neighs.setCloud(farNodesQ);

    auto cmp = [] (const JNodePtr &v0, const JNodePtr &v1)
    {
        double d0, d1;
        v0->getAttribute("Radius", d0);
        v1->getAttribute("Radius", d1);
        return d0 > d1;
    };

    size_t numnodes = nodesQ.size();
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0   = nodesQ[i];
        const Point3D  &xyz  = v0->getXYZCoords();
        const JNodePtr &v1   = neighs.getNearest(xyz);
        double r             = JNodeGeometry::getLength(v0,v1);
        v0->setAttribute("Radius", r);
    }
    boost::sort(nodesQ, cmp);

    double r0, r1;
    bool   emptyspace;
    JNodePtr nextNode;
    for( size_t j = 0; j < numnodes; j++) {
        const JNodePtr &v0 = nodesQ[j];
        v0->getAttribute("Radius", r0);
        bool emptyspace = 1;
        for( const JNodePtr &v1 : nodesSelected) {
            double l = JNodeGeometry::getLength(v0,v1);
            v1->getAttribute("Radius", r1);
            if( l < r0 + r1) {
                emptyspace = 0;
                break;
            }
        }
        if( emptyspace ) {
            nextNode = v0;
            break;
        }
    }

    if( nextNode == nullptr) {
        nodesQ.pop_front();
        return 1;
    }

    nodesSelected.push_back(nextNode);

    double radius;
    Point3D center = nextNode->getXYZCoords();
    circle.setCenter(center);
    nextNode->getAttribute("Radius", radius);
    circle.setRadius(radius);

    JNodeSequence withinNodes;
    for( size_t j = 0; j < numnodes; j++) {
        double dist = JNodeGeometry::getLength(nextNode, nodesQ[j] );
        if( dist <= radius) withinNodes.push_back(nodesQ[j] );
    }

    boost::sort( nodesQ );
    boost::sort( withinNodes );

    JNodeSequence outsideNodes;
    boost::set_difference(nodesQ, withinNodes, back_inserter(outsideNodes));
    nodesQ.clear();

    boost::copy( outsideNodes, std::back_inserter(nodesQ));

    if( covermesh == nullptr) return 0;

    newnodes.resize(72);
    newfaces.resize(72);
    JNodeSequence trinodes(3);

    Point3D xyz;
    xyz[2] = 0.0;
    double dtheta = M_PI/36.0;
    for( int j = 0; j < 72; j++) {
        xyz[0] = center[0] + radius*cos( j*dtheta );
        xyz[1] = center[1] + radius*sin( j*dtheta );
        newnodes[j] = JNode::newObject();
        newnodes[j]->setXYZCoords(xyz);
    }

    boost::push_back( farNodesQ, newnodes);

    JNodePtr apex = JNode::newObject();
    apex->setXYZCoords(center);
    for( int j = 0; j < 72; j++) {
        trinodes[0] = apex;
        trinodes[1] = newnodes[j];
        trinodes[2] = newnodes[(j+1)%72];
        newfaces[j] = JFace::newObject(trinodes);
    }
    covermesh->addObject( apex );
    covermesh->addObjects( newnodes );
    covermesh->addObjects( newfaces );

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////

vector<JCircle> JCircleCover :: getCircles()
{
    vector<JCircle>  circles;
    if( mesh == nullptr) return circles;

    farNodesQ.clear();
    mesh->getTopology()->getBoundary( farNodesQ );

    size_t numnodes = nodesQ.size();
    JCircle circle;
    for( size_t i = 0; i < numnodes; i++) {
        int err = getNewCircle( circle );
        if( !err ) circles.push_back(circle);
    }
    return circles;
}

/////////////////////////////////////////////////////////////////////////////////////////

JMeshPtr JCircleCover :: getCoverMesh()
{
    covermesh = JMesh::newObject();
    getCircles();
    return covermesh;
}

/////////////////////////////////////////////////////////////////////////////////////////
