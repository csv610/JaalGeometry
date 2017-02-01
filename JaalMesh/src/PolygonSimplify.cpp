#include "PolygonSimplify.hpp"
////////////////////////////////////////////////////////////////////////////////

int JPolygonSimplify:: getAlgorithmID( const string &str)
{
    if( str == "Douglas-Peucker")  return DOUGLAS_PEUCKER;
    if( str == "Douglas-Peucker-N") return DOUGLAS_PEUCKERN;
    if( str == "Nth Point")        return NTH_POINT;
    if( str == "Opheim")           return   OPHEIM;
    if( str == "Perpendicular Distance") return PERPENDICULAR_DISTANCE;
    if( str == "Radial Distance")  return RADIAL_DISTANCE;
    if( str == "Reumann-Witkam")   return REUMANN_WITKAM;
    if( str == "Lang")             return LANG;
    return -1;
}
////////////////////////////////////////////////////////////////////////////////

void JPolygonSimplify :: setMesh( const JMeshPtr &m)
{
    mesh = m;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JPolygonSimplify :: getSimplified( int alg_, double tol_)
{
    if( mesh == nullptr) return nullptr;
    algo =  alg_;
    tol  =  tol_;

    simplified = JMesh::newObject();
    vector<JEdgeSequence>   contours;
    mesh->getTopology()->getBoundary(contours);

    for( size_t i = 0; i < contours.size(); i++)
        getSimplified( contours[i] );

    return simplified;
}

////////////////////////////////////////////////////////////////////////////////

void  JPolygonSimplify :: getSimplified( const JEdgeSequence &acontour)
{
    vector<double> inpoints, result;

    JNodeSequence bndnodes;
    JEdgeTopology::getChainNodes( acontour, bndnodes);

    if( JEdgeGeometry::getOrientation(acontour) < 0)
        boost::reverse( bndnodes );

    size_t numnodes  = bndnodes.size();
    for( size_t i  = 0; i < numnodes; i++) {
        Point2D  xy = bndnodes[i]->getXYCoords();
        inpoints.push_back( xy[0] );
        inpoints.push_back( xy[1] );
    }

    switch(algo)
    {
    case DOUGLAS_PEUCKER:
        psimpl::simplify_douglas_peucker<2>( inpoints.begin(), inpoints.end(),
                                             tol, std::back_inserter(result));
        break;
    case DOUGLAS_PEUCKERN:
        count = 0.9*numnodes;
        psimpl::simplify_douglas_peucker_n<2>( inpoints.begin(), inpoints.end(),
                                               count, std::back_inserter(result));
        break;
    case NTH_POINT:
        psimpl::simplify_nth_point<2>( inpoints.begin(), inpoints.end(),
                                       2, std::back_inserter(result));
        break;
    case PERPENDICULAR_DISTANCE:
        psimpl::simplify_perpendicular_distance<2>( inpoints.begin(), inpoints.end(),
                tol, std::back_inserter(result));
        break;
    case REUMANN_WITKAM:
        psimpl::simplify_reumann_witkam<2>( inpoints.begin(), inpoints.end(),
                                            tol, std::back_inserter(result));
        break;
    case RADIAL_DISTANCE:
        psimpl::simplify_radial_distance<2>( inpoints.begin(), inpoints.end(),
                                             tol, std::back_inserter(result));
        break;
    case LANG:
        psimpl::simplify_lang<2>( inpoints.begin(), inpoints.end(),
                                  tol, 4, std::back_inserter(result));
        break;
    case OPHEIM:
        psimpl::simplify_lang<2>( inpoints.begin(), inpoints.end(),
                                  tol, 1.5*tol, std::back_inserter(result));
        break;
    }

    int numSamples = result.size()/2;
    Point3D xyz;
    JNodeSequence newnodes;
    for( int i = 0; i < numSamples; i++) {
        xyz[0] = result[2*i];
        xyz[1] = result[2*i+1];
        xyz[2] = 0.0;
        JNodePtr  newnode = JNode::newObject();
        newnode->setXYZCoords(xyz);
        newnodes.push_back(newnode);
    }

    JEdgeSequence newedges;
    for( int i = 0;  i < numSamples; i++) {
        const JNodePtr &v0 = newnodes[i];
        const JNodePtr &v1 = newnodes[(i+1)%numSamples];
        JEdgePtr newedge   = JEdge::newObject(v0,v1);
        newedges.push_back(newedge);
    }

    if( JEdgeGeometry::getOrientation(newedges) < 0) {
        JEdgeTopology::reverse( newedges );
        boost::reverse( newnodes );
    }

    simplified->addObjects( newnodes);
    simplified->addObjects( newedges);

    if( skipParity == 0) return;

    vector<int> minindex( newnodes.size() );
    for( size_t j = 0; j < newnodes.size(); j++) {
        double maxerror = 0.90*std::numeric_limits<double>::max();
        for( size_t i = 0; i < bndnodes.size(); i++) {
            double l = JNodeGeometry::getLength2(newnodes[j], bndnodes[i] );
            if( l < maxerror) {
                minindex[j] = i;
                maxerror    = l;
            }
        }
    }

    if( !boost::algorithm::is_increasing( minindex ) ) {
        cout << "Warning: Simpliied contour nodes are not sorted " << endl;
        return;
    }

}

