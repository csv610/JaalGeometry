#include "SwapEdges.hpp"
#include "basic_math.hpp"

#include <iostream>
#include <iomanip>

using namespace Jaal;

//#############################################################################

bool JSwapTriEdge::isSwappable(const JEdgePtr &edge, bool check567)
{
    if( mesh == nullptr|| edge == nullptr) return 0;

    if( !edge->isActive() ) return 0;

    swapedge = edge;
    nodes[0] = edge->getNodeAt(0);
    nodes[1] = edge->getNodeAt(1);

    JEdge::getRelations(edge, faceNeighs);
    if( faceNeighs.size() != 2 ) return 0;

    if( faceNeighs[0] == faceNeighs[1] ) return 0;

    nodes[2] = JTriangle::getOppositeNode(faceNeighs[0], nodes[0], nodes[1] );
    nodes[3] = JTriangle::getOppositeNode(faceNeighs[1], nodes[0], nodes[1] );
    if( nodes[2] == nullptr || nodes[3] == nullptr ) return 0;

    if( check567 ) {
        if( nodes[0]->getNumRelations(2) <= 6) return 0;
        if( nodes[1]->getNumRelations(2) <= 6) return 0;
        if( nodes[2]->getNumRelations(2) >= 6) return 0;
        if( nodes[3]->getNumRelations(2) >= 6) return 0;
        if( nodes[2]->isBoundary() ) return 0;
        if( nodes[3]->isBoundary() ) return 0;
    }

    if(  isSharp() || isDart() ) return 0;

    return 1;
}

//#############################################################################

int JSwapTriEdge::commit()
{
    JFacePtr oldt1 = faceNeighs[0];
    JFacePtr oldt2 = faceNeighs[1];

    JTrianglePtr newt1, newt2;

    const JNodePtr &v1  = nodes[0];
    const JNodePtr &v2  = nodes[1];
    const JNodePtr &ov1 = nodes[2];
    const JNodePtr &ov2 = nodes[3];

    int pos = oldt1->getPosOf(v1);
    if( oldt1->getNodeAt(pos+1) == ov1) {
        newt1 = JTriangle::newObject( v1, ov1, ov2);
        newt2 = JTriangle::newObject( v2, ov2, ov1);
    }

    if( oldt1->getNodeAt(pos+2) == ov1) {
        newt1 = JTriangle::newObject( v1, ov2, ov1);
        newt2 = JTriangle::newObject( v2, ov1, ov2);
    }

    oldt1->setStatus(JMeshEntity::REMOVE);
    oldt2->setStatus(JMeshEntity::REMOVE);
    swapedge->setStatus(JMeshEntity::REMOVE);

    int err;
    err = mesh->addObject( newt1 );
    assert(!err);
    err = mesh->addObject( newt2 );
    assert(!err);
    return 0;
}
//#############################################################################

int JSwapTriEdge::applyAt(const JEdgePtr &edge)
{
    if( !isSwappable(edge) ) return 1;

    JEdgePtr newedge = nullptr;
    double len0,len1;

    int err = 1;
    switch( fliprule ) {
    case JEdgeSwap::NO_RULE_FLIP:
        err = commit();
        break;
    case JEdgeSwap::DELAUNAY_FLIP:
        len0 = JMath::length2( nodes[0]->getXYZCoords(), nodes[1]->getXYZCoords() );
        len1 = JMath::length2( nodes[2]->getXYZCoords(), nodes[3]->getXYZCoords() );
        if( len1 < len0) err = commit();
        break;
    }

    return err;
}

///////////////////////////////////////////////////////////////////////////////

int JSwapTriEdge::execute()
{
    if( mesh == nullptr) return 1;

    num_edges_flipped = 0;
    size_t numEdges = mesh->getSize(1);

    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &oldedge = mesh->getEdgeAt(i);
        int err = applyAt(oldedge);
        if( !err ) num_edges_flipped++;
    }

    if( num_edges_flipped) mesh->pruneAll();

    return num_edges_flipped;
}

//#############################################################################

bool JSwapTriEdge::isSharp() const
{
    const JNodePtr &v1  = nodes[0];
    const JNodePtr &v2  = nodes[1];
    const JNodePtr &ov1 = nodes[2];
    const JNodePtr &ov2 = nodes[3];

    Vec3D A, B;
    A = JTriGeometry::getNormal(ov1, v1, v2);
    B = JTriGeometry::getNormal(ov2, v2, v1);

    double angle = JMath::getVecAngle(A, B, ANGLE_IN_DEGREES);

    if( angle > creaseAngle)  {
//        cout << "Sharp edge detected : flipping not done for " << angle << " Crease Angle :" << creaseAngle << endl;
        return 1;
    }

    return 0;
}

//#############################################################################

bool JSwapTriEdge::isDart() const
{
    const Point3D &p1 = nodes[0]->getXYZCoords();
    const Point3D &p2 = nodes[1]->getXYZCoords();
    const Point3D &p3 = nodes[2]->getXYZCoords();
    const Point3D &p4 = nodes[3]->getXYZCoords();

    Point3D  alpha, beta;

    JMath::getTriAngles(p1, p2, p3, alpha);
    JMath::getTriAngles(p1, p2, p4, beta );

    if( alpha[0] + beta[0]  > 179 ) return 1;
    if( alpha[1] + beta[1]  > 179 ) return 1;

    return 0;
}

//#############################################################################

#ifdef CSV
bool JSwapTriEdge::is_edge_flip_allowed(const FlipEdge &edge, int rule) const
{
    if (!edge.isValid())   return 0;
    if (edge.isConcave()) return 0;

    JNode  *v1  = edge.getNodeAt(0);
    JNode  *v2  = edge.getNodeAt(1);
    JNode *ov1 = edge.opposite_nodes[0];
    JNode *ov2 = edge.opposite_nodes[1];

    if (v1->isBoundary()  && v2->isBoundary()  ) return 0;
    if (ov1->isBoundary() || ov2->isBoundary() ) return 0;

    int d1 = v1->getNumRelations(2);
    int d2 = v2->getNumRelations(2);
    int d3 = ov1->getNumRelations(2);
    int d4 = ov2->getNumRelations(2);

    if( rule == DELAUNAY_RULE) {
        double len1 = JNode::length2(v1, v2);
        double len2 = JNode::length2(ov1, ov2);
        if (len1 <= len2) return 0;

        if (edge.isSharp(creaseAngle)) return 0;
        return 1;
    }

    if( rule == ADVANCE_FRONT_RULE ) {
        int ideal_v1 =  v1->get_ideal_face_degree(3);
        int ideal_v2 =  v2->get_ideal_face_degree(3);
        int ideal_v3 =  ov1->get_ideal_face_degree(3);
        int ideal_v4 =  ov2->get_ideal_face_degree(3);

        int l1 = 0;
        v1->getAttribute("Layer",  l1);

        int l2 = 0;
        v2->getAttribute("Layer",  l2);

        int l3 = 0;
        ov1->getAttribute("Layer", l3);

        int l4 = 0;
        ov2->getAttribute("Layer", l4);

        // Decrease the vertex degree
        if( (d1  > ideal_v1) && (l2 > l1) && (l3 > l1) && (l4 > l1) ) return 1;
        if( (d2  > ideal_v2) && (l1 > l2) && (l3 > l2) && (l4 > l2) ) return 1;

        // Increase the vertex degree ...
        if( (d3  < ideal_v3) && (l1 > l3) && (l2 > l3) && (l4 > l3) ) return 1;
        if( (d4  < ideal_v4) && (l1 > l4) && (l2 > l4) && (l3 > l4) ) return 1;

        return 0;
    }

    if( rule == DEGREE_REDUCTION_RULE) {
        int ideal_v1 =  v1->get_ideal_face_degree(3);
        int ideal_v2 =  v2->get_ideal_face_degree(3);

        if( v1->isBoundary() &&  (d1 > ideal_v1) && (d2 >3) )  return 1;
        if( v2->isBoundary() &&  (d2 > ideal_v2) && (d1 >3) )  return 1;

        int relaxation_index = d1 + d2 - d3 - d4;
        if (relaxation_index < 3) return 0;
        return 1;
    }

    return 0;
}

//#############################################################################

int JSwapTriEdge ::atomicOp(const JNodePtr &apexVertex, int rule)
{
    JFaceSequence vneighs;
    apexVertex->getRelations( vneighs );
    int numNeighs = vneighs.size();

    for( int i = 0; i < numNeighs; i++) {
        if( unchecked(vneighs[i] ) ) {
            int err = atomicOp( vneighs[i], rule);
            if( err  == 0) return 0;
        }
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////
int JSwapTriEdge :: apply_advance_front_rule()
{
    int relexist2 = mesh->buildRelations(0, 2);
    int relexist0 = mesh->buildRelations(0, 0);

    mesh->getTopology()->searchBoundary();

    size_t numNodes = mesh->getSize(0);
    JNodeSequence currlayer, vneighs;
    JNodeSet   nextlayer;

    for(size_t i = 0; i < numNodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        if( v->isBoundary() ) {
            v->setAttribute("Layer", 0);
            currlayer.push_back(v);
        } else
            v->setAttribute("Layer",INT_MAX);
    }

    Jaal::MeshOptimization mopt;
    size_t ncount, nSize, num_edges_flipped = 0;
    mesh->getTopology()->getConsistent();
    mopt.shape_optimize(mesh);

    LaplaceNoWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(100);

    int curr_layer_id = 0;
    int nirregular0, nirregular1;

    while(1) {
        nSize = currlayer.size();
        nirregular0 = 0;
        for( size_t i = 0; i < nSize; i++) {
            Vertex *v = currlayer[i];
            if( !v->isRemoved() && !isIdeal(v) ) nirregular0++;
        }

        for( int k = 0; k < 3; k++) {  // Three attempts for single layer
            while(1) {
                nSize = currlayer.size();
                ncount = 0;
                for( size_t i = 0; i < nSize; i++) {
                    Vertex *v = currlayer[i];
                    if( !v->isRemoved() &&  !isIdeal(v)) {
                        int err  = atomicOp( v, ADVANCE_FRONT_RULE );
                        if( !err ) ncount++;
                    }
                }
                if( ncount == 0) break;
                num_edges_flipped += ncount;
            }

            nSize = currlayer.size();
            nirregular1 = 0;
            for( size_t i = 0; i < nSize; i++) {
                Vertex *v = currlayer[i];
                if( !v->isRemoved() &&  !isIdeal(v) ) nirregular1++;
            }
            assert( nirregular1 <= nirregular0);
            if( nirregular1 == 0) break;
        }
        lapsmooth.execute();

//        mopt.shape_optimize(mesh);

        cout << "Layer : " << curr_layer_id << endl;
        cout << "# of Irregular nodes before swapping : " << nirregular0 << endl;
        cout << "# of Irregular nodes after swapping  : " << nirregular1 << endl;

        int lid = 0;

        nextlayer.clear();
        nSize = currlayer.size();
        for( size_t i = 0; i < nSize; i++) {
            JNodePtr v = currlayer[i];
            v->getRelations( vneighs );
            for( size_t k = 0; k < vneighs.size(); k++) {
                vneighs[k]->getAttribute( "Layer", lid);
                if( lid > curr_layer_id ) {
                    vneighs[k]->setAttribute("Layer", curr_layer_id+1 );
                    nextlayer.insert( vneighs[k] );
                }
            }
        }
        if( nextlayer.empty() ) break;

        JNodeSet::const_iterator it;
        currlayer.resize(nextlayer.size() );
        int index = 0;
        for( it = nextlayer.begin(); it != nextlayer.end(); ++it)
            currlayer[index++] = *it;
        curr_layer_id++;
    }
    cout << "# of Edges Swapped " << num_edges_flipped << endl;

    vector<int>  less_than_ideal, more_than_ideal, total_ideal;

    int numLayers = curr_layer_id;
    less_than_ideal.resize( numLayers );
    more_than_ideal.resize( numLayers );
    total_ideal.resize( numLayers );

    for( int i = 0; i < numLayers; i++) {
        less_than_ideal[i] = 0;
        more_than_ideal[i] = 0;
        total_ideal[i] = 0;
    }

    int lid = 0;

    numNodes = mesh->getSize(0);
    int final_irregular = 0;
    for( size_t i = 0; i < numNodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        if( !v->isRemoved()) {
            v->getAttribute("Layer", lid );
            int curr_degree  = v->getNumRelations(2);
            int ideal_degree = v->get_ideal_face_degree(3);
            if( curr_degree != ideal_degree ) {
                final_irregular++;
                if( curr_degree < ideal_degree) less_than_ideal[lid]++;
                if( curr_degree > ideal_degree) more_than_ideal[lid]++;
            } else
                total_ideal[lid]++;
        }
    }

    cout << " Layer   Less   More  Ideal " << endl;
    for( int i = 0; i < numLayers; i++)
        cout << i << setw(10) <<  less_than_ideal[i]
             << setw(10) <<  more_than_ideal[i]
             << setw(10) <<  total_ideal[i] << endl;
    cout << " Final # of irregular nodes : " << final_irregular << endl;

    mopt.shape_optimize(mesh);

    if (!relexist2) mesh->clearRelations(0, 2);
    if (!relexist0) mesh->clearRelations(0, 0);

    return num_edges_flipped;
}
#endif

