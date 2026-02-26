#include "SwapEdges.hpp"
#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

bool JSwapQuadEdge :: isSwappable( const JEdgePtr &e)
{
    swapedge = e;

    if( swapedge == nullptr) return 0;

    if (buildBoundary() != 0) return 0;

    int d1 = swapedge->getNodeAt(0)->getNumRelations(2);
    int d2 = swapedge->getNodeAt(1)->getNumRelations(2);

    if( d1 < 1 || d2 < 1 ) {
        cout <<"Warning: Edge swapping not performed because vertex-face relations absent" << endl;
        return 0;
    }

    int start_pos = getPosOf( swapedge->getNodeAt(0));
    assert(start_pos >= 0);

    for (int i = 0; i < 2; i++) {
        int pos = (start_pos + i + 1) % 6;
        JNodePtr v0 = boundNodes[ pos ];
        if (v0->isBoundary()) continue;
        int d3 = v0->getNumRelations(2);

        JNodePtr v1 = boundNodes[(pos + 3) % 6];
        if (v1->isBoundary()) continue;
        int d4 = v1->getNumRelations(2);

        if (JSwapQuadEdge::is_topologically_valid_swap(d1, d2, d3, d4)) return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool
JSwapQuadEdge::is_topologically_valid_swap(int d1, int d2, int d3, int d4)
{
    if (d1 < 4 || d2 < 4) return 0;
    if ((d1 > 4) && (d2 > 4) && (d3 == 3) && (d4 == 3)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 3) && (d4 == 4)) return 1;
    if ((d1 == 5) && (d2 == 5) && (d3 == 4) && (d4 == 3)) return 1;
    if (max(d1, d2) > max(d3 + 1, d4 + 1)) return 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::makeDiagonalAt(int pos, bool bound_check)
{
    //
    // Given closednodes ( exactly six nodes) in an order ( clockwise or
    // counter-clockwise) and create two quads having common diagonal between
    // (pos, pos+3).
    //
    assert(pos >= 0 && pos < 6);

    if (bound_check) {
        if (boundNodes[(pos + 0) % 6]->isBoundary()) return 1;
        if (boundNodes[(pos + 3) % 6]->isBoundary()) return 1;
    }
    newfaces.resize(2);

    newfaces[0] = JQuadrilateral::newObject();
    newfaces[1] = JQuadrilateral::newObject();

    // Change the connectivity of the two quads.
    JNodeSequence qConnect(4);
    qConnect[0] = boundNodes[(pos + 0) % 6];
    qConnect[1] = boundNodes[(pos + 1) % 6];
    qConnect[2] = boundNodes[(pos + 2) % 6];
    qConnect[3] = boundNodes[(pos + 3) % 6];
    newfaces[0]->setNodes(qConnect);

    /*
        if( !JFaceGeometry::isSimple( newfaces[0]) ) {
            cout << "Concave face " << endl;
            return GEOMETRIC_ERROR;
        }
    */

    qConnect[0] = boundNodes[(pos + 3) % 6];
    qConnect[1] = boundNodes[(pos + 4) % 6];
    qConnect[2] = boundNodes[(pos + 5) % 6];
    qConnect[3] = boundNodes[(pos + 6) % 6];
    newfaces[1]->setNodes(qConnect);
    /*
        if( !JFaceGeometry::isSimple( newfaces[1]) ) {
            cout << "Concave face " << endl;
            return GEOMETRIC_ERROR;
        }
    */

    // Old edge goes away
    swapedge->setStatus(JMeshEntity::REMOVE);

    // A new edge is formed...
    JNodePtr v0 = boundNodes[(pos + 0) % 6];
    JNodePtr v1 = boundNodes[(pos + 3) % 6];
    JEdgePtr newedge = JEdge::newObject(v0,v1);
    mesh->addObject(newedge);

    // Old faces about the edges are removed.
    edgefaces[0]->setStatus( JMeshEntity::REMOVE);
    edgefaces[1]->setStatus( JMeshEntity::REMOVE);

    // New faces are added ..
    mesh->addObject( newfaces[0] );
    mesh->addObject( newfaces[1] );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
JSwapQuadEdge::getPosOf(const JNodePtr &vertex) const
{
    for (int i = 0; i < 6; i++)
        if (boundNodes[i] == vertex) return i;
    return -1;
}

/////////////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::buildBoundary()
{
    if( swapedge == nullptr ) return 1;
    assert( mesh != nullptr );

    assert( mesh->getAdjTable(0,2) );

    // Since the degree of each node of an existing edge decreases by
    // one, by restricting the swapping to degree greater than 3 will
    // ensure that no doublet is created.
    //
    JNodePtr v0 = swapedge->getNodeAt(0);
    JNodePtr v1 = swapedge->getNodeAt(1);

    int nsize1 = v0->getNumRelations(2);
    assert( nsize1 > 0);

    if (v0->isBoundary()) {
        if (nsize1 < 3) return 1;
    } else {
        if (nsize1 < 4) return 1;
    }

    int nsize2 = v1->getNumRelations(2);
    assert( nsize2 > 0);

    if (v1->isBoundary()) {
        if (nsize2 < 3) return 1;
    } else {
        if (nsize2 < 4) return 1;
    }

    boundNodes.resize(6);
    boundNodes[0] = swapedge->getNodeAt(0);
    boundNodes[3] = swapedge->getNodeAt(1);

    JFaceSequence nghs;
    JEdge::getRelations(swapedge, nghs);

    assert( !nghs.empty() ) ;

    if (nghs.size() != 2) return 1;

    edgefaces.resize(2);
    if (firstFace == nullptr) {
        edgefaces[0] = nghs[0];
        edgefaces[1] = nghs[1];
    } else {
        if (nghs[0] == firstFace) {
            edgefaces[0] = nghs[0];
            edgefaces[1] = nghs[1];
        }

        if (nghs[1] == firstFace) {
            edgefaces[0] = nghs[1];
            edgefaces[1] = nghs[0];
        }
    }

    JFacePtr face;
    int pos;
    // No ambiguity for the second node in face1
    face = edgefaces[0];
    pos = face->getPosOf( boundNodes[0] );
    boundNodes[2] = face->getNodeAt( pos + 2 );

    if( face->getNodeAt(pos+1) ==  boundNodes[3] )
        boundNodes[1] = face->getNodeAt(pos+3);
    else
        boundNodes[1] = face->getNodeAt(pos+1);


    // No ambiguity for the second node in face 2
    face = edgefaces[1];
    pos = face->getPosOf( boundNodes[0] );
    boundNodes[4] = face->getNodeAt( pos + 2 );

    if( face->getNodeAt(pos+1) ==  boundNodes[3] )
        boundNodes[5] = face->getNodeAt(pos+3);
    else
        boundNodes[5] = face->getNodeAt(pos+1);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JSwapQuadEdge :: applyAt( const JEdgePtr &edge, JFacePtr face )
{
    swapedge = edge;
    if( swapedge->isBoundary() ) return 1;

    firstFace = face;

    int err = applyDegreeRule();
    return err;
}

////////////////////////////////////////////////////////////////////////////////

int JSwapQuadEdge :: applyAt( const JEdgePtr &edge, const JNodePtr &vtx)
{
    swapedge = edge;
    if( swapedge->isBoundary() ) return 1;

    if (buildBoundary() != 0) return 1;

    int pos = -1;
    for( int i = 0; i < 6; i++)
        if( boundNodes[i] == vtx ) pos = i;

    if( pos >= 0) {
        int err = makeDiagonalAt(pos, 0);
        if (!err) return 0;
    }

    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int JSwapQuadEdge :: applyAt( const JNodePtr &vtx)
{
    if( mesh == nullptr) return 1;

    JEdgeSequence edges;
    mesh->getTopology()->getDerivedRelations(vtx, edges);
    int numedges = edges.size();

    int nCount = 0;
    for( int i = 0; i < numedges; i++) {
        int err = applyAt( edges[i] );
        if( !err) nCount;
    }

    if( nCount > 0) return 0;

    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::applyDegreeRule()
{
    if( swapedge == nullptr ) return 1;

    if (buildBoundary() != 0) return 1;

    int d1 = swapedge->getNodeAt(0)->getNumRelations(2);
    int d2 = swapedge->getNodeAt(1)->getNumRelations(2);

    int start_pos = getPosOf( swapedge->getNodeAt(0));
    assert(start_pos >= 0);

    for (int i = 0; i < 2; i++) {
        int pos = (start_pos + i + 1) % 6;
        JNodePtr v0 = boundNodes[ pos ];
        if (v0->isBoundary()) continue;
        int d3 = v0->getNumRelations(2);

        JNodePtr v1 = boundNodes[(pos + 3) % 6];
        if (v1->isBoundary()) continue;
        int d4 = v1->getNumRelations(2);

        if (JSwapQuadEdge::is_topologically_valid_swap(d1, d2, d3, d4)) {
            int err = makeDiagonalAt(pos);
            if (!err) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::applyDeficientRule(const JNodePtr &vertex)
{

    if (vertex->isBoundary()) return 1;

#ifdef NEEDS_RECONSIDER
    assert(firstFace);

    int d1 = connect[0]->getNumRelations( 2 );
    int d2 = connect[1]->getNumRelations( 2 );

    // Note 1: Just near the deficient node, there must be atleast 5 nodes in
    // order to swap the edge,otherwise, things may not converge, on the
    // same level, the deficient node will keep changing the location.
    //
    // Note 2: The node connect[1] is opposite to the deficient node and
    // most probably on the next level. If we keep moving the deficient
    // node inside the domain, then it is fine to have its degree to 4, so
    // that the deficient node will be generated there.
    if( d2 < 4 ) return 1;

    JFaceSequence faces;
    JFacePtr f0, f1;
    if (d1 == 4) {
        vertex->getRelations( faces );
        if (faces.size() != 3) return 1;
        int layerid = connect[0]->getLayerID();
        f0 = firstFace;
        Mesh::getRelations112(connect[0], connect[1], faces);
        assert(faces.size() == 2);
        if (faces[0] == f0) f1 = faces[1];
        if (faces[1] == f0) f1 = faces[0];
        assert(f1);

        int pos = f1->getPosOf(connect[0]);
        JNodePtr vopp = f1->getNodeAt(pos + 2);
        if (connect[1]->getLayerID() > layerid && vopp->getLayerID() > layerid) {
            set_no_tags(mesh);
            connect[1]->setTag(1);
            vopp->setTag(1);
            mesh->saveAs("dbg.dat");
            Break();
            JSwapQuadEdge edge(mesh, connect[1], vopp, f1);
            int err = edge.apply_deficient_rule(connect[0]);
            if (err) return 1;
        }
    }
    exit(0);

    d1 = connect[0]->getNumRelations(2);
    d2 = connect[1]->getNumRelations(2);

    // Don't create doublet at two nodes...
    if( d1 == 3 || d2 == 3 ) return 1;

    // Having these conditions set, now build the structure.
    if (build_boundary() != 0) return 3;

    update_front();

    //Find the position of triplet in the contour ...
    int pos = this->getPosOf(vertex);
    assert(pos == 1 || pos == 5);

    if (check_fronts) {
        int layerid = vertex->getLayerID();
        if (boundNodes[ (pos + 3) % 6]->getLayerID() <= layerid) return 1;

        if( connect[0]->getLayerID() <= layerid && d1 < 5 ) return 1;
        if( connect[1]->getLayerID() <= layerid && d2 < 5 ) return 1;
    }

    // Create new quads whose common diagonal must contain the singlet.
    int err = make_new_diagonal_at(pos);

    /*
        if( err  ) {
            for (int i = 0; i < 6; i++) {
                 JNodeSequence vneighs = boundNodes[i]->getRelations0();
                 int minid = boundNodes[i]->getLayerID();
                 for( int j = 0; j < vneighs.size(); j++)
                      minid = min( minid, vneighs[j]->getLayerID() );
                assert( boundNodes[i]->getLayerID() == minid+1);
             }
        }
    */

    return err;
#endif
    return 1;

}

/////////////////////////////////////////////////////////////////////////////

#ifdef OBSOLETE
int
QuadCleanUp::swap_concave_faces()
{
    /*
         int relexist0 = mesh->buildRelations(0, 0);
         int relexist2 = mesh->buildRelations(0, 2);

         mesh->getTopology()->search_boundary();

         size_t numfaces = mesh->getSize(2);
         for (size_t i = 0; i < numfaces; i++) {
              Face *face = mesh->getFaceAt(i);
              face->setVisitBit(0);
         }

         size_t numnodes = mesh->getSize(0);
         for (size_t i = 0; i < numnodes; i++) {
              Vertex *vertex = mesh->getNodeAt(i);
              vertex->setVisitBit(0);
         }

         size_t ncount = 0;
         for (size_t iface = 0; iface < numfaces; iface++) {
              Face *face = mesh->getFaceAt(iface);
              int pos = face->concaveAt();
              if (pos >= 0) {
                   for (int i = 0; i < 2; i++) {
                        Vertex *v0 = face->getNodeAt(pos + 1 + i);
                        Vertex *v1 = face->getNodeAt(pos + 2 + i);
                        JSwapQuadEdge edge(mesh, v0, v1);
                        int err = edge.apply_concave_rule();
                        if (!err) {
                             ncount++;
                             break;
                        }
                   }
              }
         }

         if (ncount) {
              mesh->prune();
              mesh->enumerate(0);
              mesh->enumerate(2);
              cout << "Info: # of swapped edges " << ncount << endl;
         }

         if (!relexist0) mesh->clearRelations(0, 0);
         if (!relexist2) mesh->clearRelations(0, 2);

         return ncount;
    */
}
#endif

///////////////////////////////////////////////////////////////////////////////
//             OBSOLETE Functions...
///////////////////////////////////////////////////////////////////////////////


/*
int
JSwapQuadEdge::apply_advance_front_rule()
{
    if (build_boundary() != 0) return 1;

         int layerid = connect[0]->getLayerID();
         for (int i = 0; i < 3; i++) {
              Vertex *v0 = boundNodes[(i + 0) % 6];
              Vertex *v3 = boundNodes[(i + 3) % 6];
              if ((v0->getLayerID() > layerid) && (v3->getLayerID() > layerid)) {
                   int err = make_new_diagonal_at(i);
                   if (!err) return 0;
              }
         }
    return 3;
}

/////////////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::apply_singlet_rule(Vertex *singlet)
{
          // The first node must be boundary and the other node must be internal..
          assert(connect[0]->isBoundary());

          if (connect[1]->isBoundary()) return 1;

          if (build_boundary() != 0) return 1;
          ////////////////////////////////////////////////////////////////////////////
          // Objective :: Swap quads common diagonal such that the resulting diagonal
          //              contains the Singlet node.
          ////////////////////////////////////////////////////////////////////////////
          assert(QuadCleanUp::isSinglet(singlet));

          if (adjFaces[0] == nullptr || adjFaces[1] == nullptr) return 2;

          assert(QuadCleanUp::hasSinglet(adjFaces[0]));

          // Make sure that other face doesn't have a singlet
          if (QuadCleanUp::hasSinglet(adjFaces[1])) return 3;

          //Find the position of singlet in the contour ...
          int pos = this->getPosOf(singlet);
          assert(pos >= 0);

          // Create new quads whose common diagonal must contain the singlet.
          return make_new_diagonal_at(pos, 0);
}
 */

/////////////////////////////////////////////////////////////////////////////


int
JSwapQuadEdge::applyConcaveRule()
{
#ifdef CSV
    if (build_boundary() != 0) return 1;

    int start_pos = getPosOf(connect[0]);
    for (int i = 0; i < 2; i++) {
        int err = make_new_diagonal_at((start_pos + 1 + i) % 6);
        if (!err) return 0;
    }
#endif
    return 1;
}

/////////////////////////////////////////////////////////////////////

int
JSwapQuadEdge::applyBoundRule()
{
    /*
         if (!connect[0]->isBoundary()) return 1;
         if (connect[1]->isBoundary()) return 1;

         if (build_boundary() != 0) return 2;

         for (int i = 0; i < 3; i++) {
              Vertex *v0 = boundNodes[(i + 0) % 6];
              Vertex *v3 = boundNodes[(i + 3) % 6];
              if ((!v0->isBoundary()) && (!v3->isBoundary())) {
                   int err = make_new_diagonal_at(i, 0);
                   if (!err) return 0;
              }
         }
    */

    return 3;
}

