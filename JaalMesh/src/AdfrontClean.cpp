#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_singlet_rule(Vertex *singlet)
{
    int err;

    int lid;
    singlet->getAttribute("Layer", lid);
    if ( lid  != 0) return 1;
    if (!QuadCleanUp::isSinglet(singlet)) return 1;

    JFaceSequence vfaces = singlet->getRelations2();
    assert(vfaces.size() == 1);

    Face *face = vfaces[0];
    int pos = face->getPosOf(singlet);

    Vertex *v1 = face->getNodeAt(pos + 1);
    Vertex *v2 = face->getNodeAt(pos + 2);
    Vertex *v3 = face->getNodeAt(pos + 3);

    SwapQuadEdge edge1(mesh, v1, v2);
    err = edge1.apply_singlet_rule(singlet);
    if (!err) return 0;

    SwapQuadEdge edge2(mesh, v3, v2);
    err = edge2.apply_singlet_rule(singlet);
    if (!err) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_triplet_rule(Vertex *vertex)
{
    //****************************************************************************
    // Implementation ideas:
    //
    // This is the case of 5-434 transition, Where 5 degree node is opposite
    // to 3 degree node and with first swap, it becomes 4-335 or 4-534 and with
    // the second swap, it becomes 3-444 and therefore, the 3 degree node moves
    // inside the domain.
    //
    // Having degree 5 is a somewhat strict condition, because degree four is
    // sufficient, but it will create doublet, and as per our design goals, we
    // do not want to introduce any new doublet in the mesh. In future, we may
    // change the conditions, if that is useful.
    //
    //****************************************************************************
    if (vertex->isBoundary()) return 1;

    JNodeSequence vnodes = vertex->getRelations0();
    if (vnodes.size() != 3) return 1;

    int layerid = vertex->getLayerID();
    JFaceSequence vfaces = vertex->getRelations2();

    int pos;
    // Find out which face which contains the higher level node than the given node.
    Face *face = NULL;
    for (int j = 0; j < 3; j++) {
        Face *f = vfaces[j];
        pos = f->getPosOf( vertex );
        Vertex *v = f->getNodeAt(  pos+2 );
        if( v->getLayerID() > layerid ) {
            face = f;
            break;
        }
    }

    if( face == NULL ) return 1;

    pos = face->getPosOf(vertex);
    Vertex *v1 = face->getNodeAt((pos + 1) % 4);
    Vertex *v2 = face->getNodeAt((pos + 2) % 4);
    Vertex *v3 = face->getNodeAt((pos + 3) % 4);

    // If v1 and v3 are at the lower level, they are already processed.
    // v2 must always be higher.
    if (v1->getLayerID() <  layerid   ||
            v2->getLayerID() <= layerid   ||
            v3->getLayerID() <  layerid ) return 1;

    int d1 = v1->getRelations2().size();
    int d2 = v2->getRelations2().size();
    int d3 = v3->getRelations2().size();

    // Avoids creating doublet at v2
    if( d2 < 4 ) return 1;

    // Avoids creating doublet at v1 or v3
    if (max(d1, d3) < 4 ) return 1;

    SwapQuadEdge edge1(mesh, v1, v2, face);
    edge1.modify_fronts(1);
    int err1 = edge1.apply_deficient_rule(vertex);
    if (!err1) return 0;

    SwapQuadEdge edge2(mesh, v3, v2, face);
    edge2.modify_fronts(1);
    int err2 = edge2.apply_deficient_rule(vertex);
    if (!err2 ) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_bridge_rule(Vertex *v0, Vertex *v1)
{
    int layerid = v0->getLayerID();
    if (v1->getLayerID() != layerid) return 1;

    // Although it is checked once again, but it is essential.
    JFaceSequence adjFaces;
    adjFaces = v0->getRelations2();
    if (adjFaces.size() != 3) return 2;

    adjFaces = v1->getRelations2();
    if (adjFaces.size() != 3) return 3;

    // Our assumption is that all the edges are simple.
    Mesh::getRelations112(v0, v1, adjFaces);
    if (adjFaces.size() != 2) return 4;

    if (adjFaces[0]->getLayerID() == adjFaces[1]->getLayerID()) return 5;

    // Check for the internal face.
    Face *internal_face = NULL;
    internal_face = adjFaces[0]->getLayerID() > adjFaces[1]->getLayerID() ? adjFaces[0] : adjFaces[1];

    // Swap the opposite edge:
    Vertex *v2, *v3;
    Face::opposite_nodes(internal_face, v0, v1, v2, v3);

    // Try swapping at the first node.
    int err;
    SwapQuadEdge edge1(mesh, v2, v3, internal_face);
    err = edge1.apply_deficient_rule(v0);
    if (!err) return 0;

    // Try swapping at the second node.
    SwapQuadEdge edge2(mesh, v2, v3, internal_face);
    err = edge2.apply_deficient_rule(v1);
    if (!err) return 0;

    return 1;

}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::apply_advance_front_excess_rule(Vertex *vertex)
{
    if (vertex->isBoundary()) return 1;

    int layerid = vertex->getLayerID();
    JNodeSequence vneighs = vertex->getRelations0();
    int degree = vneighs.size();

    if (degree < 5) return 1;

    int ncount = 0;
    for (int k = 0; k < degree; k++)
        if (vneighs[k]->getLayerID() > layerid) ncount++;

    if (ncount < 2) return 1;

    for (int k = 0; k < degree; k++)
    {
        if (vneighs[k]->getLayerID() > layerid)
        {
            SwapQuadEdge edge(mesh, vertex, vneighs[k]);
            int err = edge.apply_advance_front_rule();
            if (!err) return 0;
        }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
QuadCleanUp::advance_front_edges_swap_once(int layerid)
{
    if (mesh->getAdjTable(0, 0))
        mesh->clearRelations(0, 0);

    if (mesh->getAdjTable(0, 2))
        mesh->clearRelations(0, 2);

    int relexist0 = mesh->buildRelations(0, 0);
    int relexist2 = mesh->buildRelations(0, 2);

    //
    // The input mesh should be doublet free and the output mesh will
    // be doublet free.
    //
    vector<Doublet> doublets = search_interior_doublets();
    assert(doublets.empty());

    // If the boundary is unknown ...
    if (layerid == 0) mesh->search_boundary();

    // We need atleast two layer to work with.
    size_t numnodes = mesh->getSize(0);

    // Check how many iiregular nodes in the present layer...
    size_t num_irregular_nodes_before = 0;
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->getLayerID() == layerid)
        {
            if (vertex->getRelations2().size() != 4)
                num_irregular_nodes_before++;
        }
    }

    size_t ncount = 0;

    //
    // If the layer is boundary, then highest priority must be given to remove
    // the singlets. Many singlets can be removed by swapping the edges. There are
    // some cases, where "Swapping" may not remove singlets. There are two such
    // scenarios, which are handled by calling "remove_boundary_singlet" method.
    //

#ifdef REMOVE_LATER
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        if (vertex->isBoundary())
        {
            err = apply_advance_front_singlet_rule(vertex);
            if (!err) ncount++;
        }
    }

    // Second priority is given to bridges in the mesh. Removal of bridges
    // create a three degree node, which is the next case.

    if (ncount == 0)
    {
        /*
           vector<Edge33> bridges = search_edges33();
           for (size_t i = 0; i < bridges.size(); i++){
                Vertex *v0 = bridges[i].getNodeAt(0);
                Vertex *v1 = bridges[i].getNodeAt(1);
                if( v0->getLayerID() == layerid  && v1->getLayerID() == layerid ) {
                    err = apply_advance_front_bridge_rule(v0, v1);
                    if( !err ) ncount++;
                    v0->setTag(3);
                    v1->setTag(3);
                    mesh->saveAs("dbg.dat");
                    exit(0);
                }
           }
         */
    }
#endif

    // Third Priority is given to degree three nodes. if the adjacent nodes are
    // regular, then two swaps are required. In the first step, the degree of the
    // adjacent node is increased by swapping from outer layer, and then swapping
    // is done to make 3 degree node to regular node.

    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);

            if (vertex->getLayerID() == layerid)
            {
                int err = apply_advance_front_triplet_rule(vertex);
                if (!err) ncount++;
            }
        }
    }

    // Final case, this is the most intuitive and general; swapping is done to
    // reduce the vertex degree.
    if (ncount == 0)
    {
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            if (vertex->getLayerID() == layerid) {
                int err = apply_advance_front_excess_rule(vertex);
                if (!err) ncount++;
            }
        }
    }

#ifdef DEBUG
    if (ncount)
        cout << "Info: Layer " << layerid << " number of layer edges swapped " << ncount << endl;

    if( ncount ) {
        set_no_tags(mesh);
        size_t num_irregular_nodes_after = 0;
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            vertex->setTag(1);
            if (vertex->getLayerID() == layerid)
            {
                if (vertex->getRelations2().size() != 4) {
                    vertex->setTag(2);
                    num_irregular_nodes_after++;
                }
                JNodeSequence vnghs = vertex->getRelations0();
                int minid = MAXINT;
                for( int j = 0; j < vnghs.size(); j++)
                    minid = min( minid, vnghs[j]->getLayerID() );
                vertex->setTag(2);
                mesh->saveAs( "dbg.dat");
                assert( vertex->getLayerID() == minid+1 );
            }
        }

        cout << "Edges Swaaped " << num_irregular_nodes_before
             << " " << num_irregular_nodes_after << endl;
        assert( num_irregular_nodes_before > num_irregular_nodes_after);
        Break();
    }
#endif

    if (!relexist0)
        mesh->clearRelations(0, 0);

    if (!relexist2)
        mesh->clearRelations(0, 2);

    mesh->collect_garbage();

    doublets = search_interior_doublets();
    assert(doublets.empty());

    return ncount;
}

////////////////////////////////////////////////////////////////////////////////

void
QuadCleanUp::advancing_front_edges_swap()
{
    //****************************************************************************
    // Implementation ideas:
    // In general, irregular nodes are scattered around in the domain ( even after
    // initial cleanup operations. With the advance front swapping, we try to
    // (1) clean the mesh starting from the boundary (2) cluster irregular nodes
    // deep inside the mesh ( far from the boundary ) with the hope that clustering
    // will provide more opportunities for cleanup operations ( for example, enable
    // diamonds, bridges etc).
    //
    // There are two side-effects of this procedure.
    // 1)  It may create very high valance nodes in the last layers, which we may
    //     not be able to clean. ( this is not a stopper, as we can call vertex
    //     reduction modules, which will again scatter the irregular nodes.
    // 2)  The position of irregular nodes may depends on the density of mesh.
    //
    //
    //***************************************************************************
    int nfronts = mesh->setWavefront(0);

#ifdef DEBUG
    cout << "Info:  Before Advance front swapping " << endl;
    map<int,size_t> nodedegree;
    map<int,size_t>::const_iterator it;
    for (int ilayer = 1; ilayer < nfronts;  ilayer++) {
        nodedegree.clear();
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            if (vertex->getLayerID() == ilayer) {
                int nd = vertex->getRelations2().size();
                nodedegree[nd]++;
            }
        }
        cout << "In Layer " <<   ilayer << " Degree Distribution" << endl;
        for( it = nodedegree.begin(); it != nodedegree.end(); ++it)
            cout << it->first <<  "  " << it->second << endl;
    }
#endif

    /////////////////////////////////////////////////////////////////////////////

    for (int ilayer = 1; ilayer < nfronts;  ilayer++)
        advance_front_edges_swap_once(ilayer);

    /////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
    nfronts = mesh->setWavefront(0);
    cout << "Info:  After Advance front swapping " << endl;
    for (int ilayer = 1; ilayer < nfronts;  ilayer++) {
        nodedegree.clear();
        for (size_t i = 0; i < numnodes; i++)
        {
            Vertex *vertex = mesh->getNodeAt(i);
            if (vertex->getLayerID() == ilayer) {
                int nd = vertex->getRelations2().size();
                nodedegree[nd]++;
            }
        }
        cout << "In Layer " <<   ilayer << " Degree Distribution" << endl;
        for( it = nodedegree.begin(); it != nodedegree.end(); ++it)
            cout << it->first <<  "  " << it->second << endl;
    }
#endif

}
////////////////////////////////////////////////////////////////////
