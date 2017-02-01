#include "QuadCleanUp.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_bridge_tag(Mesh *mesh)
{
    int relexist = mesh->buildRelations(0, 2);

    QuadCleanUp qClean(mesh);
    vector<Edge33> bridges = qClean.search_edges33();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *v = mesh->getNodeAt(i);
        v->setTag(1);
    }

    for (size_t i = 0; i < bridges.size(); i++)
    {
        Vertex *v0 = bridges[i].connect[0];
        Vertex *v1 = bridges[i].connect[1];
        if ((v0 != NULL) && (v1 != NULL))
        {
            v0->setTag(0);
            v1->setTag(0);
        }
    }

    if (!relexist)
        mesh->clearRelations(0, 2);
}

//////////////////////////////////////////////////////////////////////

vector<Edge33>
QuadCleanUp::search_edges33_in_layer(int layerid)
{
    int relexist0 = mesh->buildRelations(0, 0);

    JEdgeSequence medges = mesh->getEdges();

    mesh->search_boundary();

    vEdges33.clear();
    for (size_t iedge = 0; iedge < medges.size(); iedge++)
    {
        Edge *edge = medges[iedge];
        Vertex *v0 = edge->getNodeAt(0);
        Vertex *v1 = edge->getNodeAt(1);
        int l0 = v0->getLayerID();
        int l1 = v1->getLayerID();
        if (l0 < layerid - 2 || l0 > layerid + 2) continue;
        if (l1 < layerid - 2 || l1 > layerid + 2) continue;
        if (isEdge33(edge))
        {
            Edge33 bridge(mesh, v0, v1);
            vEdges33.push_back(bridge);
        }
    }

    // We don't need edge pointers anymore
    for (size_t iedge = 0; iedge < medges.size(); iedge++)
        delete medges[iedge];

    if (!relexist0)
        mesh->clearRelations(0, 0);

    return vEdges33;
}

/////////////////////////////////////////////////////////////////////////////

vector<Edge33>
QuadCleanUp::search_edges33()
{
    vEdges33.clear();

    int relexist0 = mesh->buildRelations(0, 0);

    JEdgeSequence medges = mesh->getEdges();

    mesh->search_boundary();

    for (size_t iedge = 0; iedge < medges.size(); iedge++)
    {
        Edge *edge = medges[iedge];
        if (isEdge33(edge))
        {
            Vertex *v0 = edge->getNodeAt(0);
            Vertex *v1 = edge->getNodeAt(1);
            Edge33 bridge(mesh, v0, v1);
            vEdges33.push_back(bridge);
        }
        delete edge; // No more needed..
    }

    if (!relexist0)
        mesh->clearRelations(0, 0);

    return vEdges33;
}

/////////////////////////////////////////////////////////////////////////////

int
Edge33::build_boundary()
{
    assert(!connect[0]->isBoundary());
    assert(!connect[1]->isBoundary());

    adjFaces = Mesh::getRelations102(connect[0], connect[1]);
    assert(adjFaces.size() == 4);

    // Create a closed chain of bounding edges ...
    vector<Edge> boundedges;
    for (size_t i = 0; i < adjFaces.size(); i++)
    {
        Face *face = adjFaces[i];
        for (int j = 0; j < 4; j++)
        {
            Vertex *ev0 = face->getNodeAt((j + 0) % 4);
            Vertex *ev1 = face->getNodeAt((j + 1) % 4);
            if (ev0 == connect[0] || ev0 == connect[1]) continue;
            if (ev1 == connect[0] || ev1 == connect[1]) continue;
            Edge edge(ev0, ev1);
            boundedges.push_back(edge);
        }
    }

    assert(boundedges.size() == 6);
    Mesh::make_chain(boundedges);

    // Create a closed chain of bounding nodess ...
    bound_nodes = Mesh::chain_nodes(boundedges);
    assert(bound_nodes.size() == 6);

    // It is important that the first node of is attached to one of
    // the two nodes of the bridge.
    assert(mesh->getAdjTable(0, 0));
    JNodeSequence vneighs = connect[0]->getRelations0();
    assert(vneighs.size() == 3);
    JNodeSequence rotated(6);
    for (int i = 0; i < 3; i++)
    {
        if (vneighs[i] != connect[1])
        {
            int pos = -1;
            Vertex *start_vertex = vneighs[i];
            for (int j = 0; j < 6; j++)
            {
                if (bound_nodes[j] == start_vertex)
                {
                    pos = j;
                    break;
                }
            }
            assert(pos >= 0);
            for (int j = 0; j < 6; j++)
                rotated[j] = bound_nodes[ (pos + j) % 6 ];
            bound_nodes = rotated;
        }
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

int
Edge33::build()
{
    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;

    if (v0->getRelations2().size() != 3) return 1;
    if (v1->getRelations2().size() != 3) return 1;

    assert(mesh->getAdjTable(0, 2));

    // Our edges are assumed to be simple i.e. shared by at the most
    // two faces.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 1;

    // Check if the bridge touches boundary.
    boundary = 0;
    if (adjFaces[0]->has_boundary_edge())
        boundary = 1;

    if (adjFaces[1]->has_boundary_edge())
        boundary = 1;

    int err;
    err = build_boundary();
    if (err) return 2;

    // If boundary nothing is done here.
    if (boundary) return 0;

    // Build the outer shell of the bridge. Must contain closed
    // contour with six line segments. There are many ways to
    // convert a hexagon into quads ( mostly 2, 3 or 4 quads).
    // For general case, it is a hard problem.
    //

    int ncount[6];
    for (int i = 0; i < 6; i++)
    {
        JFaceSequence vfaces = bound_nodes[i]->getRelations2();
        if (vfaces.size() == 3) return 3;
        ncount[i] = vfaces.size();
    }

    if (ncount[0] == 4 && ncount[3] == 4)
    {
        err = Face::hexagon_2_quads(bound_nodes, newFaces, 0);
        if (err) return 4;
    }
    else if (ncount[1] == 4 && ncount[4] == 4)
    {
        err = Face::hexagon_2_quads(bound_nodes, newFaces, 1);
        if (err) return 4;
    }

    //
    // Check for "area invariance". Before and after the decomposition
    // area must be same. When there is concave faces, such violations
    // may occur. This method just avoid overlapping cases and produce
    // simple polygons.
    //

    double area0 = 0.0;
    for (size_t i = 0; i < adjFaces.size(); i++)
        area0 += adjFaces[i]->getArea();

    // We have decomposed a hexagon into 2 quads...
    double area1 = 0.0;
    for (size_t i = 0; i < 2; i++)
        area1 += newFaces[i]->getArea();

    if (fabs(area1 - area0) < 1.0E-06) return 0;

    return 1;
}

////////////////////////////////////////////////////////////////////

int
Edge33::remove_internal_one()
{
    if (newFaces.empty()) return 1;

    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;
    if (v0->isRemoved() || v1->isRemoved()) return 1;

    // Check for double removal..
    assert(bound_nodes.size() == 6);
    for (int i = 0; i < 6; i++)
        if (bound_nodes[i]->isVisited()) return 1;

    assert(adjFaces.size() == 4);

    // In this case, no new vertex is created, only for faces are
    // deleted and two new inserted.
    assert(mesh);
    mesh->addFace(newFaces[0]);
    mesh->addFace(newFaces[1]);

    // The bridge edge goes away along with the neighboring faces..
    for (int i = 0; i < 4; i++)
        mesh->remove(adjFaces[i]);

    mesh->remove(connect[0]);
    mesh->remove(connect[1]);

    Point3D backupCoords[6];
    for (int i = 0; i < 6; i++)
        backupCoords[i] = bound_nodes[i]->getXYZCoords();

    LaplaceLengthWeight lw;
    LaplaceSmoothing lapsmooth(mesh);
    lapsmooth.setWeight(&lw);
    lapsmooth.setNumIterations(10);
    lapsmooth.localized_at(bound_nodes);

    set<Face*> faces_to_check;
    for (int i = 0; i < 6; i++)
    {
        JFaceSequence vfaces = bound_nodes[i]->getRelations2();
        for (size_t j = 0; j < vfaces.size(); j++)
            faces_to_check.insert(vfaces[j]);
    }

    bool pass = 1;
    set<Face*>::const_iterator siter;
    for (siter = faces_to_check.begin(); siter != faces_to_check.end(); ++siter)
    {
        Face *f = *siter;
        if (f->invertedAt() >= 0)
        {
            pass = 0;
            break;
        }
    }

    if (!pass)
    {
        mesh->remove(newFaces[0]);
        mesh->remove(newFaces[1]);
        for (int i = 0; i < 4; i++)
            mesh->addFace(adjFaces[i]);
        mesh->addNode(connect[0]);
        mesh->addNode(connect[1]);
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setXYZCoords(backupCoords[i]);
    }
    else
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitBit(1);
        v0->setVisitBit(1);
        v1->setVisitBit(1);
    }

    // Destructor must not deallocate newFaces since they are not kept by
    // the mesh object.
    newFaces[0] = NULL;
    newFaces[1] = NULL;

    return 0;
}

///////////////////////////////////////////////////////////////////////////

int
Edge33::remove_boundary_one()
{
    Vertex *v0 = connect[0];
    Vertex *v1 = connect[1];

    if (v0->isVisited() || v1->isVisited()) return 1;
    if (v0->isRemoved() || v1->isRemoved()) return 1;

    // Our assumption is that all the edges are simple.
    adjFaces = Mesh::getRelations112(v0, v1);
    if (adjFaces.size() != 2) return 1;

    // Check for the internal face.
    Face *internal_face = NULL;

    int ncount_boundfaces = 0;
    if (internal_face == NULL)
    {
        if (!adjFaces[0]->has_boundary_edge())
            internal_face = adjFaces[0];
        else
            ncount_boundfaces++;
    }
    if (internal_face == NULL)
    {
        if (!adjFaces[1]->has_boundary_edge())
            internal_face = adjFaces[1];
        else
            ncount_boundfaces++;
    }

    // A valid boundary bridge must have one boundary face and one internal
    // face.
    if (internal_face == NULL) return 1;
    if (ncount_boundfaces == 2) return 2;

    // Swap the opposite edge:
    Vertex *v2, *v3;
    Face::opposite_nodes(internal_face, v0, v1, v2, v3);

    // Try swapping at the first node.
    int err;
    SwapQuadEdge edge1(mesh, v2, v3, internal_face);
    err = edge1.apply_deficient_rule(connect[0]);
    if (!err)
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitBit(1);
        v0->setVisitBit(1);
        v1->setVisitBit(1);
        return 0;
    }

    // Try swapping at the second node.
    SwapQuadEdge edge2(mesh, v2, v3, internal_face);
    err = edge2.apply_deficient_rule(connect[1]);
    if (!err)
    {
        for (int i = 0; i < 6; i++)
            bound_nodes[i]->setVisitBit(1);
        v0->setVisitBit(1);
        v1->setVisitBit(1);
        return 0;
    }

    return 3;
}
///////////////////////////////////////////////////////////////////////////////

int
Edge33::commit()
{
    int err = 1;

    if (boundary)
        err = remove_boundary_one();
    else
        err = remove_internal_one();

    return err;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges_in_layer(int layerid)
{
    cout << "CSV : " << endl;
    abort();
    /*
        int rel2exist = mesh->buildRelations(0, 2);
        mesh->search_boundary();

        if (vEdges33.empty())
            search_edges33_in_layer(layerid);

        return remove_bridges_once();
    */
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges_once()
{
    int rel0exist = mesh->buildRelations(0, 0);
    int rel2exist = mesh->buildRelations(0, 2);

    mesh->search_boundary();

    int ncount = 0;

    if (vEdges33.empty())
        search_edges33();

    size_t numnodes = mesh->getSize(0);
    for (size_t i = 0; i < numnodes; i++)
    {
        Vertex *vertex = mesh->getNodeAt(i);
        vertex->setVisitBit(0);
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++)
    {
        Face *face = mesh->getFaceAt(i);
        face->setVisitBit(0);
    }

    for (size_t i = 0; i < vEdges33.size(); i++)
    {
        int err = vEdges33[i].remove();
        if (!err) ncount++;
    }

    if (ncount)
    {
        mesh->prune();
        mesh->enumerate(0);
        mesh->enumerate(2);
        cout << "Info: number of bridges removed " << ncount << endl;
        mesh->collect_garbage();
    }

    if (!rel0exist)
        mesh->clearRelations(0, 0);

    if (!rel2exist)
        mesh->clearRelations(0, 2);

    vEdges33.clear();

    return ncount;
}

////////////////////////////////////////////////////////////////////

int
QuadCleanUp::remove_bridges()
{
    if (!mesh->isPruned()) mesh->prune();

    while (1)
    {
        int ncount = remove_bridges_once();
        if (ncount == 0) break;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
