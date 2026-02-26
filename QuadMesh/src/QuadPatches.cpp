#include "Mesh.hpp"

using namespace Jaal;

//
// QTrack is similar to "Motorcycle graph" proposed by Eppestein group.
// But somehow, I don't like this term.
//

#ifdef CSV
struct QTrack {
    const static int END_AT_TERMINALS = 0;
    const static int END_AT_CROSSINGS = 1;

    JNodeSequence sequence;

    bool operator ==(const QTrack & rhs) const {
        size_t nSize = sequence.size();
        if (nSize != rhs.sequence.size()) return 0;

        Vertex *v0src = sequence.front();
        Vertex *v0dst = sequence.back();
        Vertex *v1src = rhs.sequence.front();
        Vertex *v1dst = rhs.sequence.back();

        if (v0src == v1src && v0dst == v1dst) return 1;
        if (v0src == v1dst && v0dst == v1src) return 1;

        return 0;
    }

    bool operator<(const QTrack & rhs) const {
        return sequence.size() < rhs.sequence.size();
    }

    /////////////////////////////////////////////////////////////////////////////
    // There are two ways to advance along the track.
    // (1) Incremental:  All track propagate simulateneously and
    //                   at the intersection, only one is allowed
    //                   to proceed towards other irregular node.
    // (11) Greedy   :   Start from one vertex and complete its
    //                   track.
    // It is not clear which method is better, but "greedy" may likely
    // give very high aspect ratio quad patches. Incremental on the
    // other hand may produce many small patches..
    //
    /////////////////////////////////////////////////////////////////////////////

    int advance_single_step(int endat) {
        //////////////////////////////////////////////////////////////////////////
        //                   **********************
        //                   *         *          *
        //                   *         * Next     *
        //                   *         *          *
        //           Avoid   **********************  Avoid
        //                   *         *          *
        //                   *         * Current  *
        //                   *         *          *
        //                   *         *          *
        //                   **********************
        //                            Source
        // A Source vertex and Current edge is chosen.
        // We want to avoid two edges and want to select "Next" edge.
        //////////////////////////////////////////////////////////////////////////
        Vertex *v0, *v1, *v2, *v3, *v4;

        JFaceSequence  adjFaces;
        JNodeSequence  vneighs;
        JNodeSet       vset;

        size_t index = sequence.size();
        v0 = sequence[index - 2];
        v1 = sequence[index - 1];
        v0->setVisitBit(1);

        if (endat == END_AT_CROSSINGS && v1->getVisitBit()) return 0;
        if (v1->isBoundary()) return 0;

        v1->setVisitBit(1);
        vneighs = v1->getRelations0();
        if (vneighs.size() != 4) return 0;

        adjFaces = Mesh::getRelations112(v0, v1);
        assert(adjFaces.size() == 2);
        v2 = Face::opposite_node(adjFaces[0], v0);
        v3 = Face::opposite_node(adjFaces[1], v0);

        vset.clear();
        vset.insert(vneighs[0]);
        vset.insert(vneighs[1]);
        vset.insert(vneighs[2]);
        vset.insert(vneighs[3]);
        vset.erase(v0);
        vset.erase(v2);
        vset.erase(v3);
        assert(vset.size() == 1);
        v4 = *vset.begin();
        sequence.push_back(v4);
        return 1;
    }

    void advance(int endat) {
        assert(sequence.size() == 2);

        // Starting node is always irregular ...
        vector<Face*> vfaces = sequence[0]->getRelations2();
        assert(vfaces.size() != 4);

        while (1) {
            int progress = advance_single_step(endat);
            if (!progress) break;
        }

        /*
                // The path is reversible and therefore, we will give a direction
                // from the lower source node to higher destination node.
                if (sequence.front() > sequence.back())
                    reverse(sequence.begin(), sequence.end());
                assert(sequence.front() < sequence.back());
        */
        // Checking the correctness..
        if (endat == END_AT_TERMINALS) {
            vfaces = sequence.front()->getRelations2();
            assert(vfaces.size() != 4);
            vfaces = sequence.back()->getRelations2();
            assert(vfaces.size() != 4);
            for (int i = 1; i < sequence.size() - 1; i++) {
                vfaces = sequence[i]->getRelations2();
                assert(vfaces.size() == 4);
            }

        }
    }

};

///////////////////////////////////////////////////////////////////////////////

struct StructuredMesh2D {

    StructuredMesh2D() {
        nx = ny = 0;
    }
    int nx, ny;

    JFaceSequence faces;
    JNodeSequence nodes;
    JNodeSet      cornerNodes;

    void clearAll() {
        nodes.clear();
        faces.clear();
        neighs.clear();
        cornerNodes.clear();
        nx = 0;
        ny = 0;
    }

    bool operator<(const StructuredMesh2D & rhs) const {
        return this->getSize(2) < rhs.getSize(2);
    }

    size_t getSize(int e) const {
        if (e == 0) return nodes.size();
        if (e == 2) return faces.size();
        return 0;
    }

    int myID;
    vector<int> neighs;
};
#endif

///////////////////////////////////////////////////////////////////////////////

/*
void build_submesh_topology(JStructuredMesh &smesh)
{
    smesh.neighs.clear();

    JFaceSequence adjfaces;
    set<int> nset;

    size_t numFaces = smesh.faces.size();
    for (size_t i = 0; i < numFaces; i++) {
        Face *face = smesh.faces[i];
        for (int j = 0; j < 4; j++) {
            Vertex *v0 = face->getNodeAt(j + 0);
            Vertex *v1 = face->getNodeAt(j + 1);
            Mesh::getRelations112(v0, v1, adjfaces);
            int numneighs = adjfaces.size();
            for (int k = 0; k < numneighs; k++) {
                int nid = adjfaces[k]->getTag();
                nset.insert(nid);
            }
        }
    }

    nset.erase(smesh.myID);
    set<int>::const_iterator it1;
    for (it1 = nset.begin(); it1 != nset.end(); ++it1)
        smesh.neighs.push_back(*it1);
}
*/

///////////////////////////////////////////////////////////////////////////////

JFacePtr getSeedFace(JMeshPtr mesh) {
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        if (!f->getVisitBit()) return f;
    }
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

size_t independent_components(JMeshPtr mesh) {

    size_t numfaces = mesh->getSize(2);
    size_t numneighs;

    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        f->setVisitBit(0);
        f->setAttribute( "Partition", 0);
    }

    size_t compid = 0;
    deque<JFacePtr> faceQ;
    JFaceSequence adjfaces;
    while (1) {
        JFacePtr currface = getSeedFace(mesh);
        if (currface == NULL) break;
        faceQ.push_back(currface);
        while (!faceQ.empty()) {
            currface = faceQ.front();
            faceQ.pop_front();
            if (!currface->getVisitBit()) {
                currface->setAttribute( "Partition", compid);
                currface->setVisitBit(1);
                for (int i = 0; i < 4; i++) {
                    JNodePtr v0 = currface->getNodeAt(i + 0);
                    JNodePtr v1 = currface->getNodeAt(i + 1);
                    if (v0->getVisitBit() && v1->getVisitBit()) continue;
                    JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
                    JEdge::getRelations(edge, adjfaces);
//                  Mesh::getRelations112(v0, v1, adjfaces);
                    numneighs= adjfaces.size();
                    for (size_t j = 0; j < numneighs; j++)
                        faceQ.push_back(adjfaces[j]);
                }
            }
        }
        compid++;
    }

#ifdef DEBUG
    for (size_t i = 0; i < numfaces; i++) {
        Face *f = mesh->getFaceAt(i);
        assert(f->getVisitBit());
    }
#endif

    return compid;
}

/////////////////////////////////////////////////////////////////////////////////

/*
int merge_submesh(StructuredMesh2D &amesh, JStructuredMesh &bmesh) {
    if (amesh.myID == bmesh.myID) return 1;

    if (amesh.cornerNodes.size() != 4) return 2;
    if (bmesh.cornerNodes.size() != 4) return 2;

    JNodeSequence common;
    set_intersection(amesh.cornerNodes.begin(), amesh.cornerNodes.end(),
            bmesh.cornerNodes.begin(), bmesh.cornerNodes.end(),
            back_inserter(common));

    if (common.size() != 2) return 3;

    size_t numFaces;

    JFaceSequence vfaces;
    numFaces = amesh.faces.size();
    for (size_t i = 0; i < numFaces; i++) {
        for (int j = 0; j < 4; j++) {
            JNodePtr v = amesh.faces[i]->getNodeAt(j);
            if (!v->isBoundary()) {
                if ( v->getNumRelations(2) != 4) return 4;
            }
        }
    }

    numFaces = bmesh.faces.size();
    for (size_t i = 0; i < numFaces; i++) {
        for (int j = 0; j < 4; j++) {
            JNodePtr v = bmesh.faces[i]->getNodeAt(j);
            if (!v->isBoundary()) {
                if ( v->getNumRelations(2) != 4) return 4;
            }
        }
    }

    numFaces = bmesh.faces.size();
    for (size_t i = 0; i < numFaces; i++) {
        amesh.faces.push_back(bmesh.faces[i]);
        bmesh.faces[i]->setAttribute("Partition", amesh.myID);
    }

    JNodeSet::const_iterator it;
    for (it = bmesh.cornerNodes.begin(); it != bmesh.cornerNodes.end(); ++it)
        amesh.cornerNodes.insert(*it);

    assert(amesh.cornerNodes.size() == 6);

    amesh.cornerNodes.erase(common[0]);
    amesh.cornerNodes.erase(common[1]);

    bmesh.clearAll();

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////

JStructuredMesh getQuadPatch(JMeshPtr mesh, int compid) {
    StructuredMesh2D smesh;
    smesh.myID = compid;
    smesh.nx = 1;
    smesh.ny = 1;

    size_t numfaces = mesh->getSize(2);

    JFaceSet faceSet;
    JFaceSet::const_iterator fit;
    JNodeSet nodeSet;
    JNodeSet::const_iterator nit;

    int ival = 0;
    for (size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
        face->getAttribute("Partition", ival );
        if ( ival == compid) {
            faceSet.insert(face);
            nodeSet.insert(face->getNodeAt(0));
            nodeSet.insert(face->getNodeAt(1));
            nodeSet.insert(face->getNodeAt(2));
            nodeSet.insert(face->getNodeAt(3));
        }
        }
    }
    if (faceSet.empty()) return smesh;

    JNodeSet boundNodes, cornerNodes;
    JFaceSet boundFaces, cornerFaces;

    JFaceSequence vfaces;
    int     ncount;
    size_t  numneighs;
    for( JNodePtr vertex: nodeSet) {
        vertex->getRelations( vfaces);
        ncount = 0;
        numneighs = vfaces.size();
        for (size_t j = 0; j < numneighs; j++)
            if (faceSet.find(vfaces[j]) != faceSet.end()) ncount++;

        if (ncount == 1) {
            cornerNodes.insert(vertex);
            for (size_t j = 0; j < numneighs; j++)
                if (faceSet.find(vfaces[j]) != faceSet.end()) cornerFaces.insert(vfaces[j]);
        }

        if (ncount != 4) boundNodes.insert(vertex);
    }

    for (fit = faceSet.begin(); fit != faceSet.end(); ++fit)
        smesh.faces.push_back(*fit);

    for (nit = nodeSet.begin(); nit != nodeSet.end(); ++nit)
        smesh.nodes.push_back(*nit);

    smesh.cornerNodes = cornerNodes;

    build_submesh_topology(smesh);
    return smesh;

}
*/
/////////////////////////////////////////////////////////////////////////////////

#ifdef CSV
int Mesh::search_quad_patches()
{
    int nTopo = getElementsType();
    if (nTopo != JFace::QUADRILATERAL) {
        cout << "Error: The mesh must be all Quads " << endl;
        return 1;
    }

    int relexist2 = buildRelations(0, 2);
    int relexist0 = buildRelations(0, 0);

    search_boundary();

    size_t numnodes = getSize(0);
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        vertex->setVisitBit(0);
    }

    size_t numfaces = getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = getFaceAt(i);
        face->setVisitBit(0);
    }

    vector<QTrack> qpath;
    QTrack qp;
    qp.sequence.resize(2); // As we know the starting edge

    size_t nCount = 0;
    JNodeSequence vnodes;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        if (!vertex->isBoundary()) {
            vertex->getRelations( vnodes );
            size_t numneighs = vnodes.size();
            if (numneighs  != 4)  {
                qp.sequence[0] = vertex;
                for (size_t j = 0; j < numneighs; j++) {
                    qp.sequence[1] = vnodes[j];
                    nCount++;
                    qpath.push_back(qp);
                }
            }
        }
    }

    if (qpath.empty()) {
        cout << "Info: There are no irregular nodes in the mesh" << endl;
        return 0;
    } else
        cout << "# of branches spawned " << nCount << endl;

    for (size_t j = 0; j < qpath.size(); j++)
        qpath[j].advance(0);

    sort(qpath.begin(), qpath.end());

    /*
        saveAs("b.dat");
        for (int i = 0; i < qpath.size(); i++) {
            if (qpath[i].sequence.front()->isBoundary()) continue;
            if (qpath[i].sequence.back()->isBoundary() ) continue;
            int found = 0;
            cout << " PATH SIZE " << qpath[i].sequence.size() << endl;
            for (int k = 0; k < qpath[i].sequence.size(); k++)
                   cout << qpath[i].sequence[k]->getID() << " ";
            cout << endl;
            for (int j = 0; j < qpath.size(); j++) {
                if ((i != j) && (qpath[i] == qpath[j])) {
                    found = 1;
                    cout << "Same path " << i << " " << j << endl;
                    for (int k = 0; k < qpath[j].sequence.size(); k++)
                        cout << qpath[j].sequence[k]->getID() << " ";
                    cout << endl;
                }
            }
            assert(found);
        }
    */
    exit(0);

    int numPatches = independent_components(this);

    deque<StructuredMesh2D> submesh(numPatches);

    // New we can merge some sub-meshes
    for (int i = 0; i < numPatches; i++) {
        StructuredMesh2D sm = getQuadPatch(this, i);
        submesh[i] = sm;
    }
    sort(submesh.begin(), submesh.end());

    while (1) {
        int  nSize = submesh.size();
        cout << "#Submeshes " << submesh.size() << endl;
        size_t count_merged = 0;
        for (int i = 0; i < nSize; i++) {
            for (int j = i + 1; j < nSize; j++) {
                int err = merge_submesh(submesh[i], submesh[j]);
                if (!err) count_merged++;
            }
        }
        cout << " #of Submesh merged " << count_merged << endl;
        if (count_merged == 0) break;

        //
        // Retain only those submeshes which are not empty. ( Note
        // that when two submeshes merge, one of them become empty.
        //
        for (int i = 0; i < nSize; i++) {
            StructuredMesh2D sm = submesh.front();
            submesh.pop_front();
            if (sm.getSize(2)) submesh.push_back(sm);
        }
    }

    sort(submesh.begin(), submesh.end());

    exit(0);

    /*
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *vertex = getNodeAt(i);
        vertex->setTag(vertex->isVisited());
    }
    */

    if (relexist0)
        clearRelations(0, 0);

    if (relexist2)
        clearRelations(0, 2);
    return submesh.size();
}
#endif
