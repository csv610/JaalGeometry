#pragma once

#include "Mesh.hpp"
#include "MeshTopology.hpp"

using namespace Jaal;

class JTriDecimator
{
public:
    JTriDecimator() {
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    bool isCollapsable(const JEdgePtr &edge, bool check567 = 1);
    int  collapse(const JEdgePtr &edge);
    int  collapseEdges();

    void removeLowDegreeNodes();
    int  removeHighDegreeNodes();

    const JNodeSequence &getNewNodes() const {
        return newNodes;
    }
    const JEdgeSequence &getNewEdges() const {
        return newEdges;
    }
    const JFaceSequence &getNewFaces() const {
        return newFaces;
    }

private:
    JMeshPtr mesh;
    JNodeSequence nodeneighs;
    JEdgeSequence coveredges;
    JFaceSequence faceneighs, f0neighs, f1neighs;

    JNodeSequence newNodes;
    JEdgeSequence newEdges;
    JFaceSequence newFaces;

    void remove_degree3_nodes();
    void remove_degree4_nodes();
    int  remove_below_degree5_node(const JNodePtr &v);
    int  remove_degree3_node(const JNodePtr &v);
    int  remove_degree4_node(const JNodePtr &v);

    int  remove_above_degree7_node(const JNodePtr &v);
    int  remove_above_degree8_node(const JNodePtr &v);
    int  remove_degree8_node(const JNodePtr &v);

    void updateRelations( const JFacePtr &f);
    void smooth_local_patch( const JNodePtr &v, const JNodeSequence &neighs);

};

