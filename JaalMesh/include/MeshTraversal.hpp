#pragma once

#ifndef JMESHTRAVERSE_H
#define JMESHTRAVERSE_H

#include "Mesh.hpp"

class JMeshTraversal
{
public:
    static const int  BREADTH_FIRST_ORDER   = 0;
    static const int  DEPTH_FIRST_ORDER     = 1;

    void  setMesh( const JMeshPtr &m) {
        mesh = m;
        initialized = 0;
    }

    // Start from a node of the mesh and traverse the mesh untill maxSize is reached.
    void  getToplogicalSequence( const JNodePtr &v, JNodeSequence &seq, size_t maxSize, int method = BREADTH_FIRST_ORDER);

    // Start from a node of the mesh and traverse the mesh untill maxdistance is reached.
    void  getGeometricSequence( const JNodePtr &v,  JNodeSequence &seq, double maxdistance, int method = BREADTH_FIRST_ORDER);

    // Start from a node of the mesh and traverse the mesh untill maxSize is reached.
    void  getTopologicalSequence( const JNodePtr &v, JFaceSequence &seq, size_t maxSize, int method = BREADTH_FIRST_ORDER);

    // Start from a node of the mesh and traverse the mesh untill maxdistance is reached.
    void  getGeometricSequence( const JNodePtr &v,  JFaceSequence &seq, double maxdistance, int method = BREADTH_FIRST_ORDER);

    // Start from a node of the mesh and traverse the mesh untill maxSize is reached.
    void  getToplogicalSequence( const JFacePtr &v, JFaceSequence &seq, size_t maxSize, int method = BREADTH_FIRST_ORDER);

    // Start from a node of the mesh and traverse the mesh untill maxdistance is reached.
    void  getGeometricSequence( const JFacePtr &v,  JFaceSequence &seq, double maxdistance, int method = BREADTH_FIRST_ORDER);

private:
    JMeshPtr mesh;
    bool initialized;
    void init();
    void init02();
};

#endif
