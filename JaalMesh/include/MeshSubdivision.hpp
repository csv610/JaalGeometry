#pragma once

#include "Mesh.hpp"

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Subdivision_method_3.h>
#include <iostream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#endif

class JMeshSubdivision
{
#ifdef USE_CGAL
    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef CGAL::Polyhedron_3<Kernel>     CGALPolyhedron;
#endif

public:
    static const int  LOOP   = 0;
    static const int  SQRT3  = 1;
    static const int  DOOSABIN = 3;
    static const int  BUTTERFLAY = 4;
    static const int  CATMULL_CLARK = 5;
    static const int  MODIFIED_BUTTERFLAY = 6;

    void setMesh( const JMeshPtr &m) {
        inMesh = m;
    }

    void setAlgorithm( int a ) {
        algo = a;
    }
    void setIterations( int d)  {
        numIters = d;
    }

    JMeshPtr getSubdivided();
private:
    JMeshPtr   inMesh;
    int        algo;
    int        numIters = 1;
};

/*
class CatmullClarkSubDivision : public JMeshSubDivision
{
public:
    JMeshPtr refineAll(int n);
    void  smoothSurface(int n);
    void  limitSurface();
private:
    JFaceSequence vfaces;
    JEdgeSequence vedges;
    void apply_face_rule();
    void apply_edge_rule();
    void apply_node_rule();
    void build_newmesh();
    void limitPoint(JNodePtr v);
};

class JLoopSubDivision : public JMeshSubDivision
{
public:
    void refineAll();
private:
    void apply_edge_rule();
    void build_new_faces();
};

struct JButterflySubDivision : public JMeshSubDivision {
};

struct JModifiedButterflySubDivision : public JMeshSubDivision {
};

///////////////////////////////////////////////////////////////////////////////
*/
