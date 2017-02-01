#pragma once

#include "Mesh.hpp"
#include "MeshOptimization.hpp"
#include "MeshRefine.hpp"

class JFaceClose  : public JMeshTopologyOptimization {
public:
    static const int  MID_POINT_POSITION    = 0;
    static const int  FIRST_POINT_POSITION  = 1;
    static const int  SECOND_POINT_POSITION = 2;

    explicit JFaceClose( JMeshPtr m )
    {
        mesh = m;
        face = nullptr;
        vertex0 = nullptr;
        vertex2 = nullptr;
        check_inversion = 1;
        node_placement_policy = 1;
    }

    void setNodePlacementPolicy ( int p)
    {
        node_placement_policy = p;
    }

    void checkInversion( bool b)
    {
        check_inversion = b;
    }

    int apply(JFacePtr f, int pos)
    {
        face  = f;
        replaceNode = nullptr;
        if( !face->isActive() )  return 1;
        if( pos < 0) return 2;

        // Only for Quadrilateral element..
        if( face->getTypeID() != JFace::QUADRILATERAL ) return 2;

        vertex0 = face->getNodeAt( pos );
        vertex2 = face->getNodeAt( pos + 2);

        int err = build();
        if( err ) return 3;
        return commit();
    }

private:
    JFacePtr face;
    JNodePtr vertex0, vertex2;
    int node_placement_policy;
    int check_inversion;

    bool isSafe() const;
    int  build();
    int  commit();

    JNodePtr replaceNode;
};

///////////////////////////////////////////////////////////////////////////////////
// Diamond:  An element whose at least one of the opposite vertex is surrounded by
//           three faces. In many cases, diamonds are essential in the quadrilateral
// mesh and they can not be removed, Finding the minimum number of diamonds is hard,
// and we are working towards it.
///////////////////////////////////////////////////////////////////////////////////

class JDiamond : public JMeshTopologyOptimization {
public:
    static int  addAttribute(const JMeshPtr &mesh, const string &s = "Diamond" );
    static bool isDiamond(const JFacePtr &f, int &pos, int type = 33);

    ~JDiamond()
    {
        faceclose.reset();
    }

    JFaceSequence getDiamonds(int type = 33 );

    int remove(const JFacePtr &f)
    {
        int pos;
        if( !isDiamond(f, pos) ) return 1;
        return faceclose->apply(f, pos);
    }

    int remove(const JFacePtr &f, int pos)
    {
        return faceclose->apply(f, pos);
    }

private:
    boost::scoped_ptr<JFaceClose> faceclose;
};

///////////////////////////////////////////////////////////////////////////////////
