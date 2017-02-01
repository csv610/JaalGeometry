#pragma once

#ifndef MOAB_JAALMESH_H
#define MOAB_JAALMESH_H

#include "Mesh.hpp"
#include "MeshQuality.hpp"

#include "SimpleArray.hpp"

#ifdef HAVE_MOAB

using namespace Jaal;

class JaalMoabConverter
{
public:
    void clear()
    {
        moabnode.clear();
        moabface.clear();
        jaalnode.clear();
        jaalface.clear();
    }
    //  Converts the mesh into MOAB data structures.
    int toMOAB(JMeshPtr mesh, iMesh_Instance &imesh, iBase_EntitySetHandle eset = 0);

    //  Fill the mesh from MOAB..
    JMeshPtr fromMOAB(iMesh_Instance imesh, JMeshPtr m, iBase_EntitySetHandle eset = 0);

private:
    JMeshPtr jmesh;

    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, JNodePtr v);
    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, JFacePtr f);
    iBase_EntityHandle new_MOAB_Handle(iMesh_Instance imesh, JCellPtr c);

    std::map<JNodePtr, iBase_EntityHandle> moabnode;
    std::map<JFacePtr, iBase_EntityHandle> moabface;
    std::map<JCellPtr, iBase_EntityHandle> moabcell;

    std::map<iBase_EntityHandle, JNodePtr> jaalnode;
    std::map<iBase_EntityHandle, JFacePtr> jaalface;
    std::map<iBase_EntityHandle, JCellPtr> jaalcell;
};

#endif

#endif
