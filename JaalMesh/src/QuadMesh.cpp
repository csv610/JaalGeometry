#include "meshkit/QuadMesh.hpp"
#include "meshkit/MKCore.hpp"
#include "meshkit/ModelEnt.hpp"
#include "Mesh.hpp"
#include "JaalMoabConverter.hpp"
#include "Tri2Quad.hpp"
#include <vector>
#include <math.h>

namespace MeshKit
{
moab::EntityType QuadMesher_tps[] = {moab::MBVERTEX, moab::MBQUAD, moab::MBMAXTYPE};

const moab::EntityType* QuadMesher::output_types()
{
    return QuadMesher_tps;
}

//---------------------------------------------------------------------------//
// Construction Function for Edge Mesher
QuadMesher::QuadMesher(MKCore *mk_core, const MEntVector &me_vec)
    : MeshScheme(mk_core, me_vec) //, schemeType(EQUAL)
{
}

// setup function: set up the number of intervals for edge meshing through the
// sizing function
void QuadMesher::setup_this()
{
    setup_boundary();

    // Make sure everything has a tri mesher or is meshed
    MEntVector need_trimesher;
    MEntSelection::iterator i;
    std::vector<moab::EntityType> types;
    std::vector<MeshOp*> ops;
    std::vector<MeshOp*>::iterator j;
    for (i = me_selection().begin(); i != me_selection().end(); ++i) {
        ModelEnt* ent = i->first;
        if (ent->get_meshed_state() == COMPLETE_MESH)
            continue;
        ops.clear();
        ent->get_meshops( ops );
        bool found_one = false;
        for (j = ops.begin(); j != ops.end(); ++j) {
            types.clear();
            (*j)->mesh_types( types );
            if (std::find(types.begin(), types.end(), moab::MBTRI) != types.end()) {
                found_one = true;
                break;
            }
        }

        if (!found_one)
            need_trimesher.push_back( ent );
    }

    if (need_trimesher.empty())
        return;

    // find a tri mesher
    std::vector<MeshOpProxy*> avail_ops;
    mk_core()->meshop_by_mesh_type( moab::MBTRI, avail_ops );
    if (avail_ops.empty())
        return; // throw exception or just hope that everthing already has a tri mesh?
    MeshOp* tri_mesher = mk_core()->construct_meshop( avail_ops.front(), need_trimesher );

    // tell core that the tri mesher needs to run before this
    mk_core()->insert_node( tri_mesher, this );
}

void QuadMesher::execute_this()
{
    Jaal::Tri2Quads t2quad;
    JaalMoabConverter meshconverter;

    iMesh_Instance imesh = mk_core()->imesh_instance()->instance();
    MEntSelection::iterator i;
    for (i = me_selection().begin(); i != me_selection().end(); ++i) {
        ModelEnt* ent = i->first;
        std::auto_ptr<Jaal::Mesh> trimesh(meshconverter.fromMOAB( imesh, (iBase_EntitySetHandle)ent->mesh_handle() ));
        std::auto_ptr<Jaal::Mesh> quadmesh(t2quad.getQuadMesh( trimesh.get(), 1));
        meshconverter.toMOAB( quadmesh.get(), imesh, (iBase_EntitySetHandle)ent->mesh_handle() );
//  int ierr;
//  int len =  strlen( "JVtk.vtk" );
//  iMesh_save( imesh, 0, "JVtk.vtk", NULL, &ierr, len, 0);
        ent->commit_mesh( i->second, COMPLETE_MESH );
    }
}

} // namespace MeshKit

