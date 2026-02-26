#include <iostream>
#include <string>
#include "Mesh.hpp"

#include "StopWatch.hpp"
#include "SwapTriEdge.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int advancing_front(JMeshPtr trimesh)
{
    Jaal::MeshOptimization mopt;
    mopt.shape_optimize( trimesh );

    trimesh->get_topological_statistics();

    Jaal::advancing_front_triangle_cleanup ( trimesh );
    exit(0);

    SwapTriEdge edgeswapper(trimesh);
//    edgeswapper.apply_degree_reduction_rule();
//    cout << " After Degree Reduction Sweeping  : " << endl;
//    trimesh->get_topological_statistics();

    edgeswapper.apply_advance_front_rule();
    trimesh->get_topological_statistics();

    plot_all_tri_quality_measures( trimesh );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if( argc != 3)
    {
        cout << "Usage: Executable infile outfile" << endl;
        return 1;
    }


    JMeshPtr trimesh = Mesh::newObject();
    trimesh->readFromFile( argv[1] );

    StopWatch watch;
    watch.start();

    int nGroups = trimesh->getNumOfGroups(2);

    size_t numNodes;

    JNodeSet qnodes;
    JNodeSet::const_iterator it;

    vector<Face*> qfaces;

    if( nGroups > 1)
    {
        for( int ig = 0; ig < nGroups; ig++)
        {
            JMeshPtr subtrimesh = trimesh->getPartMesh( ig );
            assert( subtrimesh );
            numNodes = subtrimesh->getSize(0);
            for( size_t j = 0; j < numNodes; j++)
            {
                JNodePtr v = subtrimesh->getNodeAt(j);
                v->setID(j);
            }
            advancing_front( subtrimesh );
        }
    }
    else
        advancing_front(trimesh );

    watch.stop();
    Jaal::set_layer_tag(trimesh);
    trimesh->saveAs( argv[2] );
    trimesh->deleteAll();
}

///////////////////////////////////////////////////////////////////////////////
