#include <iostream>
#include <string>

#include "JaalMoabConverter.hpp"
#include "StopWatch.hpp"
#include "MeshRefine.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if( argc != 3)
    {
        cout << "Usage: Executable infile outfile" << endl;
        return 1;
    }

    MeshOptimization mopt;

    Mesh *trimesh = new Mesh2D;
    trimesh->readFromFile( argv[1] );

    Mesh *quadmesh = NULL;
    StopWatch watch;
    watch.start();

    int nGroups = 1;
// int nGroups = trimesh->getNumPartitions(2);

    size_t numNodes, numFaces;

    JNodeSet qnodes;
    JNodeSet::const_iterator it;
    JFaceSequence qfaces;

    /*
        if( nGroups > 1)
        {
            for( int ig = 0; ig < nGroups; ig++)
            {
                Mesh *subtrimesh = trimesh->getPartMesh( ig );
                assert( subtrimesh );
                numNodes = subtrimesh->getSize(0);
                for( size_t j = 0; j < numNodes; j++)
                {
                    Vertex *v = subtrimesh->getNodeAt(j);
                    v->setID(j);
                }
                Mesh *subquadmesh = tri_quad_conversion( subtrimesh );
                assert( subquadmesh );
                int nQuads = subquadmesh->getSize(2);
                for( int j = 0; j < nQuads; j++)
                {
                    Face *f = subquadmesh->getFaceAt(j);
                    f->setAttribute("Partition", ig );
                }

                numNodes = subquadmesh->getSize(0);
                for( size_t j = 0; j < numNodes; j++)
                    qnodes.insert( subquadmesh->getNodeAt(j));

                numFaces = subquadmesh->getSize(2);
                for( size_t j = 0; j < numFaces; j++)
                    qfaces.push_back( subquadmesh->getFaceAt(j));

                subquadmesh->emptyAll();
                delete subquadmesh;
            }

            quadmesh = new Mesh;
            numNodes = qnodes.size();
            int index = 0;
            for( it = qnodes.begin(); it != qnodes.end(); ++it)
            {
                Vertex *v = *it;
                v->setID(index++);
                quadmesh->addNode( v );
            }

            numFaces = qfaces.size();
            for( size_t i = 0; i < numFaces; i++)
                quadmesh->addFace( qfaces[i] );
        }
        else
        {
            quadmesh = tri_quad_conversion( trimesh );
        }
    */

    quadmesh = AllQuadMeshGenerator::BinaryTreeMatching( trimesh );
//  quadmesh = AllQuadMeshGenerator::HamiltonianQuads(trimesh );

    watch.stop();

    if( quadmesh )
    {
        cout << "Info: Tri-Quad conversion time (sec) : " << watch.getSeconds() << endl;
        quadmesh->saveAs( argv[2] );
        mopt.shape_optimize( quadmesh );
        quadmesh->saveAs( "optquad.off" );
//      plot_all_quad_quality_measures( quadmesh );
        quadmesh->deleteFaces();
        trimesh->deleteFaces();
        quadmesh->deleteNodes();
        delete quadmesh;
    }

    delete trimesh;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
