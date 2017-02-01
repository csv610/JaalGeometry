#include "Mesh.hpp"
#include "MeshRefine.hpp"
#include "QuadCleanUp.hpp"

#include "DijkstraShortestPath.hpp"
#include "ObjectPool.hpp"

extern int QuadPatches(Jaal::Mesh *mesh);

using namespace Jaal;

void usage()
{
    cout << "Usage: Executable -i in_meshfile -o out_meshfile -c cleanOp " << endl;

    cout << " *****************************************************************" << endl;
    cout << " Option :   Mesh Cleanup Operation " << endl;
    cout << " *****************************************************************" << endl;
    cout << " 0     :   Report Mesh Quality " << endl;

    cout << " 1     :   Remove interior doublets  " << endl;
    cout << " 2     :   Remove boundary singlets  " << endl;
    cout << " 3     :   Remove diamonds " << endl;
    cout << " 4     :   Vertex Degree Reduction " << endl;
    cout << " 5     :   Regularization with remeshing " << endl;
    cout << " 6     :   Shape Optimization  " << endl;
    cout << " 7     :   Everything automatic " << endl;

    cout << " 8     :   Reverse Elements Connectivity  " << endl;
    cout << " 9     :   Generate Quad-Irregular(Motorcycle) Graph " << endl;
    cout << " 10    :   Refine QuadMesh ( Scheme 14 )  " << endl;
    cout << " 11    :   Refine QuadMesh ( Scheme 15 )  " << endl;
    cout << " 12    :   Generate Quad-to-Tri4 " << endl;
    cout << " 13    :   Generate Quad-to-Tri2 " << endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    Jaal::Mesh *mesh = new Jaal::Mesh2D;
    Jaal::MeshOptimization mopt;

    string infilename, outfilename;

    int iopt, numiters = 100, topo_improve = 1;
    int cleanup_op = 0;

    while( (iopt = getopt( argc, argv, "hc:i:o:") ) != -1) {
        switch(iopt) {
        case 'c':
            cleanup_op = atoi( optarg );
            break;
        case 'h':
            usage();
            break;
        case 'i':
            infilename = optarg;   // Input QuadMesh File
            break;
        case 'o':
            outfilename = optarg; // Output QuadMesh File
            break;
        default:
            cout << "Usage: Executable -i in_meshfile -o out_meshfile -c cleanOp " << endl;
            break;
        }
    }

    if( infilename.empty() ) {
        cout <<"Warning: No input file specified " << endl;
        usage();
        return 1 ;
    }

    if( outfilename.empty() ) {
        cout <<"Warning: No output file specified " << endl;
        usage();
        return 2;
    }
    mesh->readFromFile( infilename );

#ifdef CSV
//   size_t ninvert  =  mesh->count_inverted_faces();
    size_t numfaces =  mesh->getSize(2);
    size_t numBound =  mesh->getBoundarySize(0);
    size_t nireg0   =  mesh->count_irregular_nodes(4);

    cout << "# of irregular nodes before cleanup : " << nireg0 << endl;

    if( ninvert > 0.5*numfaces ) mesh->reverse();

    QuadCleanUp qClean(mesh);

    vector<QTrack>  qpath;
    vector<Vertex*> steiner;
    Mesh *q2t;
    int  algo, numiter;

    StopWatch swatch;
    swatch.start();

    switch( cleanup_op ) {
    case 0:
        qClean.report();
        break;
    case 1:
        qClean.remove_interior_doublets();
        break;
    case 2:
        qClean.remove_boundary_singlets();
        break;
    case 3:
        qClean.remove_diamonds();
        break;
    case 4:
        qClean.vertex_degree_reduction();
        break;
    case 5:
        qClean.remesh_defective_patches();
        break;
    case 6:
        /*
                  cout << "Choose algorithm : " << endl;
                  cout << "    Steepest Descent       : 0 " << endl;
                  cout << "    Quasi Newton(default)  : 1 " << endl;
                  cout << "    Trust Region           : 2 " << endl;
                  cout << "    Feasible Newton        : 3 " << endl;
                  cout << "    Laplacian              : 4 " << endl;
                  cin  >> algo;
                  cout << "Give number of iterations " << endl;
                  cin  >> numiter;
                  mopt.shape_optimize( mesh, algo, numiter );
        */
        mopt.shape_optimize( mesh );
        break;
    case 7:
        qClean.automatic();
        break;
    case 8:
        mesh->reverse();
        break;
    case 9:
        qpath = Jaal::generate_quad_partitioning(mesh);
        break;
    case 10:
        mesh->refine_quads14();
        break;
    case 11:
        mesh->refine_quads15();
        break;
    case 12:
        q2t = Jaal::quad_to_tri4( mesh, steiner);
        q2t->saveAs("tmesh.off");
        break;
    case 13:
        q2t = Jaal::quad_to_tri2( mesh );
        mopt.shape_optimize( q2t );
        q2t->saveAs("tmesh.off");
        break;
    }

    swatch.stop();
    cout << "CleanUp time : " << swatch.getSeconds() << endl;

    if( cleanup_op ) {

        cout << "# Nodes           : " << mesh->getSize(0) << endl;
        cout << "# Faces           : " << mesh->getSize(2) << endl;
        cout << "# Inverted Faces  : " << mesh->count_inverted_faces() << endl;
        cout << "# Concave Faces   : " << mesh->count_concave_faces() << endl;
        cout << "# Irregular nodes : " << mesh->count_irregular_nodes(4) << endl;
//        cout << "Mesh Consistency  : " << mesh->is_consistently_oriented() << endl;
    }

    mesh->collect_garbage();

    cout << " Saving Mesh " << outfilename << endl;
    mesh->saveAs( outfilename);

    mesh->get_topological_statistics();
    plot_all_quad_quality_measures( mesh );

    assert( numBound == mesh->getBoundarySize(0) );

    if( mesh ) mesh->deleteAll();

    delete mesh;
#endif

    return 0;
}

