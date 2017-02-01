#include "Mesh.hpp"

using namespace Jaal;

void print_histogram( const vector<double> &quality, const string &header )
{

}
///////////////////////////////////////////////////////////////////////////////

void plot_all_quad_quality_measure( Mesh *mesh )
{
#ifdef HAVE_VERDICT
    VERDICT_REAL  coords[4][3];
    Point3D xyz;
    int num_nodes = 4;

    size_t numfaces = mesh->getSize(2);

    vector<double> quality( numfaces );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] =  v_quad_aspect( num_nodes, coords);
    }
    boost::sort( quality );
    print_histogram( quality, "Aspect Ratio"  );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_skew( num_nodes, coords );
    }
    boost::sort( quality.begin() );
    print_histogram( quality, "Shewness"  );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] =  v_quad_taper( num_nodes, coords );
    }
    boost::sort( quality.begin() );
    print_histogram( quality, "Taper" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i]  = v_quad_warpage( num_nodes, coords );
    }
    boost::sort( quality.begin() );
    print_histogram( quality, "Warpage" );


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_area( num_nodes, coords );
    }
    boost::sort( quality.begin() );
    print_histogram( quality, "Area" );


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_stretch( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Stretch" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_minimum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Min Angle" );


    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_maximum_angle( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Max Angle" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_oddy( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Oddy" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_condition( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Condition Number" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_jacobian( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Jacobian" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_scaled_jacobian( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Scaled Jacobian" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shear( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shear" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shape( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shape" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_relative_size_squared( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Relative Size Squared" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shape_and_size( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shape and Size" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_shear_and_size( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Shear and Size" );

    ///////////////////////////////////////////////////////////////////////////

    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        for( int j = 0; j < 4; j++) {
            Vertex *vertex = face->getNodeAt(j);
            xyz  = vertex->getXYZCoords();
            coords[j][0] = xyz[0];
            coords[j][1] = xyz[1];
            coords[j][2] = xyz[2];
        }
        quality[i] = v_quad_distortion( num_nodes, coords );
    }
    sort( quality.begin(), quality.end() );
    print_histogram( quality, "Distortion" );

    ///////////////////////////////////////////////////////////////////////////
#endif
}



