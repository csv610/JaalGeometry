#include "Mesh.hpp"

JMeshPtr AllQuadMeshGenerator:: getStructuredMesh(int *grid_dim, double *length, double *origin)
{
    /*
        JMeshPtr mesh = JMesh::newObject();

        int nx = 2, ny = 2;
        if( grid_dim ) {
            nx = grid_dim[0];
            ny = grid_dim[1];
        }

        double xorg = 0.0, yorg = 0.0;
        if( origin  ) {
            xorg = origin[0];
            yorg = origin[1];
        }

        double xlen = 0.0, ylen = 0.0;
        if( length ) {
            xlen = length[0];
            ylen = length[1];
        }

        double  dx = 0.0, dy = 0.0;
        if( nx > 1) dx = xlen/ (double) (nx - 1);
        if( ny > 1) dy = ylen/ (double) (ny - 1);

        JNodeSequence nodes(nx * ny);

        int offset;
        Point3D xyz;
        xyz[0] = 0.0;
        xyz[1] = 0.0;
        xyz[2] = 0.0;
        int index = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    offset = k * nx * ny + j * nx + i;
                    JNodePtr v = JVertex::newObject();
                    xyz[0] = xorg + i*dx;
                    xyz[1] = yorg + j*dy;
                    v->setXYZCoords(xyz);
                    v->setID(index++);
                    nodes[ offset ] = v;
                }
            }
        }

        mesh->addObjects(nodes);
        JNodeSequence qnodes(4);
        int bmark = 1;

        JFaceSequence faces((nx-1)*(ny-1));
        index = 0;
        for (int j = 0; j < ny - 1; j++) {
                for (int i = 0; i < nx - 1; i++) {
                    offset = j * nx + i;
                    qnodes[0] = nodes[offset];
                    qnodes[1] = nodes[offset + 1];
                    qnodes[2] = nodes[offset + 1 + nx];
                    qnodes[3] = nodes[offset + nx];
                    faces[index++]  = Quadrilateral::newObject( qnodes );
                }
        }
        mesh->addObjects(faces);

        bmark = 1;
        for (int i = 0; i < nx - 1; i++) {
             JFacePtr face = mesh->getFaceAt(i);
            JEdgePtr edge = face->getEdgeAt(0);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 3;
            for (int i = 0; i < nx - 1; i++) {
                JFacePtr face = mesh->getFaceAt((nx-1)*(ny-1)- i-1);
                JEdgePtr edge = face->getEdgeAt(2);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 4;
            for (int j = 0; j < ny - 1; j++) {
                JFacePtr face = mesh->getFaceAt(j*(nx-1));
                JEdgePtr edge = face->getEdgeAt(3);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 2;
            for (int j = 0; j < ny - 1; j++) {
                JFacePtr face = mesh->getFaceAt(j*(nx-1)+(nx-2));
                JEdgePtr edge = face->getEdgeAt(1);
                edge->setAttribute("Boundary", bmark);
            }
    */

    return mesh;
}

JMeshPtr AllHexMeshGenerator::getStructuredMesh(int *grid_dim, double *length, double *origin)
{
    JMeshPtr mesh = JMesh::newObject();


    int nx = 1, ny = 1, nz = 1;

    if( grid_dim ) {
        nx = grid_dim[0];
        ny = grid_dim[1];
        nz = grid_dim[2];
    }

    double xorg = -0.0, yorg = -0.0, zorg = -0.0;
    if( origin  ) {
        xorg = origin[0];
        yorg = origin[1];
        zorg = origin[2];
    }

    /*
        double xlen = 0.0, ylen = 0.0, zlen = 0.0;
        if( length ) {
            xlen = length[0];
            ylen = length[1];
            if (space_dim == 3) zlen = length[2];
        }

        double  dx = 0.0, dy = 0.0, dz = 0.0;
        if( nx > 1) dx = xlen/ (double) (nx - 1);
        if( ny > 1) dy = ylen/ (double) (ny - 1);
        if (nz > 1) dz = zlen/ (double) (nz - 1);

        JNodeSequence nodes(nx * ny * nz);

        int offset;
        Point3D xyz;
        int index = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    offset = k * nx * ny + j * nx + i;
                    JNodePtr v = JVertex::newObject();
                    xyz[0] = xorg + i*dx;
                    xyz[1] = yorg + j*dy;
                    xyz[2] = zorg + k*dz;
                    v->setXYZCoords(xyz);
                    v->setID(index++);
                    nodes[ offset ] = v;
                    mesh->addObject(v);
                }
            }
        }
        JNodeSequence qnodes;
        int bmark = 1;

        if (space_dim == 2) {
            qnodes.resize(4);
            index = 0;
            for (int j = 0; j < ny - 1; j++) {
                for (int i = 0; i < nx - 1; i++) {
                    offset = j * nx + i;
                    qnodes[0] = nodes[offset];
                    qnodes[1] = nodes[offset + 1];
                    qnodes[2] = nodes[offset + 1 + nx];
                    qnodes[3] = nodes[offset + nx];
                    JFacePtr face = Quadrilateral::newObject( qnodes );
                    face->setID(index++);
                    mesh->addObject(face);
                }
            }

            bmark = 1;
            for (int i = 0; i < nx - 1; i++) {
                JFacePtr face = mesh->getFaceAt(i);
                JEdgePtr edge = face->getEdgeAt(0);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 3;
            for (int i = 0; i < nx - 1; i++) {
                JFacePtr face = mesh->getFaceAt((nx-1)*(ny-1)- i-1);
                JEdgePtr edge = face->getEdgeAt(2);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 4;
            for (int j = 0; j < ny - 1; j++) {
                JFacePtr face = mesh->getFaceAt(j*(nx-1));
                JEdgePtr edge = face->getEdgeAt(3);
                edge->setAttribute("Boundary", bmark);
            }

            bmark = 2;
            for (int j = 0; j < ny - 1; j++) {
                JFacePtr face = mesh->getFaceAt(j*(nx-1)+(nx-2));
                JEdgePtr edge = face->getEdgeAt(1);
                edge->setAttribute("Boundary", bmark);
            }
        }

        if (space_dim == 3) {
            qnodes.resize(8);
            index = 0;
            for (int k = 0; k < nz - 1; k++) {
                for (int j = 0; j < ny - 1; j++) {
                    for (int i = 0; i < nx - 1; i++) {
                        offset = k*nx*ny + j * nx + i;
                        qnodes[0] = nodes[offset];
                        qnodes[1] = nodes[offset + 1];
                        qnodes[2] = nodes[offset + 1 + nx];
                        qnodes[3] = nodes[offset + nx];
                        offset = (k+1)*nx*ny + j * nx + i;
                        qnodes[4] = nodes[offset];
                        qnodes[5] = nodes[offset + 1];
                        qnodes[6] = nodes[offset + 1 + nx];
                        qnodes[7] = nodes[offset + nx];
                        JCellPtr cell = Hexahedron::newObject();
                        cell->setID(index++);
                        cell->setNodes(qnodes);
                        mesh->addObject(cell);
                    }
                }
            }
        }
    */

    return mesh;
}


