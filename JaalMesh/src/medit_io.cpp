#include "MeshImporter.hpp"

JMeshPtr JMeshMeditImporter ::readFile( const string &fname)
{
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << fname << endl;
        return nullptr;
    }
    JMeshPtr mesh = JMesh::newObject();

    size_t numnodes, numfaces, numcells;
    int n0, n1, n2, n3, n4, n5, n6, n7;
    int ref;
    JCellPtr newcell;
    JFacePtr newface;
    string str;
    JNodeSequence connect;

    while(!infile.eof() ) {
        infile >> str;
        if( str == "MeshVersionFormatted" ) infile >> str;
        if( str == "Dimension" ) infile >> str;
        if( str == "Vertices" )  {
            infile >> numnodes;
            double x, y, z = 0.0;
            Point3D xyz;
            for( size_t i = 0; i < numnodes; i++)  {
                infile >> x >> y  >> z >> ref;
                xyz[0] = x;
                xyz[1] = y;
                xyz[2] = z;
                JNodePtr v = JNode::newObject();
                v->setID(i);
                v->setXYZCoords(xyz);
                mesh->addObject(v);
            }
        }

        if( str == "Triangles" ) {
            infile >> numfaces;
            connect.resize(3);
            for( size_t i = 0; i < numfaces; i++) {
                infile >> n0 >> n1 >> n2 >> ref;
                connect[0] = mesh->getNodeAt(n0-1);
                connect[1] = mesh->getNodeAt(n1-1);
                connect[2] = mesh->getNodeAt(n2-1);
                newface = JTriangle::newObject(connect);
                mesh->addObject( newface );
            }
        }

        if( str == "Quadrilatrals" ) {
            infile >> numfaces;
            connect.resize(4);
            for( size_t i = 0; i < numfaces; i++) {
                infile >> n0 >> n1 >> n2 >> n3 >> ref;
                connect[0] = mesh->getNodeAt(n0-1);
                connect[1] = mesh->getNodeAt(n1-1);
                connect[2] = mesh->getNodeAt(n2-1);
                connect[3] = mesh->getNodeAt(n3-1);
                newface = JQuadrilateral::newObject( connect );
                mesh->addObject( newface );
            }
        }

        if( str == "Tetrahedra" ) {
            infile >> numcells;
            connect.resize(4);
            for( size_t i = 0; i < numcells; i++) {
                infile >> n0 >> n1 >> n2 >> n3 >> ref;
                connect[0] = mesh->getNodeAt(n0-1);
                connect[1] = mesh->getNodeAt(n1-1);
                connect[2] = mesh->getNodeAt(n2-1);
                connect[3] = mesh->getNodeAt(n3-1);
                newcell = JTetrahedron::newObject();
                newcell->setNodes( connect );
                mesh->addObject( newcell );
            }
        }

        if( str == "Hexahedra" ) {
            infile >> numcells;
            connect.resize(8);
            for( size_t i = 0; i < numcells; i++) {
                infile >> n0 >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> ref;
                connect[0] = mesh->getNodeAt(n0-1);
                connect[1] = mesh->getNodeAt(n1-1);
                connect[2] = mesh->getNodeAt(n2-1);
                connect[3] = mesh->getNodeAt(n3-1);

                connect[4] = mesh->getNodeAt(n4-1);
                connect[5] = mesh->getNodeAt(n5-1);
                connect[6] = mesh->getNodeAt(n6-1);
                connect[7] = mesh->getNodeAt(n7-1);
                newcell = JHexahedron::newObject();
                newcell->setNodes(connect);
                mesh->addObject( newcell );
            }
        }
    }
    cout << "Read file " << endl;
    return mesh;
}

