#include "MeshImporter.hpp"

using namespace Jaal;

JLogger* JMeshImporter::logger = JLogger::getInstance();

JMeshImporterPtr JMeshImporter :: getProduct(const string &fname)
{
    JMeshImporterPtr importer;

    if (fname.rfind(".ele") != string::npos)
        importer.reset(new JMeshTRIImporter);

    if (fname.rfind(".node") != string::npos)
        importer.reset(new JMeshTRIImporter);

    if (fname.rfind(".poly") != string::npos)
        importer.reset(new JMeshTRIImporter);

    if (fname.rfind(".vtk") != string::npos)
        importer.reset(new JMeshVTKImporter);

    if (fname.rfind(".off") != string::npos)
        importer.reset(new JMeshOFFImporter);

    if (fname.rfind(".obj") != string::npos)
        importer.reset(new JMeshOBJImporter);

    if (fname.rfind(".smf") != string::npos)
        importer.reset(new JMeshOBJImporter);

    if (fname.rfind(".xml") != string::npos)
        importer.reset(new JMeshXMLImporter);

    if (fname.rfind(".msh") != string::npos)
        importer.reset(new JMeshGmshImporter);

    if (fname.rfind(".mesh") != string::npos)
        importer.reset(new JMeshMeditImporter);

    return importer;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshIO :: readFile(const string &fname)
{
    JMeshPtr newmesh;
    JMeshImporterPtr imp = JMeshImporter::getProduct(fname);
    cout << "file name " << fname << endl;

    if( imp == nullptr) {
        cout << "Warning: No mesh importer found for the extension " << endl;
        return nullptr;
    }

    newmesh = imp->readFile(fname);

    if( newmesh == nullptr) {
        cout << "Warning: The mesh file not read properly " << endl;
        return nullptr;
    }

    newmesh->setFileName(fname);

    return newmesh;
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr read_mesh_file( const string &fname)
{
    JMeshPtr mesh;
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << fname << endl;
        return mesh;
    }
    mesh = JMesh::newObject();

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
    return mesh;
}

