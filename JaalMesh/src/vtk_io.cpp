#include "MeshExporter.hpp"
#include "MeshImporter.hpp"

using namespace Jaal;

int JMeshVTKExporter :: writeFile(const JMeshPtr &mesh, const string &fname)
{
    ofstream ofile(fname.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    mesh->pruneAll();

    size_t numNodes = mesh->getSize(0);

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << " Jaal" << endl;

    ofile << "ASCII " << endl;

    ofile << "DATASET UNSTRUCTURED_GRID " << endl;
    ofile << "POINTS " << numNodes << " float " << endl;
    ofile << std::setiosflags(ios::fixed);

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        const Point3D &xyz = vertex->getXYZCoords();
        ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
    }
    int nCount1 = 0, nCount2 = 0;

    size_t numEdges = mesh->getSize(1);

    numEdges = 0;
    nCount1 = numEdges;
    nCount2 = 3*numEdges;

    size_t numFaces = mesh->getSize(2);
    nCount1 += numFaces;

    size_t numCells = mesh->getSize(3);
    nCount1 += numCells;

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        nCount2 += face->getSize(0) + 1;
    }

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        nCount2 += cell->getSize(0) + 1;
    }

    ofile << "CELLS " << nCount1 << "  " << nCount2 << endl;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        ofile << "2  ";
        for(int i = 0; i < 2; i++) {
            const JNodePtr &v = edge->getNodeAt(i);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        ofile << nsize << " ";
        for(int i = 0; i < nsize; i++) {
            JNodePtr v = face->getNodeAt(i);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        int nsize = cell->getSize(0);
        ofile << nsize << " ";
        for(int i = 0; i < nsize; i++) {
            JNodePtr v = cell->getNodeAt(i);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }

    ofile << "CELL_TYPES " << nCount1 << endl;
    for( size_t i = 0; i < numEdges; i++)
        ofile << "3 ";
    ofile << endl;

    int type = -1;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int nsize = face->getSize(0);
        type  = 7;     // Polygon
        if( nsize == 3 ) type = 5;  // Triangle
        if( nsize == 4 ) type = 9;  // Quad
        ofile << type << " ";
    }
    ofile << endl;

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        int nsize = cell->getSize(0);
        if( nsize == 4 ) type = 10;   // Tet
        if( nsize == 8 ) type = 12;   // Hex
        ofile << type << " ";
    }
    ofile << endl;

    double val;
    for( size_t i = 0; i < nodeAttribs.size(); i++) {
        string name = nodeAttribs[i];
        ofile << "POINT_DATA " << numNodes << endl;
        ofile << "SCALARS " << name << " float 1 " << endl;
        ofile << "LOOKUP_TABLE  default" << endl;
        for( size_t j = 0; j <  numNodes; j++) {
            const JNodePtr &v = mesh->getNodeAt(j);
            v->getAttribute(name, val);
            ofile << val << endl;
        }
    }

    for( size_t i = 0; i < faceAttribs.size(); i++) {
        string name = faceAttribs[i];
        ofile << "CELL_DATA " << numFaces << endl;
        ofile << "SCALARS " << name << " float 1 " << endl;
        ofile << "LOOKUP_TABLE  default" << endl;
        for( size_t j = 0; j <  numFaces; j++) {
            const JFacePtr &f = mesh->getFaceAt(j);
            f->getAttribute(name, val);
            ofile << val << endl;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshVTKImporter ::readFile( const string &fname)
{
    JMeshPtr mesh;
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << fname << endl;
        return mesh;
    }

    mesh = JMesh::newObject();

    size_t numnodes, numelems;
    double x, y, z = 0.0;
    Point3D xyz;

    int nv, numelemnodes, nsize;
    JNodeSequence connect;
    JEdgePtr newedge;
    JFacePtr newface;
    JCellPtr newcell;
    vector<int> vconn, elmtype;

    string str;
    while( !infile.eof() ) {
        infile >> str;

        if( str == "POINTS") {
            infile >> numnodes >> str;
            mesh->reserve( numnodes, 0);
            for( size_t i = 0; i < numnodes; i++)  {
                infile >> x >> y  >> z;
                xyz[0] = x;
                xyz[1] = y;
                xyz[2] = z;
                JNodePtr v = JNode::newObject();
                v->setID(i);
                v->setXYZCoords(xyz);
                mesh->addObject(v);
            }
        }

        if( str == "CELLS" || str == "POLYGONS") {
            infile >> numelems >> nsize;
            vconn.reserve(nsize);
            for( size_t i = 0; i < numelems; i++) {
                infile >> numelemnodes;
                vconn.push_back(numelemnodes);
                for( int j = 0; j < numelemnodes; j++) {
                    infile >> nv;
                    vconn.push_back(nv);
                }
            }
        }

        if( str == "CELL_TYPES") {
            infile >> nsize;
            elmtype.resize(nsize);
            for( size_t i = 0; i < numelems; i++)
                infile >> elmtype[i];
        }
    }

    size_t offset = 0;
    for( size_t i = 0; i < numelems; i++) {
        switch( elmtype[i] ) {
        case 3:
            assert( vconn[offset] == 2 );
            connect.resize(2);
            connect[0] = mesh->getNodeAt(vconn[offset+1]);
            connect[1] = mesh->getNodeAt(vconn[offset+2]);
            newedge    = JEdge::newObject( connect[0], connect[1] );
            mesh->addObject( newedge );
            offset += 3;
            break;
        case 5:
            assert( vconn[offset] == 3 );
            connect.resize(3);
            connect[0] = mesh->getNodeAt(vconn[offset+1]);
            connect[1] = mesh->getNodeAt(vconn[offset+2]);
            connect[2] = mesh->getNodeAt(vconn[offset+3]);
            newface = JTriangle::newObject( connect );
            mesh->addObject( newface );
            offset += 4;
            break;
        case 9:
            assert( vconn[offset] == 4 );
            connect.resize(4);
            connect[0] = mesh->getNodeAt(vconn[offset+1]);
            connect[1] = mesh->getNodeAt(vconn[offset+2]);
            connect[2] = mesh->getNodeAt(vconn[offset+3]);
            connect[3] = mesh->getNodeAt(vconn[offset+4]);
            newface = JQuadrilateral::newObject(connect);
            mesh->addObject( newface );
            offset += 5;
            break;
        case 10:
            assert( vconn[offset] == 4 );
            connect.resize(4);
            connect[0] = mesh->getNodeAt(vconn[offset+1]);
            connect[1] = mesh->getNodeAt(vconn[offset+2]);
            connect[2] = mesh->getNodeAt(vconn[offset+3]);
            connect[3] = mesh->getNodeAt(vconn[offset+4]);
            newcell = JTetrahedron::newObject(connect);
            mesh->addObject( newcell );
            offset += 5;
            break;
        case 12:
            assert( vconn[offset] == 8 );
            connect.resize(8);
            connect[0] = mesh->getNodeAt(vconn[offset+1]);
            connect[1] = mesh->getNodeAt(vconn[offset+2]);
            connect[2] = mesh->getNodeAt(vconn[offset+3]);
            connect[3] = mesh->getNodeAt(vconn[offset+4]);
            connect[4] = mesh->getNodeAt(vconn[offset+5]);
            connect[5] = mesh->getNodeAt(vconn[offset+6]);
            connect[6] = mesh->getNodeAt(vconn[offset+7]);
            connect[7] = mesh->getNodeAt(vconn[offset+8]);
            newcell = JHexahedron::newObject(connect);
            mesh->addObject( newcell );
            offset += 9;
            break;
        default:
            cout << "At present this cell not implemented" << endl;
            exit(0);
        }
    }

    return mesh;
}
