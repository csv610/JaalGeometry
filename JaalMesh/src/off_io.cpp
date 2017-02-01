#include "MeshImporter.hpp"
#include "MeshExporter.hpp"

using namespace Jaal;

//##############################################################################
void skip_comments(FILE *f)
{
    // Skip comments in an ASCII file (lines beginning with #)
    int c;
    bool in_comment = false;
    while (1) {
        c = fgetc(f);
        if (c == EOF) return;
        if (in_comment) {
            if (c == '\n') in_comment = false;
        } else if (c == '#') {
            in_comment = true;
        } else if (!isspace(c)) {
            break;
        }
    }
    ungetc(c, f);
}
//##############################################################################

JMeshPtr JMeshOFFImporter ::readFile(const string &fname)
{
    JMeshPtr mesh;
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() ) {
//     logger->setWarn("Can not open Off File");
        cout << "Warning: Can not open file " << fname << endl;
        return mesh;
    }

    string str;
    infile >> str;

    bool valid_header = 0;

    if( str == "OFF") valid_header = 1;
    if( str == "COFF") valid_header = 1;

    if( !valid_header ) {
        cout << "Warning: A valid Off format must start with OFF as first string " << endl;
        return mesh;
    }
    mesh = JMesh::newObject();

    logger->setInfo("Reading Off File ");

    //  The codelet is borrowed from TriMesh Software
    vector<int> facevtx;
    double  x, y, z;

    JNodeSequence vnodes, connect(3);

    size_t  numNodes, numFaces, numEdges;
    infile >> numNodes >> numFaces >> numEdges;

    mesh->reserve( numNodes, 0);
    mesh->reserve( numFaces, 2);
    string line;

    JNodeSequence newnodes = JNode::newObjects(numNodes);

    Point3D p3d;
    for( size_t i = 0; i < numNodes; i++) {
        infile >> x >> y >> z;
        getline(infile, line);
        p3d[0] = x;
        p3d[1] = y;
        p3d[2] = z;
        newnodes[i]->setXYZCoords(p3d);
        newnodes[i]->setID(i);
    }
    mesh->addObjects( newnodes );

    vector<size_t> nodeConnect;
    size_t nodeid, nodeindex;

    if( numFaces ) {
        JFaceSequence triangles, quads, polygons;
        size_t numTriangles = 0, numQuads = 0, numPolygons = 0;
        vector<short>  elemType(numFaces);
        nodeConnect.reserve(3*numNodes);
        for( size_t i = 0; i < numFaces; i++) {
            infile >> elemType[i];
            switch( elemType[i] ) {
            case 3:
                numTriangles++;
                break;
            case 4:
                numQuads++;
                break;
            default:
                numPolygons++;
                break;
            }
            for( int j = 0; j < elemType[i]; j++) {
                infile >> nodeid;
                nodeConnect.push_back(nodeid);
            }
            getline(infile, str);  // Skipping everything after the connectivity..
        }
        if( numTriangles ) triangles = JTriangle::newObjects(numTriangles);
        if( numQuads )     quads     = JQuadrilateral::newObjects(numQuads);
        if( numPolygons )  polygons  = JPolygon::newObjects(numPolygons);

        size_t triindex  = 0;
        size_t quadindex = 0;
        size_t polyindex = 0;
        nodeindex = 0;
        for( size_t i = 0; i < numFaces; i++) {
            switch( elemType[i] )
            {
            case 3:
                connect.resize(3);
                connect[0] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                connect[1] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                connect[2] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                triangles[triindex++]->setNodes( connect );
                break;
            case 4:
                connect.resize(4);
                connect[0] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                connect[1] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                connect[2] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                connect[3] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                quads[quadindex++]->setNodes( connect );
                break;
            default:
                connect.resize(elemType[i]);
                for( int j = 0; j < elemType[i]; j++)
                    connect[j] = mesh->getNodeAt( nodeConnect[nodeindex++]);
                polygons[polyindex++]->setNodes( connect );
                break;
            }
        }
        assert( triindex  == numTriangles );
        assert( quadindex == numQuads);
        assert( polyindex == numPolygons );
        if( numTriangles ) mesh->addObjects( triangles);
        if( numQuads )     mesh->addObjects( quads );
        if( numPolygons )  mesh->addObjects( polygons);
    }

    JEdgeSequence edges;

    // Edges are created by the faces, but if they are present in the file,
    // we can change their ID...
    if( numEdges && numFaces == 0)   {
        nodeConnect.clear();
        int dummy;
        for( size_t i = 0; i < numEdges; i++) {
            infile >> dummy;
            for( int j = 0; j < 2; j++) {
                infile >> nodeid;
                nodeConnect.push_back(nodeid);
            }
            getline(infile, str);  // Skipping everything after the connectivity..
        }
//        size_t edgeindex = 0;
        nodeindex = 0;
        connect.resize(2);
        JEdgePtr  edge;
        for( size_t i = 0; i < numEdges; i++) {
            const JNodePtr &v0 = mesh->getNodeAt( nodeConnect[2*i] );
            const JNodePtr &v1 = mesh->getNodeAt( nodeConnect[2*i+1] );
            edge = JSimplex::getEdgeOf(v0,v1);
            if( edge == nullptr) {
                edge = JEdge::newObject(v0,v1);
                mesh->addObject(edge);
            }
            edge->setID(i);
        }
    }

    logger->setInfo("Reading mesh file complete");

    return mesh;
}

//##############################################################################

int
JMeshOFFExporter ::writeFile(const JMeshPtr &mesh, const string &s)
{
    if( mesh == nullptr) return 1;

    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);
    size_t numfaces = mesh->getSize(2);

    logger->setError("Saving mesh file in off format");
    mesh->pruneAll();

    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() )
        return 1;

    // We do not want to write both faces and edges ...
    if( numfaces  ) numedges = 0;

    ofile << "OFF" << endl;
    ofile << numnodes << " " << numfaces <<  "  " << numedges  << endl;

    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        assert( v->isActive() );
        const Point3D &p3d = v->getXYZCoords();
        ofile << fixed << setprecision(precision );
        ofile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
        v->setID(i);
    }

    JNodeSequence oldConnect, newConnect;
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive()) {
            if (face->getSize(0) == 4) {
                oldConnect = face->getNodes();
                JQuadrilateral::quad_tessalate(oldConnect, newConnect); // Because of OpenGL
            } else
                newConnect = face->getNodes();

            int nnodes = newConnect.size();
            assert( nnodes);
            ofile << nnodes << " ";
            for (int j = 0; j < nnodes; j++) {
                size_t vid = newConnect[j]->getID();
                if (vid >= numnodes) {
                    assert(!newConnect[j]->isRemoved());
                    cout << face->getStatus() << endl;
                    cout << "Node indexing out of range " << vid << " Total : " << numnodes << endl;
                    exit(0);
                }
                ofile << vid << " ";
            }
            ofile << endl;
        }
    }

    size_t numCells = mesh->getSize(3);
    if( numCells ) {
        ofstream volfile;
        string basefile;
        size_t pos1 = filename.rfind(".off");
        if( pos1 != string::npos) basefile  = filename.substr(0, pos1);
        basefile += "_volmesh.off";
        volfile.open(basefile.c_str(), ios::out);

        volfile << "OFF" << endl;
        volfile << numnodes << " " << numCells << " 0  " << endl;

        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            assert( v->isActive() );
            const Point3D &p3d = v->getXYZCoords();
            volfile << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
        }

        for (size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive()) {
                int nnodes = cell->getSize(0);
                volfile << nnodes << " ";
                for (int j = 0; j < nnodes; j++) {
                    volfile << cell->getNodeAt(j)->getID() << " ";
                }
                volfile << endl;
            }
        }
    }

    for (size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive()) {
            ofile << "2  " << edge->getNodeAt(0)->getID() << "  " << edge->getNodeAt(1)->getID() << endl;
        }
    }


    return 0;
}
