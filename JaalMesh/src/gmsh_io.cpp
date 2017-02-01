#include "MeshImporter.hpp"

using namespace Jaal;

//##############################################################################

JMeshPtr JMeshGmshImporter ::readFile(const string &fname)
{
    JMeshPtr mesh;
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() ) {
        cout << "Warning: Cann't open file " << fname << endl;
        return mesh;
    }
    mesh = JMesh::newObject();

    string str;
    size_t  numNodes, numElems;

    infile >> str;
    assert( str == "$MeshFormat");
    int  file_type, data_size, id;
    float version;
    infile >> version >> file_type >> data_size;
    infile >> str;
    assert( str == "$EndMeshFormat");

    infile >> str;
    assert( str == "$Nodes");

    infile >> numNodes;

    mesh->reserve( numNodes, 0);
    Point3D p3d;
    double  x, y, z;
    JNodePtr vertex;

    for( size_t i = 0; i < numNodes; i++) {
        infile >> id >> x >> y >> z;
        p3d[0] = x;
        p3d[1] = y;
        p3d[2] = z;
        vertex = JNode::newObject();
        vertex->setXYZCoords(p3d);
        vertex->setID(i);
        mesh->addObject( vertex );
    }

    infile >> str;
    assert( str == "$EndNodes");
    infile >> str;
    assert( str == "$Elements");

    int elm_number, elm_type, numTags, tagval;
    vector<int> facevtx;
    JNodeSequence vnodes, connect(3);
    string line;

    infile >> numElems;
    mesh->reserve( numElems, 2);

    JFacePtr newface;
    for( size_t i = 0; i < numElems; i++) {
        infile >> elm_number >> elm_type >> numTags;
        for( int j = 0; j < numTags; j++)
            infile >> tagval;

        int valid_elem = 0;
        switch( elm_type ) {
        case 2:
            facevtx.resize(3);
            for( size_t j = 0; j < 3; j++)
                infile >> facevtx[j];
            connect.resize(3);
            connect[0] = mesh->getNodeAt(facevtx[0]-1);
            connect[1] = mesh->getNodeAt(facevtx[1]-1);
            connect[2] = mesh->getNodeAt(facevtx[2]-1);
            newface = JTriangle::newObject(connect);
            valid_elem = 1;
            break;
        case 3:
            facevtx.resize(4);
            for( size_t j = 0; j < 4; j++)
                infile >> facevtx[j];
            connect.resize(4);
            connect[0] = mesh->getNodeAt(facevtx[0]-1);
            connect[1] = mesh->getNodeAt(facevtx[1]-1);
            connect[2] = mesh->getNodeAt(facevtx[2]-1);
            connect[3] = mesh->getNodeAt(facevtx[3]-1);
            newface = JQuadrilateral::newObject(connect);
            valid_elem = 1;
            break;
        default:
            getline( infile, line );
        }

        if( valid_elem ) mesh->addObject(newface);
    }

    return mesh;
}

//##############################################################################
int JMeshGmshExporter ::writeFile( const JMeshPtr &mesh, const string &fname)
{
    if( mesh == nullptr) return 1;

    ofstream ofile( fname.c_str(), ios::out);

    JNodeSequence nodes;
    mesh->getTopology()->getBoundary( nodes );

    map<JNodePtr, int> vmap;
    int index = 1;
    ofile << "lc = 0.1;" << endl;

    for( const JNodePtr &vtx: nodes ) {
        vmap[vtx] = index;
        const Point3D &xyz = vtx->getXYZCoords();
        ofile << "Point(" << index << ")= {" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ",lc};" << endl;
        index++;
    }

    vector<JEdgeSequence> edges;
    mesh->getTopology()->getBoundary( edges );

    vector<int> loopindex;
    for( size_t i = 0; i < edges.size(); i++) {
        JEdgeTopology::getChain( edges[i] );
        int lineindex = index;
        for( const JEdgePtr &e: edges[i] ) {
            int v1 = vmap[e->getNodeAt(0)];
            int v2 = vmap[e->getNodeAt(1)];
            ofile << "Line(" << index++ << ") = {" << v1 << "," << v2 << "};" << endl;
        }
        loopindex.push_back(index);
        ofile << "Line Loop(" << index++ <<") = {";
        for( size_t j = 0; j < edges[i].size(); j++) {
            ofile << lineindex+j;
            if( j < edges[i].size()-1) ofile << ",";
        }
        ofile << "};" << endl;
    }
    ofile << "Plane Surface(" << index++ << ") = {";
    for( size_t i = 0; i < edges.size(); i++)  {
        ofile << loopindex[i];
        if( i < edges.size()-1) ofile << ",";
    }
    ofile << "};" << endl;

    return 0;


}

