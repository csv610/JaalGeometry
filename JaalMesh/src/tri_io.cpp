#include "MeshImporter.hpp"
#include "MeshExporter.hpp"
#include "MeshTopology.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshTRIImporter ::readPoly(const string &filename)
{
    JMeshPtr mesh;
    size_t pos = filename.rfind(".poly");
    if( pos ==  string::npos) {
        cout << "Warning: Non-Standard poly file not read: " <<  filename << endl;
        return mesh;
    }

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return mesh;
    }

    int id, numnodes, ndim, numattrib, boundflag, bmark;

    infile >> numnodes >> ndim >> numattrib >> boundflag;
    if( numnodes < 1) return mesh;

    mesh = JMesh::newObject();
    mesh->reserve(numnodes, 0);

    JNodeSequence newnodes = JNode::newObjects(numnodes);

    start_index = std::numeric_limits<int>::max();

    bmark = 1;
    double x, y, z = 0.0;
    Point3D xyz;
    for( int i = 0; i < numnodes; i++)  {
        infile >> id >> x >> y ;
        if( ndim == 3) infile >> z;
        if( boundflag ) infile >> bmark;
        global2local[id] = i;
        start_index = min(start_index, id);
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        newnodes[i]->setID(i);
        newnodes[i]->setXYZCoords(xyz);
        newnodes[i]->setAttribute("Boundary", bmark);
    }
    mesh->addObjects( newnodes);
    mesh->enumerate(0);

    int numSegments;
    infile >> numSegments >> boundflag;

    JEdgeSequence newedges = JEdge::newObjects(numSegments);

    int vid1, vid2;
    for( int i = 0; i < numSegments; i++) {
        infile >> id >> vid1 >> vid2;
        if( boundflag) infile >> bmark;
        const JNodePtr &v0 = mesh->getNodeAt(vid1-start_index);
        const JNodePtr &v1 = mesh->getNodeAt(vid2-start_index);
        newedges[i]->setNodes(v0,v1);
    }

    mesh->addObjects( newedges);
    mesh->enumerate(1);

    int numholes;
    infile >> numholes;
    vector<Point3D> holes;
    if( numholes > 0) {
        holes.resize( numholes);
        for( int i = 0; i < numholes; i++) {
            infile >> x >> y ;
            if( ndim == 3) infile >> z;
            xyz[0] = x;
            xyz[1] = y;
            xyz[2] = z;
            holes[i] = xyz;
        }
        mesh->setAttribute("HolesCoord", holes);
    }
    return mesh;
}
/////////////////////////////////////////////////////////////////////////////////

int JMeshTRIImporter ::readNodes( const string &filename, const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;

    size_t pos = filename.rfind(".node");
    if( pos ==  string::npos) {
        cout << "Warning: Non-Standard node file not read: " <<  filename << endl;
        return 1;
    }

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    int id, numnodes, ndim, numattrib, boundflag, bmark;
    start_index = std::numeric_limits<int>::max();

    infile >> numnodes >> ndim >> numattrib >> boundflag;

    if( numnodes < 1) return 1;

    mesh->reserve(numnodes, 0);


    JNodeSequence newnodes = JNode::newObjects(numnodes);

    double x, y, z = 0.0;
    Point3D xyz;
    for( int i = 0; i < numnodes; i++)  {
        infile >> id >> x >> y ;
        if( ndim == 3) infile >> z;
        if( boundflag ) {
            infile >> bmark;
            if( bmark ) newnodes[i]->setAttribute("Boundary", bmark);
        }
        global2local[id] = i;
        start_index = min(id, start_index);
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        newnodes[i]->setID(i);
        newnodes[i]->setXYZCoords(xyz);
    }
    mesh->addObjects( newnodes);
    mesh->enumerate(0);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int JMeshTRIImporter::readEdges( const string &filename, const JMeshPtr &mesh)
{
    if( mesh == nullptr) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() ) {
        return 1;
    }
    int id, numedges, boundflag;
    int n0, n1, bmark;

    infile >> numedges >> boundflag;
    if( numedges < 1) return 1;
    JEdgeSequence newedges = JEdge::newObjects(numedges );

    for( int i = 0; i < numedges; i++) {
        infile >> id >> n0 >> n1;
        n0 = global2local[n0];
        n1 = global2local[n1];
        const JNodePtr &v0 = mesh->getNodeAt(n0);
        const JNodePtr &v1 = mesh->getNodeAt(n1);
        JEdgePtr edge = JSimplex::getEdgeOf(v0, v1, 1);
        assert(edge);
        if( boundflag ) {
            infile >> bmark;
            if( bmark ) edge->setAttribute("Boundary", bmark);
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTRIImporter ::readFaces( const string &filename, const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() ) {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    int id, numfaces, numelemnodes, boundflag, bmark;
    infile >> numfaces >> numelemnodes >> boundflag;
    if( numfaces < 1) return 1;

    assert(numelemnodes == 3 );

    JNodeSequence connect(3);
    JFaceSequence newfaces = JTriangle::newObjects(numfaces);

    int n0, n1, n2;
    for( int i = 0; i < numfaces; i++) {
        infile >> id >> n0 >> n1 >> n2;
        if( boundflag )  infile >> bmark;
        n0 = global2local[n0];
        n1 = global2local[n1];
        n2 = global2local[n2];
        connect[0] = mesh->getNodeAt(n0);
        connect[1] = mesh->getNodeAt(n1);
        connect[2] = mesh->getNodeAt(n2);
        newfaces[i]->setNodes( connect );
    }
    mesh->enumerate(2);
    mesh->addObjects( newfaces );
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshTRIImporter ::readCells( const string &filename, const JMeshPtr &mesh)
{
    if( mesh == nullptr) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() ) {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    int id, numtets, numelemnodes, boundflag, bmark;
    infile >> numtets >> numelemnodes >> boundflag;
    if( numelemnodes == 3) {
        infile.close();
        readFaces( filename, mesh);
        return 0;
    }
    if( numtets < 1) return 1;

    assert(numelemnodes == 4 );
    JCellSequence newtets = JTetrahedron::newObjects(numtets);

    JNodeSequence connect(4);
    int n0, n1, n2, n3;
    for( int i = 0; i < numtets; i++)  {
        infile >> id >> n0 >> n1 >> n2 >> n3;
        if( boundflag )  infile >> bmark;
        n0 = global2local[n0];
        n1 = global2local[n1];
        n2 = global2local[n2];
        n3 = global2local[n3];
        connect[0] = mesh->getNodeAt(n0);
        connect[1] = mesh->getNodeAt(n1);
        connect[2] = mesh->getNodeAt(n2);
        connect[3] = mesh->getNodeAt(n3);
        newtets[i]->setNodes(connect);
    }
    mesh->enumerate(3);
    mesh->addObjects( newtets );
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshTRIImporter ::readFile( const string &fname)
{
    size_t pos = fname.rfind(".poly");
    if( pos != string::npos) {
        return readPoly(fname);
    }

    JMeshPtr mesh = JMesh::newObject();


    string basefile;
    pos = fname.rfind(".node");
    if( pos != string::npos) basefile  = fname.substr(0, pos);

    if( basefile.empty() )  {
        pos = fname.rfind(".ele");
        if( pos != string::npos) basefile  = fname.substr(0, pos);
    }

    if( !basefile.empty() ) {
        readNodes( basefile  + ".node", mesh);
        readEdges( basefile  + ".edge", mesh);
        readCells( basefile  + ".ele",  mesh);
        global2local.clear();
        return mesh;
    }

    cout << "Warning: Not a valid file suffix for triangle file " << endl;
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTRIExporter ::writeNodes(const JMeshPtr &mesh, const string &fname)
{
    if( mesh == nullptr) return 1;

    string filename = fname;
    size_t pos = fname.rfind(".node");
    if( pos ==  string::npos) filename += ".node";

    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    size_t numnodes = mesh->getSize(0);
    int numattrib  = 0;
    int boundflag  = 0;
//    int bmark;

    ofile  << numnodes <<  "  " << dim  << " " <<  numattrib << " " << boundflag << endl;

    for( size_t i = 0; i < numnodes; i++)  {
        JNodePtr v = mesh->getNodeAt(i);
        const Point3D &xyz = v->getXYZCoords();
        ofile << i << " ";
        ofile << fixed << setprecision(precision) << scientific;
        ofile <<  xyz[0] << " " << xyz[1];
        if( dim == 3 ) ofile << " " << xyz[2];
        ofile << endl;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTRIExporter ::writeFaces(const JMeshPtr &mesh, const string &filename)
{
//      string filename = fname + ".ele";
    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    size_t numfaces = mesh->getSize(2);
    int numelemnodes  = 3;
    int boundflag  = 0;

    ofile << numfaces <<  " " << numelemnodes << " " <<  boundflag << endl;

    for( size_t id = 0; id < numfaces; id++)  {
        const JFacePtr &face = mesh->getFaceAt(id);
        if( face->isActive() ) {
            int n0 = face->getNodeAt(0)->getID();
            int n1 = face->getNodeAt(1)->getID();
            int n2 = face->getNodeAt(2)->getID();
            ofile <<  id << " " << n0 << " " << n1 << " " << n2 << endl;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTRIExporter ::writeCells(const JMeshPtr &mesh, const string &filename)
{
//      string filename = fname + ".ele";
    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    size_t numcells = mesh->getSize(3);
    int numelemnodes  = 4;
    int boundflag  = 0;

    ofile << numcells <<  " " << numelemnodes << " " <<  boundflag << endl;

    for( size_t id = 0; id < numcells; id++)  {
        const JCellPtr &cell = mesh->getCellAt(id);
        if( cell->isActive() ) {
            int n0 = cell->getNodeAt(0)->getID();
            int n1 = cell->getNodeAt(1)->getID();
            int n2 = cell->getNodeAt(2)->getID();
            int n3 = cell->getNodeAt(3)->getID();
            ofile <<  id << " " << n0 << " " << n1 << " " << n2 << " " << n3 << endl;
        }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int
JMeshTRIExporter ::writeFacets(const JMeshPtr &mesh, const string &s)
{
    // Ref: http://tetgen.berlios.de/fformats.smesh.html
    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) {
        logger->setWarn("Unable to open file for writing mesh");
        return 1;
    }

    int nattrib = 0, bmarker = 0, ndim = 3;
    size_t numnodes = mesh->getSize(0);

    ofile <<  numnodes << " " << ndim << " " << nattrib << " " << bmarker << endl;

    ofile << setprecision(10);
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        assert( v->isActive() );
        const Point3D &p3d = v->getXYZCoords();
        ofile << i << " " << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
    }

    size_t numfaces = mesh->getSize(2);

    ofile << numfaces << " " << bmarker << endl;

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
            ofile << nnodes << " ";
            for (int j = 0; j < nnodes; j++) {
                size_t vid = newConnect[j]->getID();
                if (vid >= numnodes) {
                    assert(!newConnect[j]->isRemoved());
                    logger->setError("Vertex indexing out of range " );
                    exit(0);
                }
                ofile << vid << " ";
            }
            ofile << endl;
        }
    }

    int numholes = 0;
    ofile << numholes << endl;

    int numregions = 0;
    ofile << numregions << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////

int JMeshTRIExporter ::writeFile(const JMeshPtr &mesh, const string &fname)
{
    string basefile = fname;
    size_t pos = fname.rfind(".node");
    if( pos != string::npos) basefile  = fname.substr(0, pos);

    if( basefile.empty() )  {
        pos = fname.rfind(".ele");
        if( pos != string::npos) basefile  = fname.substr(0, pos);
    }

    if( basefile.empty() ) {
        cout << "Warning: basefile is empty " << endl;
        return 1;
    }

    int topDim = mesh->getTopology()->getDimension();
    writeNodes(mesh, basefile + ".node");

    if( topDim == 2 ) writeFaces(mesh, basefile + ".ele");

    if( topDim == 3 ) {
        writeFacets(mesh, basefile + ".smesh");
        writeCells(mesh, basefile + ".ele");
    }

    return 0;

}
///////////////////////////////////////////////////////////////////////////////
