#include "MeshExporter.hpp"
#include "MeshImporter.hpp"

using namespace Jaal;

int JMeshXMLImporter :: readAttributeHeader( ifstream &infile, string &name, string &type, size_t &ncount)
{
    string str;
    infile >> str;

    if( str != "<Name>") {
        cout <<"Warning: attribute name not specified at line # : " << endl;
        return 1;
    }

    infile >> name;
    infile >> str;
    if( str != "</Name>") {
        cout << "Warning: name value not closed at line # : " << __LINE__ << endl;
        return 1;
    }

    infile >> str;
    if( str != "<Datatype>" )  {
        cout << "Warning: datatype not specified at line # : " << __LINE__ << endl;
        return 1;
    }
    infile >> type;
    infile >> str;
    if( str != "</Datatype>" )  {
        cout << "Warning: datatype not closed at line # :  " << __LINE__ << endl;
        return 1;
    }

    infile >> str;
    if( str != "<Count>" )  {
        cout << "Warning: count not specified at line # : " << __LINE__ << endl;
        return 1;
    }
    infile >> ncount;
    infile >> str;
    if( str != "</Count>" )  {
        cout << "Warning: count value not closed at line # : " << __LINE__ << endl;
        return 1;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int JMeshXMLImporter :: readNodes( const string &filename, const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return 2;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    string str;
    bool hasnodes = 0;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "<Nodes>" ) {
            hasnodes = 1;
            break;
        }
    }
    if( !hasnodes ) return 1;

    int numNodes;

    infile >> str;
    assert( str == "<Count>" );

    infile >> numNodes;
    infile >> str;
    assert( str == "</Count>" );

    JNodeSequence newnodes = JNode::newObjects(numNodes);

    double x, y, z;
    Point3D p3d;
    infile >> str;
    assert(str == "<NodeCoordinates>" );
    for( int i = 0; i < numNodes; i++) {
        infile >> x >> y >> z;
        p3d[0] = x;
        p3d[1] = y;
        p3d[2] = z;
        newnodes[i]->setXYZCoords(p3d);
        newnodes[i]->setID(i);
    }
    mesh->addObjects( newnodes );
    infile >> str;
    assert(str == "</NodeCoordinates>" );

    int err;
    string name;
    string type;
    size_t ncount;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "</Nodes>" ) break;
        if( str == "<NodeAttribute>") {
            err = readAttributeHeader( infile, name, type, ncount );
            if( !err ) readNodeAttribute( name, type, ncount, infile);
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshXMLImporter :: readEdges(const string &filename, const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    string str;
    bool hasedges = 0;
    while(!infile.eof() ) {
        infile >> str;
        if( str == "<Edges>" ) {
            hasedges = 1;
            break;
        }
    }
    if( !hasedges) return 1;

    infile >> str;
    if( str != "<Count>" ) {
        cout << "Warning:  <Count> missing in edge field " << endl;
        return 1;
    }

    size_t numedges;
    infile >> numedges;
    infile >> str;
    if( str != "</Count>" ) {
        cout << "Warning:  </Count> missing in edge field " << endl;
        return 1;
    }

    if( numedges == 0) return 1;

    JEdgeSequence newedges;
    int n1, n2;
    for( size_t i = 0; i < numedges; i++) {
        infile >> n1 >> n2;
        const JNodePtr &v0 = mesh->getNodeAt(n1);
        const JNodePtr &v1 = mesh->getNodeAt(n2);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        if( edge == nullptr) {
            edge = JSimplex::getEdgeOf(v0,v1,1);
            newedges.push_back(edge);
        }
    }
    mesh->addObjects(newedges);

    int err;
    string name;
    string type;
    size_t ncount;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "</Edges>" ) break;
        if( str == "<EdgeAttribute>") {
            err = readAttributeHeader( infile, name, type, ncount );
            if( !err )
                readEdgeAttribute( name, type, ncount, infile);
        }
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshXMLImporter :: readFaces( const string &filename, const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    string str;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "<Faces>" ) break;
    }

    int err;
    string name;
    string type;

    JNodeSequence nodes;
    size_t n1, n2, n3, n4, ncount, nnodes, numfaces = 0;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "</Faces") break;

        if( str == "<Triangles>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            nodes.resize(3);
            JFaceSequence newtriangles = JTriangle::newObjects(ncount);
            for( size_t i = 0; i < ncount; i++) {
                infile >> n1 >> n2 >> n3;
                nodes[0] =  mesh->getNodeAt(n1);
                nodes[1] =  mesh->getNodeAt(n2);
                nodes[2] =  mesh->getNodeAt(n3);
                newtriangles[i]->setNodes(nodes);
            }
            mesh->addObjects(newtriangles);
            numfaces += ncount;
            infile >> str;
            assert( str == "</Triangles>" );
        }

        if( str == "<Quads>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            nodes.resize(4);
            JFaceSequence newquads = JQuadrilateral::newObjects(ncount);
            for( size_t i = 0; i < ncount; i++) {
                infile >> n1 >> n2 >> n3 >> n4;
                nodes[0] =  mesh->getNodeAt(n1);
                nodes[1] =  mesh->getNodeAt(n2);
                nodes[2] =  mesh->getNodeAt(n3);
                nodes[3] =  mesh->getNodeAt(n4);
                newquads[i]->setNodes(nodes);
            }
            mesh->addObjects( newquads );
            numfaces += ncount;
            infile >> str;
            assert( str == "</Quads>" );
        }

        if( str == "<Polygons>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            JFaceSequence newpolys = JPolygon::newObjects(ncount);
            for( size_t i = 0; i < ncount; i++) {
                infile >> nnodes ;
                nodes.resize(nnodes);
                for( size_t j = 0; j < nnodes; j++) {
                    infile >> n1;
                    nodes[j] =  mesh->getNodeAt(n1);
                }
                newpolys[i]->setNodes(nodes);
            }
            mesh->addObjects(newpolys);
            infile >> str;
            numfaces += ncount;
            assert( str == "</Polygons>" );
        }

        if( str == "<FaceAttribute>" ) {
            err = readAttributeHeader( infile, name, type, ncount );
            if( !err )
                readFaceAttribute( name, type, ncount, infile);
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshXMLImporter :: readCells( const string &filename, const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return 1;

    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return 1;
    }

    string str;
    bool hascells  = 0;

    while( !infile.eof() ) {
        infile >> str;
        if( str == "<Cells>" ) {
            hascells = 1;
            break;
        }
    }

    if( !hascells) return 1;

    JNodeSequence nodes;
    size_t n1, ncount, numcells = 0;
    int index = 0;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "</Cells") break;
        if( str == "<Tetrahedra>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            nodes.resize(4);
            JCellSequence newtets = JTetrahedron::newObjects(ncount);
            for( size_t i = 0; i < ncount; i++) {
                for( size_t j = 0; j < 4; j++) {
                    infile >> n1 ;
                    nodes[j] =  mesh->getNodeAt(n1);
                }
                newtets[i]->setNodes(nodes);
                newtets[i]->setID( index++);
            }
            mesh->addObjects(newtets);
            numcells += ncount;
            infile >> str;
            assert( str == "</Tetrahedra>" );
        }

        if( str == "<Hexahedra>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            nodes.resize(8);
            JCellSequence newhexs = JHexahedron::newObjects(ncount);
            for( size_t i = 0; i < ncount; i++) {
                for( size_t j = 0; j < 8; j++) {
                    infile >> n1 ;
                    nodes[j] =  mesh->getNodeAt(n1);
                }
                newhexs[i]->setNodes(nodes);
                newhexs[i]->setID( index++);
            }
            mesh->addObjects(newhexs);
            numcells += ncount;
            infile >> str;
            assert( str == "</Hexahedra>" );
        }

        if( str == "<TriangularPrism>" ) {
            infile >> str;
            assert( str == "<Count>" );
            infile >> ncount;
            infile >> str;
            assert( str == "</Count>" );
            nodes.resize(6);
            for( size_t i = 0; i < ncount; i++) {
                for( size_t j = 0; j < 6; j++) {
                    infile >> n1 ;
                    nodes[j] =  mesh->getNodeAt(n1);
                }
                JCellPtr c  =  JTriangularPrism::newObject();
                c->setNodes(nodes);
                c->setID( index++);
                mesh->addObject(c);
            }
            numcells += ncount;
            infile >> str;
            assert( str == "</TriangularPrism>" );
        }
    }

    int err;
    string name;
    string type;
    while( !infile.eof() ) {
        infile >> str;
        if( str == "</Nodes>" ) break;
        if( str == "<CellAttribute>") {
            err = readAttributeHeader( infile, name, type, ncount );
            if( !err )
                readCellAttribute( name, type, ncount, infile);
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshXMLImporter ::readFile(const string &fname)
{
    filename = fname;
    infile.open(fname.c_str(), ios::in);
    if( infile.fail() ) {
        logger->setError("Can not open file ");
        return mesh;
    }


    int topodim = -1;
    string str;
    while(!infile.eof() )  {
        infile >> str;
        if( str == "<TopoDimension>")  {
            infile >> topodim;
        }
    }

    if( topodim < 0 ) {
//      newlogger->setError( "Invalid topological dimension: File not read");
        return nullptr;
    }

    // Nodes must be the first entities to be read and simplices are from lower to higher
    JMeshPtr newmesh = JMesh::newObject();
    readNodes( filename, newmesh);
    readEdges( filename, newmesh);
    readFaces( filename, newmesh);
    readCells( filename, newmesh);
    return mesh;
}

//##############################################################################
int
JMeshXMLExporter ::writeNodes(ofstream &ofile)
{
    size_t numnodes = mesh->getSize(0);
    if( numnodes == 0) return 1;

    ofile << "<Nodes>" << endl;
    ofile << "<Count> " <<  mesh->getActiveSize(0) << " </Count>" << endl;
    ofile << "<NodeCoordinates>" << endl;

    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            const Point3D &p3d = v->getXYZCoords();
            ofile << fixed << p3d[0] << " " << p3d[1] << " " << p3d[2] << endl;
        }
    }
    ofile << "</NodeCoordinates>" << endl;

    vector<string> attribnames;
    mesh->getAttributeNames(attribnames, 0);

    for( size_t i = 0; i < attribnames.size(); i++) {
        string name = attribnames[i];
        string type = JNode::getAttributeTypeName( attribnames[i] );
        if( !type.empty() ) {
            if( type == "int") {
                int  ival = 0;
                writeNodeAttribute(name, ival, ofile);
            }
            if( type == "char") {
                char cval = '0';
                writeNodeAttribute(name, cval, ofile);
            }
            if( type == "uchar") {
                unsigned char ucval = '0';
                writeNodeAttribute(name, ucval, ofile);
            }
            if( type == "double") {
                double dval = '0';
                writeNodeAttribute(name, dval, ofile);
            }
            if( type == "float") {
                double fval = '0';
                writeNodeAttribute(name, fval, ofile);
            }
        }
    }
    ofile << "</Nodes>" << endl << endl;
    return 0;
}
//##############################################################################

int
JMeshXMLExporter ::writeEdges(ofstream &ofile)
{
    size_t numedges = mesh->getSize(1);
    ofile << "<Edges>" << endl;
    ofile << "<Count> " <<  numedges << " </Count>" << endl;
    for (size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            ofile << edge->getNodeAt(0)->getID() << " "
                  << edge->getNodeAt(1)->getID() << endl;
        }
    }

    vector<string> attribnames;
    mesh->getAttributeNames(attribnames, 1);
    for( size_t i = 0; i < attribnames.size(); i++) {
        string name = attribnames[i];
        string type = JEdge::getAttributeTypeName( attribnames[i] );
        if( !type.empty() ) {
            if( type == "int") {
                int  ival = 0;
                writeEdgeAttribute(name, ival, ofile);
            }
            if( type == "char") {
                char cval = '0';
                writeEdgeAttribute(name, cval, ofile);
            }
            if( type == "uchar") {
                unsigned char ucval = '0';
                writeEdgeAttribute(name, ucval, ofile);
            }
            if( type == "double") {
                double dval = '0';
                writeEdgeAttribute(name, dval, ofile);
            }
            if( type == "float") {
                double fval = '0';
                writeEdgeAttribute(name, fval, ofile);
            }
        }
    }
    ofile << "</Edges>" << endl << endl;
    return 0;
}

int
JMeshXMLExporter ::writeFaces(ofstream &ofile)
{
    size_t numfaces = mesh->getSize(2);
    if( numfaces == 0) return 1;

    size_t ncount;

    ofile << "<Faces>" << endl;
    ncount = 0;

    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 3 ) ncount++;
    }
    if( ncount > 0) {
        ofile << "<Triangles> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() && face->getSize(0) == 3 ) {
                for( int j = 0; j < 3; j++)
                    ofile << face->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</Triangles> " << endl << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) == 4 ) ncount++;
    }

    // Write down all quad  faces
    if( ncount > 0) {
        ofile << "<Quads> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() && face->getSize(0) == 4 ) {
                for( int j = 0; j < 4; j++)
                    ofile << face->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</Quads> " << endl << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && face->getSize(0) > 4 ) ncount++;
    }

    // Write down all polygonal faces
    if( ncount > 0) {
        ofile << "<Polygons> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() && face->getSize(0) > 4 ) {
                int nn = face->getSize(0);
                ofile << nn << " ";
                for( int j = 0; j < nn; j++)
                    ofile << face->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</Polygons> " << endl << endl;
    }

    vector<string> attribnames;
    mesh->getAttributeNames(attribnames, 2);
    for( size_t i = 0; i < attribnames.size(); i++) {
        string name = attribnames[i];
        string type = JFace::getAttributeTypeName( attribnames[i] );
        if( !type.empty() ) {
            if( type == "int") {
                int  ival = 0;
                writeFaceAttribute(name, ival, ofile);
            }
            if( type == "char") {
                char cval = '0';
                writeFaceAttribute(name, cval, ofile);
            }
            if( type == "uchar") {
                unsigned char ucval = '0';
                writeFaceAttribute(name, ucval, ofile);
            }
            if( type == "double") {
                double dval = '0';
                writeFaceAttribute(name, dval, ofile);
            }
            if( type == "float") {
                double fval = '0';
                writeNodeAttribute(name, fval, ofile);
            }
        }
    }

    ofile << "</Faces>" << endl;
    return 0;
}

int
JMeshXMLExporter ::writeCells(ofstream &ofile)
{
    size_t numcells = mesh->getSize(3);
    if( numcells == 0) return 1;

    size_t ncount;
    ofile << "<Cells>" << endl;

    ncount = 0;
    for (size_t i = 0; i < numcells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() && cell->getSize(0) ==  4 ) ncount++;
    }

    // Write down all polygonal faces
    if( ncount > 0) {
        ofile << "<Tetrahedra> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;

        for (size_t i = 0; i < numcells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            if( cell->isActive() && cell->getSize(0) == 4 ) {
                for( int j = 0; j < 4; j++)
                    ofile << cell->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</Tetrahedra> " << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numcells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() && cell->getSize(0) ==  8 ) ncount++;
    }

    // Write down all polygonal faces
    if( ncount > 0) {
        ofile << "<Hexahedra> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;

        for (size_t i = 0; i < numcells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            if( cell->isActive() && cell->getSize(0) == 8 ) {
                for( int j = 0; j < 8; j++)
                    ofile << cell->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</Hexahedra> " << endl;
    }

    ncount = 0;
    for (size_t i = 0; i < numcells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() && cell->getSize(0) ==  6 ) ncount++;
    }

    // Write down all polygonal faces
    if( ncount > 0) {
        ofile << "<TriangularPrism> " << endl;
        ofile << "<Count> " <<  ncount <<  " </Count> " <<  endl;

        for (size_t i = 0; i < numcells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            if( cell->isActive() && cell->getSize(0) == 6 ) {
                for( int j = 0; j < 6; j++)
                    ofile << cell->getNodeAt(j)->getID() << " ";
                ofile << endl;
            }
        }
        ofile << "</TriangularPrism> " << endl;
    }

    vector<string> attribnames;
    mesh->getAttributeNames(attribnames, 3);
    for( size_t i = 0; i < attribnames.size(); i++) {
        string name = attribnames[i];
        string type = JCell::getAttributeTypeName( attribnames[i] );
        if( !type.empty() ) {
            if( type == "int") {
                int  ival = 0;
                writeCellAttribute(name, ival, ofile);
            }
            if( type == "char") {
                char cval = '0';
                writeCellAttribute(name, cval, ofile);
            }
            if( type == "uchar") {
                unsigned char ucval = '0';
                writeCellAttribute(name, ucval, ofile);
            }
            if( type == "double") {
                double dval = '0';
                writeCellAttribute(name, dval, ofile);
            }
            if( type == "float") {
                double fval = '0';
                writeCellAttribute(name, fval, ofile);
            }
        }
    }
    ofile << "</Cells>" << endl;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
JMeshXMLExporter ::writeFile(const JMeshPtr &m, const string &s)
{
    mesh = m;
    if( mesh == nullptr) return 1;

    string filename = s;
    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) {
        logger->setWarn("Could not open file for writing mesh ");
        return 1;
    }

    mesh->pruneAll();
    mesh->enumerate(0);

    logger->setInfo("Writing mesh into xml format" );

    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);
    size_t numfaces = mesh->getSize(2);
    size_t numcells = mesh->getSize(3);

    int topodim = -1;
    if( numnodes) topodim = 0;
    if( numedges) topodim = 1;
    if( numfaces) topodim = 2;
    if( numcells) topodim = 3;

    ofile << "<Mesh> " << endl;
    ofile << "<Name> " << mesh->getName() << " </Name> " << endl;
    ofile << "<Format> ASCII </Format> " << endl;
    ofile << "<TopoDimension> " << topodim << " </TopoDimension>" << endl;
    ofile << "<GeomDimension> 3 </GeomDimension>" << endl;

    writeNodes(ofile);
    writeEdges(ofile);
    writeFaces(ofile);
    writeCells(ofile);

    logger->setInfo("Writing mesh info the file completed");

    return 0;
}

