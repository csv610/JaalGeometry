#ifndef MESH_IMPORTER_H
#define MESH_IMPORTER_H

#include "Mesh.hpp"
using namespace Jaal;

class JMeshImporter;
typedef boost::shared_ptr<JMeshImporter> JMeshImporterPtr;

class JMeshImporter
{
public:

    static const int XML_FORMAT      = 0;   // Default format ...
    static const int OFF_FORMAT      = 1;
    static const int OBJ_FORMAT      = 2;
    static const int VTK_FORMAT      = 3;
    static const int TRIANGLE_FORMAT = 4;
    static const int CUBIT_FORMAT    = 5;
    static const int GMSH_FORMAT     = 6;
    static const int SIMPLE_FORMAT   = 7;
    static const int SMF_FORMAT      = 8;
    static const int MESH_FORMAT     = 9;

    static JMeshImporterPtr getProduct( const string &fname);

    virtual ~JMeshImporter() {
        infile.close();
    }

    virtual JMeshPtr readFile( const string &s) = 0;

protected:
    static JLogger *logger;
    std::map<int, int> global2local;
    string  filename;
    ifstream infile;

    int reopen_file() {
        if( infile.is_open() ) infile.close();
        infile.open( filename.c_str(), ios::in);
        if( infile.fail() ) return 1;
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

struct JMeshXMLImporter : public JMeshImporter {
    ~JMeshXMLImporter() {}
    JMeshPtr readFile( const string &s);

    int readNodes( const string &s, const JMeshPtr &m);
    int readEdges( const string &s, const JMeshPtr &m);
    int readFaces( const string &s, const JMeshPtr &m);
    int readCells( const string &s, const JMeshPtr &m);

private:
    JMeshPtr mesh;
    int readAttributeHeader( ifstream &infile, string &name, string &type, size_t &ncount);

    template<class T>
    int readNodeAttribute(const string &name , const T &, size_t ncount, ifstream &infile)
    {
        T  val;
        size_t numnodes = mesh->getSize(0);
        if( ncount == numnodes) {
            for( size_t i = 0; i < numnodes; i++) {
                infile >> val;
                const JNodePtr &vtx = mesh->getNodeAt(i);
                vtx->setAttribute(name, val);
            }
            return 0;
        }

        size_t id;
        for( size_t i = 0; i < ncount; i++) {
            infile >> id >> val;
            const JNodePtr &vtx = mesh->getNodeAt(id);
            vtx->setAttribute(name, val);
        }
        return 0;
    }

    template<class T>
    int readEdgeAttribute(const string &name , const T &, size_t ncount, ifstream &infile)
    {
        T  val;
        size_t numedges = mesh->getSize(1);
        if( ncount == numedges) {
            for( size_t i = 0; i < numedges; i++) {
                infile >> val;
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                edge->setAttribute(name, val);
            }
            return 0;
        }

        size_t id;
        for( size_t i = 0; i < ncount; i++) {
            infile >> id >> val;
            const JEdgePtr &edge = mesh->getEdgeAt(id);
            edge->setAttribute(name, val);
        }
        return 0;
    }

    template<class T>
    int readFaceAttribute(const string &name , const T &, size_t ncount, ifstream &infile)
    {
        cout << "reading face attrib " << endl;
        T  val;
        size_t numfaces = mesh->getSize(2);
        if( ncount == numfaces) {
            for( size_t i = 0; i < numfaces; i++) {
                infile >> val;
                const JFacePtr &face = mesh->getFaceAt(i);
                face->setAttribute(name, val);
            }
            return 0;
        }

        size_t id;
        for( size_t i = 0; i < ncount; i++) {
            infile >> id >> val;
            const JFacePtr &face = mesh->getFaceAt(id);
            face->setAttribute(name, val);
        }
        return 0;
    }

    template<class T>
    int readCellAttribute(const string &name , const T &, size_t ncount, ifstream &infile)
    {
        T  val;
        size_t numcells = mesh->getSize(3);
        if( ncount == numcells) {
            for( size_t i = 0; i < numcells; i++) {
                infile >> val;
                const JCellPtr &cell = mesh->getCellAt(i);
                cell->setAttribute(name, val);
            }
            return 0;
        }

        size_t id;
        for( size_t i = 0; i < ncount; i++) {
            infile >> id >> val;
            const JCellPtr &cell = mesh->getCellAt(id);
            cell->setAttribute(name, val);
        }
        return 0;
    }
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshTRIImporter : public JMeshImporter {
    ~JMeshTRIImporter() {}

    JMeshPtr readFile( const string &s);
    JMeshPtr readPoly( const string &s);

    int readNodes(const string & s, const JMeshPtr &m);
    int readEdges(const string & s, const JMeshPtr &m);
    int readFaces(const string & s, const JMeshPtr &m);
    int readCells(const string & s, const JMeshPtr &m);
private:
    int start_index;
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshVTKImporter : public JMeshImporter {
    ~JMeshVTKImporter() {}
    JMeshPtr readFile( const string &s);
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshOFFImporter : public JMeshImporter {
    ~JMeshOFFImporter() {}
    JMeshPtr readFile( const string &s);
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshOBJImporter : public JMeshImporter {
    ~JMeshOBJImporter() {}

    JMeshPtr readFile( const string &s);

private:
    vector<double> nodeCoords, uvCoords;
    vector<size_t> triConnect, quadConnect;
    vector<vector<size_t>> polyConnect;
    vector<size_t> connect, uvConnect;

    vector<size_t> uvTriConnect, uvQuadConnect, uvPolyConnect;

    bool has_texture;
    void readFaceLine(const string &s );
    JMeshPtr getMesh();
    JMeshPtr getTextureMesh();
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshCubitImporter : public JMeshImporter {
    ~JMeshCubitImporter() {}
    JMeshPtr readFile( const string &s);
};
///////////////////////////////////////////////////////////////////////////////

struct JMeshMeditImporter : public JMeshImporter {
    ~JMeshMeditImporter() {}
    JMeshPtr readFile( const string &s);
};

struct JMeshGmshImporter : public JMeshImporter {
    ~JMeshGmshImporter() {}
    JMeshPtr readFile( const string &s);
};
///////////////////////////////////////////////////////////////////////////////

#endif
