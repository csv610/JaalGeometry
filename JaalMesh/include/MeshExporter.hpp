#ifndef MESH_EXPORTER_H
#define MESH_EXPORTER_H

#include "Mesh.hpp"
#include <iomanip>

using namespace Jaal;

class JMeshExporter
{
public:
    static const int XML_FORMAT    = 0;   // Most Versatile format( Native format )
    static const int SIMPLE_FORMAT = 1;   // Non-Standard and Primiitive.
    static const int OFF_FORMAT    = 2;
    static const int OBJ_FORMAT    = 3;   //
    static const int SMF_FORMAT    = 4;   // SMF is a truncated OBJ file. We don't need SMF extensions
    static const int VTK_FORMAT    = 5;   // For ParaView Rendering
    static const int TRIANGLE_FORMAT = 6; // Jonathan Shewchuk's Triangle ...
    static const int POV_FORMAT   = 7;    // PovRay for visualization( Surface only)
    static const int RIB_FORMAT   = 8;    // Renderman for visualization( Surface Only)

    static JMeshExporter *getProduct( const string &name);
    JMeshExporter() {
        precision = 6;
    }

    void clearAttributes()
    {
        nodeAttribs.clear();
        edgeAttribs.clear();
        faceAttribs.clear();
        cellAttribs.clear();
    }

    void addNodeAttribute( const string &s) {
        if( find( nodeAttribs.begin(), nodeAttribs.end(), s) == nodeAttribs.end() )
            nodeAttribs.push_back(s);
    }

    void addEdgeAttribute( const string &s) {
        if( find( edgeAttribs.begin(), edgeAttribs.end(), s) == edgeAttribs.end() )
            edgeAttribs.push_back(s);
    }

    void addFaceAttribute( const string &s) {
        if( find( faceAttribs.begin(), faceAttribs.end(), s) == faceAttribs.end() )
            faceAttribs.push_back(s);
    }
    void addCellAttribute( const string &s) {
        if( find( cellAttribs.begin(), cellAttribs.end(), s) == cellAttribs.end() )
            cellAttribs.push_back(s);
    }

    virtual ~JMeshExporter() {}

    virtual int writeFile(const JMeshPtr &, const string &) {
        return 1;
    }
    void setPrecision( int p) {
        precision = p;
    }

    vector<string> nodeAttribs;
    vector<string> edgeAttribs;
    vector<string> faceAttribs;
    vector<string> cellAttribs;
protected:
    static JLogger *logger;
    int precision;
};

struct JMeshOFFExporter : public JMeshExporter {
    ~JMeshOFFExporter() {}
    int writeFile(const JMeshPtr &m, const string & s);
};

struct JMeshOBJExporter : public JMeshExporter {
    ~JMeshOBJExporter() {}
    int writeFile(const JMeshPtr &m, const string & s);
};

struct JMeshPOVExporter : public JMeshExporter {
    ~JMeshPOVExporter() {}

    void setCamera( const Point3D &p) {
         cameraPos = p; 
    }
    void setLookAt( const Point3D &p) {
         lookAtPos = p;
    }

    void addLight(const Point3D &p) {
         lights.push_back(p);
    }

    int writeFile(const JMeshPtr &m, const string & s);
    
private:
    Point3D  cameraPos;
    Point3D  lookAtPos;
    vector<Point3D> lights;
    void writeheader(const JMeshPtr &m, ofstream & s);
    void writemesh(const JMeshPtr &m,   ofstream & s);
    void writemesh2(const JMeshPtr &m,   ofstream & s);
};

struct JMeshXMLExporter : public JMeshExporter {

    ~JMeshXMLExporter() {}
    int writeFile(const JMeshPtr &m, const string & s);
private:
    JMeshPtr mesh;

    template<class T>
    int writeNodeAttribute(const string &name , const T &dummy, ofstream &ofile)
    {
        size_t numActive = 0, nCount = 0;
        size_t numnodes = mesh->getSize(0);
        T val;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) numActive++;
            int err = vtx->getAttribute(name, val);
            if( !err) nCount++;
        }
        if( nCount == 0) return 1;

        ofile << "<NodeAttribute>" << endl;
        ofile << "<Name>\t"     <<  name << "\t</Name>" << endl;
        ofile << "<DataType>\t" <<  get_pod_name(dummy) << "\t</DataType>" << endl;
        ofile << "<Count>\t"    << nCount <<  "\t</Count> " << endl;

        if( numActive == nCount) {
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &vtx = mesh->getNodeAt(i);
                int err = vtx->getAttribute(name, val);
                if( !err) ofile << val << endl;
            }
        } else {
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &vtx = mesh->getNodeAt(i);
                int err = vtx->getAttribute(name, val);
                if( !err) ofile << vtx->getID() << " " <<  val << endl;
            }
        }
        ofile << "</NodeAttribute>" << endl;
        return 0;
    }

    template<class T>
    int writeEdgeAttribute(const string &name , const T &dummy, ofstream &ofile)
    {
        T val;
        size_t numActive= 0, nCount = 0;
        size_t numedges = mesh->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) numActive++;
            int err = edge->getAttribute(name, val);
            if( !err) nCount++;
        }

        ofile << "<EdgeAttribute>" << endl;
        ofile << "<Name>\t" << name << "\t</Name>" << endl;
        ofile << "<DataType>\t" << get_pod_name(dummy) << "\t</DataType>" << endl;
        ofile << "<Count>\t" << nCount <<  "\t</Count> " << endl;

        if( numActive == nCount) {
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                int err = edge->getAttribute(name, val);
                if( !err) ofile << val << endl;
            }
        } else {
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                int err = edge->getAttribute(name, val);
                if( !err) ofile << edge->getID() << "  " << val << endl;
            }
        }

        ofile << "</EdgeAttribute>" << endl;
        return 0;
    }

    template<class T>
    int writeFaceAttribute(const string &name , const T &dummy, ofstream &ofile)
    {
        T val;
        size_t numActive =0, nCount = 0;
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) numActive++;
            int err = face->getAttribute(name, val);
            if( !err) nCount++;
        }


        ofile << "<FaceAttribute>" << endl;
        ofile << "<Name>\t" <<  name << "\t</Name>" << endl;
        ofile << "<DataType>\t" <<  get_pod_name(dummy) << "\t</DataType>" << endl;
        ofile << "<Count>\t" << nCount <<  "\t</Count> " << endl;

        if( numActive == nCount) {
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                int err = face->getAttribute(name, val);
                if( !err ) ofile << val << endl;
            }
        } else {
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                int err = face->getAttribute(name, val);
                if( !err ) ofile << face->getID() << " " << val << endl;
            }
        }
        ofile << "</FaceAttribute>" << endl;
        return 0;
    }

    template<class T>
    int writeCellAttribute(const string &name , const T &dummy, ofstream &ofile)
    {
        ofile << "<CellAttribute>" << endl;
        ofile << "<Name>\t" <<  name << "\t</Name>" << endl;
        ofile << "<DataType>\t" <<  get_pod_name(dummy) << "\t</DataType>" << endl;

        size_t numcells = mesh->getSize(3);
        T val;
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &face = mesh->getCellAt(i);
            int err = face->getAttribute(name, val);
            if( !err ) ofile << val << endl;
        }
        ofile << "</CellAttribute>" << endl;
        return 0;
    }

    int writeNodes(ofstream &s);
    int writeEdges(ofstream &s);
    int writeFaces(ofstream &s);
    int writeCells(ofstream &s);
};

struct JMeshTRIExporter : public JMeshExporter {
    JMeshTRIExporter() {
        dim = 3;
    }
    ~JMeshTRIExporter() {}

    void setDimension(int d) {
        dim = d;
    }
    int writeFile(const JMeshPtr &m, const string & s);

    int writeNodes( const JMeshPtr &m, const string &s);
    int writeEdges( const JMeshPtr &m, const string &s);
    int writeFaces( const JMeshPtr &m, const string &s);
    int writeCells( const JMeshPtr &m, const string &s);
    int writeFacets(const JMeshPtr &m, const string & s);
private:
    int dim;
};

struct JMeshVTKExporter : public JMeshExporter {
    ~JMeshVTKExporter() {}
    int writeFile(const JMeshPtr &m, const string & s);
};

struct JMeshGmshExporter : public JMeshExporter {
    JMeshGmshExporter() {
        lc = 1.0;
    }
    ~JMeshGmshExporter() {}
    void setCharacteristicLength( double l)  {
        lc = l;
    }
    int writeFile(const JMeshPtr &m, const string & s);
    double lc;
};

#endif

