#include "MeshExporter.hpp"

using namespace Jaal;

JLogger* JMeshExporter :: logger = JLogger::getInstance();

JMeshExporter* JMeshExporter :: getProduct(const string & fname)
{
    JMeshExporter *exporter = NULL;

    if (fname.rfind(".off") != string::npos)
        exporter = new JMeshOFFExporter;

    if (fname.rfind(".obj") != string::npos)
        exporter = new JMeshOBJExporter;

    if (fname.rfind(".pov") != string::npos)
        exporter = new JMeshPOVExporter;

    if (fname.rfind(".smf") != string::npos)
        exporter = new JMeshOBJExporter;

    if (fname.rfind(".xml") != string::npos)
        exporter = new JMeshXMLExporter;

    if (fname.rfind(".node") != string::npos)
        exporter = new JMeshTRIExporter;

    if (fname.rfind(".ele") != string::npos)
        exporter = new JMeshTRIExporter;

    if (fname.rfind(".vtk") != string::npos)
        exporter = new JMeshVTKExporter;

    if (fname.rfind(".geo") != string::npos)
        exporter = new JMeshGmshExporter;

    /*
         if (fname.rfind(".dat") != string::npos)
              err = write_simple_file(m, fname);
    */

    return exporter;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshIO :: saveAs(const JMeshPtr &m, const string &fname)
{
    boost::scoped_ptr<JMeshExporter> exporter(JMeshExporter::getProduct(fname));

    if( exporter == nullptr) return 1;
    int err = exporter->writeFile(m, fname);

    return err;
}

////////////////////////////////////////////////////////////////////////////////

