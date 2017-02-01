#include "MeshImporter.hpp"
#include "MeshExporter.hpp"
#include "MeshTopology.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"
#include "MeshQuality.hpp"
#include <iostream>

using namespace Jaal;


///////////////////////////////////////////////////////////////////////////////
bool AllTetMeshGenerator :: isDelaunay(const JMeshPtr &mesh)
{
}
///////////////////////////////////////////////////////////////////////////////

void AllTetMeshGenerator :: replaceNodes( const JNodeSequence &oldnodes, JNodeSequence &newnodes)
{
    int numNodes = min( oldnodes.size(), newnodes.size() );

    for( int i = 0;  i < numNodes; i++) {
        const JNodePtr &oldnode = oldnodes[i];
        const JNodePtr &newnode = newnodes[i];
        const Point3D  &pa = oldnode->getXYZCoords();
        const Point3D  &pb = newnode->getXYZCoords();
        double d  = JMath::length2(pa,pb);
        if( d < 1.0E-06) newnodes[i] = oldnode;
    }
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: genMesh(const JMeshPtr &mesh, const string &cmd)
{
    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2 ) {
        cout << "Error: Tetgen require triangle surface mesh " << endl;
        return nullptr;
    }


    JMeshTRIExporter mexp;
    mexp.writeFacets(mesh, "tmp.smesh");

    // Calling Tetgen  ....
    string cmd1 = cmd +  "  tmp.smesh";
    system( cmd1.c_str() );

    JMeshTRIImporter mimp;
    JMeshPtr tetmesh = mimp.readFile( "tmp.1.ele");

    return tetmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getConvexHull(const JNodeSequence &nodes)
{
    string filename = "tmp.node";

    ofstream ofile( filename.c_str(), ios::out);
    if( ofile.fail() )  {
        cout << "Warning: cann't open node file " << filename << endl;
        return NULL;
    }
    size_t numnodes, ndim = 3, numattrib = 0, boundflag = 0, bmark = 0;

    numnodes = nodes.size();

    ofile << numnodes << " " << ndim << " " << numattrib << " " << boundflag << endl;

    for( size_t i = 0; i < numnodes; i++)  {
        const Point3D &p3d = nodes[i]->getXYZCoords();
        nodes[i]->setID(i);
        ofile <<  i << " " <<  p3d[0] << " " <<  p3d[1]  << " " << p3d[2];
        if( boundflag ) ofile <<  bmark;
        ofile << endl;
    }
    ofile.close();

    int err = system( "tetgen -CO tmp.node");
    if( err < 0) return nullptr;

    JMeshTRIImporter mimp;
    JMeshPtr tetmesh = mimp.readFile( "tmp.1.ele");

    tetmesh->getTopology()->searchBoundary();
    size_t numfaces = tetmesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = tetmesh->getFaceAt(i);
        if( face->getNumRelations(3) == 2)
            face->setStatus(JMeshEntity::REMOVE);
    }

    return tetmesh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTetMeshGenerator :: getConvexHull(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return nullptr;
    JNodeSequence nodes = mesh->getNodes();
    return getConvexHull(nodes);
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getConstrainedMesh(const JMeshPtr &msh)
{
    mesh = msh;
    if( mesh == nullptr) return nullptr;
    if( mesh->getSize(3) ) mesh->deleteVolumeMesh();

//  preprocess();
//  tetrahedralize("pYS0", &inmesh, &outmesh);
//  postprocess();

    JNodeSequence oldnodes = mesh->getNodes();
    string cmd  = "tetgen -pYS0 tmp.smesh";
    JMeshPtr surfmesh = mesh;

    if( surfmesh == nullptr) return nullptr;

    JMeshTRIExporter mexp;
    mexp.writeFacets(mesh, "tmp.smesh");

    int err = system( cmd.c_str() );
    if( err < 0) return nullptr;

    JMeshTRIImporter mimp;

    JMeshPtr tetmesh = JMesh::newObject();
    mimp.readNodes( "tmp.1.node", tetmesh);

    JNodeSequence newnodes = tetmesh->getNodes();
    tetmesh->clearAll();

    replaceNodes( oldnodes, newnodes);
    tetmesh->addObjects( newnodes);

    mimp.readCells( "tmp.1.ele", tetmesh);

    if( tetmesh == nullptr) return nullptr;

    for( size_t i = 0; i < mesh->getSize(0); i++) {
        if( tetmesh->getNodeAt(i) != mesh->getNodeAt(i)  ) {
            cout << "Warning: Boundary nodes not matching " << endl;
            break;
        }
        if( !mesh->getNodeAt(i)->isBoundary() )  {
            cout << "Warning: input node is not on boundary" << endl;
            break;
        }
        if( !tetmesh->getNodeAt(i)->isBoundary() )  {
            cout << "Warning: output node is not on boundary" << endl;
            break;
        }
    }

    for( size_t i = 0; i < mesh->getSize(2); i++) {
        if( !mesh->getFaceAt(i)->isBoundary() ) {
            cout << "Warning: input face not on the boundary " << endl;
            break;
        }
        if( !tetmesh->getFaceAt(i)->isBoundary() ) {
            cout << "Warning: input face not on the boundary " << endl;
            break;
        }
    }

    if( tetmesh->getSize(0) != surfmesh->getSize(0) ) {
        cout << "Warning : new nodes introduced the contrained mesh " << endl;
    }

    return tetmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getReMesh(const JMeshPtr &mesh, const vector<Point3D> &addPoints)
{
    if(mesh == nullptr) return nullptr;

    JMeshTRIExporter mexp;
    ostringstream  oss;
    oss << "tetgen ";

//  mexp.writeNodes(mesh, "tmp.node");
    mexp.writeFacets(mesh, "tmp.smesh");

    if( mesh->getSize(3) == 0)
        oss << "-pQ";
    else {
        oss << "-CQYdzOq1.4r";
        mexp.writeCells(mesh, "tmp.ele");
    }

    if( !addPoints.empty() ) {
        oss << "i";
        ofstream ofile("tmp.a.node", ios::out);
        ofile <<   addPoints.size() << "  3  0  0  " << endl;
        int index = 0;
        for( const Point3D &xyz : addPoints)
            ofile << index++ << "  " << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
        ofile.close();
    }

    oss << " tmp.smesh";

    string cmd = oss.str();

    int err = system( cmd.c_str() );
    if( err < 0) return nullptr;

    JMeshTRIImporter mimp;
    JMeshPtr tetmesh = mimp.readFile( "tmp.1.ele");
    if( tetmesh == nullptr) return tetmesh;

    /*
        tetmesh->getTopology()->search_boundary();

        size_t numnodes = tetmesh->getTopology()->getBoundarySize(0);
        if( numnodes != mesh->getSize(0) )
            cout << "Warning: new nodes on the boundary introduced " << endl;

        size_t numedges = tetmesh->getTopology()->getBoundarySize(1);
        if( numedges != mesh->getSize(1) )
            cout << "Warning: new edges on the boundary introduced " << endl;

        size_t numfaces = tetmesh->getTopology()->getBoundarySize(2);
        if( numfaces != mesh->getSize(2) )
            cout << "Warning: new faces on the boundary introduced " << endl;
    */

    return tetmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getQualityMesh(const JMeshPtr &mesh, bool modifySurface)
{
    ostringstream  oss;
//  oss << "tetgen -CzO10/7q1.5/10.0V";
//  if( !modifySurface) oss << "Y";
//  if( maxVolume > 0.0) oss << "a" << maxVolume;

    JNodeSequence oldnodes = mesh->getNodes();

    JMeshTRIExporter mexp;
    mexp.writeNodes(mesh, "tmp.node");
    mexp.writeFacets(mesh, "tmp.smesh");

    oss << "tetgen -" << options <<  " tmp.smesh";

    string cmd = oss.str();
    cout << "cmd " << cmd << endl;
    int err = system( cmd.c_str() );
    if( err < 0) return nullptr;

    JMeshPtr tetmesh = JMesh::newObject();

    JMeshTRIImporter mimp;
    mimp.readNodes( "tmp.1.node", tetmesh);

    JNodeSequence newnodes = tetmesh->getNodes();

    tetmesh->clearAll();
    replaceNodes( oldnodes, newnodes);

    tetmesh->addObjects(newnodes);
    mimp.readCells( "tmp.1.ele", tetmesh);

    size_t nsize;
    if( !modifySurface) {
        tetmesh->getTopology()->collectFaces();
        tetmesh->getTopology()->collectEdges();

        nsize = mesh->getSize(0);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getNodeAt(i) ) ) {
                cout << "Warning: Input node not found in the tetmesh " << endl;
                exit(0);
            }
        }

        nsize = mesh->getSize(1);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getEdgeAt(i) ) ) {
                cout << "Warning: Input edge not found in the tetmesh " << endl;
                exit(0);
            }
        }

        nsize = mesh->getSize(2);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getFaceAt(i) ) ) {
                cout << "Warning: Input face not found in the tetmesh " << endl;
                exit(0);
            }
        }
    }

    return tetmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getQualityMesh(const JMeshPtr &mesh, const vector<Point3D> &addPoints, bool modifySurface)
{
    ostringstream  oss;
//  oss << "tetgen -CzO10/7q1.5/10.0i";
    oss << "tetgen -CziQ";
    if( !modifySurface) oss << "Y";

    /*
        JMeshQuality mq;
        mq.setMesh(mesh);
        double len = mq.getMeanEdgeLength();
        double vol = TetGeometry::getRegularTetrahedronVolume(len);
        oss << "a" << 1.5*vol;
    */

    JMeshTRIExporter mexp;
    mexp.writeNodes(mesh, "tmp.node");
    mexp.writeFacets(mesh, "tmp.smesh");

    if( !addPoints.empty() ) {
        oss << "i";
        ofstream ofile("tmp.a.node", ios::out);
        ofile <<   addPoints.size() << "  3  0  0  " << endl;
        int index = 0;
        for( const Point3D &xyz : addPoints)
            ofile << index++ << "  " << xyz[0] << "  " << xyz[1] << "  " << xyz[2] << endl;
        ofile.close();
    }

    oss << " tmp.smesh";

    string cmd = oss.str();
    int err = system( cmd.c_str() );
    if( err < 0) return nullptr;

    JMeshPtr tetmesh = JMesh::newObject();

    JMeshTRIImporter mimp;
    mimp.readNodes( "tmp.1.node", tetmesh);

    JNodeSequence newnodes = tetmesh->getNodes();

    tetmesh->clearAll();

    int numNodes = mesh->getSize(0);
    for( int i = 0;  i < numNodes; i++) {
        const JNodePtr &oldnode = mesh->getNodeAt(i);
        const JNodePtr &newnode = newnodes[i];
        const Point3D  &pa = oldnode->getXYZCoords();
        const Point3D  &pb = newnode->getXYZCoords();
        double d  = JMath::length2(pa,pb);
        if( d < 1.0E-06)
            newnodes[i] = oldnode;
        else {
            cout << "Warning: Input mesh nodes not matched in the tetmesh " << d << endl;
            cout << "OLD " << pa[0] << " " << pa[1] << " " << pa[2] << endl;
            cout << "NEW " << pb[0] << " " << pb[1] << " " << pb[2] << endl;
        }
    }
    tetmesh->addObjects(newnodes);
    mimp.readCells( "tmp.1.ele", tetmesh);

    size_t nsize;
    if( !modifySurface) {
        tetmesh->getTopology()->collectFaces();
        tetmesh->getTopology()->collectEdges();

        nsize = mesh->getSize(0);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getNodeAt(i) ) ) {
                cout << "Warning: Input node not found in the tetmesh " << endl;
                exit(0);
            }
        }

        nsize = mesh->getSize(1);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getEdgeAt(i) ) ) {
                cout << "Warning: Input edge not found in the tetmesh " << endl;
                exit(0);
            }
        }

        nsize = mesh->getSize(2);
        for( size_t i = 0; i < nsize; i++)  {
            if (!tetmesh->contains( mesh->getFaceAt(i) ) ) {
                cout << "Warning: Input face not found in the tetmesh " << endl;
                exit(0);
            }
        }
    }

    return tetmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: fromHexMesh(const JMeshPtr &hexmesh)
{
    JMeshPtr tetmesh = JMesh::newObject();

    size_t numnodes = hexmesh->getSize(0);
    size_t numcells = hexmesh->getSize(3);

    tetmesh->reserve( numnodes, 0);
    tetmesh->reserve( 5*numcells, 3);

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = hexmesh->getNodeAt(i);
        if( vtx->isActive() ) tetmesh->addObject(vtx);
    }

    vector<JCellPtr> celltets;
    size_t nCount = 0;
    for( size_t i = 0; i < numcells; i++) {
        JHexahedronPtr hex = JHexahedron::down_cast(hexmesh->getCellAt(i));
        if( hex )
            if( hex->isActive() ) {
                hex->getTetrahedra( celltets );
                tetmesh->addObjects(celltets);
                hex->setStatus(JMeshEntity::REMOVE);
                nCount++;
            }
    }
    assert( tetmesh->getSize(3) == 6*nCount );

    return tetmesh;
}
///////////////////////////////////////////////////////////////////////////////
void AllTetMeshGenerator :: preprocess()
{
    if( mesh == nullptr) return;

    inmesh.firstnumber = 0;
    size_t numnodes = mesh->getSize(0);
    inmesh.pointlist    = new double[3*numnodes];
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D  &xyz = vtx->getXYZCoords();
        inmesh.pointlist[3*i+0] = xyz[0];
        inmesh.pointlist[3*i+1] = xyz[1];
        inmesh.pointlist[3*i+2] = xyz[2];
    }

    size_t   numfaces = mesh->getSize(2);
    inmesh.numberoffacets = numfaces;
    inmesh.facetlist = new tetgenio::facet[numfaces];
    inmesh.facetmarkerlist = new int[numfaces];

    tetgenio::facet *f;
    tetgenio::polygon *p;
    for( size_t i = 0;  i < numfaces; i++) {
        const JFacePtr &tri = mesh->getFaceAt(i);
        f = &inmesh.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = tri->getNodeAt(0)->getID();
        p->vertexlist[1] = tri->getNodeAt(1)->getID();
        p->vertexlist[2] = tri->getNodeAt(2)->getID();
    }
}
///////////////////////////////////////////////////////////////////////////////

void AllTetMeshGenerator :: postprocess()
{
    int numPoints = outmesh.numberofpoints;
    int numSurfPoints = inmesh.numberofpoints;

    JNodeSequence newnodes;
    newnodes.reserve( numPoints);

    for( int i = 0; i < numSurfPoints; i++) {
        double x0 = inmesh.pointlist[3*i+0];
        double y0 = inmesh.pointlist[3*i+1];
        double z0 = inmesh.pointlist[3*i+2];

        double x1 = outmesh.pointlist[3*i+0];
        double y1 = outmesh.pointlist[3*i+1];
        double z1 = outmesh.pointlist[3*i+2];
        if( x0 == x1 && y0 == y1 && z0 == z1 )
            newnodes.push_back( mesh->getNodeAt(i) );
        else {
            JNodePtr vtx = JNode::newObject();
            vtx->setXYZCoords(x1,y1,z1);
            newnodes.push_back( vtx ) ;
        }
    }

    for( int i = numSurfPoints; i < numPoints; i++) {
        double x1 = outmesh.pointlist[3*i+0];
        double y1 = outmesh.pointlist[3*i+1];
        double z1 = outmesh.pointlist[3*i+2];
        JNodePtr vtx = JNode::newObject();
        vtx->setXYZCoords(x1,y1,z1);
        newnodes.push_back( vtx ) ;
    }

    int numCells = outmesh.numberoftetrahedra;
    JNodeSequence connect(4);

    JCellSequence newtets;
    newtets.resize(numCells);
    for( int i = 0; i < numCells; i++) {
        connect[0] = newnodes[ outmesh.tetrahedronlist[4*i+0] ];
        connect[1] = newnodes[ outmesh.tetrahedronlist[4*i+1] ];
        connect[2] = newnodes[ outmesh.tetrahedronlist[4*i+2] ];
        connect[3] = newnodes[ outmesh.tetrahedronlist[4*i+3] ];
        newtets[i] = JTetrahedron::newObject( connect );
    }

    newtetmesh = JMesh::newObject();
    newtetmesh->addObjects( newnodes);
    newtetmesh->addObjects( newtets );
}
///////////////////////////////////////////////////////////////////////////////
