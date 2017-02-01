#include "AllTriMeshGenerator.hpp"
#include "MeshAffineTransforms.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTriMeshGenerator :: getCylinder( const Point3D &pbase, const Point3D &ptip, double radius, int nR)
{
    double height = JMath::length( pbase, ptip);
    int    nH     = 0.25*nR*height/(2.0*M_PI*radius);

    ostringstream oss;
    oss << "mesh_make ccyl " << nR << " " << nH << "  " << radius << " cyl.off";

    string cmd = oss.str();
    system( cmd.c_str() );

    JMeshOFFImporter mimp;
    JMeshPtr  cylmesh = mimp.readFile( "cyl.off");

    JMeshAffineTransform affine;
    affine.setMesh(cylmesh);
    affine.translate(0.0, 0.0, 1.0);
    Point3D p0, p1;
    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;
    JNodePtr vbase = cylmesh->getGeometry()->getNearest(p0);
    p1[0] = 0.0;
    p1[1] = 0.0;
    p1[2] = 2.0;
    JNodePtr vtip  = cylmesh->getGeometry()->getNearest(p1);
    affine.scale(1.0, 1.0, 0.5*height);

    p0 = vbase->getXYZCoords();
    affine.translate(0.0, 0.0, -p0[2]);

    p0 = vbase->getXYZCoords();
    p1 = vtip->getXYZCoords();
    assert(fabs(JMath::length(p0,p1) - height) < 1.0E-10);

    Vec3D vsrc, vdst;
    vsrc[0] = 0.0;
    vsrc[1] = 0.0;
    vsrc[2] = 1.0;
    double dx = ptip[0] - pbase[0];
    double dy = ptip[1] - pbase[1];
    double dz = ptip[2] - pbase[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);
    vdst[0]   = dx/dl;
    vdst[1]   = dy/dl;
    vdst[2]   = dz/dl;

    affine.alignAlong( vsrc, vdst);
    affine.translate( pbase[0], pbase[1], pbase[2] );

    double dist;
    p0 = vbase->getXYZCoords();
    dist = JMath::length(pbase, p0);
    if( dist > 1.0E-06)
        cout << "Warning: Cylinder base does not match " << dist << endl;

    p0 = vtip->getXYZCoords();
    dist = JMath::length(ptip, p0);
    if( dist > 1.0E-06) {
        cout << "Warning: Cylinder tip does not match " << dist << endl;
        cout << "Base " << pbase[0] << " " << pbase[1] << " " << pbase[2] << endl;
        cout << "tip  " << ptip[0]  << " " << ptip[1]  << " " << ptip[2] << endl;
    }

    return cylmesh;
}
////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTriMeshGenerator :: getFromQuad2Tri( int *dim,  double *len, double *org, bool texCoord)
{
    if( dim == nullptr ) return nullptr;

    int nx = dim[0];
    int ny = dim[1];
    if( nx < 2 || ny < 2 ) return nullptr;


    double xlen = 1.0, ylen = 1.0;
    if( len ) {
        xlen = len[0];
        ylen = len[1];
    }

    double xorg = 0.0, yorg = 0.0;
    if( org ) {
        xorg = org[0];
        yorg = org[1];
    }

    double dx = xlen/(double)(nx - 1);
    double dy = ylen/(double)(ny - 1);

    double du = 1.0/(double)(nx-1);
    double dv = 1.0/(double)(ny-1);

    Point3D xyz;
    Point2D uv;


    JMeshPtr mesh = JMesh::newObject();
    JNodeSequence nodes(nx*ny);

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = xorg + i * dx;
            xyz[1] = yorg + j * dy;
            xyz[2] = 0.0;
            JNodePtr vnew = JNode::newObject();
            vnew->setXYZCoords(xyz);
            if( texCoord) {
                uv[0] = i*du;
                uv[1] = j*dv;
                vnew->setAttribute("UVCoords", uv);
            }
            nodes[index] = vnew;
            vnew->setID(index);
            index++;
        }
    }
    mesh->addObjects(nodes);

    JFaceSequence faces(2*(nx-1)*(ny-1));

    index = 0;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            const JNodePtr &v0 = mesh->getNodeAt(n0);
            const JNodePtr &v1 = mesh->getNodeAt(n1);
            const JNodePtr &v2 = mesh->getNodeAt(n2);
            const JNodePtr &v3 = mesh->getNodeAt(n3);
            JFacePtr t1 =  JTriangle::newObject( v0, v1, v2);
            JFacePtr t2 =  JTriangle::newObject( v0, v2, v3);
            faces[index++] = t1;
            faces[index++] = t2;
        }
    }

    mesh->addObjects(faces);

    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTriMeshGenerator :: getFromQuad4Tri( int *dim,  double *len, double *org, bool texCoord)
{
    if( dim == nullptr ) return nullptr;

    int nx = dim[0];
    int ny = dim[1];
    if( nx < 2 || ny < 2 ) return nullptr;


    double xlen = 1.0, ylen = 1.0;
    if( len ) {
        xlen = len[0];
        ylen = len[1];
    }

    double xorg = 0.0, yorg = 0.0;
    if( org ) {
        xorg = org[0];
        yorg = org[1];
    }

    double dx = xlen/(double)(nx - 1);
    double dy = ylen/(double)(ny - 1);

    double du = 1.0/(double)(nx-1);
    double dv = 1.0/(double)(ny-1);

    Point3D xyz;
    Point2D uv;

    JMeshPtr mesh = JMesh::newObject();
    JNodeSequence nodes(nx*ny);

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = xorg + i * dx;
            xyz[1] = yorg + j * dy;
            xyz[2] = 0.0;
            JNodePtr vnew = JNode::newObject();
            vnew->setXYZCoords(xyz);
            if( texCoord) {
                uv[0] = i*du;
                uv[1] = j*dv;
                vnew->setAttribute("UVCoords", uv);
            }
            nodes[index] = vnew;
            vnew->setID(index);
            index++;
        }
    }
    mesh->addObjects(nodes);

    JFaceSequence faces(4*(nx-1)*(ny-1));

    index = 0;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            const JNodePtr &v1 = mesh->getNodeAt(n0);
            const JNodePtr &v2 = mesh->getNodeAt(n1);
            const JNodePtr &v3 = mesh->getNodeAt(n2);
            const JNodePtr &v4 = mesh->getNodeAt(n3);
            JNodePtr v0 = JFaceGeometry::getCentroid(v1,v2,v3,v4);
            mesh->addObject(v0);
            faces[index++] = JTriangle::newObject(v0, v1, v2);
            faces[index++] = JTriangle::newObject(v0, v2, v3);
            faces[index++] = JTriangle::newObject(v0, v3, v4);
            faces[index++] = JTriangle::newObject(v0, v4, v1);
        }
    }
    mesh->addObjects(faces);

    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTriMeshGenerator :: getFromQuadMesh( const JMeshPtr &quadmesh, JNodeSequence &steiner, 
int trisPerQuad, bool randomize)
{
    JMeshPtr trimesh;

    if( quadmesh == nullptr ) return trimesh;

    quadmesh->pruneNodes();
    quadmesh->pruneFaces();

    size_t numfaces = quadmesh->getSize(2);

    if( numfaces == 0) {
        cout << "Warning: the mesh is empty" << endl;
        return NULL;
    }

    srand (time(NULL)); // Initialize necessary

    size_t numnodes = quadmesh->getSize(0);

    trimesh = JMesh::newObject();

    trimesh->reserve( numnodes + numfaces, 0);

    JNodeSequence oldnodes = quadmesh->getNodes();
    trimesh->addObjects( oldnodes);

    JEdgeSequence oldedges = quadmesh->getEdges();
    trimesh->addObjects(oldedges);

    steiner.clear();

    Point3D p3d;
    JFacePtr  tface;
    JNodeSequence nodes(3);

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = quadmesh->getFaceAt(i);
        if( f->isActive() && f->getSize(0) == 3)
            trimesh->addObject(f);
    }


    if( trisPerQuad  == 4 ) {
        trimesh->reserve( 4*numfaces, 2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = quadmesh->getFaceAt(i);
            if( f->isActive()  && f->getSize(0) == 4) {
                int pos = JFaceGeometry::reflexAngleAt(f);
                if( pos >= 0) {
                    p3d = JNodeGeometry::getMidPoint( f->getNodeAt(pos), f->getNodeAt(pos+2));
                } else {
                    f->getAvgXYZ( p3d );
                }
                JNodePtr v0 = JNode::newObject();
                v0->setXYZCoords( p3d );
                steiner.push_back(v0);
                trimesh->addObject(v0);
                nodes[0] = v0;
                for( int j = 0; j < 4; j++) {
                    nodes[1] = f->getNodeAt(j);
                    nodes[2] = f->getNodeAt(j+1);
                    tface  = JTriangle::newObject(nodes);
                    trimesh->addObject( tface );
                }
            }
        }
    }

    int pos;
    if( trisPerQuad == 2 ) {
        trimesh->reserve( 2*numfaces, 2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = quadmesh->getFaceAt(i);
            if( f->isActive() && f->getSize(0) == 4 ) {
                pos = JFaceGeometry::reflexAngleAt(f);
                if( pos >= 0) {
                    nodes[0] = f->getNodeAt(pos);
                    nodes[1] = f->getNodeAt(pos+1);
                    nodes[2] = f->getNodeAt(pos+2);
                    tface  = JTriangle::newObject( nodes );
                    trimesh->addObject( tface );

                    nodes[0] = f->getNodeAt(pos);
                    nodes[1] = f->getNodeAt(pos+2);
                    nodes[2] = f->getNodeAt(pos+3);
                    tface    = JTriangle::newObject( nodes );
                    trimesh->addObject( tface );
                } else {
                    pos = 0;
                    if( randomize ) {
                        int val = rand();
                        if( val%2 > 0 ) pos = 1;
                    }
                    nodes[0] = f->getNodeAt(pos+0);
                    nodes[1] = f->getNodeAt(pos+1);
                    nodes[2] = f->getNodeAt(pos+2);
                    tface    = JTriangle::newObject( nodes );
                    trimesh->addObject( tface );

                    nodes[0] = f->getNodeAt(pos+0);
                    nodes[1] = f->getNodeAt(pos+2);
                    nodes[2] = f->getNodeAt(pos+3);
                    tface    = JTriangle::newObject( nodes );
                    trimesh->addObject( tface );
                }
            }
        }
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = quadmesh->getFaceAt(i);
        if( f->getSize(0) == 4) f->setStatus(JMeshEntity::REMOVE);
    }

    trimesh->enumerate(0);
    trimesh->enumerate(2);

    return trimesh;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTriMeshGenerator :: getFromQuadMesh( const JMeshPtr &quadmesh, int trisPerQuad, bool randomize)
{
    JNodeSequence steiner;
    return AllTriMeshGenerator::getFromQuadMesh(quadmesh, steiner, trisPerQuad, randomize);
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTriMeshGenerator :: getFromPolyMesh( const JMeshPtr &polymesh, JNodeSequence &steiner)
{
    JMeshPtr trimesh;
    if( polymesh == nullptr) return trimesh;

    size_t numnodes = polymesh->getSize(0);
    size_t numfaces = polymesh->getSize(2);

    if( numfaces == 0) {
        cout << "Warning: the mesh is empty" << endl;
        return NULL;
    }

    trimesh = JMesh::newObject();
    trimesh->reserve( numnodes, 0);
    trimesh->reserve( numfaces, 2);

    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr v = polymesh->getNodeAt(i);
        if( v->isActive() ) trimesh->addObject(v);
    }
    steiner.clear();

    Point3D p3d;
    JNodeSequence nodes(3);
    JFacePtr  tface;

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = polymesh->getFaceAt(i);
        if( face->isActive() ) {
            int nnodes = face->getSize(0);
            if( nnodes == 3) {
                trimesh->addObject( face );
            } else {
                face->getAvgXYZ( p3d );
                JNodePtr v0 = JNode::newObject();
                v0->setXYZCoords( p3d );
                steiner.push_back(v0);
                trimesh->addObject(v0);
                nodes[0] = v0;
                for( int j = 0; j < nnodes; j++) {
                    nodes[1] = face->getNodeAt(j);
                    nodes[2] = face->getNodeAt((j+1)%nnodes);
                    tface  = JTriangle::newObject(nodes);
                    trimesh->addObject( tface );
                }
            }
        }
    }
    return trimesh;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr AllTriMeshGenerator :: getFromPolyMesh( const JMeshPtr &polymesh)
{
    JNodeSequence steiner;
    return AllTriMeshGenerator :: getFromPolyMesh(polymesh, steiner);
}
///////////////////////////////////////////////////////////////////////////////

void  AllTriMeshGenerator :: getCongruentMesh( const JMeshPtr &mesh)
{
    /*
        if( mesh == nullptr) return;

        mesh->setVisitBits(0, 0);
        size_t numfaces = mesh->getSize(2);

        size_t numfixed = 0;
        Point3D xyz;
        vector<double> faceangles;

        mesh->delete_node_attribute("Constraint");

        bool fixed  = 1;
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            JNodePtr vtx = mesh->getNodeAt(i);
            if( vtx->isActive() && vtx->isBoundary() ) vtx->setAttribute("Constraint", fixed);
        }

        JMeshNonlinearOptimization mopt;

        int numAdjusted = 0;
        for( int j = 0; j < 5; j++) {
            double angle = 60 + 5*j;
            double lowangle  = M_PI*angle/180.0;
            double highangle = M_PI*(angle+5)/180.0;

            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr face = mesh->getFaceAt(i);
                FaceGeometry::getAngles(face, faceangles, ANGLE_IN_RADIANS);
                double maxval = 0;
                int    maxpos = 0;
                for( int k = 0; k < 3; k++) {
                    if( faceangles[k] > maxval) {
                        maxval = faceangles[k];
                        maxpos = k;
                    }
                }
                JNodePtr vtx = face->getNodeAt(maxpos);
                bool visit  = vtx->getVisitBit();

                if( !vtx->isBoundary() && visit == 0 && maxval >= lowangle && maxval < highangle) {
                    xyz = TriGeometry::getIdealPosition(face, vtx, lowangle);
                    vtx->setXYZCoords(xyz);
                    face->getNodeAt(0)->setVisitBit(1);
                    face->getNodeAt(1)->setVisitBit(1);
                    face->getNodeAt(2)->setVisitBit(1);
                    numAdjusted++;
                }
            }

            for( size_t i = 0; i < numnodes; i++) {
                JNodePtr vtx = mesh->getNodeAt(i);
                if( vtx->isActive() && vtx->getVisitBit() ) vtx->setAttribute("Constraint", fixed);
            }

            mopt.setNumIterations(10);
            mopt.setMesh(mesh);
            mopt.improveShapes();
        }
    */

}
///////////////////////////////////////////////////////////////////////////////

/*
Mesh *
AllTriMeshGenerator:: getStructuredMesh(int nx, int ny)
{
    Mesh *trimesh = Mesh::newObject();

    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);

    Point3D xyz;

    int index = 0;
    for (int j = 0; j < ny; j++) {
	for (int i = 0; i < nx; i++) {
	    xyz[0] = -1.0 + i * dx;
	    xyz[1] = -1.0 + j * dy;
	    xyz[2] = 0.0;
	    Vertex *vnew = Vertex::newObject();
	    vnew->setID(index++);
	    vnew->setXYZCoords(xyz);
	    trimesh->addObject(vnew);
	}
    }

    JNodeSequence nodes(3);
    index = 0;
    Face *newtri;
    for (int j = 0; j < ny - 1; j++) {
	for (int i = 0; i < nx - 1; i++) {
	    int n0 = j * nx + i;
	    int n1 = n0 + 1;
	    int n2 = n1 + nx;
	    int n3 = n0 + nx;
	    nodes[0] = trimesh->getNodeAt(n0);
	    nodes[1] = trimesh->getNodeAt(n1);
	    nodes[2] = trimesh->getNodeAt(n2);
	    newtri = Triangle::newObject( nodes);
	    trimesh->addObject(newtri);

	    nodes[0] = trimesh->getNodeAt(n0);
	    nodes[1] = trimesh->getNodeAt(n2);
	    nodes[2] = trimesh->getNodeAt(n3);
	    newtri = Triangle::newObject( nodes );
	    trimesh->addObject(newtri);
	}
    }
    return trimesh;
}
*/

