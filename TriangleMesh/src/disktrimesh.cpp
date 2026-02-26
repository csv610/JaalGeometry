#include "Mesh.hpp"
#include "MeshRefine2D.hpp"

using namespace Jaal;

int create_disk_triangle_mesh(const JEdgeSequence &eseq, int nlevel, JNodeSequence &newnodes, JFaceSequence &newfaces)
{
    /*
       map<Vertex*, JNodeSequence>   vmap;
       int nSize = eseq.size();
       for( int i = 0; i <  nSegments; i++) {
            Vertex *v0 = eseq[i]->getNodeAt(0);
            Vertex *v1 = eseq[i]->getNodeAt(1);
            vmap[v0].push_back(v1);
            vmap[v1].push_back(v0);
       }

       map<Vertex*, JNodeSequence> :: const_iterator it;
       vector<Vertex*> edgenodes;
       for( it = vmap.begin(); it != vmap.end(); ++it) {
            size_t n = it->second.size();
            if( n != 2 ) {
                cout << "Warning: Edge Sequence is not simple, triangulation not done" << endl;
                return 1;
            }
            edgenodes.push_back( it->first );
       }

       Point3D pCenter;
       pCenter[0] = 0.0;
       pCenter[1] = 0.0;
       pCenter[2] = 0.0;
       nSize = edgenodes.size();

       for( int i = 0; i < nSize; i++) {
            const Point3D &p3d = edgenodes[i]->getXYZCoords();
    	pCenter[0] += p3d[0];
            pCenter[1] += p3d[1];
            pCenter[2] += p3d[2];
       }
       pCenter[0] /= ( double) nSize;
       pCenter[1] /= ( double) nSize;
       pCenter[2] /= ( double) nSize;

       Point3D prel, pref;
       double radius, theta;
       for( int i = 0; i < nSize; i++) {
            Vertex *vertex = edgenodes[i];
            const Point3D  &p3d = vertex->getXYZCoords();
            prel[0]  = p3d[0] - pCenter[0];
            prel[1]  = p3d[1] - pCenter[1];
            prel[2]  = p3d[2] - pCenter[2];
            radius   = JMath::length( prel );
            if( i = 0 ) {
               pref[0]  = p3d[0] - pCenter[0];
               pref[1]  = p3d[1] - pCenter[1];
               pref[2]  = p3d[2] - pCenter[2];
               theta    = 0.0;
            } else {



            }

       }
    */

}


int main(int argc, char **argv)
{
    Jaal::Mesh *mesh = new Jaal::Mesh;
    Jaal::MeshOptimization mopt;

    int ntheta = 100;
    double dtheta = 2*M_PI/(ntheta);
    double radius = 1.0;

    Point3D p3d;
    Vertex *vertex = Vertex::newObject();
    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    vertex->setXYZCoords(p3d);
    vertex->setID(0);
    mesh->addNode(vertex);

    for( int i = 0; i < ntheta; i++) {
        vertex = Vertex::newObject();
        p3d[0] = radius*cos( i*dtheta);
        p3d[1] = radius*sin( i*dtheta);
        p3d[2] = 0.0;
        vertex->setXYZCoords(p3d);
        vertex->setID(i+1);
        mesh->addNode(vertex);
    }

    JNodeSequence conn(3);
    for( int i = 0; i < ntheta; i++) {
        conn[0] = mesh->getNodeAt(0);
        conn[1] = mesh->getNodeAt((i)%ntheta + 1);
        conn[2] = mesh->getNodeAt((i+1)%ntheta+ 1);
        Face *f = Face::newObject();
        f->setNodes(conn);
        mesh->addFace(f);
    }

    Sqrt3Refiner2D refiner;
    refiner.refineAll(mesh, 1);

    mesh->saveAs("tmp.off");

    return 0;
}

