#include "MeshRefine2D.hpp"

using namespace Jaal;

#ifdef CSV

void LoopSubDivision :: refineAll( Mesh *m)
{
    mesh = m;
    mesh->build_edges();

    mesh->buildRelations(1,2);
    apply_edge_rule();
}

void LoopSubDivision :: apply_edge_rule()
{
    size_t numedges = mesh->getSize(1);

    Point3D p3d;
    JFaceSequence efaces;
    Vertex *vf0 = NULL, *vf1 = NULL;

    for( size_t i = 0; i < numedges; i++)  {
        Edge  *e = mesh->getEdgeAt(i);
        if( e->isActive() ) {
            e->getRelations(efaces);
            Vertex *ve0 = e->getNodeAt(0);
            Vertex *ve1 = e->getNodeAt(1);
            const Point3D &p0 = ve0->getXYZCoords();
            const Point3D &p1 = ve1->getXYZCoords();
            p3d[0] = 0.50*( p0[0] + p1[0] );
            p3d[1] = 0.50*( p0[1] + p1[1] );
            p3d[2] = 0.50*( p0[2] + p1[2] );
            /*
                           if( efaces.size() == 2 ) {
                                efaces[0]->getAttribute("Steiner", vf0);
                                efaces[1]->getAttribute("Steiner", vf1);
                                const Point3D &p2  = vf0->getXYZCoords();
                                const Point3D &p3  = vf1->getXYZCoords();
                                p3d[0] = 0.25*( p0[0] + p1[0] + p2[0] + p3[0] );
                                p3d[1] = 0.25*( p0[1] + p1[1] + p2[1] + p3[1] );
                                p3d[2] = 0.25*( p0[2] + p1[2] + p2[2] + p3[2] );
                           } else {
                                p3d[0] = 0.50*( p0[0] + p1[0] );
                                p3d[1] = 0.50*( p0[1] + p1[1] );
                                p3d[2] = 0.50*( p0[2] + p1[2] );
                           }
            */
            Vertex *edgenode = Vertex::newObject();
            edgenode->setXYZCoords(p3d);
            e->setAttribute("Steiner", edgenode);
            mesh->addNode(edgenode);
        }
    }
    cout << "Edge points generated : " << endl;
}

void LoopSubDivision :: build_new_faces()
{
    int err;
    JNodeSequence tnodes(3), edgenodes(3);

#ifdef CSV

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        const JNodeSequence &fnodes =  face->getNodes();
        /*
                   for( int j = 0; j < face->getSize(1); j++) {
                        Edge *edge = face->getEdgeAt(j);
                        edge->getAttribute("Steiner", edgenode[j]);
                   }
        */
        Face *tri0 = Face::newObject();
        tnodes[0] = fnodes[0];
        tnodes[1] = edgenodes[2];
        tnodes[2] = edgenodes[1];
        err = mesh->addFace(tri0);

        Face *tri1 = Face::newObject();
        tnodes[0] = fnodes[1];
        tnodes[1] = edgenodes[0];
        tnodes[2] = edgenodes[2];
        err = mesh->addFace(tri1);

        Face *tri2 = Face::newObject();
        tnodes[0] = fnodes[2];
        tnodes[1] = edgenodes[1];
        tnodes[2] = edgenodes[0];
        err = mesh->addFace(tri2);

        Face *tri3 = Face::newObject();
        tnodes[0] = edgenodes[0];
        tnodes[1] = edgenodes[1];
        tnodes[2] = edgenodes[2];
        err = mesh->addFace(tri3);

        mesh->remove(face);
    }
#endif

}
#endif
