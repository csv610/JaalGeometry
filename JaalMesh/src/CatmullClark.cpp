#include "MeshRefine.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void CatmullClarkSubDivision :: apply_face_rule()
{
    Point3D p3d;
    size_t numfaces = inmesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr f = inmesh->getFaceAt(i);
        if( f->isActive() ) {
            f->getAvgXYZ(p3d);
            JNodePtr facenode = JVertex::newObject();
            facenode->setXYZCoords(p3d);
            f->setAttribute("Steiner", facenode);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void CatmullClarkSubDivision :: apply_edge_rule()
{

    size_t numedges = inmesh->getSize(1);
    Point3D p3d;
    JFaceSequence efaces;
    JNodePtr vf0, vf1;

    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr e = inmesh->getEdgeAt(i);
        if( e->isActive() ) {
            JEdge::getRelations(e, efaces);
            JNodePtr ve0 = e->getNodeAt(0);
            JNodePtr ve1 = e->getNodeAt(1);
            const Point3D &p0 = ve0->getXYZCoords();
            const Point3D &p1 = ve1->getXYZCoords();
            if( efaces.size() == 2 ) {
                efaces[0]->getAttribute("Steiner", vf0);
                efaces[1]->getAttribute("Steiner", vf1);
                const Point3D &p2 = vf0->getXYZCoords();
                const Point3D &p3 = vf1->getXYZCoords();
                p3d[0] = 0.25*( p0[0] + p1[0] + p2[0] + p3[0] );
                p3d[1] = 0.25*( p0[1] + p1[1] + p2[1] + p3[1] );
                p3d[2] = 0.25*( p0[2] + p1[2] + p2[2] + p3[2] );
            } else {
                p3d[0] = 0.50*( p0[0] + p1[0] );
                p3d[1] = 0.50*( p0[1] + p1[1] );
                p3d[2] = 0.50*( p0[2] + p1[2] );
            }
            JNodePtr edgenode = JVertex::newObject();
            edgenode->setXYZCoords(p3d);
            e->setAttribute("Steiner", edgenode);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void  CatmullClarkSubDivision :: apply_node_rule()
{
    size_t numnodes = inmesh->getSize(0);
    vfaces.clear();
    vedges.clear();
    JNodePtr vf;

    Point3D  Q, R, newxyz, p3d, p0, p1;
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = inmesh->getNodeAt(i);
        if( vtx->isActive() ) {
            JVertex::getRelations(vtx, vfaces);
            int  n = vfaces.size();
            Q[0] = 0.0;
            Q[1] = 0.0;
            Q[2] = 0.0;
            for( int j = 0; j < n; j++) {
                vfaces[j]->getAttribute("Steiner", vf);
                p3d = vf->getXYZCoords();
                Q[0] += p3d[0];
                Q[1] += p3d[1];
                Q[2] += p3d[2];
            }
            Q[0] /= (double)n;
            Q[1] /= (double)n;
            Q[2] /= (double)n;


            JVertex::getRelations(vtx, vedges);
            R[0] = 0.0;
            R[1] = 0.0;
            R[2] = 0.0;
            assert( vedges.size() == (size_t) n);
            for( int  j = 0; j < n; j++) {
                p0 = vedges[j]->getNodeAt(0)->getXYZCoords();
                p1 = vedges[j]->getNodeAt(1)->getXYZCoords();
                R[0] += 0.5*(p0[0] + p1[0]);
                R[1] += 0.5*(p0[1] + p1[1]);
                R[2] += 0.5*(p0[2] + p1[2]);
            }
            R[0] /= (double)n;
            R[1] /= (double)n;
            R[2] /= (double)n;
            const Point3D &vp = vtx->getXYZCoords();
            newxyz[0] = (Q[0] + 2.0*R[0] + (n-3)*vp[0])/(double)n;
            newxyz[1] = (Q[1] + 2.0*R[1] + (n-3)*vp[1])/(double)n;
            newxyz[2] = (Q[2] + 2.0*R[2] + (n-3)*vp[2])/(double)n;
            vtx->setXYZCoords(newxyz);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void  CatmullClarkSubDivision :: limitPoint(JNodePtr vertex )
{
    if( !vertex->isActive() ) return;

    Point3D  Q, R, newxyz, p3d, p0, p1;
    JVertex::getRelations(vertex, vfaces);
    int  n = vfaces.size();
    Q[0] = 0.0;
    Q[1] = 0.0;
    Q[2] = 0.0;
    JNodePtr vf;
    for( int j = 0; j < n; j++) {
        vfaces[j]->getAttribute("Steiner", vf);
        p3d = vf->getXYZCoords();
        Q[0] += p3d[0];
        Q[1] += p3d[1];
        Q[2] += p3d[2];
    }
    JVertex::getRelations(vertex, vedges);
    R[0] = 0.0;
    R[1] = 0.0;
    R[2] = 0.0;
    assert( vedges.size() == (size_t) n);
    for( int  j = 0; j < n; j++) {
        p0 = vedges[j]->getNodeAt(0)->getXYZCoords();
        p1 = vedges[j]->getNodeAt(1)->getXYZCoords();
        R[0] += 0.5*(p0[0] + p1[0]);
        R[1] += 0.5*(p0[1] + p1[1]);
        R[2] += 0.5*(p0[2] + p1[2]);
    }
    const Point3D &vp = vertex->getXYZCoords();
    newxyz[0] = (n*n*vp[0] + 4*R[0] + Q[0])/(double)(n*(n+5));
    newxyz[1] = (n*n*vp[1] + 4*R[1] + Q[1])/(double)(n*(n+5));
    newxyz[2] = (n*n*vp[2] + 4*R[2] + Q[2])/(double)(n*(n+5));
    vertex->setXYZCoords(newxyz);
}


///////////////////////////////////////////////////////////////////////////////
void  CatmullClarkSubDivision :: build_newmesh()
{
    /*
        JNodePtr vertex;
        JNodeSequence newnodes;
        JFaceSequence newfaces;

        size_t numnodes = inmesh->getSize(0);
        size_t numedges = inmesh->getSize(1);
        size_t numfaces = inmesh->getSize(2);

        outmesh->reserve(numnodes+numedges+numfaces, 0);
        outmesh->reserve(4*numfaces, 2);

        for( size_t i = 0; i < numedges; i++) {
            JEdgePtr edge = mesh->getEdgeAt(i);
            edge->getAttribute("Steiner", vertex);
            outmesh->addObject(vertex);
        }

        for( size_t i = 0; i < numfaces; i++) {
            JFacePtr face = mesh->getFaceAt(i);
            face->getAttribute("Steiner", vertex);
            outmesh->addObject(vertex);
        }

        TriRefiner  triRefiner;
        QuadRefiner quadRefiner;

        triRefiner.setMesh(mesh);
        quadRefiner.setMesh(mesh);

        for( size_t iface = 0; iface < numfaces; iface++) {
            JFacePtr face = mesh->getFaceAt(iface);
            int numnodes = face->getSize(0);
            switch(numnodes) {
            case 3:
               triRefiner.tri2quads(face);
                break;
            case 4:
                quadRefiner.refine4(face);
                break;
            default:
                cout << "Error: Poly->Quad refiner not implemented" << endl;
                exit(0);
            }
        }

        for( size_t i = 0; i < numedges; i++) {
            JEdgePtr edge = mesh->getEdgeAt(i);
            edge->setStatus(JMeshEntity::REMOVE);
        }

        mesh->delete_face_attribute( "Steiner");
        mesh->delete_edge_attribute( "Steiner");

        mesh->clearRelations(0,1);
        mesh->clearRelations(1,2);
        mesh->clearRelations(0,2);

    //  mesh->collect_garbage_faces();
    //  mesh->collect_garbage_edges();

        assert(mesh->getSize(1) == 0);
    */
}

///////////////////////////////////////////////////////////////////////////////

void CatmullClarkSubDivision :: limitSurface()
{
    if( inmesh == NULL ) return;
    /*
        mesh = refineAll(m, 1);

        mesh->buildRelations(0,1);
        mesh->buildRelations(1,2);
        mesh->buildRelations(0,2);

        apply_face_rule();
        apply_edge_rule();

        // Limit points only for the original nodes..
        size_t numNodes = orgmesh->getSize(0);
        for( size_t i  = 0; i < numNodes; i++) {
            Vertex *src = orgmesh->getNodeAt(i);
            Vertex *dst = mesh->getNodeAt(i);
            limitPoint( dst );
    	src->setXYZCoords( dst->getXYZCoords() );
        }
        mesh->deleteAll();
        delete mesh;
    */
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr CatmullClarkSubDivision :: refineAll(int numIter)
{

    if( inmesh == NULL ) return NULL;

    /*
        orgmesh = mesh;
        mesh = orgmesh->deepCopy();

        for( int i = 0; i < numIter; i++)
        {
            mesh->getTopology()->collect_edges();
            mesh->buildRelations(0,1);
            mesh->buildRelations(1,2);
            mesh->buildRelations(0,2);
            apply_face_rule();
            apply_edge_rule();
            apply_node_rule();
            build_newmesh();
        }
        mesh->enumerate(0);
        return mesh;
    */
}

///////////////////////////////////////////////////////////////////////////////

void CatmullClarkSubDivision :: smoothSurface(int numIter)
{
    if( inmesh == NULL ) return;

    /*
        orgmesh = NULL;
        size_t numNodes = m->getSize(0);
        size_t numEdges = m->getSize(1);
        size_t numFaces = m->getSize(2);

        refineAll(m, numIter);

        for( size_t i  = 0; i < numNodes; i++) {
            Vertex *src = orgmesh->getNodeAt(i);
            Vertex *dst = mesh->getNodeAt(i);
    	src->setXYZCoords( dst->getXYZCoords() );
        }

        mesh->deleteFaces();
        mesh->deleteEdges();
        mesh->deleteNodes();

        assert( m->getSize(0) == numNodes);
        assert( m->getSize(1) == numEdges);
        assert( m->getSize(2) == numFaces);
    */
}

///////////////////////////////////////////////////////////////////////////////
