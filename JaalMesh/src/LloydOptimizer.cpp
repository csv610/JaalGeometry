#include "LloydOptimizer.hpp"

///////////////////////////////////////////////////////////////////////////////

int JLloydMeshOptimizer :: getVoroCentroid(const JNodePtr &vertex, Point3D &pcenter)
{
    if( vertex == nullptr ) return 1;
    if( !vertex->isActive() ) return 2;
    if( vertex->isBoundary()) return 3;
    if( vertex->hasAttribute("Constraint") ) return 4;

    JNode::getRelations(vertex, neighedges);

    int numedges =  neighedges.size();

    Point3D p0, p1, p2, xyz;

    p0 = vertex->getXYZCoords();

    pcenter[0] =  0.0;
    pcenter[1] =  0.0;
    pcenter[2] =  0.0;
    double area, sumarea = 0.0;

    JEdgePtr dualedge = nullptr;
    for( int i = 0; i < numedges; i++) {
        int err = neighedges[i]->getAttribute("DualEdge", dualedge);
        assert(!err);
        JNodePtr v1 = dualedge->getNodeAt(0);
        JNodePtr v2 = dualedge->getNodeAt(1);
        p1 = v1->getXYZCoords();
        p2 = v2->getXYZCoords();
        area = fabs(JTriGeometry::getArea(p0, p1, p2));
        JFaceGeometry::getCentroid( vertex, v1, v2, xyz);
        pcenter[0] += area*xyz[0];
        pcenter[1] += area*xyz[1];
        pcenter[2] += area*xyz[2];
        sumarea    += area;
    }

    pcenter[0] /= sumarea;
    pcenter[1] /= sumarea;
    pcenter[2] /= sumarea;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JLloydMeshOptimizer :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    JMeshDualGraph  jd;
    jd.setMesh(mesh);
    jd.setBoundaryNodes(1);
    dualGraph  = jd.getGraph();
}

///////////////////////////////////////////////////////////////////////////////

int JLloydMeshOptimizer :: updateDual()
{
    JNodePtr dnode;
    Point3D   xyz;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr  &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("DualNode", dnode);
            if( !err) {
                JFaceGeometry::getDualPosition(face,xyz);
                dnode->setXYZCoords(xyz);
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JLloydMeshOptimizer :: updatePrimal()
{
    JFaceSequence faceneighs;
    JNodePtr  dualnode;

    Point3D  xyz, newCoord;
    maxResidue = 0;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &primalnode = mesh->getNodeAt(i);
        if( preserve_boundary == 1 && primalnode->isBoundary() ) {
            continue;
        }
        JNode::getRelations(primalnode, faceneighs);
        xyz[0] = 0.0;
        xyz[1] = 0.0;
        xyz[2] = 0.0;
        for( const JFacePtr &face : faceneighs) {
            face->getAttribute("DualNode", dualnode);
            const Point3D &pj = dualnode->getXYZCoords();
            xyz[0] += pj[0];
            xyz[1] += pj[1];
            xyz[2] += pj[2];
        }
        newCoord[0] = xyz[0]/(double)faceneighs.size();
        newCoord[1] = xyz[1]/(double)faceneighs.size();
        newCoord[2] = xyz[2]/(double)faceneighs.size();
        maxResidue    = max( maxResidue, fabs(newCoord[0] - xyz[0]));
        maxResidue    = max( maxResidue, fabs(newCoord[1] - xyz[1]));
        maxResidue    = max( maxResidue, fabs(newCoord[2] - xyz[2]));
        xyz[0]      = newCoord[0];
        xyz[1]      = newCoord[1];
        xyz[2]      = newCoord[2];
        primalnode->setXYZCoords(xyz);
    }

    if( maxResidue < tolerance) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JLloydMeshOptimizer :: smoothAll()
{
    if( mesh == nullptr) return 1;

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    int nUsed = 0;
    for( int i = 0; i < numIters; i++) {
        nUsed++;
        updateDual();
        int status = updatePrimal();
        if( status == 0) break;
    }
    numIters = nUsed;

    return 2;
}

///////////////////////////////////////////////////////////////////////////////
