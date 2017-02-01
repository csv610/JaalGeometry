#pragma once

#include "Mesh.hpp"

class JSimpleLaplacianSkeleton
{
    struct MyNode
    {
        double lambda;
        double area;
        double maxAngle;
        double distance;
        Point3D  xyz;
        Point3F  normal;
        vector<size_t> adjNodes;
        vector<size_t> adjFaces;
    };

    struct MyFace
    {
        double  area;
        Point3F normal;
        size_t  nodes[3];
    };

    struct MyMesh
    {
        vector<MyFace> faces;
        vector<MyNode> nodes;
    };

    MyMesh  mesh;
    vector<Point3D> newCoords;
    vector<double>  normArea;

    void setFaceArea(size_t id) {
        size_t nodes[3];
        nodes[0] = mesh.faces[id].nodes[0];
        nodes[1] = mesh.faces[id].nodes[1];
        nodes[2] = mesh.faces[id].nodes[2];

        double abc[3];
        for( int i = 0; i < 3; i++) {
            size_t v1 = nodes[(i+1)%3];
            size_t v2 = nodes[(i+2)%3];
            double dx = mesh.nodes[v2].xyz[0] - mesh.nodes[v1].xyz[0];
            double dy = mesh.nodes[v2].xyz[1] - mesh.nodes[v1].xyz[1];
            double dz = mesh.nodes[v2].xyz[2] - mesh.nodes[v1].xyz[2];
            abc[i] = sqrt(dx*dx + dy*dy + dz*dz);
        }
        double a = abc[0];
        double b = abc[1];
        double c = abc[2];
        double s = 0.5*(a+b+c);

        mesh.faces[id].area = sqrt(s*(s-a)*(s-b)*(s-c));
    }

    void setNodeArea( size_t vi)
    {
        int nSize = mesh.nodes[vi].adjFaces.size();
        double sumArea = 0.0;
        for( int j = 0; j < nSize; j++) {
            int fd = mesh.nodes[vi].adjFaces[j];
            sumArea += mesh.faces[fd].area;
        }
        mesh.nodes[vi].area = sumArea/(double)nSize;
    }

    void setFaceNormal( size_t id) {
        size_t n0 = mesh.faces[id].nodes[0];
        size_t n1 = mesh.faces[id].nodes[1];
        size_t n2 = mesh.faces[id].nodes[2];

        double a[3];
        a[0] = mesh.nodes[n1].xyz[0] - mesh.nodes[n0].xyz[0];
        a[1] = mesh.nodes[n1].xyz[1] - mesh.nodes[n0].xyz[1];
        a[2] = mesh.nodes[n1].xyz[2] - mesh.nodes[n0].xyz[2];

        double b[3];
        b[0] = mesh.nodes[n2].xyz[0] - mesh.nodes[n0].xyz[0];
        b[1] = mesh.nodes[n2].xyz[1] - mesh.nodes[n0].xyz[1];
        b[2] = mesh.nodes[n2].xyz[2] - mesh.nodes[n0].xyz[2];

        double nx = a[1]*b[2] - a[2]*b[1];
        double ny = a[2]*b[0] - a[0]*b[2];
        double nz = a[0]*b[1] - a[1]*b[0];
        double nl = sqrt(nx*nx + ny*ny + nz*nz);
        mesh.faces[id].normal[0] = nx/nl;
        mesh.faces[id].normal[1] = ny/nl;
        mesh.faces[id].normal[2] = nz/nl;
    }

    void setNodeNormal(size_t vi)
    {
        int nSize = mesh.nodes[vi].adjFaces.size();
        double nx = 0.0;
        double ny = 0.0;
        double nz = 0.0;
        for( int j = 0; j < nSize; j++) {
            int fd = mesh.nodes[vi].adjFaces[j];
            nx += mesh.faces[fd].normal[0];
            ny += mesh.faces[fd].normal[1];
            nz += mesh.faces[fd].normal[2];
        }

        double nl = sqrt(nx*nx + ny*ny + nz*nz);
        mesh.nodes[vi].normal[0] = nx/nl;
        mesh.nodes[vi].normal[1] = ny/nl;
        mesh.nodes[vi].normal[2] = nz/nl;
    }

    void setAngles( size_t id )
    {
        size_t nodes[3];
        nodes[0] = mesh.faces[id].nodes[0];
        nodes[1] = mesh.faces[id].nodes[1];
        nodes[2] = mesh.faces[id].nodes[2];

        double abc[3];
        for( int i = 0; i < 3; i++) {
            size_t v1 = nodes[(i+1)%3];
            size_t v2 = nodes[(i+2)%3];
            double dx = mesh.nodes[v2].xyz[0] - mesh.nodes[v1].xyz[0];
            double dy = mesh.nodes[v2].xyz[1] - mesh.nodes[v1].xyz[1];
            double dz = mesh.nodes[v2].xyz[2] - mesh.nodes[v1].xyz[2];
            abc[i] = sqrt(dx*dx + dy*dy + dz*dz);
        }
        double a = abc[0];
        double b = abc[1];
        double c = abc[2];

        double cosA = (b*b + c*c - a*a)/(2.0*b*c);
        double cosB = (a*a + c*c - b*b)/(2.0*a*c);
        double cosC = (a*a + b*b - c*c)/(2.0*a*b);

        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        double A = 180.0*acos(cosA)/M_PI;
        mesh.nodes[0].maxAngle = max(mesh.nodes[0].maxAngle, A);

        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        double B = 180.0*acos(cosB)/M_PI;
        mesh.nodes[1].maxAngle = max(mesh.nodes[1].maxAngle, B);

        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        double C = 180.0*acos(cosC)/M_PI;
        mesh.nodes[2].maxAngle = max(mesh.nodes[2].maxAngle, C);
    }

    void setNodesAngle()
    {
      size_t numNodes = mesh.nodes.size();
      for( size_t i = 0; i < numNodes; i++) 
           mesh.nodes[i].maxAngle = 0.0;

      size_t numFaces = mesh.faces.size();
      for( size_t i = 0; i < numFaces; i++) 
           setAngles(i);
    }


public:
    void setMesh (const JMeshPtr &m);
    JMeshPtr getWorkingMesh() {
        return jmesh;
    }

    void applyOneStep();

private:
    JMeshPtr jmesh;
    void     setNodesNormal();
    void     setNodesArea();
    void     atomicOp( size_t v);
};
