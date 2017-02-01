#pragma once

#include "Mesh.hpp"
#include "JaalViewer.hpp"

using namespace Jaal;

class JPickingTexture
{
public:
    void   setViewManager( JaalViewer *v) {
        viewer = v;
    }
    void   genTexture( const JMeshPtr &m);

    void   setPickingShape( int p );
    void   setPickingRegionSize( double p);

    JFaceSequence getPickedFaces(const Point2I &p);

private:
    JaalViewer *viewer;
    JMeshPtr   mesh;
    int brushShape;
    int brushSize;
    vector<unsigned char> pixelData;

    // Initialize the texture buffer ....
    void    init();

    // get an unique color of each face of the mesh.
    void    getColor( size_t  id, Point3I &p);

    // from the  picked color, get the ID ....
    size_t  getID( const Point3I &data) const;

    size_t readPixel( const Point2I &p);
    void   getBuffer();
};

