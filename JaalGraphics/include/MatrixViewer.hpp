#pragma once

#include "Mesh.hpp"

class JMatrixViewer : public JViewComponent
{
public:
    static const int  PRIMAL_GRAPH = 0;
    static const int  DUAL_GRAPH   = 1;

    static JViewComponentPtr registerComponent(JaalViewer *v);
 
    explicit JMatrixViewer( JaalViewer *p);

    // Set the primal mesh...
    void setMesh( JMeshPtr m )
    {
        mesh = m;
    }

    // What kind of mesh it is: Primal or Dual graph. Default value
    // is PRIMAL_MESH.
    void setGraphType( int t )
    {
        graphType = t;
    }

    // In the Dual mesh, do you wish to include boundary anchors. Remember
    // that all internal simplicies have at the most 2 lower entities as
    // neighbours. Since entities near the boundary may have only one
    // neighbour, so in many applications, we can dual edge from the
    // boundary simplex to the simplex on the boundary ( which we
    // call "Anchor" as they act like constrains when there are
    // geometric transformations ( for example: Laplacian smoothing).
    //
    // Default: No anchors.
    void setBoundaryAnchors( bool b )
    {
        bound_anchors = b;
    }

    // Do you wish to create a grid on the Matrix viewer. Sometimes
    // useful for small graphs.
    // Default: No grid...
    void setMatrixGrid( bool f )
    {
        matrix_grid = f;
    }

    // Set the font scale, to display the nodes iD on the top and
    // to the right of the matrix..
    // Default: 1.0
    void setFontsScale( double fs )
    {
        fontsScale = fs;
    }

    // For large graphs, displaying all the node Ids on the
    // top and to the rights is not essentials and may be
    // inefficient. We can skip nodes to display the IDs.
    // For example display:1, 10, 20, 30 .... ).
    // Default : 1.
    // Acceptable value : Any n >= 1.
    void setStepLabels(int  n)
    {
        stepLabels = n;
    }

    void setGlyph( int g )
    {
        glyph = g;
    }

    void setPointSize( double s)
    {
        pointSize = s;
    }

    // Virtual function: Draw the matrix.
    void draw();

private:
    JMeshPtr mesh;
    GLUquadricObj *disk;
    int    glyph = 0;
    bool   matrix_grid = 0;
    double dl;
    double fontsScale = 1.0;
    double pointSize = 1.0;

    bool  bound_anchors = 0;  // In a dual graph, we can put anchors on the boundary entities.
    int   graphType = 0;; // 0: Primary mesh  1: Dual Mesh.
    int   stepLabels = 1;

    void  drawPrimal();
    void  drawDual();
    void  drawEdge(const JEdgePtr &e);
    void  fillIJ(double px, double py);
    void  fillMatrix();
    void  drawIDs();
};

typedef boost::shared_ptr<JMatrixViewer> JMatrixViewerPtr;
