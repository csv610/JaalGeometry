#include "MeshViewer.hpp"

JNodePickColor ::  JNodePickColor()
{
    highlightColor[0]  = 1.0;
    highlightColor[1]  = 0.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = alpha;
}

///////////////////////////////////////////////////////////////////////////////
int JNodePickColor :: assign(const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;

    JNodeRenderPtr nAttrib;
    int err = vertex->getAttribute("Render", nAttrib);
    if( err ) return 1;
    nAttrib->color = highlightColor;
    if( nAttrib->glyph == 0)
        nAttrib->scale = 2.0;
    else
        nAttrib->scale = 1.1;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JEdgePickColor ::  JEdgePickColor()
{
    highlightColor[0]  = 1.0;
    highlightColor[1]  = 0.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = alpha;
}

///////////////////////////////////////////////////////////////////////////////
int JEdgePickColor :: assign( const JEdgePtr &edge)
{
    if( edge == nullptr ) return 1;

    JEdgeRenderPtr eAttrib;
    int err = edge->getAttribute("Render", eAttrib);
    if( err ) return 1;
    eAttrib->color = highlightColor;
    eAttrib->scale = 2.0;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JFacePickColor ::  JFacePickColor()
{
    highlightColor[0]  = 1.0;
    highlightColor[1]  = 0.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = alpha;
}
///////////////////////////////////////////////////////////////////////////////

int  JFacePickColor :: assign( const JFacePtr &face)
{
    JFaceRenderPtr fAttrib;
    int err = face->getAttribute("Render", fAttrib);
    if( err ) return 1;
    fAttrib->color = highlightColor;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JCellPickColor :: JCellPickColor()
{
    highlightColor[0]  = 1.0;
    highlightColor[1]  = 0.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = 0.0;
//  faceColor = new JFacePickColor();
}
///////////////////////////////////////////////////////////////////////////////

int JCellPickColor :: assign(const JCellPtr &cell)
{
    if( cell ) {
        cell->setAttribute("Color", highlightColor);
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

JMeshEntityPicker :: JMeshEntityPicker()
{
    pickOp      =  0;
    active      =  0;
    selection_mode   = SINGLE_SELECTION;
    pickNodeColor.reset(new JNodePickColor);
    pickEdgeColor.reset(new JEdgePickColor);
    pickFaceColor.reset(new JFacePickColor);
    pickCellColor.reset(new JCellPickColor);
    meshViewer = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
void JMeshEntityPicker :: setViewManager (JaalViewer *v)
{
    viewManager = v;
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    clearAll();
}
////////////////////////////////////////////////////////////////////////////////

void JMeshEntityPicker :: setHighlightColor( const JColor &color, int e)
{

    switch(e)
    {
    case 0:
        assert( pickNodeColor );
        pickNodeColor->setHighlightColor(color);
        break;
    case 1:
        assert( pickEdgeColor );
        pickEdgeColor->setHighlightColor(color);
        break;
    case 2:
        assert( pickFaceColor );
        pickFaceColor->setHighlightColor(color);
        break;
    case 3:
        assert( pickCellColor );
        pickCellColor->setHighlightColor(color);
        break;
    }
    viewManager->refreshDisplay();
}
////////////////////////////////////////////////////////////////////////////////

void JMeshEntityPicker :: clearAll()
{
    picked_nodes.clear();
    picked_edges.clear();
    picked_faces.clear();
    picked_cells.clear();

    if( viewManager) {
        viewManager->setSelectedName(-1);
        viewManager->refreshDisplay();
    }
}

////////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshEntityPicker :: getPickedNodes() const
{
    JNodeSequence  seq;
    if( mesh == nullptr) {
        cout << "Warning: No mesh attached for picking nodes " << endl;
        return seq;
    }

    for( size_t vid : picked_nodes)
    {
        const JNodePtr &node = mesh->getNodeAt( vid );
        if( node->isActive() ) seq.push_back(node);
    }
    return seq;
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMeshEntityPicker :: getPickedEdges() const
{
    JEdgeSequence  seq;
    if( mesh == nullptr) {
        cout << "Warning: No mesh attached for picking edges " << endl;
        return seq;
    }

    int edgeID;
    for( size_t edgeID: picked_edges)
    {
        const JEdgePtr &edge = mesh->getEdgeAt( edgeID );
        if( edge->isActive() ) seq.push_back(edge);
    }
    return seq;

}
///////////////////////////////////////////////////////////////////////////////

JFaceSequence JMeshEntityPicker :: getPickedFaces() const
{
    JFaceSequence  seq;
    if( mesh == nullptr) {
        cout << "Warning: No mesh attached for picking faces " << endl;
        return seq;
    }

    for( size_t faceID : picked_faces)
    {
        const JFacePtr &face = mesh->getFaceAt( faceID );
        if( face->isActive() ) seq.push_back(face);
    }
    return seq;
}
///////////////////////////////////////////////////////////////////////////////

JCellSequence JMeshEntityPicker :: getPickedCells() const
{
    JCellSequence  seq;
    if( mesh == nullptr)  {
        cout << "Warning: No mesh attached for picking faces " << endl;
        return seq;
    }

    for( size_t cellID: picked_cells)
    {
        const JCellPtr &cell = mesh->getCellAt( cellID );
        if( cell->isActive() ) seq.push_back(cell);
    }
    return seq;
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshEntityPicker ::drawPickRegion()
{
    if( !active ) return;

    viewManager->startScreenCoordinatesSystem();
    int width  = viewManager->selectRegionWidth();
    int height = viewManager->selectRegionHeight();

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);

    Point2I currpos =  viewManager->getMouseCurrentPixelPosition();
    /*
      glColor4f(0.0, 0.0, 0.3f, 0.3f);
      glBegin(GL_QUADS);
      glVertex2i(rectangle_.left(),  rectangle_.top());
      glVertex2i(rectangle_.right(), rectangle_.top());
      glVertex2i(rectangle_.right(), rectangle_.bottom());
      glVertex2i(rectangle_.left(),  rectangle_.bottom());
      glEnd();
    */

    glLineWidth(2.0);
    glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
    glBegin(GL_LINE_LOOP);
    glVertex2i( currpos[0]-0.5*width,  currpos[1]-0.5*height);
    glVertex2i( currpos[0]+0.5*width,  currpos[1]-0.5*height);
    glVertex2i( currpos[0]+0.5*width,  currpos[1]+0.5*height);
    glVertex2i( currpos[0]-0.5*width,  currpos[1]+0.5*height);
    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    viewManager->stopScreenCoordinatesSystem();
}

////////////////////////////////////////////////////////////////////////////////////

void
JMeshEntityPicker ::draw()
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    set<int>::const_iterator siter;
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    int pick_entity =  mrender->pickableEntity;

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked nodes ...
    //////////////////////////////////////////////////////////////////////////////
    if( pick_entity == 0)
    {
        JNodeDraw  *drawnode  = meshViewer->getNodeDraw();
        for (siter = picked_nodes.begin(); siter != picked_nodes.end(); ++siter)
            drawnode->draw(mesh->getNodeAt(*siter));
    }

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked edges ...
    //////////////////////////////////////////////////////////////////////////////

    if( pick_entity == 1 )
    {
        JEdgeDraw  *drawedge  = meshViewer->getEdgeDraw();
        for (siter = picked_edges.begin(); siter != picked_edges.end(); ++siter)
            drawedge->draw(mesh->getEdgeAt(*siter));
    }

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked faces ...
    //////////////////////////////////////////////////////////////////////////////

    if( pick_entity == 2 )
    {
        JFaceDraw  *drawface  = meshViewer->getFaceDraw();
        for (siter = picked_faces.begin(); siter != picked_faces.end(); ++siter)
            drawface->draw(mesh->getFaceAt(*siter));
    }

    if( pick_entity == 3 )
    {
        JCellDraw  *drawcell  = meshViewer->getCellDraw();
        for (siter = picked_cells.begin(); siter != picked_cells.end(); ++siter)
            drawcell->draw(mesh->getCellAt(*siter));
    }

    drawPickRegion();
}

////////////////////////////////////////////////////////////////////////////////

int JMeshEntityPicker::select_entity()
{
    if( meshViewer == nullptr || mesh == nullptr ) return 2;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    int pick_entity =  mrender->pickableEntity;

    int id = meshViewer->getViewManager()->selectedName();

    if( id < 0 ) return 3;

    if (pick_entity == 0)
    {
        if( id < (int)mesh->getSize(0) )
        {
            JNodeRenderPtr nAttrib;
            int nodeid = id;
            if( pickOp == 0)
            {
                if( selection_mode == SINGLE_SELECTION )
                {
                    picked_nodes.clear();
                }
                picked_nodes.insert(nodeid);
                currPickedNode = mesh->getNodeAt( nodeid );
                pickNodeColor->assign(currPickedNode);
            }
            else
            {
                picked_nodes.erase(nodeid);
            }
        }
        return 0;
    }

    if (pick_entity == 1)
    {
        if( id < (int) mesh->getSize(1) )
        {
            JEdgeRenderPtr eAttrib;
            int edgeid = id;
            if( pickOp == 0 )
            {
                if( selection_mode == SINGLE_SELECTION )
                {
                    picked_edges.clear();
                }
                picked_edges.insert(edgeid);
                currPickedEdge = mesh->getEdgeAt( edgeid );
                pickEdgeColor->assign(currPickedEdge);
            }
            else
            {
                picked_edges.erase(edgeid);
            }
        }
        return 0;
    }

    if (pick_entity == 2)
    {
        if( id < (int) mesh->getSize(2) )
        {
            int faceid = id;
            if( pickOp == 0 )
            {
                if( selection_mode == SINGLE_SELECTION )
                {
                    picked_faces.clear();
                }
                picked_faces.insert(faceid);
                currPickedFace = mesh->getFaceAt( faceid );
                pickFaceColor->assign(currPickedFace);
            }
            else
            {
                picked_faces.erase(faceid);
            }
        }
        return 0;
    }

    if (pick_entity == 3)
    {
        if( id < (int) mesh->getSize(3) )
        {
            int cellid = id;
            if( pickOp == 0 )
            {
                if( selection_mode == SINGLE_SELECTION )
                {
                    picked_cells.clear();
                }
                picked_cells.insert(cellid);
                currPickedCell = mesh->getCellAt( cellid );
                pickCellColor->assign(currPickedCell);
            }
            else
            {
                picked_cells.erase(cellid);
            }
        }
        return 0;
    }
    return 6;
}

///////////////////////////////////////////////////////////////////////////////
