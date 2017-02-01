#include "MeshViewer.hpp"


JMeshViewer :: JMeshViewer( JaalViewer *vm)
{
    viewManager = vm;
    init();
}

/////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::init()
{
    nodeDraw.reset( new JNodeDraw() );
    edgeDraw.reset( new JEdgeDraw() );
    faceDraw.reset( new JFaceDraw() );
    cellDraw.reset( new JCellDraw() );

    nodeDraw->setViewManager( viewManager );
    edgeDraw->setViewManager( viewManager );
    faceDraw->setViewManager( viewManager );
    cellDraw->setViewManager( viewManager );

    currCounter  = 0;
    display_boundary   = 0;
    display_dual_graph = 0;

    display_enclosure  = 0;
    enclosure_type     = AXIS_ALIGNED_BOUNDING_BOX;

    entityPicker.reset( new JMeshEntityPicker());
    assert( entityPicker );
}

///////////////////////////////////////////////////////////////////////////////
JMeshViewer :: ~JMeshViewer()
{
}

///////////////////////////////////////////////////////////////////////////////

JViewComponentPtr JMeshViewer :: registerComponent(JaalViewer *viewer)
{
    boost::shared_ptr<JMeshViewer> obj;
    obj.reset(new JMeshViewer(viewer));
    obj->setName("MeshViewer");
    viewer->attach(obj);
    return obj;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshViewer::  addObject( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return 1;

    if( find(meshdb.begin(), meshdb.end(), mesh) == meshdb.end() ) {
        initMesh(mesh);
        meshdb.push_back(mesh);
    }
    if( meshdb.size() == 1) lookAt(mesh);

    setCurrentMesh(mesh);
    viewManager->updateGL();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshViewer:: removeObject( const JMeshPtr &msh)
{
    vector<JMeshPtr>::iterator it;
    it = std::remove(meshdb.begin(), meshdb.end(), msh);
    meshdb.erase(it, meshdb.end());

    currMesh = nullptr;
    if( !meshdb.empty() ) currMesh = meshdb.back();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer:: lookAt( const JMeshPtr &msh)
{
    if( msh == nullptr) return;

    viewManager->resetView( JViewDirection:: FRONT_VIEW);

    JBoundingBox box = msh->getGeometry()->getBoundingBox();
    Point3D  pC   = box.getCenter();
    double   len  = box.getMaxLength();

    qglviewer::Vec  vec;
    vec[0] = pC[0];
    vec[1] = pC[1];
    vec[2] = pC[2] + 2.0*len;
    viewManager->camera()->setPosition(vec);

    vec[0] = pC[0];
    vec[1] = pC[1];
    vec[2] = pC[2];
    viewManager->camera()->lookAt(vec);
    viewManager->camera()->setSceneCenter(vec);

    vec[0] = 0.0;
    vec[1] = 1;
    vec[2] = 0;
    viewManager->camera()->setUpVector(vec);

    refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer::displayEnclosure(bool val, int t)
{
    display_enclosure = val;
    enclosure_type = t;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: refreshDisplay()
{
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshViewer :: loadMesh( const string &filename)
{
    static int nCounter = 1;
    if (filename.empty()) return nullptr;

    /*
        QtWaitingSpinnerWidget  spinner;
        spinner.setInnerRadius(100);
        spinner.start();
    */

    JMeshPtr mesh = readData(filename);
    size_t pos0 = filename.rfind("/");
    size_t pos1 = filename.rfind('.');
    string meshname = filename.substr(pos0+1, pos1-pos0-1);
    mesh->setName(meshname);
    nCounter++;

//  spinner.stop();
    return mesh;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: animate()
{
    currCounter++;
//  QGLViewer :: animate();
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr
JMeshViewer::readData(const string &fname)
{
    JMeshPtr mesh = JMeshIO::readFile( fname);
    initMesh( mesh );
    meshdb.push_back(mesh);
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

size_t JMeshViewer :: getNumVisible( const JMeshPtr &mesh, int entity )
{
    size_t nCount = 0;
    size_t nSize;

    JNodeRenderPtr vAttrib;
    JEdgeRenderPtr eAttrib;
    JFaceRenderPtr fAttrib;
    JCellRenderPtr cAttrib;

    switch( entity ) {
    case 0:
        nSize = mesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                int err = vertex->getAttribute("Render", vAttrib) ;
                if( !err && vAttrib->display ) nCount++;
            }
        }
        break;
    case 1:
        nSize = mesh->getSize(1);
        for( size_t i = 0; i < nSize; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                int err = edge->getAttribute("Render", eAttrib) ;
                if( !err && eAttrib->display ) nCount++;
            }
        }
        break;
    case 2:
        nSize = mesh->getSize(2);
        for( size_t i = 0; i < nSize; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                int err = face->getAttribute("Render", fAttrib) ;
                if( !err && fAttrib->display ) nCount++;
            }
        }
        break;
    case 3:
        nSize = mesh->getSize(3);
        for( size_t i = 0; i < nSize; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                int err = cell->getAttribute("Render", cAttrib) ;
                if( !err && cAttrib->display ) nCount++;
            }
        }
        break;
    }
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

size_t JMeshViewer :: getNumVisible( int entity )
{
    if( meshdb.empty() ) return 0;

    size_t nCount = 0;
    for( const JMeshPtr &mesh: meshdb)
        nCount += getNumVisible(mesh, entity);
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::drawWithNames()
{
    if( !isActive()  ) return;

    for( const JMeshPtr &mesh: meshdb ) {
        glPushMatrix();
        nodeDraw->withName(mesh);
        edgeDraw->withName(mesh);
        faceDraw->withName(mesh);
//      cellDraw->withName(mesh);
        glPopMatrix();
    }
}

//////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: actionMouseEvent(int id)
{
    if( entityPicker ) {
        if( id == 3) {    // Only for the mouse release event, an entity will be picked.
            entityPicker->select_entity();
        }
    }

    /*
        if( entityPicker ) {
            if( painter ) {
                entityPicker->select_entity();
    //          entityPicker->setMode(2);
                refreshDisplay();
            }
            if(id == 1) {    // Only for the mouse release event, an entity will be picked.
                entityPicker->select_entity();
            }
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

JBoundingBox JMeshViewer :: getAxisAlignedBox()
{
    JBoundingBox box;

    Point3D  minC, maxC;
    if( meshdb.empty() ) {
        minC[0] = -1.0;
        minC[1] = -1.0;
        minC[2] = -1.0;
        maxC[0] =  1.0;
        maxC[1] =  1.0;
        maxC[2] =  1.0;
        box.setPoints(minC, maxC);
        return box;
    }

    minC[0] = 0.0;
    minC[1] = 0.0;
    minC[2] = 0.0;

    maxC[0] = 0.0;
    maxC[1] = 0.0;
    maxC[2] = 0.0;

    box.setPoints(minC, maxC);
    for( const JMeshPtr &mesh : meshdb) {
        if( mesh->isActive() ) {
            JBoundingBox b = mesh->getGeometry()->getBoundingBox();
            box.setUnion(b);
        }
    }
    return box;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: drawEnclosure( const JMeshPtr &mesh)
{
    if( !display_enclosure) return;

    JBoundingBox box;
    JHexahedronPtr  minBox;
    glPushMatrix();

    glDisable(GL_LIGHTING);
    switch( enclosure_type ) {
    case AXIS_ALIGNED_BOUNDING_BOX:
        mesh->getAttribute("AxisBoundingBox", box);
        JaalViewer::drawBox(&box);
        break;
    case MINIMUM_BOUNDING_BOX:
        mesh->getAttribute("MinBoundingBox", minBox);
        JaalViewer::drawBox(minBox);
        break;
    }
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshViewer:: setNodesGlyph( int glyph, double r )
{
    if( r == 0) {
        if(glyph == JNodeRender::NODE_AS_SPHERE)
            r = nodeDraw->getBallRadius();
        else
            r = nodeDraw->getPointSize();
    }

    JNodeRenderPtr attrib;

    for( const JMeshPtr &mesh: meshdb) {
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i <  numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx  ) {
                int err = vtx->getAttribute("Render", attrib);
                assert( !err );
                if(glyph == JNodeRender::NODE_AS_SPHERE)
                    attrib->ballRadius = r;
                else
                    attrib->pointSize  = r;

                attrib->glyph  = glyph;

            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void JMeshViewer:: setEdgesGlyph(int glyph, double r)
{
    if( r == 0) {
        if( glyph == JEdgeDraw ::EDGE_AS_CYLINDER )
            r = edgeDraw->getCylinderRadius();
        else
            r = edgeDraw->getLineWidth();
    }

    JEdgeRenderPtr attrib;

    for( const JMeshPtr &mesh: meshdb) {
        size_t numedges = mesh->getSize(1);
        for( size_t i = 0; i <  numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge  ) {
                int err = edge->getAttribute("Render", attrib);
                assert( !err );
                if( glyph == JEdgeDraw::EDGE_AS_CYLINDER )
                    attrib->cylinderRadius = r;
                else
                    attrib->lineWidth = r;
                attrib->glyph  = glyph;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::initMesh( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    JMeshRenderPtr mrender(new JMeshRender);
    mesh->setAttribute("Render", mrender);
    mrender->pickableEntity = -1;

    mesh->getTopology()->searchBoundary();
    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);
    mesh->enumerate(3);

    JMeshQuality mq;
    mq.setMesh(mesh);
    JStatisticalInfo stat = mq.getEdgesLength();
    double avglen = 0.10*stat.getAverage();

    nodeDraw->setNormalsLength( avglen ) ;
    faceDraw->setNormalsLength( avglen ) ;
    nodeDraw->setBallRadius( avglen );
    edgeDraw->setCylinderRadius( avglen );

    initBuffers(mesh);

    if( entityPicker ) {
        entityPicker->setViewManager(viewManager);
        entityPicker->clearAll();
    }

    JBoundingBox box;
    mesh->getAttribute("AxisBoundingBox", box);
    viewManager->addBox(box);
/*
    JTransformNode tnode;
    mesh->setAttribute("TransformNode", tnode);
*/
}

////////////////////////////////////////////////////////////////////////////////
void JMeshViewer :: initBuffers( const JMeshPtr &mesh)
{
    updateGeometryBuffers( mesh );
    updateTopologyBuffers( mesh );
}
////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: updateGeometryBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    nodeDraw->initRenderAttributes(mesh);
    edgeDraw->initRenderAttributes(mesh);
    faceDraw->initRenderAttributes(mesh);
    cellDraw->initRenderAttributes(mesh);

    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    mesh->setAttribute("AxisBoundingBox", box);

    faceDraw->updateGeometryBuffers(mesh);
    edgeDraw->updateGeometryBuffers(mesh);
    nodeDraw->updateBuffers(mesh);

    refreshDisplay();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: updateTopologyBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    nodeDraw->initRenderAttributes(mesh);
    edgeDraw->initRenderAttributes(mesh);
    faceDraw->initRenderAttributes(mesh);

    faceDraw->updateTopologyBuffers(mesh);
    edgeDraw->updateTopologyBuffers(mesh);
    refreshDisplay();
}
////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::updateBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    cellDraw->initRenderAttributes(mesh);
    nodeDraw->initRenderAttributes(mesh);
    edgeDraw->initRenderAttributes(mesh);
    faceDraw->initRenderAttributes(mesh);

    cellDraw->updateBuffers(mesh);
    faceDraw->updateBuffers(mesh);
    edgeDraw->updateBuffers(mesh);
    nodeDraw->updateBuffers(mesh);

    refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::updateBuffers()
{
    for( JMeshPtr mesh: meshdb) updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::fastDraw()
{
    nodeDraw->setDisplay(1);
    edgeDraw->setDisplay(0);
    faceDraw->setDisplay(0);
    cellDraw->setDisplay(0);
    draw();
    nodeDraw->setDisplay(1);
    edgeDraw->setDisplay(1);
    faceDraw->setDisplay(1);
    cellDraw->setDisplay(1);
}

////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw()
{
    if( !isActive() )  return;

    glPushMatrix();

    if( entityPicker) entityPicker->draw();

    JMeshRenderPtr mrender;
    for( const JMeshPtr &mesh: meshdb) {
        mesh->getAttribute("Render", mrender);
        if(mrender->display) {
            drawEnclosure(mesh);
            if( mrender->displayEntity[3]) cellDraw->draw(mesh);
            if( mrender->displayEntity[2]) faceDraw->draw(mesh);
            if( mrender->displayEntity[1]) edgeDraw->draw(mesh);
            if( mrender->displayEntity[0]) nodeDraw->draw(mesh);
        }
    }
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: displayAll( const JMeshPtr &mesh, int entity, bool val)
{

    if( entity == 0) {
        JNodeRenderPtr nodeAttrib;
        size_t numNodes = mesh->getSize(0);
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                vertex->getAttribute("Render", nodeAttrib) ;
                if( nodeAttrib ) nodeAttrib->display = val;
            }
        }
        nodeDraw->updateBuffers(mesh);
    }

    if( entity == 1 ) {
        JEdgeRenderPtr eAttrib;
        size_t numEdges = mesh->getSize(1);
        for( size_t i = 0; i < numEdges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = val;
        }
        edgeDraw->updateBuffers(mesh);
    }

    if( entity == 2 ) {
        JFaceRenderPtr fAttrib;
        size_t numFaces = mesh->getSize(2);
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            face->getAttribute("Render", fAttrib);
            fAttrib->display = val;
        }
        faceDraw->updateBuffers(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshViewer :: displayAll( int entity, bool val)
{
    for( const JMeshPtr &mesh: meshdb)
        displayAll( mesh, entity, val);
}
*/
///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: lookAt( const JNodePtr &vtx)
{
    qglviewer::Vec vec;
    const Point3D xyz = vtx->getXYZCoords();
    vec[0] = xyz[0];
    vec[1] = xyz[1];
    vec[2] = xyz[2];
    viewManager->camera()->setSceneCenter(vec);

    /*
        Vec3F srcVec, dstVec, perpAxis;

        Vec3F normal;
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 1.0;
        vtx->getAttribute("Normal", normal);

        Vec vec;
        vec = viewManager->camera()->position();
        double dist = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
        vec[0] = dist*normal[0];
        vec[1] = dist*normal[1];
        vec[2] = dist*normal[2];

        viewManager->camera()->setPosition(vec);

        const Point3D &p3d = vtx->getXYZCoords();
        vec[0] = p3d[0];
        vec[1] = p3d[1];
        vec[2] = p3d[2];

        viewManager->camera()->lookAt( vec );
    */
    refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer:: alignAlong( const JMeshPtr &mesh, const Vec3F &srcVec, const Vec3F &dstVec)
{
    Vec3F  perpAxis;
    JMath::cross_product( dstVec, srcVec, perpAxis);
    double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
    qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );
    qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);
    qglviewer::Vec prot;

    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        Point3D p3d = v->getXYZCoords();
        prot[0] = p3d[0];
        prot[1] = p3d[1];
        prot[2] = p3d[2];
        prot    = quaternion.rotate(prot);
        p3d[0]  = prot[0];
        p3d[1]  = prot[1];
        p3d[2]  = prot[2];
        v->setXYZCoords(p3d);
    }
}

/////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: alignAlong( const JMeshPtr &mesh, const JEdgePtr &currEdge, int along)
{
    if( currEdge == nullptr ) return;

    JMeshAffineTransform af;
    af.setMesh(mesh);

    const Point3D &p1 = currEdge->getNodeAt(0)->getXYZCoords();
    const Point3D &p2 = currEdge->getNodeAt(1)->getXYZCoords();
    af.translate(-p1[0], -p1[1], -p1[2] );

    // Specify the destination vector ...
    Vec3F dstVec;
    dstVec[0] = 0.0;
    dstVec[1] = 0.0;
    dstVec[2] = 0.0;
    dstVec[along] = 1.0;

    // Where is the vector now ...
    Vec3F currVec;
    currVec[0] = p2[0] - p1[0];
    currVec[1] = p2[1] - p1[1];
    currVec[2] = p2[2] - p1[2];

    // Let the Quaternion rotate the "Current Vector" to "Destination Vector"
    alignAlong(mesh, currVec, dstVec);
    refreshDisplay();
}

/////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: alignAlong( const JMeshPtr &mesh, const JFacePtr &currFace, int along )
{
    if( currFace == nullptr ) return;

    Point3D p3d;
    JFaceGeometry::getCentroid(currFace, p3d);

    qglviewer::Vec pos(p3d[0], p3d[1], p3d[2] );
    viewManager->camera()->setSceneCenter(pos);

/*
    Vec3F normal;
    currFace->getAttribute("Normal", normal);
    qglviewer::Vec viewdir(-normal[0], -normal[1], -normal[2] );
    viewManager->camera()->setViewDirection(viewdir);


    const Point3D &p1 = currFace->getNodeAt(0)->getXYZCoords();
    const Point3D &p2 = currFace->getNodeAt(1)->getXYZCoords();
    const Point3D &p3 = currFace->getNodeAt(2)->getXYZCoords();

    JMeshAffineTransform af(mesh);
    af.translate(-p1[0], -p1[1], -p1[2] );

    // Specify the destination vector ...
    Vec3F dstVec;
    dstVec[0] = 0.0;
    dstVec[1] = 0.0;
    dstVec[2] = 0.0;
    dstVec[along] = 1.0;

    // Where is the vector now ...
    Vec3F vec1, vec2, currVec;

    vec1[0] = p2[0] - p1[0];
    vec1[1] = p2[1] - p1[1];
    vec1[2] = p2[2] - p1[2];

    vec2[0] = p3[0] - p1[0];
    vec2[1] = p3[1] - p1[1];
    vec2[2] = p3[2] - p1[2];

    JMath::cross_product( vec1, vec2, currVec);

    // Let the Quaternion rotate the "Current Vector" to "Destination Vector"
    alignAlong( mesh, currVec, dstVec);

*/
    refreshDisplay();
}

void JMeshViewer :: setRenderMode( int m)
{
    /*
        switch( m ) {
        case JRenderMode::POINTCLOUD:
            displayAll( 0, 1);
            displayAll( 1, 0);
            displayAll( 2, 0);
            displayAll( 3, 0);
            setNodesGlyph( JNodeDraw::NODE_AS_POINT, 1.0);
            break;
        case JRenderMode::WIREFRAME:
            displayAll( 0, 0);
            displayAll( 1, 1);
            displayAll( 2, 0);
            displayAll( 3, 0);
            break;
        case JRenderMode::FLAT_SHADE:
            faceDraw->setRenderMode(m);
            displayAll( 0, 0);
            displayAll( 1, 0);
            displayAll( 2, 1);
            displayAll( 3, 0);
            break;
        case JRenderMode::SMOOTH_SHADE:
            faceDraw->setRenderMode(m);
            displayAll( 0, 0);
            displayAll( 1, 0);
            displayAll( 2, 1);
            displayAll( 3, 0);
            break;
        case JRenderMode::HIDDENLINES:
            faceDraw->setRenderMode(m);
            displayAll( 0, 0);
            displayAll( 1, 0);
            displayAll( 2, 1);
            displayAll( 3, 0);
            break;
        case JRenderMode::SURFACE_AND_WIREFRAME:
            faceDraw->setRenderMode(m);
            displayAll( 0, 0);
            displayAll( 1, 1);
            displayAll( 2, 1);
            displayAll( 3, 0);
            break;
        }
    */
}



