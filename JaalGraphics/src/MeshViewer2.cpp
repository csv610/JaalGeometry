#include "MeshViewer2.h"

size_t Mesh :: globalID = 1;

MeshViewer :: MeshViewer( QWidget *parent) : QGLViewer(parent)
{ }

///////////////////////////////////////////////////////////////////////////////
void MeshViewer:: draw( Edge *edge)
{
    Vertex *v0 = edge->connect[0];
    Vertex *v1 = edge->connect[1];
    glVertex3f( v0->coords[0], v0->coords[1], v0->coords[2] );
    glVertex3f( v1->coords[0], v1->coords[1], v1->coords[2] );
}
///////////////////////////////////////////////////////////////////////////////

void MeshViewer:: draw( Face *face)
{
    int numNodes = face->connect.size();
    for( int i = 0; i < numNodes; i++) {
        Vertex *v = face->connect[i];
        glVertex3f( v->coords[0], v->coords[1], v->coords[2] );
    }
}
///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::drawMeshNodes(Mesh *mesh)
{
    glColor3fv( &nodesColor[0] );

    if( node_attrib_color ) {
        for( int i = 0; i < mesh->nodes.size(); i++) {
            Vertex *v = mesh->nodes[i];
            v->rgba[0] = colormap[v->groupID][0];
            v->rgba[1] = colormap[v->groupID][1];
            v->rgba[2] = colormap[v->groupID][2];
            v->rgba[3] = alpha;
        }
    }

    if( low_complexity) {
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

        glEnable( GL_POINT_SMOOTH );

        glPointSize(pointSize);

        glBegin( GL_POINTS);
        for( int i = 0; i < mesh->nodes.size(); i++) {
            Vertex *v = mesh->nodes[i];
            if( node_attrib_color) glColor4fv( &v->rgba[0] );
            glVertex3f( v->coords[0], v->coords[1], v->coords[2] );
        }
        glEnd();

        glDisable( GL_POINT_SMOOTH );
        glDisable( GL_BLEND );
        return;
    }

    glEnable( GL_LIGHTING );
    gluQuadricDrawStyle(sphere, GLU_FILL);
    for( int i = 0; i < mesh->nodes.size(); i++) {
        Vertex *v = mesh->nodes[i];
        glPushMatrix();
        glTranslatef( v->coords[0], v->coords[1], v->coords[2] );
        gluSphere( sphere, 0.005, 20, 20);
        glPopMatrix();
    }
    glDisable( GL_LIGHTING );
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::drawMeshFaces( Mesh *mesh, bool filled )
{
    if( filled )
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    else
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

    if( face_attrib_color ) {
        for( int i = 0; i < mesh->faces.size(); i++) {
            Face *f = mesh->faces[i];
            f->rgba[0] = colormap[f->groupID][0];
            f->rgba[1] = colormap[f->groupID][1];
            f->rgba[2] = colormap[f->groupID][2];
            f->rgba[3] = alpha;
        }
    }

    glBegin( GL_TRIANGLES);
    for( int i = 0; i < mesh->faces.size(); i++) {
        Face *f = mesh->faces[i];
        int numNodes = f->connect.size();
        if( face_attrib_color) glColor4fv( &f->rgba[0] );
        if( numNodes == 3 ) draw( f  );
    }
    glEnd();

    glBegin( GL_QUADS);
    for( int i = 0; i < mesh->faces.size(); i++) {
        Face *f = mesh->faces[i];
        int numNodes = f->connect.size();
        if( face_attrib_color) glColor4fv( &f->rgba[0] );
        if( numNodes == 4 ) draw( f );
    }
    glEnd();

    glBegin( GL_POLYGON);
    for( int i = 0; i < mesh->faces.size(); i++) {
        Face *f = mesh->faces[i];
        int numNodes = f->connect.size();
        if( face_attrib_color) glColor4fv( &f->rgba[0] );
        if( numNodes > 4 ) draw( f );
    }
    glEnd();

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::drawMeshEdges( Mesh *mesh )
{
    if( mesh->edges.empty() ) mesh->build_edges();

    glBegin( GL_LINES);
    for( int i = 0; i < mesh->edges.size(); i++) {
        draw( mesh->edges[i] );
    }
    glEnd();
}
///////////////////////////////////////////////////////////////////////////////

void MeshViewer :: normalize()
{
    if( meshes.empty() ) return;

    Vertex *v = meshes[0]->nodes[0];

    double xmin = v->coords[0];
    double xmax = v->coords[0];
    double ymin = v->coords[1];
    double ymax = v->coords[1];
    double zmin = v->coords[2];
    double zmax = v->coords[2];

    for( int im = 0; im < meshes.size(); im++) {
        int numNodes = meshes[im]->nodes.size();
        for( int j = 0; j < numNodes; j++) {
            v = meshes[im]->nodes[j];
            xmin = min( xmin, v->coords[0] );
            xmax = max( xmax, v->coords[0] );

            ymin = min( ymin, v->coords[1] );
            ymax = max( ymax, v->coords[1] );

            zmin = min( zmin, v->coords[2] );
            zmax = max( zmax, v->coords[2] );
        }
    }

    double xlen = fabs( xmax - xmin );
    double ylen = fabs( ymax - ymin );
    double zlen = fabs( zmax - zmin );

    double maxlen = max( xlen, max( ylen, zlen ));

    for( int im = 0; im < meshes.size(); im++) {
        int numNodes = meshes[im]->nodes.size();
        for( int j = 0; j < numNodes; j++) {
            v = meshes[im]->nodes[j];
            v->coords[0] /= maxlen;
            v->coords[1] /= maxlen;
            v->coords[2] /= maxlen;
        }
    }

    xmin = v->coords[0];
    xmax = v->coords[0];
    ymin = v->coords[1];
    ymax = v->coords[1];
    zmin = v->coords[2];
    zmax = v->coords[2];
    for( int im = 0; im < meshes.size(); im++) {
        int numNodes = meshes[im]->nodes.size();
        for( int j = 0; j < numNodes; j++) {
            v = meshes[im]->nodes[j];
            xmin = min( xmin, v->coords[0] );
            xmax = max( xmax, v->coords[0] );

            ymin = min( ymin, v->coords[1] );
            ymax = max( ymax, v->coords[1] );

            zmin = min( zmin, v->coords[2] );
            zmax = max( zmax, v->coords[2] );
        }
    }

    box = BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
}

///////////////////////////////////////////////////////////////////////////////

void MeshViewer :: drawIDs()
{
    char number[128];
    double x, y, z;
    Point3D xyz;

    if (display_face_ids ) {
        glColor3f(0.0, 0.0, 1.0);
        for( int k = 0; k < meshes.size(); k++) {
            Mesh *mesh = meshes[k];
            for (size_t i = 0; i < mesh->getSize(2); i++) {
                x = 0.0;
                y = 0.0;
                z = 0.0;
                Face *face = mesh->getFaceAt(i);
                int nsize = face->getSize(0);
                for (int j = 0; j < nsize; j++) {
                    Vertex *n0 = face->getNodeAt(j);
                    xyz = n0->getXYZCoords();
                    x += xyz[0];
                    y += xyz[1];
                    z += xyz[2];
                }
                x /= (double) nsize;
                y /= (double) nsize;
                z /= (double) nsize;
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glPushMatrix();
                glTranslatef(x, y, z + 0.00001);
                glScalef(fontscale, fontscale, fontscale);
                sprintf(number, "%ld", i);
                font->Render(number);
                glPopMatrix();
            }
        }
    }

    if (display_node_ids ) {
        glColor3f(0.0, 1.0, 0.0);
        for( int k = 0; k < meshes.size(); k++) {
            Mesh *mesh = meshes[k];
            for (int i = 0; i < mesh->getSize(0); i++) {
                Vertex *vtx = mesh->getNodeAt(i);
                xyz = vtx->getXYZCoords();
                x = xyz[0];
                y = xyz[1];
                z = xyz[2];
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glPushMatrix();
                glTranslatef(x, y, z + 0.001);
                glScalef(fontscale, fontscale, fontscale);
                sprintf(number, "%d", i);
                font->Render(number);
                glPopMatrix();
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::drawNormals()
{ }

///////////////////////////////////////////////////////////////////////////////

Vertex* MeshViewer::search_node_picked( size_t id)
{
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->nodes.size(); i++) {
            Vertex *v = mesh->nodes[i];
            if( v->globalID == id ) return v;
        }
    }
    return nullptr;
}
///////////////////////////////////////////////////////////////////////////////

Edge* MeshViewer::search_edge_picked( size_t id)
{
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->edges.size(); i++) {
            Edge *e = mesh->edges[i];
            if( e->globalID == id ) return e;
        }
    }
    return nullptr;
}
///////////////////////////////////////////////////////////////////////////////

Face* MeshViewer::search_face_picked( size_t id)
{
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->faces.size(); i++) {
            Face *f = mesh->faces[i];
            if( f->globalID == id ) return f;
        }
    }
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

void MeshViewer :: drawGroupFaces( int gid )
{
    vector<Face*> gfaces;
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        mesh->getGroup( gid, gfaces);

        glBegin( GL_TRIANGLES);
        for( int i = 0; i < gfaces.size(); i++) {
            int numNodes = gfaces[i]->connect.size();
            if( numNodes == 3 ) draw( gfaces[i] );
        }
        glEnd();

        glBegin( GL_QUADS);
        for( int i = 0; i < gfaces.size(); i++) {
            int numNodes = gfaces[i]->connect.size();
            if( numNodes == 4 ) draw( gfaces[i] );
        }
        glEnd();

        glBegin( GL_POLYGON);
        for( int i = 0; i < gfaces.size(); i++) {
            int numNodes = gfaces[i]->connect.size();
            if( numNodes > 4 ) draw( gfaces[i] );
        }
        glEnd();
    }
}

void
MeshViewer::drawWithNames()
{
    glColor3fv( &nodesColor[0] );
    set<int>::const_iterator it;

    if( pickEntity == 0) {
        glPointSize(pointSize);
        for( int k = 0; k < meshes.size(); k++) {
            Mesh *mesh = meshes[k];
            for( int i = 0; i < mesh->nodes.size(); i++) {
                Vertex *v = mesh->nodes[i];
                glPushMatrix();
                glPushName(v->globalID);
                glBegin( GL_POINTS);
                glVertex3f( v->coords[0], v->coords[1], v->coords[2] );
                glEnd();
                glPopName();
                glPopMatrix();
            }
        }
    }


    if( pickEntity == 1 ) {

        glColor3fv( &edgesColor[0] );

        for( int k = 0; k < meshes.size(); k++) {
            Mesh *mesh = meshes[k];
            for( int i = 0; i < mesh->edges.size(); i++) {
                Edge *e = mesh->edges[i];
                glPushMatrix();
                glPushName(e->globalID);
                glBegin( GL_LINES);
                draw( e );
                glEnd();
                glPopName();
                glPopMatrix();
            }
        }
    }

    if( pickEntity == 2 ) {

        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        glColor3fv( &facesColor[0] );

        if( pickGeometry ) {
            for( it = faceGroups.begin(); it != faceGroups.end(); ++it) {
                glPushMatrix();
                glPushName(*it);
                drawGroupFaces( *it );
                glPopName();
                glPopMatrix();
            }
        }

        if( pickMesh) {
            for( int k = 0; k < meshes.size(); k++) {
                Mesh *mesh = meshes[k];
                for( int i = 0; i < mesh->faces.size(); i++) {
                    Face *f = mesh->faces[i];
                    int numNodes = f->connect.size();
                    if( numNodes == 3 ) {
                        glPushMatrix();
                        glPushName(f->globalID);
                        glBegin( GL_TRIANGLES);
                        draw( f );
                        glEnd();
                        glPopName();
                        glPopMatrix();
                    }
                }
            }

            for( int k = 0; k < meshes.size(); k++) {
                Mesh *mesh = meshes[k];
                for( int i = 0; i < mesh->faces.size(); i++) {
                    Face *f = mesh->faces[i];
                    int numNodes = f->connect.size();
                    if( numNodes == 4 ) {
                        glPushMatrix();
                        glPushName(f->globalID);
                        glBegin( GL_QUADS);
                        draw( f );
                        glEnd();
                        glPopName();
                        glPopMatrix();
                    }
                }
            }

            for( int k = 0; k < meshes.size(); k++) {
                Mesh *mesh = meshes[k];
                for( int i = 0; i < mesh->faces.size(); i++) {
                    Face *f = mesh->faces[i];
                    int numNodes = f->connect.size();
                    if( numNodes > 4 ) {
                        glPushMatrix();
                        glPushName(f->globalID);
                        glBegin( GL_POLYGON);
                        draw( f );
                        glEnd();
                        glPopName();
                        glPopMatrix();
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void MeshViewer :: drawPickedObjects()
{
    glColor3fv( &pickedColor[0] );

    int id = selectedName();
    if( id >= 0) {
        if( pickMesh ) {
            if( pickEntity == 0) {
                Vertex *vertex_picked = search_node_picked(id);
                if( vertex_picked ) {
                    glPointSize(10.0);
                    glBegin( GL_POINTS);
                    glVertex3f( vertex_picked->coords[0],
                                vertex_picked->coords[1],
                                vertex_picked->coords[2] );
                    glEnd();
                    glPointSize(1.0);
                }
            }

            if( pickEntity == 1) {
                Edge *edge_picked = search_edge_picked(id);
                if( edge_picked ) {
                    glLineWidth( 2.0);
                    glBegin(GL_LINES);
                    draw( edge_picked );
                    glEnd();
                    glLineWidth( 1.0);
                }
            }

            if( pickEntity == 2) {
                Face *face_picked = search_face_picked(id);
                if( face_picked ) {
                    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
                    glBegin(GL_POLYGON);
                    draw( face_picked );
                    glEnd();
                }
            }
        }

        if( pickGeometry ) {
            if( pickEntity == 2) drawGroupFaces(id);
        }
    }

}

void
MeshViewer::drawObjects()
{
    Point3D xyz;
    qglviewer::Vec wc, sc;

    glClearColor(backGroundColor[0], backGroundColor[1], backGroundColor[2], 0.0);

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    drawPickedObjects();


    if( display_box ) drawBoundingBox();

    if( !hidden_line_removal )  {
        if( display_surface ) {
            glColor3fv( &facesColor[0] );
            glEnable (GL_BLEND);
            glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            for( int i = 0; i < meshes.size(); i++)
                drawMeshFaces( meshes[i], 1);
            glDisable (GL_BLEND);
        }

        if( display_edges )  {
            glColor3fv( &edgesColor[0] );
            glLineWidth( lineWidth );
            for( int i = 0; i < meshes.size(); i++)
                drawMeshEdges( meshes[i] );
        }
    } else {
        for( int i = 0; i < meshes.size(); i++) {
            glEnable(GL_DEPTH_TEST);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glColor3fv( &foreGroundColor[0] );
            drawMeshFaces( meshes[i], 0);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);
            glColor3fv( &backGroundColor[0] );
            drawMeshFaces( meshes[i], 1);
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }

//  if( display_ids  ) drawIDs();
    if( display_normals ) drawNormals();

    if (display_nodes) {
        for( int i = 0; i < meshes.size(); i++)
            drawMeshNodes(meshes[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////

void MeshViewer :: draw()
{
    drawObjects();
}

///////////////////////////////////////////////////////////////////////////////

Mesh* MeshViewer::readOffData( const string &fname)
{
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() ) {
        cout << "Warning: Cann't load the file " << fname << endl;
        return nullptr;
    }
    string str;
    infile >> str;
    if( str != "OFF") {
        cout << "Warning: File is not in off format " << endl;
        return nullptr;
    }
    int numNodes, numFaces, numEdges;
    infile >> numNodes >> numFaces >> numEdges;

    Mesh *mesh = new Mesh;
    mesh->nodes.reserve( numNodes );
    mesh->faces.reserve( numFaces );

    Point3D p;
    for( size_t i = 0; i < numNodes; i++) {
        infile >> p[0] >> p[1] >> p[2];
        Vertex *v = new Vertex;
        v->globalID = Mesh::globalID++;
        v->coords = p;
        mesh->nodes.push_back(v);
    }

    int nnodes, vid;
    vector<Vertex*> fconn;
    for( size_t i = 0; i < numFaces; i++) {
        infile >> nnodes;
        fconn.resize(nnodes);
        for( int j = 0; j < nnodes; j++) {
            infile >> vid;
            fconn[j] =  mesh->nodes[vid];
        }
        Face *f = new Face;
        f->globalID = Mesh::globalID++;
        f->connect = fconn;
        mesh->faces.push_back(f);
    }

    infile >> str;

    if( str == "#NODE_GROUP") {
        for( size_t i = 0; i < numNodes; i++)
            infile >> mesh->nodes[i]->groupID;
        infile >> str;
    }

    if( str == "#FACE_GROUP") {
        for( size_t i = 0; i < numFaces; i++)
            infile >> mesh->faces[i]->groupID;
        infile >> str;
    }

    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

int MeshViewer :: init_mesh()
{
    normalize();

    nodeGroups.clear();
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->nodes.size(); i++) {
            Vertex *v = mesh->nodes[i];
            nodeGroups.insert(v->groupID );
        }
    }
    int numNodeGroups = nodeGroups.size();

    edgeGroups.clear();
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->edges.size(); i++) {
            Edge *e = mesh->edges[i];
            edgeGroups.insert(e->groupID );
        }
    }
    int numEdgeGroups = edgeGroups.size();

    faceGroups.clear();
    for( int k = 0; k < meshes.size(); k++) {
        Mesh *mesh = meshes[k];
        for( int i = 0; i < mesh->faces.size(); i++) {
            Face *f = mesh->faces[i];
            faceGroups.insert(f->groupID );
        }
    }
    int numFaceGroups = faceGroups.size();

    int n = max( numNodeGroups, max( numEdgeGroups, numFaceGroups ));
    create_random_colormap( n );

}

///////////////////////////////////////////////////////////////////////////////

int MeshViewer::readMeshFile( const string &fname)
{
    Mesh *m = readOffData(fname);
    if(m) meshes.push_back(m);

    init_mesh();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::init()
{
    int type = 3;
    if (type < 3) {
        camera()->setPosition(Vec((type == 0) ? 1.0 : 0.0,
                                  (type == 1) ? 1.0 : 0.0,
                                  (type == 2) ? 1.0 : 0.0));
        camera()->lookAt(sceneCenter());
        camera()->showEntireScene();

        // Forbid rotation
        WorldConstraint* constraint = new WorldConstraint();
        constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
        camera()->frame()->setConstraint(constraint);
    }

    font = new FTPolygonFont(FONT_FILE);
    assert(!font->Error());
    font->FaceSize(1);
    font->CharMap(ft_encoding_unicode);
    fontscale = 0.0010;

    backGroundColor[0] = WhiteColor[0]/255.0;
    backGroundColor[1] = WhiteColor[1]/255.0;
    backGroundColor[2] = WhiteColor[2]/255.0;

    foreGroundColor[0] = BlackColor[0]/255.0;
    foreGroundColor[1] = BlackColor[1]/255.0;
    foreGroundColor[2] = BlackColor[2]/255.0;

    nodesColor[0] = RedColor[0]/255.0;
    nodesColor[1] = RedColor[1]/255.0;
    nodesColor[2] = RedColor[2]/255.0;

    edgesColor[0] = BlueColor[0]/255.0;
    edgesColor[1] = BlueColor[1]/255.0;
    edgesColor[2] = BlueColor[2]/255.0;

    facesColor[0] = YellowColor[0]/255.0;
    facesColor[1] = YellowColor[1]/255.0;
    facesColor[2] = YellowColor[2]/255.0;

    pickedColor[0] = GreenColor[0]/255.0;
    pickedColor[1] = GreenColor[1]/255.0;
    pickedColor[2] = GreenColor[2]/255.0;

    glDisable( GL_LIGHTING );

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

    pointSize = 5.0;
    lineWidth = 1.0;
    display_nodes = 0;
    enable_lighting = 1;
    enumerate = 1;
    pickEntity = 2;
    pickGeometry = 1;
    pickMesh     = 0;
    display_surface = 1;
    display_edges   = 1;
    hidden_line_removal = 0;
    low_complexity = 1;
    node_attrib_color = 1;
    edge_attrib_color = 1;
    face_attrib_color = 1;
    alpha = 0.90;

    sphere = gluNewQuadric();
    restoreStateFromFile ();

#ifdef CSV
    wire_frame = 0;
    bounding_box = 0;
    faces_normal = 0;
    normal_sign = 1;
    enumerate = 0;
    use_color = 1;
    background = 0;
    pick_entity = 0;
    modify_node_position = 0;
    flip_edge = 0;
    swapedge = nullptr;
    display_nodes = 0;
    mesh_backup = nullptr;
    use_material = 1;
    color_filled_model = 1;
    delete_last_selection = 0;
    defect_remesh_step = 0;
    patch = nullptr;

//  camera()->setType(Camera::ORTHOGRAPHIC);

    /*
      camera ()->setZNearCoefficient (0.001);
      camera ()->setSceneRadius (0.01);
      restoreStateFromFile ();
     */


    face_numbering = 0;
    node_numbering = 0;

    if (!mesh_filename.empty())
        readData(mesh_filename);


    /*
        if (!bsurf_filename.empty())
        {
            ifstream infile(bsurf_filename.c_str(), ios::in);
            if (!infile.fail())
                infile >> bsurf;
        }
    */

    init_mesh();


    float ambient[] = {0.19125, 0.0735, 0.0225, 1.0};
    float diffuse[] = {0.7038, 0.27048, 0.0828, 1.0};
    float specular[] = {1.0, 1.0, 1.0, 1.0};
    float shininess = 120;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    glEnable(GL_NORMALIZE);
    if( enable_lighting );
    glEnable(GL_LIGHTING);


    /*
    glEnable( GL_DEPTH_TEST );

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        GLfloat light0_position[] = { 1.0, 0.2, 1.0, 0.0 };
        GLfloat light1_position[] = { 0.0, 0.0, 0.0, 1.0 };
        GLfloat light1_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat light1_specular[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat lm_ambient[] = { 0.2, 0.2, 0.2, 1.0 };

        glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
        glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
        glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_diffuse);
        glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
        glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.2);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient);

        glEnable( GL_LIGHT0 );
        glEnable(GL_LIGHT1);

        GLfloat localview = 0.0;
        glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, localview);
     */

    help();
#endif
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::drawBoundingBox()
{
    static GLuint dlistID = 0;

    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glLineWidth(2.0);
    glColor3f(1.0, 1.0, 1.0);

    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    }

    dlistID = glGenLists(1) + 1;
    glNewList( dlistID, GL_COMPILE_AND_EXECUTE);

    const Point3D &lower = box.getLowerLeftCorner();
    const Point3D &upper = box.getUpperRightCorner();

    float xmin = lower[0];
    float ymin = lower[1];
    float zmin = lower[2];

    float xmax = upper[0];
    float ymax = upper[1];
    float zmax = upper[2];

    glBegin(GL_LINES);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmax, ymin, zmin);

    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymax, zmin);

    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmin, ymax, zmin);

    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymin, zmin);

    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmax, ymin, zmax);

    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymax, zmax);

    glVertex3f(xmax, ymax, zmax);
    glVertex3f(xmin, ymax, zmax);

    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmin, ymin, zmax);

    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymin, zmax);

    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymin, zmax);

    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymax, zmax);

    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymax, zmax);
    glEnd();

    /*
        if( enable_lighting )
            glEnable(GL_LIGHTING);
    */
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );
    glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::create_random_colormap(size_t n)
{
    colormap.clear();

    float r, g, b;
    Point3F  rgb;

    rgb[0] = 1.0;
    rgb[1] = 0.0;
    rgb[2] = 0.0;
    colormap[0] = rgb;

    rgb[0] = 0.0;
    rgb[1] = 1.0;
    rgb[2] = 0.0;
    colormap[1] = rgb;

    rgb[0] = 0.0;
    rgb[1] = 0.0;
    rgb[2] = 1.0;
    colormap[2] = rgb;

    rgb[0] = 1.0;
    rgb[1] = 0.0;
    rgb[2] = 1.0;
    colormap[3] = rgb;

    for (size_t i = 4; i < n; i++) {
        while (1) {
            r = drand48();
            if (r > 0.2) break;
        }
        while (1) {
            g = drand48();
            if (g > 0.2) break;
        }

        while (1) {
            b = drand48();
            if (b > 0.2) break;
        }

        rgb[0] = r;
        rgb[1] = g;
        rgb[2] = b;
        colormap[i] = rgb;
    }
}

#ifdef CSV
///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::draw_faces_normal()
{
    size_t numfaces = mesh->getSize(2);
    double x, y, z;

    glColor3f(0.0, 0.0, 1.0);

    glLineWidth(2.0);

    Vec3D normal;

    glBegin(GL_LINES);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        Point3D centroid = face->getCentroid();

        x = centroid[0];
        y = centroid[1];
        z = centroid[2];
        glVertex3f(x, y, z);

        normal = face->getNormal();
        x += 0.01 * normal_sign * normal[0];
        y += 0.01 * normal_sign * normal[1];
        z += 0.01 * normal_sign * normal[2];
        glVertex3f(x, y, z);
    }
    glEnd();

    glPointSize(2);
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        Point3D centroid = face->getCentroid();

        x = centroid[0];
        y = centroid[1];
        z = centroid[2];
        glVertex3f(x, y, z);
    }
    glEnd();
    glPointSize(1);
    glLineWidth(1.0);

}

///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::draw_faces(bool fillcolor)
{
    int itag = 0;
    size_t numfaces = mesh->getSize(2);
    Point3D xyz;

    glShadeModel(GL_FLAT);

    if (!use_material)
        glDisable(GL_LIGHTING);

    if (color_filled_model) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        if( enable_lighting) {
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
        }
    } else
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Vec3D normal;

    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        if (nnodes == 3) {
            if (fillcolor) {

                if (use_color) {
//                  itag = face->getTag();
                    itag = face->getGroupID();
                    Color clr = colormap[itag];
                    glColor3f(clr.rgb[0], clr.rgb[1], clr.rgb[2]);
                } else {
                    glColor3f(0.5, 0.5, 0.5);
                }
            }

            if (use_material) {
                normal = face->getNormal();
                float nx = 1.0 * normal_sign * normal[0];
                float ny = 1.0 * normal_sign * normal[1];
                float nz = 1.0 * normal_sign * normal[2];
                glNormal3f(nx, ny, nz);
            }

            for (int j = 0; j < 3; j++) {
                Vertex *vtx = face->getNodeAt(j);
                xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
        }

    }
    glEnd();

    glBegin(GL_QUADS);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        if (nnodes == 4) {
            if (fillcolor) {
                if (use_color) {
//                    itag = face->getTag();
                    itag = face->getGroupID();
                    assert(itag >= 0);
                    Color clr = colormap[itag];
                    glColor3f(clr.rgb[0], clr.rgb[1], clr.rgb[2]);
                } else
                    glColor3f(0.5, 0.5, 0.5);
            }

            for (int j = 0; j < 4; j++) {
                Vertex *vtx = face->getNodeAt((j+1)%4);  // Workaround for concave quads.
                xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
        }
    }
    glEnd();

    glBegin(GL_POLYGON);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        if (nnodes > 4) {
            for (int j = 0; j < nnodes; j++) {
                Vertex *vtx = face->getNodeAt(j);
                xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
        }
    }
    glEnd();
}


///////////////////////////////////////////////////////////////////////////////

void
MeshViewer::draw_picked_entities()
{
    Point3D xyz;
    vector<int> ::const_iterator siter;
    glDisable(GL_LIGHTING);
    glColor3f(0.0, 0.0, 1.0);

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked nodes ...
    //////////////////////////////////////////////////////////////////////////////
    glPointSize(5);
    glBegin(GL_POINTS);
    for (siter = picked_nodes.begin(); siter != picked_nodes.end(); ++siter) {
        Vertex *vtx = mesh->getNodeAt(*siter);
        xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
    glPointSize(1);

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked edges ...
    //////////////////////////////////////////////////////////////////////////////

    glColor3f(0.0, 1.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glLineWidth(2);
    glBegin(GL_LINES);
    for (siter = picked_edges.begin(); siter != picked_edges.end(); ++siter) {
        Edge *edge = mesh->getEdgeAt(*siter);
        Vertex *n0 = edge->getNodeAt(0);
        xyz = n0->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        Vertex *n1 = edge->getNodeAt(1);
        xyz = n1->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked faces ...
    //////////////////////////////////////////////////////////////////////////////

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for (siter = picked_faces.begin(); siter != picked_faces.end(); ++siter) {
        glBegin(GL_POLYGON);
        Face *face = mesh->getFaceAt(*siter);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *vtx = face->getNodeAt(j);
            xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
    }
    if( enable_lighting )
        glEnable(GL_LIGHTING);
}

void MeshViewer::draw_geodesics()
{
    Point3D xyz;
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 1.0, 1.0);

    glLineWidth(2.0);

    for (int i = 0; i < geodesics.size(); i++) {
        if (most_recent_geodesic_picked == i)
            glColor3f(0.0, 1.0, 0.0);
        else
            glColor3f(0.5, 1.0, 5.0);
        glBegin(GL_LINES);
        for (size_t j = 0; j < geodesics[i].size() - 1; j++) {
            xyz = geodesics[i][j]->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
            xyz = geodesics[i][j + 1]->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
    }
    glLineWidth(1.0);

    glColor3f(0.0, 0.0, 1.0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < geodesics.size(); i++) {
        xyz = geodesics[i].front()->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        xyz = geodesics[i].back()->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
    glPointSize(1.0);

    if( enable_lighting )
        glEnable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////

void MeshViewer::select_entity()
{
    static int currstat = 0;
    static int curr_pick_entity_id = -1;

    most_recent_vertex_picked = -1;
    most_recent_edge_picked = -1;
    most_recent_face_picked = -1;
    most_recent_geodesic_picked = -1;

    int id = selectedName();

    if (id >= 0) {
        if (currstat == 2) currstat = 0;
        if (pick_entity == 0) {
            if (delete_last_selection) {
                if (!picked_nodes.empty()) picked_nodes.pop_back();
                delete_last_selection = 0;
            } else {
                if (currstat == 0) {
                    int nodeid = id;
                    if (std::find(picked_nodes.begin(), picked_nodes.end(), nodeid) == picked_nodes.end())
                        picked_nodes.push_back(nodeid);
                }
            }
        }

        if (pick_entity == 1) {
            if (id >= mesh->getSize(0) && id < mesh->getSize(0) + mesh->getSize(1)) {
                if (delete_last_selection) {
                    if (!picked_edges.empty()) picked_edges.pop_back();
                    delete_last_selection = 0;
                } else {
                    if (currstat == 0) {
                        int edgeid = id - mesh->getSize(0);
                        if (std::find(picked_edges.begin(), picked_edges.end(), edgeid) == picked_edges.end())
                            picked_edges.push_back(edgeid);
                    }
                }
            }
        }

        if (pick_entity == 2) {
            if (id >= mesh->getSize(0) + mesh->getSize(1)) {
                if (delete_last_selection) {
                    if (!picked_faces.empty()) picked_faces.pop_back();
                    delete_last_selection = 0;
                } else {
                    if (currstat == 0) {
                        int faceid = id - mesh->getSize(0) - mesh->getSize(1);
                        if (std::find(picked_faces.begin(), picked_faces.end(), faceid) == picked_faces.end())
                            picked_faces.push_back(faceid);
                    }
                }
            }
        }

        if (pick_entity == 3) {
            if (id >= mesh->getSize(0) + mesh->getSize(1) + mesh->getSize(2)) {
                int curr_geodesic_id = id - mesh->getSize(0) - mesh->getSize(1) - mesh->getSize(2);
                if (delete_last_selection) {
                    geodesics.pop_back();
                    delete_last_selection = 0;
                }
            }
        }
        currstat++;
    }
    draw_picked_entities();
}
////////////////////////////////////////////////////////////////////////////////

void
MeshViewer::init()
{
    curr_layer_id = 0;
    wire_frame = 0;
    bounding_box = 0;
    faces_normal = 0;
    normal_sign = 1;
    enumerate = 0;
    use_color = 1;
    background = 0;
    pick_entity = 0;
    modify_node_position = 0;
    flip_edge = 0;
    swapedge = nullptr;
    display_nodes = 0;
    mesh_backup = nullptr;
    use_material = 1;
    color_filled_model = 1;
    delete_last_selection = 0;
    enable_lighting = 1;
    fontscale = 0.010;
    defect_remesh_step = 0;
    patch = nullptr;

//  camera()->setType(Camera::ORTHOGRAPHIC);

    int type = 3;
    if (type < 3) {
        // Move camera according to viewer type (on X, Y or Z axis)
        camera()->setPosition(Vec((type == 0) ? 1.0 : 0.0, (type == 1) ? 1.0 : 0.0, (type == 2) ? 1.0 : 0.0));
        camera()->lookAt(sceneCenter());

        camera()->showEntireScene();

        // Forbid rotation
        WorldConstraint* constraint = new WorldConstraint();
        constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
        camera()->frame()->setConstraint(constraint);
    }

    /*
      camera ()->setZNearCoefficient (0.001);
      camera ()->setSceneRadius (0.01);
      restoreStateFromFile ();
     */

    font = new FTPolygonFont(FONT_FILE);
    assert(!font->Error());
    font->FaceSize(1);
    font->CharMap(ft_encoding_unicode);

    face_numbering = 0;
    node_numbering = 0;

    if (!mesh_filename.empty())
        readData(mesh_filename);


    /*
        if (!bsurf_filename.empty())
        {
            ifstream infile(bsurf_filename.c_str(), ios::in);
            if (!infile.fail())
                infile >> bsurf;
        }
    */

    init_mesh();

    glPolygonOffset(1.0, 2.0);
    glEnable(GL_POLYGON_OFFSET_FILL);

    float ambient[] = {0.19125, 0.0735, 0.0225, 1.0};
    float diffuse[] = {0.7038, 0.27048, 0.0828, 1.0};
    float specular[] = {1.0, 1.0, 1.0, 1.0};
    float shininess = 120;

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

    glEnable(GL_NORMALIZE);
    if( enable_lighting );
    glEnable(GL_LIGHTING);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset ( 1.0, 1.0);

    /*
    glEnable( GL_DEPTH_TEST );

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        GLfloat light0_position[] = { 1.0, 0.2, 1.0, 0.0 };
        GLfloat light1_position[] = { 0.0, 0.0, 0.0, 1.0 };
        GLfloat light1_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat light1_specular[] = { 1.0, 1.0, 1.0, 1.0 };
        GLfloat lm_ambient[] = { 0.2, 0.2, 0.2, 1.0 };

        glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
        glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
        glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_diffuse);
        glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
        glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.2);
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient);

        glEnable( GL_LIGHT0 );
        glEnable(GL_LIGHT1);

        GLfloat localview = 0.0;
        glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, localview);
     */

    help();

}

void
MeshViewer::init_mesh()
{
    if (mesh) {

        normalize();
        mesh->setFacesNormal();
        mesh->search_boundary();
        create_colormap();

        lapsmooth = new Jaal::LaplaceSmoothing(mesh, 10);
        LaplaceWeight *lapwght = new Jaal::LaplaceLengthWeight;
        lapsmooth->setWeight(lapwght);

        mesh->buildRelations(0, 0);
        mesh->buildRelations(0, 2);
        cleanup.setMesh(mesh);

        /*
           build_dijkstra_mesh();
           if (!geodesic_file.empty())
                    read_geodesic_from_file();
        */
    }
}

void MeshViewer::build_dijkstra_mesh()
{
    if (mesh == nullptr) return;
    string tmpfile = "/tmp/model.off";
    mesh->saveAs(tmpfile);

    std::vector<double> points;
    std::vector<unsigned> faces;
    geodesic::read_mesh_off_file(tmpfile.c_str(), points, faces);

    geomesh = new geodesic::Mesh();
    geomesh->initialize_mesh_data(points, faces);

    dijkstra = new geodesic::GeodesicAlgorithmDijkstra(geomesh);
}

////////////////////////////////////////////////////////////////////////////////

vector<unsigned> MeshViewer::get_shortest_path(int src, int dst)
{

    std::vector<geodesic::SurfacePoint> sources;
    sources.push_back(geodesic::SurfacePoint(&geomesh->vertices()[src])); //one source is located at vertex zero

    std::vector<geodesic::SurfacePoint> targets; //same thing with targets
    targets.push_back(geodesic::SurfacePoint(&geomesh->vertices()[dst]));

    dijkstra->propagate(sources); //cover the whole mesh

    std::vector<geodesic::SurfacePoint> path;
    dijkstra->trace_back(targets[0], path);

    vector<unsigned> spath;
    for (int i = 0; i < path.size(); i++)
        spath.push_back(path[i].base_element()->id());

    return spath;
}

QString
MeshViewer::helpString() const
{
    QString text("<h2>S i m p l e V i e w e r</h2>");
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
}

#endif


int
main(int argc, char** argv)
{


    assert(argc >= 2);
    // Read command lines arguments.
    QApplication application(argc, argv);

    // Instantiate the viewer.
    MeshViewer viewer;
    viewer.readMeshFile(argv[1]);

    viewer.setWindowTitle("simpleMeshViewer");

    // Make the viewer window visible on screen.
    viewer.show();

    // Run main loop.
    return application.exec();
}

