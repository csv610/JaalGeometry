#include "MeshCutter.h"
#include <assert.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sstream>

using namespace Jaal;
using namespace std;

#ifdef CSV


GLfloat whiteColor[] = {1.0, 1.0, 1.0};
GLfloat bgColor[] = {1.0, 1.0, 1.0};

/////////////////////////////////////////////////////////////////////////////////////

void
MeshCutter::save_features()
{
    cout << "Info: Saving Features in file " << endl;
    ofstream ofile("mesh_feature.dat", ios::out);
    if (ofile.fail()) {
        return;
    }
    vector<int>::const_iterator siter;
    ofile << picked_nodes.size() << " " << picked_edges.size() << " " << picked_faces.size() << endl;
    for (siter = picked_nodes.begin(); siter != picked_nodes.end(); ++siter) {
        ofile << *siter << endl;
    }

    for (siter = picked_edges.begin(); siter != picked_edges.end(); ++siter) {
        Edge *edge = mesh->getEdgeAt(*siter);
        Vertex *n0 = edge->getNodeAt(0);
        Vertex *n1 = edge->getNodeAt(1);
        ofile << n0->getID() << " " << n1->getID() << endl;
    }

    for (siter = picked_faces.begin(); siter != picked_faces.end(); ++siter) {
        Face *face = mesh->getFaceAt(*siter);
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v = face->getNodeAt(j);
            ofile << v->getID() << " ";
        }
        ofile << endl;
    }
    ofile.close();
}

/////////////////////////////////////////////////////////////////////////////////////

void
MeshCutter::keyPressEvent(QKeyEvent *event)
{
    static bool camera_type = 0;

    if (event->key() == Qt::Key_B) {
        bounding_box = !bounding_box;
        updateGL();
    }

    if (event->key() == Qt::Key_E) {
        enumerate = !enumerate;
        updateGL();
    }

    if (event->key() == Qt::Key_F) {
        display_faces = !display_faces;
        updateGL();
    }

    if (event->key() == Qt::Key_F2) {
        fontscale += 0.001;
        updateGL();
    }

    if (event->key() == Qt::Key_F3) {
        fontscale -= 0.001;
        updateGL();
    }

    if (event->key() == Qt::Key_F4) {
        if( camera_type == 0) {
            camera()->setType(Camera::PERSPECTIVE);
            camera_type = 1;
            updateGL();
            return;
        }

        if( camera_type == 1) {
            camera()->setType(Camera::ORTHOGRAPHIC);
            camera_type = 0;
            updateGL();
            return;
        }
    }

    if (event->key() == Qt::Key_G) {
        int nSize = picked_nodes.size();
        if(  nSize >= 2 )  {
            int src = picked_nodes[nSize-2];
            int dst = picked_nodes[nSize-1];
            vector<unsigned> spath;
            get_shortest_path(src, dst, spath );
            geodesics.push_back(spath);
            picked_nodes.clear();
            updateGL();
            return;
        }
    }


    if (event->key() == Qt::Key_N) {
        display_normals = !display_normals;
        updateGL();
    }

    if (event->key() == Qt::Key_P) {
        pick_entity += 1;
        pick_entity = pick_entity % 3;
        switch (pick_entity) {
        case 0:
            cout << "Vertex Picking On " << endl;
            break;
        case 1:
            cout << "Edge Picking On " << endl;
            break;
        case 2:
            cout << "Face Picking On " << endl;
            break;
        }
        updateGL();
    }

    if (event->key() == Qt::Key_Period) {
        display_nodes = !display_nodes;
        updateGL();
    }

    if (event->key() == Qt::Key_R) {
        normal_sign *= -1;
        updateGL();
    }

    if (event->key() == Qt::Key_S) {
        save_features();
        return;
    }

    if (event->key() == Qt::Key_W) {
        display_edges = !display_edges;
        updateGL();
    }

    if (event->key() == Qt::Key_X) {
        updateGL();
    }

    QGLViewer::keyPressEvent(event);
}

///////////////////////////////////////////////////////////////////////////////

void MeshCutter::build_dijkstra_mesh()
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

void MeshCutter::get_shortest_path(int src, int dst, vector<unsigned> &spath )
{
    spath.clear();

    std::vector<geodesic::SurfacePoint> sources;
    sources.push_back(geodesic::SurfacePoint(&geomesh->vertices()[src])); //one source is located at vertex zero

    std::vector<geodesic::SurfacePoint> targets; //same thing with targets
    targets.push_back(geodesic::SurfacePoint(&geomesh->vertices()[dst]));

    dijkstra->propagate(sources); //cover the whole mesh

    std::vector<geodesic::SurfacePoint> path;
    dijkstra->trace_back(targets[0], path);

    for (int i = 0; i < path.size(); i++)
        spath.push_back(path[i].base_element()->id());
}

////////////////////////////////////////////////////////////////////////////////

void
MeshCutter::normalize()
{
    mesh->normalize();

    BoundingBox mbox = mesh->getBoundingBox();
    Point3D lower = mbox.getLowerLeftCorner();
    Point3D upper = mbox.getUpperRightCorner();

    double xlen = mbox.getLength(0);
    double ylen = mbox.getLength(1);
    double zlen = mbox.getLength(2);

    lower[0] -= 0.010 * xlen;
    lower[1] -= 0.010 * ylen;
    lower[2] -= 0.010 * zlen;

    upper[0] += 0.010 * xlen;
    upper[1] += 0.010 * ylen;
    upper[2] += 0.010 * zlen;

    box.setLowerLeftCorner(lower);
    box.setUpperRightCorner(upper);
}

///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_bounding_box()
{
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(2.0);
    glColor3f(1.0, 1.0, 1.0);

    static GLuint dlistID = 0;
    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    }
    dlistID = glGenLists(1) + 1;
    glNewList( dlistID, GL_COMPILE_AND_EXECUTE);


    Point3D lower = box.getLowerLeftCorner();
    Point3D upper = box.getUpperRightCorner();

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

    glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_nodes()
{
    glDisable( GL_LIGHTING );
    glEnable( GL_SMOOTH );

    glPointSize(2.0);
    Point3D xyz;
    glColor3f(0.0, 1.0, 0.0);

    static GLuint dlistID = 0;
    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    }
    dlistID = glGenLists(1) + 1;
    glNewList( dlistID, GL_COMPILE_AND_EXECUTE);

    size_t numnodes = mesh->getSize(0);

    glBegin(GL_POINTS);
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v  = mesh->getNodeAt(i);
        const Point3D &xyz = v->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
    glDisable( GL_SMOOTH );
    glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_edges()
{
    glDisable( GL_LIGHTING );

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glLineWidth( 1.0);
    Point3D xyz;
    glColor3f(0.1, 0.1, 0.1);

    static GLuint dlistID = 0;
    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    }
    dlistID = glGenLists(1) + 1;
    glNewList( dlistID, GL_COMPILE_AND_EXECUTE);

    size_t numedges = mesh->getSize(1);

    glBegin(GL_LINES);
    for (size_t i = 0; i < numedges; i++) {
        Edge *edge = mesh->getEdgeAt(i);
        Vertex *n0 = edge->getNodeAt(0);
        Vertex *n1 = edge->getNodeAt(1);
        xyz = n0->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        xyz = n1->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
    glLineWidth(1.0);

    glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_geodesics()
{
    glDisable( GL_LIGHTING );

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glLineWidth( 5.0);
    Point3D xyz;
    glColor3f(0.0, 1.0, 1.0);

    glBegin(GL_LINES);
    int nSize = geodesics.size();
    for (size_t i = 0; i < nSize; i++) {
        int mSize = geodesics[i].size();
        for (size_t j = 0; j < mSize-1; j++) {
            int v0 =  geodesics[i][j];
            int v1 =  geodesics[i][j+1];
            Vertex *n0 = mesh->getNodeAt(v0);
            Vertex *n1 = mesh->getNodeAt(v1);
            xyz = n0->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
            xyz = n1->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
    }
    glEnd();
    glLineWidth(1.0);
}


///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_normals()
{
    double x, y, z;

    glColor3f(0.0, 0.0, 1.0);

    glLineWidth(2.0);

    static GLuint dlistID = 0;
    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    }
    dlistID = glGenLists(1) + 1;
    glNewList( dlistID, GL_COMPILE_AND_EXECUTE);

    Vec3D normal;
    Point3D centroid;

    size_t numfaces = mesh->getSize(2);
    glBegin(GL_LINES);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        face->getAvgPos( centroid );

        x = centroid[0];
        y = centroid[1];
        z = centroid[2];
        glVertex3f(x, y, z);

        face->getAttribute("Normal", normal);
        x += 0.01 * normal_sign * normal[0];
        y += 0.01 * normal_sign * normal[1];
        z += 0.01 * normal_sign * normal[2];
        glVertex3f(x, y, z);
    }
    glEnd();

    /*
         size_t numnodes = mesh->getSize(0);
         glBegin(GL_LINES);
         for (size_t i = 0; i < numnodes; i++) {
              Vertex *vertex = mesh->getNodeAt(i);
              const Point3D &xyz = vertex->getXYZCoords();
              x = xyz[0];
              y = xyz[1];
              z = xyz[2];
              glVertex3f(x, y, z);
              vertex->getAttribute("Normal", normal);
              x += 0.01 * normal_sign * normal[0];
              y += 0.01 * normal_sign * normal[1];
              z += 0.01 * normal_sign * normal[2];
              glVertex3f(x, y, z);
         }
         glEnd();
    */

    glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_faces()
{
    size_t numfaces = mesh->getSize(2);
    glColor3f( 1.0, 0.0, 0.0);

    glShadeModel(GL_FLAT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    Vec3D normal;

    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        if (nnodes == 3) {
            face->getAttribute( "Normal", normal );
            float nx = 1.0 * normal_sign * normal[0];
            float ny = 1.0 * normal_sign * normal[1];
            float nz = 1.0 * normal_sign * normal[2];
            glNormal3f(nx, ny, nz);

            for (int j = 0; j < 3; j++) {
                Vertex *vtx = face->getNodeAt(j);
                const Point3D &xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
        }
    }
    glEnd();

}


///////////////////////////////////////////////////////////////////////////////

void
MeshCutter::drawWithNames()
{
    glDisable(GL_LIGHTING);
    glPointSize(2.0);

    size_t numnodes = mesh->getSize(0);
    Point3D xyz;

    int offset  = 0;
    if (pick_entity == 0) {
        glColor3f(0.0, 1.0, 0.0);
        glPointSize(5.0);
        for (size_t i = 0; i < numnodes; i++) {
            Vertex *vtx = mesh->getNodeAt(i);
            glPushMatrix();
            glPushName(offset + i);
            glBegin(GL_POINTS);
            xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
            glEnd();
            glPopName();
            glPopMatrix();
        }
    }

    size_t numedges = mesh->getSize(1);

    offset = numnodes;
    if (pick_entity == 1) {
        glLineWidth(5.0);
        offset = mesh->getSize(0);
        glColor3f(0.0, 1.0, 0.0);
        for (int i = 0; i < numedges; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            glPushMatrix();
            glPushName(offset + i);

            glBegin(GL_LINES);
            Vertex *n0 = edge->getNodeAt(0);
            xyz = n0->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);

            Vertex *n1 = edge->getNodeAt(1);
            xyz = n1->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
            glEnd();

            glPopName();
            glPopMatrix();
        }
    }

    size_t numfaces = mesh->getSize(2);
    offset = numnodes + numedges;

    if (pick_entity == 2) {
        offset = numnodes + numedges;
        glColor3f(0.0, 1.0, 0.0);
        for (int i = 0; i < numfaces; i++) {
            Face *face = mesh->getFaceAt(i);
            glPushMatrix();
            glPushName(offset + i);

            glBegin(GL_POLYGON);
            int nnodes = face->getSize(0);
            for (int j = 0; j < nnodes; j++) {
                Vertex *vtx = face->getNodeAt(j);
                xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
            glEnd();
            glPopName();
            glPopMatrix();
        }
    }

    glEnable(GL_LIGHTING);

    glPointSize(1.0);
}

//////////////////////////////////////////////////////////////////////////////////

void
MeshCutter::draw_picked_entities()
{
    Point3D xyz;
    vector<int> ::const_iterator siter;
    glColor3f(1.0, 1.0, 0.0);

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked nodes ...
    //////////////////////////////////////////////////////////////////////////////
    if( pick_entity == 0) {
        glEnable(GL_LIGHTING);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        gluQuadricDrawStyle(sphere, GLU_FILL);
        for (siter = picked_nodes.begin(); siter != picked_nodes.end(); ++siter) {
            Vertex *vtx = mesh->getNodeAt(*siter);
            glPushMatrix();
            xyz = vtx->getXYZCoords();
            glTranslatef( xyz[0], xyz[1], xyz[2] );
            gluSphere( sphere, 0.01, 20, 20);
            glPopMatrix();
        }
    }

    glDisable(GL_LIGHTING);

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked edges ...
    //////////////////////////////////////////////////////////////////////////////

    if( pick_entity == 1 ) {
        /*
                  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                  gleDouble endPoints[3][3];
                  float  endColor[3][3];
                  endColor[0][0] = 0.0;
                  endColor[0][1] = 1.0;
                  endColor[0][2] = 0.0;

                  endColor[1][0] = 0.0;
                  endColor[1][1] = 1.0;
                  endColor[1][2] = 0.0;

                  endColor[2][0] = 0.0;
                  endColor[2][1] = 1.0;
                  endColor[2][2] = 0.0;

                  for (siter = picked_edges.begin(); siter != picked_edges.end(); ++siter) {
                       Edge *edge = mesh->getEdgeAt(*siter);
                       Vertex *n0 = edge->getNodeAt(0);
                       xyz = n0->getXYZCoords();
                       endPoints[0][0] = xyz[0];
                       endPoints[0][1] = xyz[1];
                       endPoints[0][2] = xyz[2];

                       Vertex *n1 = edge->getNodeAt(1);
                       xyz = n1->getXYZCoords();
                       endPoints[2][0] = xyz[0];
                       endPoints[2][1] = xyz[1];
                       endPoints[2][2] = xyz[2];

                       endPoints[1][0] = 0.5*( endPoints[0][0] + endPoints[2][0] );
                       endPoints[1][1] = 0.5*( endPoints[0][1] + endPoints[2][1] );
                       endPoints[1][2] = 0.5*( endPoints[0][2] + endPoints[2][2] );

                       glePolyCylinder (3, endPoints, endColor, 1.01);
                  }
        */

        glLineWidth(5);
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
    }

    //////////////////////////////////////////////////////////////////////////////
    // Draw picked faces ...
    //////////////////////////////////////////////////////////////////////////////

    if( pick_entity == 2 ) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for (siter = picked_faces.begin(); siter != picked_faces.end(); ++siter) {
            Face *face = mesh->getFaceAt(*siter);
            glBegin(GL_POLYGON);
            int nnodes = face->getSize(0);
            for (int j = 0; j < nnodes; j++) {
                Vertex *vtx = face->getNodeAt(j);
                xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
            glEnd();
        }
    }

    glEnable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////

void MeshCutter::select_entity()
{
    int id = selectedName();
    if( id < 0 ) return;

    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);
    size_t numfaces = mesh->getSize(2);

    if (pick_entity == 0) {
        picked_edges.clear();
        picked_faces.clear();
        if( id < numnodes ) {
            int nodeid = id;
            if( !multi_selection) picked_nodes.clear();
            if( picked_nodes.empty() ) {
                picked_nodes.push_back(nodeid);
            } else {
                if( picked_nodes.back() != nodeid )
                    picked_nodes.push_back(nodeid);
            }
        }
    }

    if (pick_entity == 1) {
        picked_nodes.clear();
        picked_faces.clear();
        if( id < numnodes + numedges  ) {
            size_t edgeid = id - numnodes;
            if( edgeid >= 0 && edgeid < numedges ) {
                if( !multi_selection) picked_edges.clear();
                if( picked_edges.empty() ) {
                    picked_edges.push_back( edgeid );
                } else {
                    if( picked_edges.back() != edgeid )
                        picked_edges.push_back( edgeid );
                }
            }
        }
    }

    if (pick_entity == 2) {
        picked_nodes.clear();
        picked_edges.clear();
        if( id < numnodes + numedges + numfaces ) {
            size_t faceid = id - (numnodes + numedges ) ;
            if( faceid >= 0 && faceid < numfaces ) {
                if( !multi_selection) picked_faces.clear();
                if( picked_faces.empty() ) {
                    picked_faces.push_back( faceid);
                } else {
                    if( picked_faces.back() != faceid )
                        picked_faces.push_back( faceid);
                }
            }
        }
    }

    draw_picked_entities();
}

////////////////////////////////////////////////////////////////////////////////
void
MeshCutter::draw()
{
    qglviewer::Vec wc, sc;

    glClearColor(bgColor[0], bgColor[1], bgColor[2], 1.0);

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glEnable(GL_LIGHTING);

    select_entity();

    if (display_nodes) draw_nodes();
    if (display_edges) draw_edges();
    if( display_faces) draw_faces();

    if (display_normals) draw_normals();
    if (bounding_box) draw_bounding_box();

    draw_geodesics();

    char number[128];
    glColor3f(1.0, 1.0, 0.0);
    glLineWidth(1.0);

    if (enumerate) {
        Vec3D normal;
        glColor3f(0.0, 1.0, 0.0);
        for (int i = 0; i < mesh->getSize(0); i++) {
            Vertex *vtx = mesh->getNodeAt(i);
            const Point3D &xyz = vtx->getXYZCoords();
            vtx->getAttribute("Normal", normal);
            wc.x = xyz[0] + 0.005*normal[0];
            wc.y = xyz[1] + 0.005*normal[1];
            wc.z = xyz[2] + 0.005*normal[2];
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glPushMatrix();
            glTranslatef(wc.x, wc.y, wc.z);
            glScalef(fontscale, fontscale, fontscale);
            sprintf(number, "%d", i);
            font->Render(number);
            glPopMatrix();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void
MeshCutter::init()
{
    display_nodes = 0;
    display_edges = 1;
    display_faces = 1;
    display_normals = 0;

    bounding_box  = 0;
    normal_sign = 1;
    enumerate = 0;
    pick_entity = 0;
    fontscale = 0.0010;
    multi_selection = 1;

    sphere = gluNewQuadric();

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
    */

    font = new FTPolygonFont(FONT_FILE);
    assert(!font->Error());
    font->FaceSize(1);
    font->CharMap(ft_encoding_unicode);

    face_numbering = 0;
    node_numbering = 0;

    if (!mesh_filename.empty()) {
        mesh = new Jaal::Mesh();
        cout << "Reading Mesh file ... " << endl;
        mesh->readFromFile( mesh_filename );
    }


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
    glEnable(GL_LIGHTING);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset ( 1.0, 1.0);

    glEnable( GL_DEPTH_TEST);
    glCullFace( GL_BACK );

//   restoreStateFromFile ();
}

void
MeshCutter::init_mesh()
{
    if (mesh) {
        normalize();
        cout << "Setting Normals " << endl;
        mesh->setNormals();

        cout << "Building edges ... " << endl;
        mesh->build_edges();

        cout << "Building Dijkstra ... " << endl;
        build_dijkstra_mesh();
    }
}

int
main(int argc, char** argv)
{
    assert(argc >= 2);
    // Read command lines arguments.
    QApplication application(argc, argv);

    // Instantiate the viewer.
    MeshCutter viewer;
    viewer.setMeshFile(argv[1]);

    viewer.setWindowTitle("MeshCutter");

    // Make the viewer window visible on screen.
    viewer.show();

    // Run main loop.
    return application.exec();
}
#endif

