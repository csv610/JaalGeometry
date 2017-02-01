#include "JaalViewer.hpp"
#include "MatrixViewer.hpp"

////////////////////////////////////////////////////////////////////////
JViewComponentPtr JMatrixViewer :: registerComponent(JaalViewer *viewer)
{
    boost::shared_ptr<JMatrixViewer> obj;
    obj.reset(new JMatrixViewer(viewer));
    obj->setName("MatrixViewer");
    viewer->attach(obj);
    return obj;
}
////////////////////////////////////////////////////////////////////////

JMatrixViewer :: JMatrixViewer( JaalViewer *vm)
{
   viewManager = vm;
   disk = gluNewQuadric();
   gluQuadricDrawStyle( disk, GLU_FILL );
   name  = "MatrixViewer";
}
////////////////////////////////////////////////////////////////////////

void JMatrixViewer :: draw()
{
    switch( graphType ) {
    case PRIMAL_GRAPH:
        drawPrimal();
        break;
    case DUAL_GRAPH:
        drawDual();
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMatrixViewer :: fillIJ( double px, double py)
{
    if( glyph == 0) {
        glBegin( GL_POINTS);
        glVertex2f( px, py);
        glEnd();
    } else {
        glPushMatrix();
        glTranslatef( px, py, 0.0 );
        gluDisk( disk, 0.0, pointSize, 10, 10);
        glPopMatrix();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMatrixViewer :: fillMatrix()
{
    int topDim  = mesh->getTopology()->getDimension();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f( 0.0, 0.0, 0.0);
    glBegin( GL_QUADS);
        glVertex2f( -0.5, -0.5);
        glVertex2f(  0.5, -0.5);
        glVertex2f(  0.5,  0.5);
        glVertex2f( -0.5,  0.5);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f( 1.0, 1.0, 1.0);
    glBegin( GL_QUADS);
        glVertex2f( -0.5, -0.5);
        glVertex2f(  0.5, -0.5);
        glVertex2f(  0.5,  0.5);
        glVertex2f( -0.5,  0.5);
    glEnd();

    size_t numNodes = mesh->getSize(0);
    dl = 1.0/(double)(numNodes);

    glColor3f( 1.0, 0.0, 0.0);

    glPointSize(pointSize);
    // Fill the diagonal....
    for( size_t i = 0; i < numNodes; i++) {
        double px  = -0.5 + (i+0.5)*dl;
        double py  =  0.5 - (i+0.5)*dl;
        fillIJ(px, py);
    }

    // Off Diagonal entries...
    glColor3f( 0.0, 0.0, 0.0);
    if( topDim == 2 ) {
        size_t numFaces = mesh->getSize(2);
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                int nedges = face->getSize(1);
                for( int j = 0; j < nedges; j++) {
                    const JEdgePtr &edge = face->getEdgeAt(j);
                    drawEdge(edge);
                }
            }
        }
    }

    if( topDim == 3 ) {
        size_t numCells = mesh->getSize(3);
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                int nedges = cell->getSize(1);
                for( int j = 0; j < nedges; j++) {
                    const JEdgePtr &edge = cell->getEdgeAt(j);
                    drawEdge(edge);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMatrixViewer ::drawIDs()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    /*
         static GLuint dlistID = 0;
         ViewComponent::displayList(distID);
    */

    FTFont *font = FontsManager::Instance().getFont();

    float fS = 0.001*fontsScale;

    if( font == nullptr) return;
    char number[128];

    size_t numNodes = mesh->getSize(0);

    double x, y, z = 0.0;
    glColor3f(1.0, 1.0, 0.0);

    for (size_t i = 0; i < numNodes; i+=stepLabels) {
        x = -0.5 + (i+0.5)*dl;
        y = 0.51;
        glPushMatrix();
        glTranslatef(x, y, z);
        glRotatef( 90, 0.0, 0.0, 1.0);
        glScalef(fS, fS, fS);
        sprintf(number, "%ld", i);
        font->Render(number);
        glPopMatrix();

        x = 0.51;
        y = 0.5 - (i+0.5)*dl;
        glPushMatrix();
        glTranslatef(x, y, z);
        glScalef(fS, fS, fS);
        sprintf(number, "%ld", i);
        font->Render(number);
        glPopMatrix();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMatrixViewer :: drawEdge( const JEdgePtr &edge)
{
    int n0  = edge->getNodeAt(0)->getID();
    int n1  = edge->getNodeAt(1)->getID();

    double px  = -0.5 + (n0+0.5)*dl;
    double py  =  0.5 - (n1+0.5)*dl;
    fillIJ(px,py);

    px  = -0.5 + (n1+0.5)*dl;
    py  =  0.5 - (n0+0.5)*dl;
    fillIJ(px,py);
}

///////////////////////////////////////////////////////////////////////////////

void JMatrixViewer :: drawPrimal()
{
    if( mesh == nullptr)return;

    size_t numNodes = mesh->getSize(0);

    if( numNodes < 1 ) {
        cout << "Info: There are no nodes in the mesh " << endl;
        return;
    }

    glColor3f( 1.0, 1.0, 0.0);
    glDisable( GL_LIGHTING );
    fillMatrix();
    drawIDs();
}


///////////////////////////////////////////////////////////////////////////////
void JMatrixViewer :: drawDual()
{
}
