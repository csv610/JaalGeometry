#include "PickingTexture.hpp"

void JPickingTexture:: getColor( size_t  id, Point3I &rgb)
{
    rgb[0] = (id & 0x000000FF) >>  0;
    rgb[1] = (id & 0x0000FF00) >>  8;
    rgb[2] = (id & 0x00FF0000) >> 16;
}

/////////////////////////////////////////////////////////////////////////

size_t  JPickingTexture :: getID( const Point3I &data) const
{
    return  data[0] + data[1]*256 + data[2]*256*256;
}

/////////////////////////////////////////////////////////////////////////

void JPickingTexture:: genTexture( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);
    mesh->enumerate(2);

    vector<Point3I> faceColor;

    faceColor.reserve(numfaces);
    Point3I rgb;
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            size_t id = face->getID();
            getColor(id, rgb);
            faceColor.push_back(rgb);

        }
    }

    glDisable( GL_LIGHTING);
    glEnable( GL_DEPTH_TEST);
    glDisable( GL_BLEND);

    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int id = face->getID();
            int numnodes = face->getSize(0);
            Point3I rgb = faceColor[id];
            glColor3f( rgb[0], rgb[1], rgb[2] );
            for( int j = 0; j <  numnodes; j++) {
                const Point3D &p3d = face->getNodeAt(j)->getXYZCoords();
                glBegin(GL_POLYGON);
                glVertex3d(p3d[0], p3d[1], p3d[2] );
                glEnd();
            }
        }
    }
    glFlush();
    glFinish();
    getBuffer();
}

/////////////////////////////////////////////////////////////////////////
void JPickingTexture:: getBuffer()
{
    int width  = viewer->width();
    int height = viewer->height();
    int numPixels = width*height;
    pixelData.resize(3*numPixels);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &pixelData[0]);
}
/////////////////////////////////////////////////////////////////////////

JFaceSequence JPickingTexture:: getPickedFaces( const Point2I &pos)
{
    int width  = viewer->width();
    int height = viewer->height();

    int xorg = pos[0] - 0.5*brushSize;
    int yorg = pos[1] - 0.5*brushSize;

    if( xorg < 0) xorg = 0;
    if( xorg > width) xorg = width -1;

    if( yorg < 0) yorg = 0;
    if( yorg > height) yorg = height -1;

    JFaceSet faceSet;
    Point3I  color;

    for( int i = 0; i < width; i++) {
        for( int j = 0; j < height; j++) {
            int ii = xorg + i;
            int jj = yorg + j;
            size_t offset = jj*width + ii;
            color[0] = pixelData[3*offset+0];
            color[1] = pixelData[3*offset+1];
            color[2] = pixelData[3*offset+2];
            int id = getID( color);
            faceSet.insert( mesh->getFaceAt(id) );
        }
    }

    JFaceSequence faces;
    int numfaces = faceSet.size();
    if( numfaces == 0) return faces;
    faces.resize( numfaces);
    copy( faceSet.begin(), faceSet.end(), faces.begin() );
    return faces;

}

/////////////////////////////////////////////////////////////////////////
void JPickingTexture::init()
{
    /*
        // Create the FBO
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);

        // Create the texture object for the primitive information buffer
        glGenTextures(1, &pickingTexture);
        glBindTexture(GL_TEXTURE_2D, pickingTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight,
                    0, GL_RGB, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                    pickingTexture, 0);

        // Create the texture object for the depth buffer
        glGenTextures(1, &depthTexture);
        glBindTexture(GL_TEXTURE_2D, depthTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, windowWidth, windowHeight,
                    0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                     depthTexture, 0);

        // Disable reading to avoid problems with older GPUs
        glReadBuffer(GL_NONE);

        glDrawBuffer(GL_COLOR_ATTACHMENT0);

        // Verify that the FBO is correct
        GLenum Status = glCheckFramebufferStatus(GL_FRAMEBUFFER);

        if (Status != GL_FRAMEBUFFER_COMPLETE) {
            printf("FB error, status: 0x%x\n", Status);
            return false;
        }

        // Restore the default framebuffer
        glBindTexture(GL_TEXTURE_2D, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return GLCheckError();
    */
}

//////////////////////////////////////////////////////////////////////////////////////
size_t JPickingTexture :: readPixel(const Point2I &p)
{
    /*
        glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
        glReadBuffer(GL_COLOR_ATTACHMENT0);

        Point3D pixel;
        glReadPixels(x, y, 1, 1, GL_RGB, GL_FLOAT, &pixel[0]);

        glReadBuffer(GL_NONE);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
        return pixel;
    */
}

//////////////////////////////////////////////////////////////////////////////////////
