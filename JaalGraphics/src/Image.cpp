#include "Image.hpp"


void JImage :: readFrom( const string &s)
{
    fileName = s;
    if (fileName.empty()) return;
    image = new Magick::Image();
    image->read(fileName.c_str() );
    image->flip();
    size_t pos0 = fileName.rfind("/");
    size_t pos1 = fileName.rfind('.');
    name = fileName.substr(pos0+1, pos1-pos0-1);
    genTexture();
}
///////////////////////////////////////////////////////////////////////////////
void
JImage:: plainImage()
{
    if( image == nullptr ) return;

    double zc = -0.0001;
    double len[2];
    double h =  1.0*getHeight()/(double)getWidth();
    if( getHeight() > getWidth() ) {
        len[0] = 1.0;
        len[1] = h;
    } else  {
        len[0] = 1/h;
        len[1] = 1.0;
    }
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(-0.5*len[0], -0.5*len[1], zc);

    glTexCoord2f(0.0, 1.0);
    glVertex3f(-0.5*len[0],  0.5*len[1], zc);

    glTexCoord2f(1.0, 1.0);
    glVertex3f(0.5*len[0],   0.5*len[1], zc);

    glTexCoord2f(1.0, 0.0);
    glVertex3f(0.5*len[0],  -0.5*len[1], zc);
    glEnd();
}

//////////////////////////////////////////////////////////////////////////////////

void
JImage:: texturedMesh()
{
    /*
            if( image == nullptr || mesh == nullptr ) return;

            double zc = -0.0001;

            size_t numfaces = mesh->getSize(2);

            Point3D xyz;
            Point2D uv;

            glBegin(GL_TRIANGLES);
            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr face = mesh->getFaceAt(i);
                if( face->isActive() ) {
                    int nnodes = face->getSize(0);
                    if( nnodes == 3 ) {
                        for( int j = 0; j < 3; j++) {
                            const JNodePtr &vtx = face->getNodeAt(j);
                            xyz = vtx->getXYZCoords();
                            int err = vtx->getAttribute("UVCoords", uv);
                            assert( !err );
                            glTexCoord2f( uv[0], uv[1] );
                            glVertex3f(xyz[0], xyz[1], zc );
                        }
                    }
                }
            }
            glEnd();

            glBegin(GL_QUADS);
            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr face = mesh->getFaceAt(i);
                if( face->isActive() ) {
                    int nnodes = face->getSize(0);
                    if( nnodes == 4 ) {
                        for( int j = 0; j < 4; j++) {
                            const JNodePtr &vtx = face->getNodeAt(j);
                            xyz = vtx->getXYZCoords();
                            int err = vtx->getAttribute("UVCoords", uv);
                            assert( !err );
                            glTexCoord2f( uv[0], uv[1] );
                            glVertex3f(xyz[0], xyz[1], zc );
                        }
                    }
                }
            }
            glEnd();
    */
}
//////////////////////////////////////////////////////////////////////////////////

void
JImage::draw()
{
    if( activeBit == 0) return;

    glDisable( GL_LIGHTING );
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable( GL_BLEND );
    glEnable( GL_TEXTURE_2D );
    glBindTexture(GL_TEXTURE_2D, textureID);

    /*
        if( mesh )
            texturedMesh();
        else
    */
    plainImage();

    glDisable( GL_TEXTURE_2D );
    glEnable( GL_LIGHTING );
}

/////////////////////////////////////////////////////////////////////////////////////////////////

GLuint JImage :: genTexture()
{
    // allocate a texture name
    glGenTextures(1, &textureID);

    // select our current texture
    glBindTexture(GL_TEXTURE_2D, textureID);

    // select modulate to mix texture with color for shading
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_DECAL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_DECAL);

    // when texture area is small, bilinear filter the closest mipmap
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
    // when texture area is large, bilinear filter the first mipmap
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // texture should tile
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glTexImage2D(GL_TEXTURE_2D,     // Type of texture
                 0,                 // Pyramid level (for mip-mapping) - 0 is the top level
                 GL_RGB,            // Internal colour format to convert to
                 getWidth(),        // Image width  i.e. 640 for Kinect in standard mode
                 getHeight(),       // Image height i.e. 480 for Kinect in standard mode
                 0,                 // Border width in pixels (can either be 1 or 0)
                 GL_RGBA,            // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
                 GL_UNSIGNED_BYTE,  // Image data type
                 getData() );        // The actual image data itself

    // build our texture mipmaps
    gluBuild2DMipmaps(GL_TEXTURE_2D, 4, getWidth(), getHeight(), GL_RGBA,
                      GL_UNSIGNED_BYTE, getData() );

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    return textureID;
}

