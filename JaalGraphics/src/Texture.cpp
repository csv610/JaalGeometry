#include "Texture.hpp"

#include <iostream>
#include <sstream>

/////////////////////////////////////////////////////////////////////////////////////

GLuint JTexture :: fromImage( const string &s)
{
    JImage img;
    img.readFrom(s);

    // allocate a texture name
    glGenTextures(1, &texID);

    // select our current texture
    glBindTexture(GL_TEXTURE_2D, texID);

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
                 img.getWidth(),        // Image width  i.e. 640 for Kinect in standard mode
                 img.getHeight(),       // Image height i.e. 480 for Kinect in standard mode
                 0,                 // Border width in pixels (can either be 1 or 0)
                 GL_RGBA,            // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
                 GL_UNSIGNED_BYTE,  // Image data type
                 img.getData() );        // The actual image data itself

    // build our texture mipmaps
    gluBuild2DMipmaps(GL_TEXTURE_2D, 4, img.getWidth(), img.getHeight(), GL_RGBA,
                      GL_UNSIGNED_BYTE, img.getData() );

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    return texID;
}

/////////////////////////////////////////////////////////////////////////////////////

GLuint JTexture :: getCheckerBoard(int height, int width, int xcells, int ycells)
{
    ostringstream oss;
// oss << "convert -size " << width <<"x"<<height << " xc:wheat  floor.gif";
    oss << "convert -size " << width <<"x"<<height << " pattern:checkerboard -auto-level floor.gif";

    string cmd = oss.str();
    system(cmd.c_str());
    return fromImage( "floor.gif");
}


