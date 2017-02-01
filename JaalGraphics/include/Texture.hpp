#pragma once

#include <string>
#include "Image.hpp"

class JTexture
{
public:
    JTexture() {
        texID = 0;
    }
    GLuint  getID() const {
        return texID;
    }
    GLuint  fromImage( const string &s);
    GLuint  getCheckerBoard(int height, int width, int xcells, int ycells);
private:
    GLuint texID;
};
typedef boost::shared_ptr<JTexture>  JTexturePtr;

