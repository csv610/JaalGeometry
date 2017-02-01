#pragma once

#include <string>
#include <vector>

#ifdef HAVE_LOG4
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/helpers/exception.h>
#endif

#include <GL/gl.h>
#include <GL/gle.h>
#include <GL/glu.h>

#include <boost/smart_ptr.hpp>

#include <Magick++.h>
//#include <opencv/cv.h>

using namespace std;

//using namespace Magick;
//using namespace cv;

class JImage
{
    typedef Magick::Image Image_t;
public:
    JImage() {
        image  = nullptr;
        activeBit = 1;
    }

    void setCheckerBoardSize(int n, int m)
    {
        checkerWidth = n, checkerHeight = m;
    }

    void setActiveBit( bool v) {
        activeBit = v;
    }
    bool isActive() const {
        return activeBit;
    }

    void readFrom( const string &s);

    string getName() const {
        return name;
    }

    string getFileName() const {
        return fileName;
    }

    Image_t*  getImage() {
        return image;
    }

    int  getWidth() const {
        return image->columns();
    }

    int  getHeight() const {
        return image->rows();
    }

    const void*  getData() {
        image->magick("RGBA");
        image->write(&blob);
        return blob.data();
    }

    void    draw();
    GLuint  genTexture();

    void saveAs(const string &s);
private:
    Magick::Blob blob;
    string    fileName, name;
    GLuint    textureID;
    Image_t  *image;
    bool     activeBit;
    void     plainImage();
    void     texturedMesh();
    int      checkerWidth, checkerHeight;
};

typedef boost::shared_ptr<JImage> JImagePtr;

