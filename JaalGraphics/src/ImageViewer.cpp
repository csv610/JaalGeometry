#include "ImageViewer.hpp"

using namespace Jaal;
using namespace std;

///////////////////////////////////////////////////////////////////////////////

#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/gle.h>

///////////////////////////////////////////////////////////////////////////////

JViewComponentPtr JImageViewer :: registerComponent(JaalViewer *viewer)
{
    boost::shared_ptr<JImageViewer> obj;
    obj.reset(new JImageViewer(viewer));
    obj->setName("ImageViewer");
    viewer->attach(obj);
    return obj;
}

///////////////////////////////////////////////////////////////////////////////

JImageViewer :: JImageViewer( JaalViewer *vm)
{
    viewManager = vm;
    Magick::InitializeMagick(nullptr);
}
///////////////////////////////////////////////////////////////////////////////

JImageViewer :: ~JImageViewer()
{
}

///////////////////////////////////////////////////////////////////////////////

void
JImageViewer::draw()
{
    if( !active ) return;
    for(const auto &img : imagedb) img->draw();
}

///////////////////////////////////////////////////////////////////////////////

void JImageViewer :: refreshDisplay()
{
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageViewer :: addImage( const string &str)
{
    JImagePtr img = getImage(str);
    if( img == nullptr ) {
        img = boost::shared_ptr<JImage>(new JImage);
        img->readFrom(str);
        imagedb.push_back(img);
    }
}
///////////////////////////////////////////////////////////////////////////////
int JImageViewer:: remove(const JImagePtr &img)
{
    boost::remove_erase(imagedb, img);
    /*
        vector<JImagePtr>::iterator it;
        it = std::remove(imagedb.begin(), imagedb.end(), img);
        imagedb.erase(it, imagedb.end());
    */
    return 0;
}
