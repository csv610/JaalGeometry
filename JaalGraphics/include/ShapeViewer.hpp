#pragma once

#include "JaalViewer.hpp"
#include "MeshViewer.hpp"
#include <boost/smart_ptr.hpp>

using namespace std;
using namespace qglviewer;
using namespace Jaal;

#include "BasicShapes.hpp"

struct JSphereRender
{
    JColor color;
    int   numSlices;
    int   numStacks;
};
typedef boost::shared_ptr<JSphereRender> JSphereRenderPtr;

struct JCylinderRender
{
    bool   cap = 1;
    int    numSlices, numStacks;
    JColor color;
};
typedef boost::shared_ptr<JCylinderRender> JCylinderRenderPtr;

struct JConeRender
{
    bool   cap;
    JColor  color;
    int    numSlices, numStacks;
};
typedef boost::shared_ptr<JConeRender> JConeRenderPtr;

struct JBoxRender
{
    JColor color;
};
typedef boost::shared_ptr<JBoxRender> JBoxRenderPtr;

struct JTubeRender
{
    JColor color;
};
typedef boost::shared_ptr<JTubeRender> JTubeRenderPtr;

class JShapeViewer : public JViewComponent {
public:
    static JViewComponentPtr registerComponent(JaalViewer *root);

    explicit JShapeViewer( JaalViewer *p);
    ~JShapeViewer();

    void clearAll() {
        spheres.clear();
        boxes.clear();
        cones.clear();
        cylinders.clear();
        tubes.clear();
    }

    void addObject( const JSpherePtr &obj) {
        spheres.push_back(obj);
    }
    void addObject( const JCylinderPtr &obj) {
        cylinders.push_back(obj);
    }
    void addObject( const JBoxPtr &obj) {
        boxes.push_back( obj);
    }
    void addObject( const JConePtr &obj) {
        cones.push_back(obj);
    }
    void addObject( const JTubePtr &obj) {
        tubes.push_back(obj);
    }

    void draw();
    void drawWithNames();

    void fastDraw();
    void refreshDisplay();

private:
    GLUquadric    *qObj;

    vector<JSpherePtr>   spheres;
    vector<JBoxPtr>      boxes;
    vector<JConePtr>     cones;
    vector<JCylinderPtr> cylinders;
    vector<JTubePtr>     tubes;

    void drawCylinder( const float *p0, const float *p1, float radius) const;

    void drawObject( const JSpherePtr &s);
    void drawObject( const JBoxPtr &s);
    void drawObject( const JConePtr &s);
    void drawObject( const JCylinderPtr &s);
    void drawObject( const JTubePtr &s);
};
typedef boost::shared_ptr<JShapeViewer> JShapeViewerPtr;
