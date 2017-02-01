#ifndef JAAL_DUAL_MESH
#define JAAL_DUAL_MESH

#include "MeshViewer.hpp"
#include "../src/MeshDual.hpp"

using namespace qglviewer;
using namespace Jaal;

class JMeshQuadEditorViewer : public JMeshViewer {
public:
    JMeshQuadEditorViewer() {
        mesh = NULL;
    }

    JMeshQuadEditorViewer( JMeshPtr m ) {
        mesh = m;
    }

    ~JMeshQuadEditorViewer() { }

protected :
    virtual void draw();
    virtual void init();
    void   animate();
    virtual void keyPressEvent(QKeyEvent *event);

private:

};

#endif
