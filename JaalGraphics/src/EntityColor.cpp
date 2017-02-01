#include "EntityColor.hpp"
#include "DrawMeshEntity.hpp"

using namespace Jaal;

//////////////////////////////////////////////////////////////////////////////

JColor JEntityColor :: getRandomColor()
{
    JColor clr;
    double r = drand48();
    r  = max(0.2, r);
    r  = min(0.9, r);

    double g = drand48();
    g  = max(0.2, g);
    g  = min(0.9, g);

    double b = drand48();
    b  = max(0.2, b);
    b  = min(0.9, b);

    clr[0] = r;
    clr[1] = g;
    clr[2] = b;
    clr[3] = 1.0;

    return clr;
}

void JColorMap :: jet( const vector<double> &darray, vector<JColor> &colors)
{
    double minVal = *boost::min_element(darray );
    double maxVal = *boost::max_element(darray );

    size_t nsize = darray.size();
    colors.resize(nsize);
    double r, g, b, val, p;

    JColor clr;
    for( size_t i = 0; i < nsize; i++) {
        val = darray[i];
        p = (val-minVal)/(maxVal-minVal);
        if( p < 0.5) {
            r = 1.0 - 2.0*p;
            g = 2.0*p;
            b = 0.0;
        } else {
            p = p - 0.5;
            r = 0.0;
            g = 1.0 - 2.0*p;
            b = 2.0*p;
        }
        clr[0] = r;
        clr[1] = g;
        clr[2] = b;
        clr[3] = 1.0;
        colors[i] = clr;
    }
}

//////////////////////////////////////////////////////////////////////////////

JColor JEntityColor :: getMinColor( int id )
{
    switch( id )  {
    case 0:
        return getColor( "Red");
        break;
    case 1:
        return getColor( "Green");
        break;
    case 2:
        return getColor( "Blue");
        break;
    case 3:
        return getColor( "Yellow");
        break;
    case 4:
        return getColor( "Purple");
        break;
    case 5:
        return getColor( "Lime");
        break;
    case 6:
        return getColor( "Moroon");
        break;
    case 7:
        return getColor( "Navy");
        break;
    case 8:
        return getColor( "Aqua");
        break;
    case 9:
        return getColor( "Olive");
        break;
    case 10:
        return getColor( "Gray");
        break;
    case 11:
        return getColor( "Pink");
        break;
    }
    cout << "No Min Color Provided " << endl;
}
//////////////////////////////////////////////////////////////////////////////

JColor JEntityColor ::getColor( const string &name)
{
    JColor color;

    color[0] = 0.5;
    color[1] = 0.5;
    color[2] = 0.5;
    color[3] = 1.0;

    if( name == "Aqua") {
        color[0] = Aqua[0];
        color[1] = Aqua[1];
        color[2] = Aqua[2];
        return color;
    }

    if( name == "Black") {
        color[0] = Black[0];
        color[1] = Black[1];
        color[2] = Black[2];
        return color;
    }

    if( name == "Blue") {
        color[0] = Blue[0];
        color[1] = Blue[1];
        color[2] = Blue[2];
        return color;
    }

    if( name == "Gray") {
        color[0] = Gray[0];
        color[1] = Gray[1];
        color[2] = Gray[2];
        return color;
    }

    if( name == "Green") {
        color[0] = Green[0];
        color[1] = Green[1];
        color[2] = Green[2];
        return color;
    }

    if( name == "Lime") {
        color[0] = Lime[0];
        color[1] = Lime[1];
        color[2] = Lime[2];
        return color;
    }

    if( name == "Magenta") {
        color[0] = Fuchsia[0];
        color[1] = Fuchsia[1];
        color[2] = Fuchsia[2];
        return color;
    }

    if( name == "Moroon") {
        color[0] = Moroon[0];
        color[1] = Moroon[1];
        color[2] = Moroon[2];
        return color;
    }

    if( name == "Navy") {
        color[0] = Navy[0];
        color[1] = Navy[1];
        color[2] = Navy[2];
        return color;
    }

    if( name == "Olive") {
        color[0] = Olive[0];
        color[1] = Olive[1];
        color[2] = Olive[2];
        return color;
    }

    if( name == "Orange") {
        color[0] = 1.0;
        color[1] = 0.6470;
        color[2] = 0.0;
        return color;
    }

    if( name == "Pink") {
        color[0] = 1.0;
        color[1] = 0.7529;
        color[2] = 0.7960;
        return color;
    }

    if( name == "Purple") {
        color[0] = Purple[0];
        color[1] = Purple[1];
        color[2] = Purple[2];
        return color;
    }

    if( name == "Red") {
        color[0] = Red[0];
        color[1] = Red[1];
        color[2] = Red[2];
        return color;
    }

    if( name == "Silver") {
        color[0] = Silver[0];
        color[1] = Silver[1];
        color[2] = Silver[2];
        return color;
    }

    if( name == "Teal") {
        color[0] = Teal[0];
        color[1] = Teal[1];
        color[2] = Teal[2];
        return color;
    }

    if( name == "Yellow") {
        color[0] = Yellow[0];
        color[1] = Yellow[1];
        color[2] = Yellow[2];
        return color;
    }

    if( name == "White") {
        color[0] = White[0];
        color[1] = White[1];
        color[2] = White[2];
        return color;
    }

    if( name == "LightBlue") {
        color[0] = 173.0/255.0;
        color[1] = 216.0/255.0;
        color[2] = 230.0/255.0;
        return color;
    }

    if( name == "LightGreen") {
        color[0] = 127.0/255.0;
        color[1] = 255.0/255.0;
        color[2] = 0.0;
        return color;
    }

    if( name == "LightGray") {
        color[0] = 119.0/255.0;
        color[1] = 136.0/255.0;
        color[2] = 153.0/255.0;
        return color;
    }

    return color;
}

//////////////////////////////////////////////////////////////////////////////////

int JNodeColor :: setMesh(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return 1;

    size_t nSize = mesh->getSize(0);
    for(size_t i = 0; i < nSize; i++) assign(mesh->getNodeAt(i) );
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////


int JEdgeColor :: setMesh(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;
    size_t nSize = mesh->getSize(1);
    for(size_t i = 0; i < nSize; i++)
        assign(mesh->getEdgeAt(i) );
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JFaceColor :: setMesh(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;
    size_t nSize = mesh->getSize(2);
    for(size_t i = 0; i < nSize; i++)
        assign(mesh->getFaceAt(i) );
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JCellColor :: setMesh(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;
    size_t nSize = mesh->getSize(3);
    for(size_t i = 0; i < nSize; i++)
        assign(mesh->getCellAt(i) );
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
AngleDefectColor ::  AngleDefectColor()
{
    highlightColor[0]  = 1.0;
    highlightColor[1]  = 0.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = alpha;
}

///////////////////////////////////////////////////////////////////////////////

int AngleDefectColor :: assign(const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;

    JNodeRenderPtr nAttrib = nullptr;

    float angle = 0.0;
    if( vertex->hasAttribute("AngleDefect")) {
        vertex->getAttribute("AngleDefect", angle);
        vertex->getAttribute("Render", nAttrib);
        if( angle > 0) {
            highlightColor[0]  = 1.0;
            highlightColor[1]  = 0.0;
            highlightColor[2]  = 0.0;
            nAttrib->color = highlightColor;
        }
        if( angle < 0) {
            highlightColor[0]  = 0.0;
            highlightColor[1]  = 1.0;
            highlightColor[2]  = 0.0;
            nAttrib->color = highlightColor;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void NodeHeightColor :: setMesh( const JMeshPtr &mesh, int dir_)
{
    dir = dir_;
    size_t nSize =  mesh->getSize(0);
    if( nSize == 0) return;

    vector<double> darray;
    darray.reserve(nSize);
    for( size_t i = 0; i < nSize; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point3D  &xyz  = vtx->getXYZCoords();
            darray.push_back(xyz[dir]);
        }
    }
    JNodeColor::assign(mesh, darray);
}

int NodeHeightColor :: assign( const JNodePtr & )
{
    cout << "Not implemented" << endl;
    return 1;
}

int ManifoldEdgeColor :: assign( const JEdgePtr &edge)
{
    if( edge == nullptr ) return 1;
    if( !edge->isActive() )  return 2;

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);

    color[3] = alpha;

    int nsize = edge->getNumRelations(2);
    if( nsize == 0) return 1;

    if( nsize > 2) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
    } else {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
    }

    eAttrib->color = color;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshCellColor :: JMeshCellColor()
{
    defaultColor[0] = 0.0;
    defaultColor[1] = 0.0;
    defaultColor[2] = 1.0;
    defaultColor[3] = alpha;
}

int JMeshCellColor :: assign(const JCellPtr &cell)
{
    if( cell ) {
        cell->setAttribute("Color", defaultColor);
        return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

JMeshNodeColor :: JMeshNodeColor()
{
    alpha    = 1.0;
    defaultColor[0] = 0.0;
    defaultColor[1] = 0.0;
    defaultColor[2] = 6.0;
    defaultColor[3] = alpha;

    boundColor[0]   = 1.0;
    boundColor[1]   = 0.0;
    boundColor[2]   = 0.6;
    boundColor[3]   = alpha;

    display_boundary  = 1;
    display_internal  = 1;
    display_interface = 0;
    display_features  = 0;
}

void JMeshNodeColor :: setInternalColor( float *rgb )
{
    defaultColor[0] = rgb[0];
    defaultColor[1] = rgb[1];
    defaultColor[2] = rgb[2];
    defaultColor[3] = alpha;
}

void JMeshNodeColor :: setBoundaryColor( float *rgb)
{
    boundColor[0] =  rgb[0];
    boundColor[1] =  rgb[1];
    boundColor[2] =  rgb[2];
    boundColor[3] =  alpha;
}

void JMeshNodeColor :: setFeatureColor( float *rgb)
{
    featureColor[0] =  rgb[0];
    featureColor[1] =  rgb[1];
    featureColor[2] =  rgb[2];
    featureColor[3] =  alpha;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNodeColor :: setInterfaceColor( float *rgb)
{
    interfaceColor[0] =  rgb[0];
    interfaceColor[1] =  rgb[1];
    interfaceColor[2] =  rgb[2];
    interfaceColor[3] =  alpha;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshNodeColor :: assign( const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;
    if( !vertex->isActive() ) return 1;

    JNodeRenderPtr nodeAttrib;
    vertex->getAttribute("Render", nodeAttrib);
    if( nodeAttrib == nullptr ) return 1;

    /*
        if( display_features ) {
          if( vertex->isFeature() ) {
                nodeAttrib->color = featureColor;
                return 0;
            }
        }
    */

    if( display_interface ) {
        if( vertex->hasAttribute("Interface") ) {
            nodeAttrib->color = interfaceColor;
            return 0;
        }
    }

    if( display_boundary &&  vertex->isBoundary() ) {
        nodeAttrib->color = boundColor;
        return 0;
    }

    nodeAttrib->color = defaultColor;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

IrregularNodeColor :: IrregularNodeColor()
{
    defaultColor[0] = 0.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 0.0;
    defaultColor[3] = alpha;
}
///////////////////////////////////////////////////////////////////////////////
int IrregularNodeColor :: buildRelations(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return 1;

    mesh->buildRelations(0,2);
    mesh->buildRelations(0,3);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int IrregularNodeColor :: assign(const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;
    if( !vertex->isActive() ) return 2;

    JNodeRenderPtr nAttrib;
    vertex->getAttribute("Render", nAttrib);

    nAttrib->display = 1;
    if( vertex->isBoundary() ) {
        nAttrib->color = defaultColor;
        return 0;
    }

    int etype = 0, degree = 0;

    getFaceType(vertex, etype, degree);
    int ideal;
    switch( etype ) {
    case JFace::QUADRILATERAL:
        ideal = 4;
        break;
    case JFace::TRIANGLE:
        ideal = 6;
        break;
    }
    nAttrib->display = 1;

    color[3] = alpha;
    if( degree < ideal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        nAttrib->color =  color;
        return 0;
    }

    if( degree == ideal ) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
        nAttrib->color =  color;
        return 0;
    }

    if( degree > ideal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        nAttrib->color =  color;
        return 0;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

JMeshEdgeColor :: JMeshEdgeColor()
{
    defaultColor[0] = 0.2;
    defaultColor[1] = 100.0/255.0;
    defaultColor[2] = 0.2;
    defaultColor[3] = alpha;

    boundColor[0] = 1.0;
    boundColor[1] = 0.5;
    boundColor[2] = 0.5;
    boundColor[3] = alpha;

    interfaceColor[0] = 1.0;
    interfaceColor[1] = 1.0;
    interfaceColor[2] = 1.0;
    interfaceColor[3] = alpha;

    featureColor[0] = 1.0;
    featureColor[1] = 0.0;
    featureColor[2] = 0.0;
    featureColor[3] = alpha;

    display_boundary  = 1;
    display_internal  = 1;
    display_interface = 1;
    display_features  = 0;
}

void JMeshEdgeColor :: setInternalColor( float *rgb )
{
    defaultColor[0] = rgb[0];
    defaultColor[1] = rgb[1];
    defaultColor[2] = rgb[2];
    defaultColor[3] = alpha;
}

void JMeshEdgeColor :: setBoundaryColor( float *rgb)
{
    boundColor[0] = rgb[0];
    boundColor[1] = rgb[1];
    boundColor[2] = rgb[2];
    boundColor[3] = alpha;
}

void JMeshEdgeColor :: setFeatureColor( float *rgb)
{
    featureColor[0] =  rgb[0];
    featureColor[1] =  rgb[1];
    featureColor[2] =  rgb[2];
    featureColor[3] =  alpha;
}

void JMeshEdgeColor :: setInterfaceColor( float *rgb)
{
    interfaceColor[0] =  rgb[0];
    interfaceColor[1] =  rgb[1];
    interfaceColor[2] =  rgb[2];
    interfaceColor[3] =  alpha;
}

int JMeshEdgeColor :: assign(const JEdgePtr &edge)
{
    if( edge == nullptr ) return 1;
    if( !edge->isActive()  ) return 1;

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);
    if( eAttrib == nullptr ) return 1;

    if( display_interface ) {
        if( edge->hasAttribute("Interface") ) {
            eAttrib->color = interfaceColor;
            return 0;
        }
    }

    if( display_boundary ) {
        if( edge->isBoundary() ) {
            eAttrib->color = boundColor;
            return 0;
        }
    }

    eAttrib->color = defaultColor;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

QuadParamLinesColor ::  QuadParamLinesColor()
{
    xcolor[0] = 1.0;
    xcolor[1] = 0.0;
    xcolor[2] = 0.0;
    xcolor[3] = alpha;

    ycolor[0] = 0.0;
    ycolor[1] = 1.0;
    ycolor[2] = 0.0;
    ycolor[3] = alpha;

    defaultcolor[0] = 1.0;
    defaultcolor[1] = 0.5;
    defaultcolor[2] = 0.5;
    defaultcolor[3] = alpha;
}

int QuadParamLinesColor :: assign(const JEdgePtr &edge)
{
    bool param;
    int err = edge->getAttribute("ParametricLine", param);
    if( !err ) {
        if( param == 0) {
            color[0] = xcolor[0];
            color[1] = xcolor[1];
            color[2] = xcolor[2];
        } else {
            color[0] = ycolor[0];
            color[1] = ycolor[1];
            color[2] = ycolor[2];
        }
        color[3] = alpha;
        edge->setAttribute("Color", color);
        return 0;
    }
    return 1;
}


////////////////////////////////////////////////////////////////////////////////

JMeshFaceColor :: JMeshFaceColor()
{
    alpha    = 0.75;
    defaultColor[0] = 0.0;
    defaultColor[1] = 175.0/255.0;
    defaultColor[2] = 1.0;
    defaultColor[3] = alpha;

    boundColor[0] = 0.3;
    boundColor[1] = 0.4;
    boundColor[2] = 0.5;
    boundColor[3] = alpha;

    display_boundary  = 1;
    display_internal  = 1;
    display_interface = 1;

    cell = nullptr;
}

void JMeshFaceColor :: setFrontColor( float *rgb)
{
    frontColor[0] = rgb[0];
    frontColor[1] = rgb[1];
    frontColor[2] = rgb[2];
    frontColor[3] = alpha;
}

void JMeshFaceColor :: setBackColor( float *rgb)
{
    backColor[0] = rgb[0];
    backColor[1] = rgb[1];
    backColor[2] = rgb[2];
    backColor[2] = alpha;
}

int JMeshFaceColor :: assign(const JFacePtr &face)
{
    if( face == nullptr ) return 1;
    if( !face->isActive()  ) return 1;

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);

    if( face->isBoundary() ) {
        fAttrib->color = boundColor;
        return 0;
    }
    fAttrib->color = defaultColor;
    return 0;
}

void JMeshFaceColor :: setInternalColor( float *rgb )
{
    defaultColor[0] = rgb[0];
    defaultColor[1] = rgb[1];
    defaultColor[2] = rgb[2];
    defaultColor[3] = alpha;
}

void JMeshFaceColor :: setBoundaryColor( float *rgb)
{
    boundColor[0] = rgb[0];
    boundColor[1] = rgb[1];
    boundColor[2] = rgb[2];
    boundColor[3] = alpha;
}

void JMeshFaceColor :: setInterfaceColor( float *rgb)
{
    interfaceColor[0] =  rgb[0];
    interfaceColor[1] =  rgb[1];
    interfaceColor[2] =  rgb[2];
    interfaceColor[3] =  alpha;
}

////////////////////////////////////////////////////////////////////////////////

/*
MeshCellFaceColor :: MeshCellFaceColor()
{
     colors.resize(8);
     colors[0][0] = 1.0;
     colors[0][1] = 0.0;
     colors[0][2] = 0.0;
     colors[0][3] = alpha;

     colors[1][0] = 0.0;
     colors[1][1] = 1.0;
     colors[1][2] = 0.0;
     colors[1][3] = alpha;

     colors[2][0] = 0.0;
     colors[2][1] = 0.0;
     colors[2][2] = 1.0;
     colors[2][3] = alpha;

     colors[3][0] = 1.0;
     colors[3][1] = 1.0;
     colors[3][2] = 0.0;
     colors[3][3] = alpha;

     colors[4][0] = 0.0;
     colors[4][1] = 1.0;
     colors[4][2] = 1.0;
     colors[4][3] = alpha;

     colors[5][0] = 1.0;
     colors[5][1] = 0.0;
     colors[5][2] = 1.0;
     colors[5][3] = alpha;

     colors[6][0] = 0.5;
     colors[6][1] = 0.5;
     colors[6][2] = 1.0;
     colors[6][3] = alpha;

     colors[7][0] = 1.0;
     colors[7][1] = 0.5;
     colors[7][2] = 0.5;
     colors[7][3] = alpha;
}
*/

/*
int MeshCellFaceColor :: getColor( const Face *face, Color &clr)
{
     assert( cell != nullptr );
     int id =  cell->getPosOf( face );
     if( id < 0 || id >= 8) return 1;
     clr = colors[id];
     return SUCCESS;
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
int  EdgeQualityColor :: assign(const JEdgePtr &edge )
{
    edge->getAttribute(attribname, val);

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);

    if( val < lowerVal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = alpha;
        eAttrib->color = color;
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        color[3] = alpha;
        eAttrib->color = color;
        return 0;
    }

    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
    color[3] = alpha;
    eAttrib->color = color;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int  FaceQualityColor :: assign(const JFacePtr &face )
{
    face->getAttribute(attribname, val);

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);

    if( val < lowerVal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = alpha;
        fAttrib->color = color;
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        color[3] = alpha;
        fAttrib->color = color;
        return 0;
    }

    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
    color[3] = alpha;
    fAttrib->color = color;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int CellQualityColor :: assign( const JCellPtr &cell)
{
    cell->getAttribute(attribname, val);

    color[3] = 0.0;
    if( val < lowerVal ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        cell->setAttribute("Color", color);
        return 0;
    }

    if( val > upperVal ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        cell->setAttribute("Color", color);
        return 0;
    }
    color[0] = 0.0;
    color[1] = 1.0;
    color[2] = 0.0;
    cell->setAttribute("Color", color);
    return 0;
}
*/
///////////////////////////////////////////////////////////////////////////////

int QuadDefectColor :: assign( const JFacePtr &face)
{
    assert( face );
    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);
    fAttrib->color = color;
    return 0;
}


///////////////////////////////////////////////////////////////////////////////

void QuadDefectColor :: assignColors( const QDefectivePatch *newpatch, const JMeshPtr &mesh )
{
    JColor c;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            c[0] = 0.0;
            c[1] = 1.0;
            c[2] = 0.0;
            c[3] = alpha;
            face->setAttribute("Color", c);
        }
    }

    if( newpatch == nullptr ) return;

    JFaceSequence faces = newpatch->getFaces();
    numfaces = faces.size();
    for( size_t i = 0; i < numfaces; i++) {
        if( faces[i]->isActive() ) {
            c[0] = 1.0;
            c[1] = 0.0;
            c[2] = 0.0;
            c[3] = alpha;
            faces[i]->setAttribute("Color", c);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int PenroseTileColor :: assign( const JFacePtr &face)
{
    int val = 0;
    face->getAttribute("Tile", val);
    if( val == 0) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        face->setAttribute("Color", color);
        return 0;
    }
    if( val == 1) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
        face->setAttribute("Color", color);
        return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
#endif
