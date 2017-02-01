#include "MeshQuality.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////
void JStatisticalInfo :: setData( vector<double> &values)
{
    numSamples = values.size();
    if( numSamples == 0) return;

    boost::sort( values );
    minVal = *boost::min_element( values );
    maxVal = *boost::max_element( values );
    sum     = std::accumulate(values.begin(), values.end(), 0.0);
    avgVal  = sum/(double)values.size();
    if (numSamples % 2 == 0)
        medianVal = (values[numSamples/ 2 - 1] + values[numSamples/ 2]) / 2;
    else
        medianVal = values[numSamples/ 2];
}

////////////////////////////////////////////////////////////////////////////////
JLogger* JMeshQuality::logger = JLogger::getInstance();

double JMeshQuality :: getEdgeLength( const JEdgePtr &edge)
{
    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);

    double val = JMath::length( v0->getXYZCoords(), v1->getXYZCoords() );

    if( !attribname.empty() ) edge->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

JStatisticalInfo JMeshQuality :: getEdgesLength() const
{
    JStatisticalInfo jstat;
    if( mesh == nullptr) return jstat;

    logger->setInfo("Collecting minimum edge length");

    vector<double> val;

    size_t numedges = mesh->getSize(1);
    for( size_t i  = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            double l2  = JEdgeGeometry::getLength(edge);
            val.push_back(l2);
        }
    }
    jstat.setData(val);
    return jstat;
}

//////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getEdgeWeight( const JEdgePtr &edge, int dir)
{
    double val;
    if( dir == 1 ) {
        if( edge->hasAttribute("EdgeWeight12"))
            edge->getAttribute("EdgeWeight12", val);
    }

    if( dir == -1) {
        if( edge->hasAttribute("EdgeWeight21") )
            edge->getAttribute("EdgeWeight21", val);
    }
    return val;
}

////////////////////////////////////////////////////////////////////////////////
void JMeshQuality :: setRange()
{
}

////////////////////////////////////////////////////////////////////////////////
double JMeshQuality :: getAspectRatio( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_aspect( num_nodes, coords );
        break;
    case 4:
        val = v_quad_aspect( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////


double JMeshQuality :: getAspectBeta( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);

    if( num_nodes != 4) {
        cout << "Warning: Aspect beta quality only for tets " << endl;
        return val;
    }

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_tet_aspect_beta( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getAspectGamma( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    if( num_nodes != 4) {
        cout << "Warning: Aspect beta quality only for tets " << endl;
        return val;
    }

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_tet_aspect_gamma( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getAspectRatio( const JCellPtr &cell)
{
    double val = 0.0;

    int num_nodes = cell->getSize(0);
    if( num_nodes != 8) return val;

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_hex_aspect( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getArea( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_area( num_nodes, coords );
        break;
    case 4:
        val = v_quad_area( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getConditionNumber( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_condition( num_nodes, coords );
        break;
    case 4:
        val = v_quad_condition( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getConditionNumber( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_condition( num_nodes, coords );
        break;
    case 8:
        val = v_hex_condition( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getDiagonalRatio( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);

    if( num_nodes != 8)  {
        cout << "Warning: Diagonal ratio only for hex element" << endl;
        return val;
    }

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_hex_diagonal( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getOddy( const JFacePtr &face)
{
    double val;
    int num_nodes = face->getSize(0);

    if( num_nodes != 4)  {
        cout << "Warning: Oddy only for quad elements" << endl;
        return val;
    }

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_quad_oddy( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getOddy( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);

    if( num_nodes != 8)  {
        cout << "Warning: Oddy is only for hex element" << endl;
        return val;
    }

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_hex_oddy( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getMinAngle( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_minimum_angle( num_nodes, coords );
        break;
    case 4:
        val = v_quad_minimum_angle( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getMaxAngle( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_maximum_angle( num_nodes, coords );
        break;
    case 4:
        val = v_quad_maximum_angle( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getDistortion( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_distortion( num_nodes, coords );
        break;
    case 4:
        val = v_quad_distortion( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////
double JMeshQuality :: getDistortion( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_distortion( num_nodes, coords );
        break;
    case 8:
        val = v_hex_distortion( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getJacobian( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch( num_nodes) {
    case 3:
        val = v_tri_scaled_jacobian( num_nodes, coords );
        break;
    case 4:
        val = v_quad_jacobian( num_nodes, coords );
        break;
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getJacobian( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_jacobian( num_nodes, coords );
        break;
    case 8:
        val = v_hex_jacobian( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getScaledJacobian( const JFacePtr &face)
{
    double val;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_scaled_jacobian( num_nodes, coords );
        break;
    case 4:
        val = v_quad_scaled_jacobian( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getScaledJacobian( const JCellPtr &cell)
{
    double val;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_scaled_jacobian( num_nodes, coords );
        break;
    case 8:
        val = v_hex_scaled_jacobian( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShear( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 4) {
        cout << "Warning: Shear is only for quad element " << endl;
        return val;
    }

    val = v_quad_shear( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShear( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 8) {
        cout << "Warning: Shear is only for hex element " << endl;
        return val;
    }
    val = v_hex_shear( num_nodes, coords );
    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;

}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShearSize( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 4 ) {
        cout << "Warning: Shear Size only for Quad element" << endl;
        return val;
    }

    val = v_quad_shear_and_size( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);

    return val;

}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShearSize( const JCellPtr &cell)
{
    double val = 0.0;

    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes  != 8) {
        cout << "Warning: ShearSize only for hex element " << endl;
        return val;
    }
    val = v_hex_shear_and_size( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);

    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getRelativeSize2( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_relative_size_squared( num_nodes, coords );
        break;
    case 4:
        val = v_quad_relative_size_squared( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getRelativeSize2( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_relative_size_squared( num_nodes, coords );
        break;
    case 8:
        val = v_hex_relative_size_squared( num_nodes, coords );
        break;
    default:
        return 1;
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShape( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_shape( num_nodes, coords );
        break;
    case 4:
        val = v_quad_shape( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////
double JMeshQuality :: getShape( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_shape( num_nodes, coords );
        break;
    case 8:
        val = v_hex_shape( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShapeSize( const JFacePtr &face)
{
    double val = 0.0;

    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 3:
        val = v_tri_shape_and_size( num_nodes, coords );
        break;
    case 4:
        val = v_quad_shape_and_size( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;

}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getShapeSize( const JCellPtr &cell )
{
    double val = 0.0;

    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch(num_nodes) {
    case 4:
        val = v_tet_shape_and_size( num_nodes, coords );
        break;
    case 8:
        val = v_hex_shape_and_size( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getSkew( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 4 ) return val;

    val = v_quad_skew( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getSkew( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 8) return val;
    val = v_hex_skew( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getStretch( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 4 ) return val;
    val = v_quad_stretch( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getStretch( const JCellPtr &cell)
{
    double val = 0.0;

    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 8) return val;
    val = v_hex_stretch( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getTaper( const JFacePtr &face)
{
    double val = 0.0;
    int num_nodes = face->getSize(0);
    if( num_nodes != 4 ) return val;

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_quad_taper( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}
////////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getTaper( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    if( num_nodes != 8) return val;

    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    val = v_hex_taper( num_nodes, coords );

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}

///////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getVolume( const JCellPtr &cell)
{
    double val = 0.0;
    int num_nodes = cell->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        const JNodePtr &vertex = cell->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    switch( num_nodes ) {
    case 4:
        val = v_tet_volume( num_nodes, coords );
        break;
    case 8:
        val = v_hex_volume( num_nodes, coords );
        break;
    default:
        JNoImpl();
    }

    if( !attribname.empty() ) cell->setAttribute(attribname, val);
    return val;
}

///////////////////////////////////////////////////////////////////////////////

double JMeshQuality :: getWarpage( const JFacePtr &face)
{
    double val = 0.0;

    int num_nodes = face->getSize(0);
    for( int j = 0; j < num_nodes; j++) {
        JNodePtr vertex = face->getNodeAt(j);
        xyz  = vertex->getXYZCoords();
        coords[j][0] = xyz[0];
        coords[j][1] = xyz[1];
        coords[j][2] = xyz[2];
    }

    if( num_nodes != 4 ) {
        cout << "Warning: Warpage is only for quad element" << endl;
        return val;
    }
    val = v_quad_warpage( num_nodes, coords );

    if( !attribname.empty() ) face->setAttribute(attribname, val);
    return val;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshQuality :: registerQuality( const JEdgePtr &edge, double val, int pos)
{
    switch( pos )
    {
    case JMeshEntity::ANY_ENTITY:
        quality.push_back(val);
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        if( edge->isBoundary() ) {
            quality.push_back(val);
        }
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        if( !edge->isBoundary() ) {
            quality.push_back(val);
        }
        break;
    }
}
///////////////////////////////////////////////////////////////////////////////

vector<double> JMeshQuality :: getEdgesQuality( int name, int pos, int sorted)
{
    quality.clear();

    int stat;
    double val;

    if( name == EDGE_LENGTH ) {
        size_t numedges = mesh->getSize(1);
        quality.reserve(numedges);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                registerQuality( edge, getEdgeLength(edge), pos);
            }
        }
        if( sorted ) boost::sort(quality );
    }

    if( name == EDGE_WEIGHT ) {
        size_t numedges = mesh->getSize(1);
        quality.reserve(numedges);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                registerQuality( edge, getEdgeWeight(edge, 1), pos);
                registerQuality( edge, getEdgeWeight(edge,-1), pos);
            }
        }
        if( sorted ) boost::sort( quality );
    }
    return quality;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshQuality :: registerQuality( const JFacePtr &face, double val, int pos)
{
    switch( pos )
    {
    case JMeshEntity::ANY_ENTITY:
        quality.push_back(val);
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        if( face->isBoundary() ) {
            quality.push_back(val);
        }
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        if( !face->isBoundary() ) {
            quality.push_back(val);
        }
        break;
    }
}
///////////////////////////////////////////////////////////////////////////////

vector<double> JMeshQuality :: getFacesQuality( int name, int pos, int sorted)
{
    quality.clear();

    JFaceSequence faces;
    switch( pos )
    {
    case JMeshEntity::ANY_ENTITY:
        faces = mesh->getFaces();
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        mesh->getTopology()->getBoundary(faces);
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        mesh->getTopology()->getInternal(faces);
        break;
    }

    if( faces.empty() ) return quality;

    size_t numfaces = faces.size();

    quality.resize(numfaces);

    switch( name )
    {
    case AREA:
        for( size_t i = 0; i < numfaces; i++)
            quality[i]  = getArea( faces[i] );
        break;
    case ASPECT_RATIO:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getAspectRatio( faces[i] );
        break;
    case CONDITION_NUMBER:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getConditionNumber(faces[i]);
        break;
    case DISTORTION:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getDistortion(faces[i]);
        break;
    case MIN_ANGLE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getMinAngle(faces[i]);
        break;
    case MAX_ANGLE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getMaxAngle(faces[i]);
        break;
    case JACOBIAN:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getJacobian(faces[i]);
        break;
    case  SCALED_JACOBIAN:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getScaledJacobian(faces[i]);
        break;
    case ODDY:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getOddy(faces[i]);
        break;
    case RELATIVE_SIZE_SQUARED:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getRelativeSize2(faces[i]);
        break;
    case SHAPE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getShape(faces[i]);
        break;
    case SHAPE_AND_SIZE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getShapeSize( faces[i] );
        break;
    case SHEAR:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] =  getShear(faces[i]);
        break;
    case SHEAR_AND_SIZE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getShearSize(faces[i] );
        break;
    case SKEW:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getSkew(faces[i]);
        break;
    case STRETCH:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getStretch( faces[i]);
        break;
    case TAPER:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getTaper( faces[i]);
        break;
    case  WARPAGE:
        for( size_t i = 0; i < numfaces; i++)
            quality[i] = getWarpage(faces[i]);
        break;
    }

    if( sorted ) boost::sort( quality );

    return quality;
}

///////////////////////////////////////////////////////////////////////////

vector<double> JMeshQuality :: getCellsQuality( int name, int pos, int sorted)
{
    quality.clear();

    size_t numcells = mesh->getSize(3);

    quality.reserve(numcells);

    int stat;
    double val;

    if( name == ASPECT_BETA ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getAspectBeta( cell );
                quality.push_back(val);
            }
        }
    }

    if( name == ASPECT_GAMMA ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getAspectGamma( cell );
                quality.push_back(val);
            }
        }
    }


    if( name == ASPECT_RATIO ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getAspectRatio( cell );
                quality.push_back(val);
            }
        }
    }

    if( name == CONDITION_NUMBER ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getConditionNumber(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == DIAGONAL_RATIO ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getDiagonalRatio(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == DISTORTION ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getDistortion(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == JACOBIAN ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getJacobian(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == ODDY ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getOddy(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == RELATIVE_SIZE_SQUARED ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getRelativeSize2(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == SCALED_JACOBIAN ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getScaledJacobian(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == SHAPE ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getShape(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == SHAPE_AND_SIZE ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getShapeSize(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == SHEAR ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getShear(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == SHEAR_AND_SIZE ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getShearSize(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == SKEW ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val  = getSkew(cell);
                quality.push_back(val);
            }
        }
    }

    if( name == STRETCH ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive()) {
                val  = getStretch(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == TAPER ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive()) {
                val  = getTaper(cell );
                quality.push_back(val);
            }
        }
    }

    if( name == VOLUME ) {
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                val = getVolume(cell );
                quality.push_back(val);
            }
        }
    }
    if( sorted ) boost::sort( quality );
    return quality;
}
