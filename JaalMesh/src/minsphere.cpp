#include <iterator>
#include "Miniball.hpp"

Sphere MeshGeometry :: getMinimumSphere() const
{
    double seed = 0;
    std::srand (seed);

    std::list<std::vector<double> > lp;
    size_t numnodes = mesh->getSize(0);
    for (size_t i=0; i< numnodes; ++i) {
        std::vector<double> p(3);
        JVertexPtr vtx = mesh->getNodeAt(i);
        const Point3D &xyz = vtx->getXYZCoords();
        p[0] =  xyz[0];
        p[1] =  xyz[1];
        p[2] =  xyz[2];
        lp.push_back(p);
    }

    // define the types of iterators through the points and their coordinates
    // ----------------------------------------------------------------------
    typedef std::list<std::vector<double> >::const_iterator PointIterator;
    typedef std::vector<double>::const_iterator CoordIterator;

    // create an instance of Miniball
    // ------------------------------
    typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> >
    MB;
    MB mb (3, lp.begin(), lp.end());

    Sphere sph;
    const double* center = mb.center();

    spr.setCenter( centr[0], center[1], center[2]);
    spr.setRadius( mb.radius() );

    return spr;
}
