#pragma once

#include <vector>
#include <iosfwd>
#include <string>

namespace DDG
{
class JOptimalVectorField
{
    typedef std::pair<int,double> Singularity;
    typedef std::pair<int,double> Generator;
public:
    void setMesh( const JMeshPtr &m);
    void solve( void ) const;      // solve the problem

protected:
    double fieldAngle;
    std::vector<Singularity> singularities;
    std::vector<Generator> generators;


};
}
