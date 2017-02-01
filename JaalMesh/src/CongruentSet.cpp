#include "CongruentSet.hpp"
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////

void CongruentSet :: lowpass( vector<JFaceSequence> &newset)
{
    newset.clear();

    if( mesh == NULL ) return;

    size_t numfaces = mesh->getSize(2);
    map<int, JFaceSequence> fmap;
    for(size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( face->getSize(0) == 3 ) {
                int id = 1000*JFaceGeometry::getMinAngle(face, ANGLE_IN_RADIANS);
                fmap[id].push_back(face);
            }
        }
    }

    map<int,JFaceSequence>::const_iterator it;
    for( it = fmap.begin(); it != fmap.end(); ++it)
        newset.push_back( it->second );
}
////////////////////////////////////////////////////////////////////////////////

void CongruentSet :: highpass( JFaceSequence &faces, vector<JFaceSequence> &newset)
{
    newset.clear();
    if( mesh == NULL ) return;

    std::map<int,JFaceSequence>  fmap;

    size_t numfaces = faces.size();
    for(size_t i = 0; i < numfaces; i++) {
        int id  = 1000*JFaceGeometry::getMaxAngle(faces[i], ANGLE_IN_RADIANS);
        fmap[id].push_back( faces[i] );
    }

    map<int,JFaceSequence>::const_iterator it;
    for( it = fmap.begin(); it != fmap.end(); ++it)
        newset.push_back( it->second );
}

///////////////////////////////////////////////////////////////////////////////

void CongruentSet :: apply()
{
    vector<JFaceSequence> aset, bset;

    lowpass( aset );

    congruentSet.clear();
    for(size_t i = 0; i < aset.size(); i++) {
        highpass( aset[i], bset);
        for( size_t j = 0; j < bset.size(); j++)
            congruentSet.push_back( bset[j] );
    }
    sort( congruentSet.begin(), congruentSet.end(), CompareSetSize() );
}

///////////////////////////////////////////////////////////////////////////////
