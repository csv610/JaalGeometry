#include "MeshPartitioner.hpp"
#include "MeshTopology.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

void JMeanShiftSegmentation :: getNeighbours( const JFacePtr &seed, JFaceSequence &seq)
{
    deque<JFacePtr> faceQ;
    faceQ.push_back(seed);

    JFaceSet faceSet;
    faceSet.insert(seed) ;

    JFaceSequence neighs;
    while(!faceQ.empty()) {
        const JFacePtr  currface = faceQ.front();
        faceQ.pop_front();
        JFace::getRelations02(currface, neighs);
        for( JFacePtr &nextface: neighs) {
            if( faceSet.find(nextface) == faceSet.end() ) {
                faceSet.insert(nextface);
                faceQ.push_back(nextface);
            }
        }
        if( faceSet.size() > numSamples) break;
    }
    seq.clear();
    seq.resize( faceSet.size() );
    copy( faceSet.begin(), faceSet.end(), back_inserter(seq));
}

///////////////////////////////////////////////////////////////////////////////

int JMeanShiftSegmentation :: getMeanValue(const JFacePtr &seed)
{
    JFaceSequence neighs;
    getNeighbours( seed, neighs);

    int n = neighs.size();
    vector<int>  faceValues(n);
    int index = 0;
    for( const JFacePtr &face: neighs)
        face->getAttribute(attribname, faceValues[index++]);

    sort( faceValues.begin(), faceValues.end() );
    if( n%2 ) return faceValues[n/2];

    return faceValues[(n+1)/2];
}

///////////////////////////////////////////////////////////////////////////////

void JMeanShiftSegmentation :: segment( const JMeshPtr &m, const string &s, int numIterations)
{
    mesh = m;
    attribname = s;
    if( mesh == nullptr) return;
    numSamples = 20;

    int val1, val2, valueChanged;
    size_t numfaces = mesh->getSize(2);
    for( int i = 0; i < numIterations; i++) {
        valueChanged = 0;
        for( size_t j = 0; j < numfaces; j++) {
            const JFacePtr &face = mesh->getFaceAt(j);
            face->getAttribute(attribname, val1);
            val2 = getMeanValue(face);
            if( val1 != val2) {
                face->setAttribute(attribname, val2);
                valueChanged = 1;
            }
        }
        if( !valueChanged ) break;
    }
}
///////////////////////////////////////////////////////////////////////////////
