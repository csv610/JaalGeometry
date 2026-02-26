#include "QuadCheckerBoard.hpp"

void JQuadCheckerBoard :: genPattern()
{
    size_t numfaces = mesh->getSize(2);
    char value = 'U';
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->setAttribute("CheckerBoard", value);
    }

    deque<JFacePtr> faceQ;

    value  = 'W';
    const JFacePtr &face0 = mesh->getFaceAt(0);
    face0->setAttribute("CheckerBoard", value);
    faceQ.push_back(face0);

    JFaceSequence faceneighs;

    char neighval, thisval;
    while(!faceQ.empty() ) {
        const JFacePtr &currface = faceQ.front();
        faceQ.pop_front();
        currface->getAttribute("CheckerBoard", thisval);
        JFace::getRelations12( currface, faceneighs);
        for( const JFacePtr &neigh : faceneighs) {
            neigh->getAttribute("CheckerBoard", neighval);
            if( neighval == 'U') {
                faceQ.push_back(neigh);
                if( thisval == 'W')
                    value = 'B';
                else
                    value = 'W';
                neigh->setAttribute("CheckerBoard", value);
            }
        }
    }
}
