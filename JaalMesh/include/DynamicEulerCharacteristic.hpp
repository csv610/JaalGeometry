#pragma once

#ifndef JDYNEULER_H
#define JDYNEULER_H

#include "Mesh.hpp"
#include <boost/unordered_set.hpp>

class JDynamicEulerCharacteristic
{
public:
    void reset();
    int  addObject(const JMeshPtr &m);

    int addObject( const JCellPtr &);
    int addObject( const JFacePtr &);
    int addObject( const JEdgePtr &);
    int addObject( const JNodePtr &);

    int removeObject( const JCellPtr &);
    int removeObject( const JFacePtr &);
    int removeObject( const JEdgePtr &);
    int removeObject( const JNodePtr &);

    vector<size_t>  getFVector() const;

    int  getEulerCharacteristic() const
    {
        return nodeSet.size()-edgeSet.size()+ faceSet.size() - cellSet.size();
    }

private:
    long  birth_time;
    boost::unordered_set<JNodePtr>  nodeSet;
    boost::unordered_set<JEdgePtr>  edgeSet;
    boost::unordered_set<JFacePtr>  faceSet;
    boost::unordered_set<JCellPtr>  cellSet;
};

#endif
