#pragma once

#ifndef CONGRUENT_H
#define CONGRUENT_H

#include "Mesh.hpp"

using namespace Jaal;

class CongruentSet
{
public:
    CongruentSet() {
        tol  = 0.0;
    }

    void setMesh( JMeshPtr m) {
        mesh = m;
    }

    void setTolerance( double t) {
        tol = t;
    }

    size_t  getSetSize() const {
        return congruentSet.size();
    }
    JFaceSequence getSet(size_t id) {
        if( id < congruentSet.size() ) return congruentSet[id];
        JFaceSequence empty;
        return empty;
    }

    void apply();

private:
    JMeshPtr mesh;
    double tol;

    vector<JFaceSequence> congruentSet;

    struct CompareSetSize {
        bool operator() (const JFaceSequence &lhs, const JFaceSequence &rhs) const {
            return lhs.size() >  rhs.size();
        }
    };


    void lowpass( vector<JFaceSequence> &newset);
    void highpass( JFaceSequence &s, vector<JFaceSequence> &newset);
};

#endif
