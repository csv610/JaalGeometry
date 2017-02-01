#pragma once

#include "Mesh.hpp"
#include "MeshOptimization.hpp"


///////////////////////////////////////////////////////////////////////////////
class JDoublet : public JMeshTopologyOptimization
{
    static JLogger *logger;
public:
    static bool isDoublet(const JNodePtr &v);
    static void addAttribute(const JMeshPtr &m, const string &s = "Doublet" );

    int  removeAll();
    int  getSize();

    JNodeSequence getDoublets();
    JNodeSequence getDoublets( const JNodeSequence &v);

    int  remove(const JNodePtr &v);

    JNodePtr insert( const JFacePtr &f);
    JNodePtr insert(const JFacePtr &face, const JNodePtr &v0, const JNodePtr &v2);

private:
    bool searched = 0;
    JNodeSequence doublets;
    JFaceSequence faceneighs;
    void searchDoublets();
};
///////////////////////////////////////////////////////////////////////////////

