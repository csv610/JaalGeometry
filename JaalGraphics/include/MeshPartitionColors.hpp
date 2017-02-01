#pragma once

#include "EntityColor.hpp"
#include "MeshEntity.hpp"

class JNodePartitionColor : public JNodeColor {
public:

    string getName() const {
        return "NodePartition";
    }

    bool isPerNode() const {
        return 1;
    }

    int assign(const JNodePtr &vtx);

private:
    map<int,JColor> colormap;
};

class JEdgePartitionColor : public JEdgeColor {
public:

    string getName() const {
        return "EdgePartition";
    }

    bool isPerEdge() const {
        return 1;
    }

    int assign(const JEdgePtr &e);

private:
    map<int,JColor> colormap;
};

class JFacePartitionColor : public JFaceColor {
public:

    string getName() const {
        return "FacePartition";
    }

    bool isPerFace() const {
        return 1;
    }

    int assign(const JFacePtr &face);
    void clear() { colormap.clear(); }

private:
    map<int,JColor> colormap;
};
///////////////////////////////////////////////////////////////////////////////
class JCellPartitionColor : public JCellColor {
public:

    string getName() const {
        return "CellPartition";
    }

    bool isPerCell() const {
        return 1;
    }

    int assign( const JCellPtr &cell);

    int  operator() (JCellPtr c) {
        return assign(c);
    }

private:
    map<int,JColor> colormap;
};

