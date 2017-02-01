#pragma once

#ifndef JGRAPH_H
#define JGRAPH_H

template<class NodeType>
class Graph
{
public:
    void addNode( NodeType &v);
    void removeNode( NodeType &v);

    void addEdge(NodeType &v, NodeType &b);
    void removeEdge(NodeType &v, NodeType &b);

    size_t getSize(int e);
private:
    typedef vector<NodeType>  Row;
    vector<Row>   nodes;
};

#endif


