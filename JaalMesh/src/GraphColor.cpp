#include "GraphColor.hpp"

void JGraphColor :: setColor( const JNodePtr &vertex)
{
    if( !vertex->isActive() ) return;

    JNode::getRelations(vertex, vneighs);

    int nSize = vneighs.size();
    int cid, err;

    usedColors.clear();

    err = vertex->getAttribute(name, cid);
    if( !err ) usedColors.push_back(cid);
    for( int i = 0; i < nSize; i++) {
        err = vneighs[i]->getAttribute(name, cid);
        if( !err ) usedColors.push_back(cid);
    }

    unusedColors.clear();
    for( int i = 0; i < nSize; i++) {
        if( find(usedColors.begin(), usedColors.end(), i) == usedColors.end() )
            unusedColors.push_back(i);
    }

    int index = 0;
    for( int i = 0; i < nSize; i++) {
        err = vneighs[i]->getAttribute(name, cid);
        if( err ) vneighs[i]->setAttribute(name, unusedColors[index++] );
    }

    err = vertex->getAttribute(name, cid);
    if( err ) vertex->setAttribute(name, nSize+1);
}

//////////////////////////////////////////////////////////////////////////////////////////
int JGraphColor:: verify()
{
    size_t numnodes = mesh->getSize(0);
    maxIndex = 0;

    int clr1, clr2;

    for( size_t i = 0; i <  numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            vertex->getAttribute(name, clr1);
            maxIndex = max( clr1, maxIndex);
            JNode::getRelations(vertex, vneighs);
            int nSize = vneighs.size();
            for( int j = 0; j < nSize; j++) {
                vneighs[j]->getAttribute(name, clr2);
                assert( clr1 != clr2);
                return 1;
                maxIndex = max( clr2, maxIndex);
            }
        }
    }
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////

int JGraphColor:: setColors( const JMeshPtr &m)
{
    mesh = m;
    maxIndex = 0;
    if( mesh == NULL ) return 1;

    name = "ColorIndex";

    mesh->buildRelations(0,0);

    mesh->deleteNodeAttribute(name);

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i <  numnodes; i++)
        setColor( mesh->getNodeAt(i));

    return verify();
}
//////////////////////////////////////////////////////////////////////////////////////////

/*
#include <iostream>
#include <list>
using namespace std;

// A class that represents an undirected graph
class Graph
{
    int V;    // No. of vertices
    list<int> *adj;    // A dynamic array of adjacency lists
public:
    // Constructor and destructor
    Graph(int V)   { this->V = V; adj = new list<int>[V]; }
    ~Graph()       { delete [] adj; }

    // function to add an edge to graph
    void addEdge(int v, int w);

    // Prints greedy coloring of the vertices
    void greedyColoring();
};

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);  // Note: the graph is undirected
}

// Assigns colors (starting from 0) to all vertices and prints
// the assignment of colors
void Graph::greedyColoring()
{
    int result[V];

    // Assign the first color to first vertex
    result[0]  = 0;

    // Initialize remaining V-1 vertices as unassigned
    for (int u = 1; u < V; u++)
        result[u] = -1;  // no color is assigned to u

    // A temporary array to store the available colors. True
    // value of available[cr] would mean that the color cr is
    // assigned to one of its adjacent vertices
    bool available[V];
    for (int cr = 0; cr < V; cr++)
        available[cr] = false;

    // Assign colors to remaining V-1 vertices
    for (int u = 1; u < V; u++)
    {
        // Process all adjacent vertices and flag their colors
        // as unavailable
        list<int>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (result[*i] != -1)
                available[result[*i]] = true;

        // Find the first available color
        int cr;
        for (cr = 0; cr < V; cr++)
            if (available[cr] == false)
                break;

        result[u] = cr; // Assign the found color

        // Reset the values back to false for the next iteration
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (result[*i] != -1)
                available[result[*i]] = false;
    }

    // print the result
    for (int u = 0; u < V; u++)
        cout << "Vertex " << u << " --->  Color "
             <<  result[u] << endl;
*/
