#include "MeshPartitioner.hpp"

//////////////////////////////////////////////////////////////////////////////////

int JMetisPartitioner :: getPartitions(int n)
{
    if( mesh == nullptr) return 1;

    clear();

    int topDim = mesh->getTopology()->getDimension();

    if( topDim == 2) getFacePartitions(n);
    if( topDim == 3) getCellPartitions(n);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

vector<JFaceSequence> JMetisPartitioner :: getPartitions( const JFaceSequence &faces, int nparts)
{
    vector<JFaceSequence> faceParts;

    int numTris   = 0;
    int numQuads = 0;
    int numPolys = 0;
    for( const JFacePtr &f : faces) {
        int nv = f->getSize(0);
        switch(nv)
        {
        case 3:
            numTris++;
            break;
        case 4:
            numQuads++;
            break;
        default:
            if( nv > 4) numPolys++;
            break;
        }
    }

    if( numTris + numQuads + numPolys == 0) return faceParts;
    if( numPolys > 0) {
        cout << "Warning: Partitions for polygons not supported yet" << endl;
        return faceParts;
    }

    vector<idx_t> eind, eptr, epart, npart;

    JNodeSequence nodes;
    JMeshTopology::getEntitySet(faces, nodes);

    idx_t nn = nodes.size();
    npart.resize(nn);
    int index = 0;
    for( const JNodePtr &vtx : nodes)
        vtx->setAttribute("LocalID", index++);

    int id;
    for( const JFacePtr &f :  faces) {
        int nv = f->getSize(0);
        for( int i = 0; i < nv; i++) {
            const JNodePtr &vtx = f->getNodeAt(i);
            vtx->getAttribute("LocalID", id);
            eind.push_back(id);
        }
    }

    idx_t ne = faces.size();
    epart.resize(ne);
    eptr.resize(ne+1);

    eptr[0]  = 0;
    index    = 1;
    int nsum = 0;
    for( idx_t i = 0; i < ne; i++)  {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() ) {
            nsum += f->getSize(0);
            eptr[index++] = nsum;
        }
    }

    idx_t ncommon = 2;
    idx_t objval;
    METIS_PartMeshDual(&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval,
                       &epart[0], &npart[0]);

    faceParts.resize( nparts);
    index = 0;
    for( const JFacePtr &f : faces) {
        int id = epart[index++];
        faceParts[id].push_back(f);
    }

    return faceParts;
}

///////////////////////////////////////////////////////////////////////////////
int JMetisPartitioner :: getFacePartitions(int nparts)
{
    if( mesh == nullptr) return 1;

    JFaceSequence faces = mesh->getFaces();

    vector<JFaceSequence> faceParts;
    faceParts = getPartitions(faces,nparts);

    sort( faceParts.begin(), faceParts.end(), []( const JFaceSequence &a, const JFaceSequence &b)
    {
        return a.size() < b.size();
    });

    int nParts = faceParts.size();
    for( int i = 0; i < nParts; i++) {
         for( const JFacePtr &face : faceParts[i])
              face->setAttribute("Partition", (int)i);
    }

    searchInterfaces();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMetisPartitioner :: getTopologicalDisks()
{
    if( mesh == nullptr) return 1;

    int euler = mesh->getTopology()->getEulerCharacteristic();
    if( euler == 1) return 0;

    int nParts = this->getNumPartitions();
    vector<JFaceSequence> faceParts, newParts;

    if( nParts == 0) {
        JFaceSequence faces = mesh->getFaces();
        faceParts = getPartitions(faces, 2);
    } else {
        faceParts.resize(nParts);
        for( int i = 0; i < nParts; i++)
            JMeshPartitioner::getPartition( i, faceParts[i] );
    }

    std::deque<JFaceSequence> partQ;
    for( size_t i = 0; i < faceParts.size(); i++)
        partQ.push_back(faceParts[i] );

    int currIndex = 0;
    while( !partQ.empty() ) {
        JFaceSequence currFaces = partQ.front();
        partQ.pop_front();
        int euler = JMeshTopology::getEulerCharacteristic(currFaces);
        if( euler == 1) {
            for( const JFacePtr &face : currFaces)
                face->setAttribute("Partition", currIndex);
            currIndex++;
        } else {
            cout << currIndex << " Split submesh with Euler Characteristic of " << euler << endl;
            newParts = getPartitions(currFaces, 2);
            cout << "CurrFace Size " << currFaces.size() << endl;
            cout << "#NewParts " << newParts.size() << endl;
            for( size_t i = 0; i < newParts.size(); i++) {
                partQ.push_back( newParts[i] );
                cout << newParts[i].size() << endl;
            }
        }
    }

    searchInterfaces();

    return  0;
}

////////////////////////////////////////////////////////////////////////////////

int JMetisPartitioner :: getCellPartitions(int n)
{
    if( mesh == nullptr) return 1;
    numParts = n;

    clear();

    mesh->enumerate(0);

    vector<size_t> eConnect;
    vector<idx_t> eind, eptr, epart, npart;

    idx_t nn = mesh->getSize(0);
    npart.resize(nn);

    int elmtype = mesh->getTopology()->getElementsType(3);

    vector<int> topo; // Not used ...
    mesh->getTopology()->getNodesArray( eConnect, topo);

    eind.resize(eConnect.size());

    boost::copy( eConnect, eind.begin() );

    idx_t ne = mesh->getSize(3);
    epart.resize(ne);
    eptr.resize(ne+1);

    int index;
    int etype = 0;
    switch(elmtype) {
    case JCell::TETRAHEDRON:
        etype = 2;
        break;
    case JCell::HEXAHEDRON:
        etype = 3;
        break;
    }
    eptr[0]  = 0;
    index    = 1;
    int nsum = 0;
    for( idx_t i = 0; i < ne; i++)  {
        const JCellPtr &c = mesh->getCellAt(i);
        if( c->isActive() ) {
            nsum += c->getSize(0);
            eptr[index++] = nsum;
        }
    }

    if( etype  == 0) return 1;

    idx_t ncommon = 2;
    idx_t objval;
    idx_t nparts = numParts;
    METIS_PartMeshDual(&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval,
                       &epart[0], &npart[0]);

    index = 0;
    for( int  i = 0; i < ne; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) cell->setAttribute("Partition", (int)epart[index++] );
    }

    index = 0;
    for( int i = 0; i < nn; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) vertex->setAttribute("Partition", (int)npart[index++] );
    }

    searchInterfaces();
    searchCorners();

    return 0;
}


