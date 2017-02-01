#include "MeshSamples.hpp"

////////////////////////////////////////////////////////////////////////////////////

JNodePtr JMeshSamples :: getRandomNode( bool fromWhere) const
{
    if( mesh == nullptr) return nullptr;
    size_t minindex = 0;
    size_t numnodes = mesh->getSize(0);
    size_t index    = JMath::random_value(minindex,numnodes-1);
    JNodePtr  nullPtr;

    if( fromWhere == 0) {
        for ( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt((index+i)%numnodes);
            if( v->isActive() && !v->isBoundary() ) return v;
        }
    } else {
        JNodeSequence vb;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt( (index+i)%numnodes);
            if( vertex->isActive() && vertex->isBoundary() ) {
                vb.push_back(vertex);
            }
        }
        if( vb.empty() ) return nullPtr;
        index = JMath::random_value(minindex, vb.size()-1);
        return vb[index];
    }

    return nullPtr;
}
///////////////////////////////////////////////////////////////////////////////

JEdgePtr JMeshSamples :: getRandomEdge( bool fromWhere ) const
{
    if( mesh == nullptr) return nullptr;
    size_t minindex = 0;
    size_t numedges = mesh->getSize(1);
    size_t index = JMath::random_value(minindex,numedges-1);
    JEdgePtr edge, nullPtr;

    // Return any active edge from the set.
    if( fromWhere == 0) {
        for( size_t i = 0; i < numedges; i++) {
            edge = mesh->getEdgeAt( (index+i)%numedges);
            if( edge->isActive() ) return edge;
        }
    } else {
        // Return any active edge which is on the boundary.
        JEdgeSequence eb;
        for( size_t i = 0; i < numedges; i++) {
            edge = mesh->getEdgeAt( (index+i)%numedges);
            if( edge->isActive() && edge->isBoundary() ) {
                eb.push_back(edge);
            }
        }
        if( eb.empty() ) return nullPtr;
        index = JMath::random_value(minindex, eb.size()-1);
        return eb[index];
    }

    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JMeshSamples :: getRandomFace( bool fromWhere) const
{
    if( mesh == nullptr) return nullptr;
    size_t minindex = 0;
    size_t numfaces = mesh->getSize(2);
    size_t index = JMath::random_value(minindex,numfaces-1);

    JFacePtr face, nullPtr;
    if( fromWhere == 0) {
        for( size_t i = 0; i < numfaces; i++) {
            face = mesh->getFaceAt( (index+i)%numfaces);
            if( face->isActive() ) return face;
        }
    } else {
        vector<size_t> boundIndex;
        for( size_t i = 0; i < numfaces; i++) {
            face = mesh->getFaceAt(i);
            if( face->isActive() && face->isBoundary() )
                boundIndex.push_back(i);
        }
        if( boundIndex.empty() ) {
            return nullptr;
        } else {
            index = JMath::random_value(minindex, boundIndex.size()-1);
            face  = mesh->getFaceAt( boundIndex[index] ) ;
            assert( face->isBoundary() );
            return face;
        }
    }

    return nullPtr;
}
////////////////////////////////////////////////////////////////////////////////////

JCellPtr JMeshSamples :: getRandomCell() const
{
    if( mesh == nullptr) return nullptr;
    size_t minindex = 0;
    size_t maxindex = mesh->getSize(3);
    size_t index    = JMath::random_value(minindex,maxindex-1);
    return mesh->getCellAt(index);
}

////////////////////////////////////////////////////////////////////////////////////

JNodePtr JMeshSamples :: getRandomNode(const JMeshPtr &submesh)
{
    int numnodes = submesh->getSize(0);
    int nr = JMath::random_value( 0, numnodes-1);
    return submesh->getNodeAt(nr);
}

////////////////////////////////////////////////////////////////////////////////////

void JMeshSamples :: getMetisSamples( int n, JNodeSequence &sampled)
{
    JMetisPartitioner mp;
    mp.setMesh(mesh);
    mp.getPartitions(n);

    int np = mp.getNumPartitions();
    for( int i = 0; i < np; i++) {
        JMeshPtr submesh = mp.getSubMesh(i);
        const JNodePtr &vtx = getRandomNode(submesh );
        sampled.push_back(vtx);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void JMeshSamples :: getRandomSamples( int n, JNodeSequence &sampled)
{
    sampled.clear();

    set<int> aset;

    srand( time(0));
    size_t numnodes = mesh->getSize(0);
    for( int i = 0; i < n; i++) {
        int id = JMath::random_value((size_t)0, numnodes);
        if(  aset.find(id) == aset.end() )  {
            const JNodePtr &vtx = mesh->getNodeAt(id);
            if( vtx->isActive() ) {
                sampled.push_back(vtx);
                aset.insert(id);
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshSamples :: getNodeSamples(int ncount)
{
    JNodeSequence samples;

    switch(algorithm)
    {
    case RANDOM:
        getRandomSamples(ncount, samples);
        break;
    case METIS:
        getMetisSamples(ncount, samples);
        break;
    case GEODESIC:
        getGeodesicSamples(ncount, samples);
        break;
    case BIHARMONIC:
        getIGLSamples(ncount, samples);
        break;
    case REGULAR:
        getRegularSamples(ncount, samples);
        break;
    }

    return samples;
}

////////////////////////////////////////////////////////////////////////////////////

void JMeshSamples :: getGeodesicSamples( int numSamples, JNodeSequence &sampleNodes)
{
    sampleNodes.clear();

    if( mesh == nullptr) return;

    std::vector<double>    points;

    size_t index = 0;
    size_t numnodes = mesh->getSize(0);
    points.reserve(3*numnodes);
    for( size_t i = 0; i <  numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->setID(index);
            const Point3D &xyz = vtx->getXYZCoords();
            points.push_back( xyz[0] );
            points.push_back( xyz[1] );
            points.push_back( xyz[2] );
            index++;
        }
    }

    std::vector<unsigned>  faces;
    size_t numfaces = mesh->getSize(2);
    faces.reserve(3*numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            faces.push_back( face->getNodeAt(0)->getID() );
            faces.push_back( face->getNodeAt(1)->getID() );
            faces.push_back( face->getNodeAt(2)->getID() );
        }
    }

    geodesic::Mesh gmesh;
    gmesh.initialize_mesh_data(points, faces);

// geodesic::GeodesicAlgorithmExact algorithm(&gmesh);
    geodesic::GeodesicAlgorithmDijkstra algorithm(&gmesh);
    std::vector<geodesic::SurfacePoint> all_sources;

    std::vector<unsigned>  sampled;
    int new_index  = 1;
    double distance, maxdist = 0;
    for( int j = 0; j < numSamples; j++) {
        sampled.push_back(new_index);
        geodesic::SurfacePoint newsource(&gmesh.vertices()[new_index]);		//create source
        all_sources.push_back(newsource);
        algorithm.propagate(all_sources);	//cover the whole mesh
        maxdist = 0.0;
        for(unsigned i=0; i< gmesh.vertices().size(); ++i)
        {
            geodesic::SurfacePoint p(&gmesh.vertices()[i]);
            unsigned best_source = algorithm.best_source(p,distance);
            if( distance > maxdist) {
                maxdist = distance;
                new_index = i;
            }
        }
    }

    sampleNodes.resize(numSamples);
    for( int j = 0; j < numSamples; j++)  {
        sampleNodes[j] = mesh->getNodeAt( sampled[j] );
        cout << "Geode Sample " << j << "  ID " << sampled[j] << endl;
    }

    /*
        JTriMeshGeodesics  geodes;
        geodes.setAlgorithm( JTriMeshGeodesics::EXACT_DIJKSTRA);
        geodes.setMesh(mesh);

        sampled.clear();

        size_t numNodes =  mesh->getSize(0);
        for( size_t i = 0; i < mesh->getSize(0); i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) {
                sampled.push_back(vtx);
                break;
            }
        }

        vector<double> dist;
        dist = geodes.getDistanceField(sampled);

        for( int i = 0; i < n-1; i++) {
            double maxdist = *boost::max_element(dist);
            for( size_t j = 0; j < numNodes; j++) {
                if( dist[j] == maxdist) {
                    sampled.push_back(mesh->getNodeAt(j));
                    break;
                }
            }
            dist = geodes.getDistanceField(sampled);
        }
    */
}

////////////////////////////////////////////////////////////////////////////////////

void JMeshSamples :: getIGLSamples( int n, JNodeSequence &sampled)
{
#ifdef USE_IGL
    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    Engine *engine;
    igl::matlab::mlinit(&engine);

    // Send Laplacian matrix to matlab
    igl::matlab::mlsetmatrix(&engine,"V",V);
    igl::matlab::mlsetmatrix(&engine,"F",F);
    igl::matlab::mlsetscalar(&engine,"n",n);

    ostringstream oss;

    switch( algorithm )
    {
    case GEODESIC:
        oss << "[~,ids] = farthest_points(V,n,F, 'distance', 'geodesic')";
        break;
    case BIHARMONIC:
        oss << "[~,ids] = farthest_points(V,n,F, 'distance', 'biharmonic')";
        break;
    }

    string cmd = oss.str();
    if( !cmd.empty() ) {
        igl::matlab::mleval(&engine, cmd);
        Eigen::MatrixXd pid;
        igl::matlab::mlgetmatrix(&engine,"ids",pid);
        for( int i = 0; i < n; i++)
            sampled[i] = mesh->getNodeAt(pid.coeff(i,0));
    }

    igl::matlab::mlclose(&engine);
#endif
}
////////////////////////////////////////////////////////////////////////////////////

void JMeshSamples :: getRegularSamples( int n, JNodeSequence &sampled)
{

}
