#include "MeshFacesClustering.hpp"

////////////////////////////////////////////////////////////////////////////////

std::pair<int,int> JMeshFacesClustering::max_frequency(const vector<int>  &items)
{
   map<int, unsigned>  freqCount;
   for( int v: items) {
        if( freqCount.find(v) == freqCount.end() ) freqCount[v] = 0;
        freqCount[v]++;
   }
   unsigned maxval = 0;
   for( auto it = freqCount.begin(); it != freqCount.end(); ++it)
        maxval = max( maxval, it->second );

   for( auto it = freqCount.begin(); it != freqCount.end(); ++it) {
        if( it->second == maxval) return std::make_pair(it->first, maxval);
   }
}
////////////////////////////////////////////////////////////////////////////////

void JMeshFacesClustering :: getPatchCenters()
{
   JMeshPartitioner mpart;
   mpart.setMesh(mesh);

   int nparts = mpart.getNumPartitions();
   for( int i = 0; i < nparts; i++) {
        JMeshPtr submesh = mpart.getSubMesh(i);
        JFacePtr fcenter = mesh->getTopology()->getCenter( submesh->getFaces() );
        if( fcenter ) seeds.push_back(fcenter);
   }
}

/////////////////////////////////////////////////////////////////////////////////
void JMeshFacesClustering :: expand()
{
    if( seeds.empty() ) getPatchCenters();

    if( seeds.empty() ) return;

    JFaceSequence faces = mesh->getFaces();

    int partid = -1;
    int levelid = std::numeric_limits<int>::max();
    for( const JFacePtr &f : faces) {
        f->setAttribute("Partition", partid);
        f->setAttribute("Level",     levelid);
    }

    partid  = 0;
    levelid = 0;
    deque<JFacePtr> faceQ;
    for( const JFacePtr &f : seeds) {
        f->setAttribute("Partition", partid++);
        f->setAttribute("Level",     levelid);
        faceQ.push_back(f);
    }

    int ipart, jpart;
    int ilevel, jlevel;

    JFaceSequence faceneighs;
    while( !faceQ.empty() ) {
        JFacePtr currface = faceQ.front();
        faceQ.pop_front();
        currface->getAttribute("Partition", ipart);
        currface->getAttribute("Level", ilevel);
        assert( ipart != -1);
        JFace::getRelations12(currface, faceneighs);
        for( const JFacePtr &nextface: faceneighs) {
            nextface->getAttribute("Level", jlevel);
            if( jlevel > ilevel+1) {
            nextface->setAttribute("Partition", ipart);
            nextface->setAttribute("Level", ilevel+1);
            faceQ.push_back(nextface);
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////
void JMeshFacesClustering :: removeZigZagInterfaces()
{
    JFaceSequence faces = mesh->getFaces();

    // The boundary may be zig-zag. The following may improve it.
    vector<int> facepartid;
    bool valueChanged = 0;
    JFaceSequence faceneighs;

    for( const JFacePtr &currface : faces) {
        JFace::getRelations12(currface, faceneighs);
        int numneighs = faceneighs.size();
        facepartid.resize( numneighs);
        for( int i = 0; i < numneighs; i++) 
             faceneighs[i]->getAttribute("Partition", facepartid[i] );
        pair<int,int> keyval = max_frequency( facepartid );
        if( keyval.second > numneighs/2) {
             currface->setAttribute("Partition", keyval.first);
             valueChanged = 1;
        }
     }

}
/////////////////////////////////////////////////////////////////////////////////
void JMeshFacesClustering :: getSpectralClusters( int numClusters)
{
   if( mesh == nullptr) return;

   cout << "Spectral " << endl;
 
   int numfaces = mesh->getSize(2);
   Eigen::MatrixXd  spmat(numfaces, numfaces);

   size_t numedges  = mesh->getSize(1);
   JFaceSequence faceneighs;
   for( size_t i = 0; i  < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        JEdge::getRelations(edge, faceneighs);
        if( faceneighs.size() == 2) {
            int i = faceneighs[0]->getID();
            int j = faceneighs[1]->getID();
            spmat.coeffRef(i,j) = 1;
            spmat.coeffRef(j,i) = 1;
        }
    }
   cout << "Spectral " << endl;
    SpectralClustering* c = new SpectralClustering(spmat, numClusters);
   cout << "Spectral " << endl;

    bool autotune = true;

    std::vector<std::vector<int> > clusters;
    if (autotune) {
        // auto-tuning clustering
        clusters = c->clusterRotate();
    } else {
        clusters = c->clusterKmeans(numClusters);
    }
   cout << "#Clusters " << clusters.size() << endl;

    // output clustered items
    // items are ordered according to distance from cluster centre
    for (unsigned int i=0; i < clusters.size(); i++) {
         for( auto j : clusters[i]) {
              const JFacePtr &f = mesh->getFaceAt(j);
              f->setAttribute("Partition", (int)i);
         }
    }
    exit(0);
}
/////////////////////////////////////////////////////////////////////////////////

