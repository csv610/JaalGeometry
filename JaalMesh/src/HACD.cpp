#include "HACD.hpp"
#include <iostream>

using namespace std;

JApproxConvexDecomposition :: JApproxConvexDecomposition()
{

    maxConcavity = 100;
    ccConnectDist = 30;
    addExtraDistPoints = 1;
    addFacesPoints = 1;

}
int JApproxConvexDecomposition :: getPartitions()
{
    int err = preprocess();
    if( err ) return 1;

    HACD::HeapManager * heapManager = HACD::createHeapManager(65536*(1000));

    HACD::HACD * const myHACD = HACD::CreateHACD(heapManager);
    myHACD->SetPoints(&points[0]);
    myHACD->SetNPoints(points.size());
    myHACD->SetTriangles(&triangles[0]);
    myHACD->SetNTriangles(triangles.size());
    myHACD->SetNClusters(nClusters);      // minimum number of clusters
    myHACD->SetCompacityWeight(0.0001);
    myHACD->SetVolumeWeight(0.0);
    myHACD->SetConnectDist(ccConnectDist);        // if two connected components are seperated by distance < ccConnectDist
    // then create a virtual edge between them so the can be merged during
    myHACD->SetNVerticesPerCH(100);       // max of 100 vertices per convex-hull
    myHACD->SetConcavity(maxConcavity);      // maximum concavity
    myHACD->SetSmallClusterThreshold(0.25); // threshold to detect small clusters
//  myHACD->SetNTargetTrianglesDecimatedMesh(targetNTrianglesDecimatedMesh); // # triangles in the decimated mesh
//  myHACD->SetCallBack(&CallBack);
    myHACD->SetAddExtraDistPoints(addExtraDistPoints);
    myHACD->SetAddFacesPoints(addFacesPoints);

    myHACD->SetNTargetTrianglesDecimatedMesh(1000); // # triangles in the decimated mesh

    clock_t start, end;
    start = clock();
    myHACD->Compute();
    end = clock();
    double elapsed = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Time " << elapsed << " s"<< std::endl;

    nClusters = myHACD->GetNClusters();

    JNodePtr v0,v1,v2;
    for(size_t c = 0; c < nClusters; ++c)
    {
        std::cout << std::endl << "Convex-Hull " << c << std::endl;
        size_t nPoints = myHACD->GetNPointsCH(c);
        size_t nTriangles = myHACD->GetNTrianglesCH(c);
        cout << "Debug " << nTriangles << endl;
        HACD::Vec3<HACD::Real> * pointsCH = new HACD::Vec3<HACD::Real>[nPoints];
        HACD::Vec3<long> * trianglesCH = new HACD::Vec3<long>[nTriangles];
        myHACD->GetCH(c, pointsCH, trianglesCH);
        /*
                std::cout << "Points " << nPoints << std::endl;
                for(size_t v = 0; v < nPoints; ++v)
                {
                    std::cout << v << "\t"
                              << pointsCH[v].X() << "\t"
                              << pointsCH[v].Y() << "\t"
                              << pointsCH[v].Z() << std::endl;
                }
                std::cout << "Triangles " << nTriangles << std::endl;
        */
        for(size_t f = 0; f < nTriangles; ++f)
        {
            int i0 = trianglesCH[f].X();
            int i1 = trianglesCH[f].Y();
            int i2 = trianglesCH[f].Z();
            v0  = mesh->getNodeAt(i0);
            v1  = mesh->getNodeAt(i1);
            v2  = mesh->getNodeAt(i2);
            JFacePtr face = Simplex::getFaceOf(v0,v1,v2);
            if( face ) {
                face->setAttribute("Partition", (int)c);
            }
        }
        delete [] pointsCH;
        delete [] trianglesCH;
    }
//    std::cout << "warninhg: memory not released " << std::endl;
//    HACD::DestroyHACD(myHACD);
//    HACD::releaseHeapManager(heapManager);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JApproxConvexDecomposition :: preprocess()
{
    if( mesh == nullptr) return 1;
    mesh->pruneAll();
    mesh->enumerate(0);
    mesh->enumerate(2);

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    points.resize(numnodes);
    for (long i = 0; i < numnodes ; i++)
    {
        const Point3D &xyz = mesh->getNodeAt(i)->getXYZCoords();
        points[i].X() = xyz[0];
        points[i].Y() = xyz[1];
        points[i].Z() = xyz[2];
    }

    triangles.resize(numfaces);
    for (long i = 0; i < numfaces ; ++i) {
        const JFacePtr &f = mesh->getFaceAt(i);
        triangles[i].X() = f->getNodeAt(0)->getID();
        triangles[i].Y() = f->getNodeAt(1)->getID();
        triangles[i].Z() = f->getNodeAt(2)->getID();
    }
    return 0;
}

