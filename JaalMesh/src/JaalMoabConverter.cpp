#include "JaalMoabConverter.hpp"

#ifdef HAVE_MOAB
iBase_EntityHandle
JaalMoabConverter :: new_MOAB_Handle(iMesh_Instance imesh, Jaal::Vertex *vertex)
{
    int err;
    iBase_EntityHandle newHandle;

    std::map<NodeHandle, iBase_EntityHandle> ::iterator it;
    it = moabnode.find(vertex);

    if( it != moabnode.end() ) {
        newHandle = it->second;
        iMesh_deleteEnt(imesh, newHandle, &err);
        assert( !err );
    }

    Point3D p = vertex->getXYZCoords();
    iMesh_createVtx(imesh, p[0], p[1], p[2], &newHandle, &err);
    assert( !err );

    moabnode[vertex]    = newHandle;
    jaalnode[newHandle] = vertex;
    return newHandle;
}


///////////////////////////////////////////////////////////////////////////////

iBase_EntityHandle
JaalMoabConverter ::new_MOAB_Handle(iMesh_Instance imesh, const JFacePtr &face)
{
    int status, err;

    iBase_EntityHandle newHandle;

    std::map<FaceHandle, iBase_EntityHandle> ::iterator fit;
    fit = moabface.find(face);

    if( fit != moabface.end() ) {
        newHandle = fit->second;
        iMesh_deleteEnt(imesh, newHandle, &err);
        assert(!err);
    }

    vector<iBase_EntityHandle> connect;

    int nnodes = face->getSize(0);
    connect.resize(nnodes);

    std::map<NodeHandle, iBase_EntityHandle> ::iterator nit;
    for (int j = 0; j < nnodes; j++) {
        Vertex *v  = face->getNodeAt(j);
        nit = moabnode.find( v );
        assert( nit != moabnode.end() );
        connect[j] = nit->second;
    }

    switch (nnodes) {
    case 3:
        iMesh_createEnt(imesh, iMesh_TRIANGLE, &connect[0], nnodes, &newHandle,
                        &status, &err);
        assert(!err);
        break;
    case 4:
        iMesh_createEnt(imesh, iMesh_QUADRILATERAL, &connect[0], nnodes,
                        &newHandle, &status, &err);
        assert(!err);
        break;
    default:
        iMesh_createEnt(imesh, iMesh_POLYGON, &connect[0], nnodes, &newHandle,
                        &status, &err);
        assert(!err);
        break;
    }

    moabface[face]  = newHandle;
    jaalface[newHandle] = face;

    return newHandle;
}

///////////////////////////////////////////////////////////////////////////////

int
JaalMoabConverter::toMOAB(const JMeshPtr &jmesh, iMesh_Instance &imesh, iBase_EntitySetHandle entitySet)
{
    int err;

    assert( jmesh->isPruned() );

    if (imesh == 0) {
        int optlen = 0;
        char *options = NULL;
        iMesh_newMesh(options, &imesh, &err, optlen);
        assert(!err);
    }

    moabnode.clear();
    jaalnode.clear();

    moabface.clear();
    jaalface.clear();

    iBase_EntityHandle newHandle;
    std::map<NodeHandle, iBase_EntityHandle> ::const_iterator niter;

    size_t numnodes = jmesh->getSize(0);

    size_t ncount0 = 0, ncount1 = 0;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = jmesh->getNodeAt(i);
        if( v->isActive() ) {
            niter = moabnode.find(v);
            if( niter == moabnode.end() ) {
                ncount0++;
                newHandle = new_MOAB_Handle(imesh, v);
                if (entitySet) {
                    iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
                    assert( !err );
                }
            } else {
                ncount1++;
                newHandle = niter->second;
                const Point3D &p3d = v->getXYZCoords();
                iMesh_setVtxCoord(imesh, newHandle, p3d[0], p3d[1], p3d[2], &err);
                assert( !err );
            }
        }
    }

    size_t numfaces = jmesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = jmesh->getFaceAt(i);
        if( f->isActive() ) {
            newHandle = new_MOAB_Handle(imesh, f);
            if (entitySet) {
                iMesh_addEntToSet(imesh, newHandle, entitySet, &err);
                assert( !err );
            }
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

JMeshPtr
JaalMoabConverter ::fromMOAB(iMesh_Instance imesh, const JMeshPtr &m, iBase_EntitySetHandle entitySet)
{

    int err;
    if (entitySet == 0)
        iMesh_getRootSet(imesh, &entitySet, &err);

    SimpleArray<iBase_EntityHandle> triHandles, quadHandles, tetHandles, hexHandles;

    iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_TRIANGLE, ARRAY_INOUT(triHandles), &err);
    iMesh_getEntities(imesh, entitySet, iBase_FACE, iMesh_QUADRILATERAL, ARRAY_INOUT(quadHandles), &err);
    iMesh_getEntities(imesh, entitySet, iBase_REGION, iMesh_TETRAHEDRON, ARRAY_INOUT(tetHandles), &err);
    iMesh_getEntities(imesh, entitySet, iBase_REGION, iMesh_HEXAHEDRON , ARRAY_INOUT(hexHandles), &err);

    size_t numTris= triHandles.size();
    size_t numQuads= quadHandles.size();
    size_t numTets = tetHandles.size();
    size_t numHexs = hexHandles.size();

    jmesh = m;
    if( jmesh == NULL  &&  numTets + numHexs > 0 )
        jmesh = JMesh::newObject();

    if( jmesh == NULL  &&  numTris + numQuads > 0 )
        jmesh = JMeshPtr::newObject();


    SimpleArray<iBase_EntityHandle> facenodes, cellnodes;
    jaalnode.clear();

    for( size_t i = 0; i < numTris; i++) {
        iMesh_getEntAdj(imesh, triHandles[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
        for (int j = 0; j < 3; j++)
            jaalnode[ facenodes[j] ] = NULL;
    }

    for( size_t i = 0; i < numQuads; i++) {
        iMesh_getEntAdj(imesh, quadHandles[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
        for (int j = 0; j < 4; j++)
            jaalnode[ facenodes[j] ] = NULL;
    }

    for( size_t i = 0; i < numTets; i++) {
        iMesh_getEntAdj(imesh, tetHandles[i], iBase_VERTEX, ARRAY_INOUT(cellnodes), &err);
        for (int j = 0; j < 4; j++)
            jaalnode[ cellnodes[j] ] = NULL;
    }

    for( size_t i = 0; i < numHexs; i++) {
        iMesh_getEntAdj(imesh, hexHandles[i], iBase_VERTEX, ARRAY_INOUT(cellnodes), &err);
        for (int j = 0; j < 8; j++)
            jaalnode[ cellnodes[j] ] = NULL;
    }

    size_t numNodes = jaalnode.size();
    jmesh->reserve(numNodes, 0);

    Point3D p3d;
    double x, y, z;

    std::map<iBase_EntityHandle, NodeHandle> ::const_iterator niter;
    Vertex *jnode;

    size_t id = 0;
    for( niter = jaalnode.begin(); niter != jaalnode.end(); ++niter) {
        iBase_EntityHandle nhandle = niter->first;
        iMesh_getVtxCoord(imesh, nhandle, &x, &y, &z, &err);
        p3d[0] = x;
        p3d[1] = y;
        p3d[2] = z;
        jnode  = Jaal::Vertex::newObject();
        jnode->setXYZCoords(p3d);
        jnode->setID(id++);
        jmesh->addObject(jnode);
        moabnode[jnode] = nhandle;
        jaalnode[ nhandle ] = jnode;
    }


    jmesh->reserve(numTris+numQuads+numTets + numHexs, 2);
    JNodeSequence connect(3);
    for (size_t i = 0; i < numTris; i++) {
        SimpleArray<iBase_EntityHandle> facenodes;
        iMesh_getEntAdj(imesh, triHandles[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
        for (int j = 0; j < 3; j++)
            connect[j] = jaalnode[facenodes[j]];
        Face *face = Triangle::newObject( connect );;
        jmesh->addObject(face);
        jaalface[triHandles[i] ] = face;
        moabface[face] = triHandles[i];
    }

    connect.resize(4);
    for (size_t i = 0; i < numQuads; i++) {
        SimpleArray<iBase_EntityHandle> facenodes;
        iMesh_getEntAdj(imesh, quadHandles[i], iBase_VERTEX, ARRAY_INOUT(facenodes), &err);
        for (int j = 0; j < 4; j++)
            connect[j] = jaalnode[facenodes[j]];
        Face *face = Quadrilateral::newObject( connect );;
        jmesh->addObject(face);
        jaalface[quadHandles[i] ] = face;
        moabface[face] = quadHandles[i];
    }

    connect.resize(4);
    cout << "#Tets " << numTets << endl;
    for (size_t i = 0; i < numTets; i++) {
        SimpleArray<iBase_EntityHandle> cellnodes;
        iMesh_getEntAdj(imesh, tetHandles[i], iBase_VERTEX, ARRAY_INOUT(cellnodes), &err);
        for (int j = 0; j < 4; j++)
            connect[j] = jaalnode[cellnodes[j]];
        Cell *cell = Tetrahedron::newObject();;
        cell->setNodes(connect);
        jmesh->addObject(cell);
        jaalcell[tetHandles[i] ] = cell;
        moabcell[cell] = tetHandles[i];
    }

    connect.resize(8);
    for (size_t i = 0; i < numHexs; i++) {
        SimpleArray<iBase_EntityHandle> cellnodes;
        iMesh_getEntAdj(imesh, hexHandles[i], iBase_VERTEX, ARRAY_INOUT(cellnodes), &err);
        for (int j = 0; j < 8; j++)
            connect[j] = jaalnode[cellnodes[j]];
        Cell *cell = Hexahedron::newObject();;
        cell->setNodes(connect);
        jmesh->addObject(cell);
        jaalcell[hexHandles[i] ] = cell;
        moabcell[cell] = hexHandles[i];
    }

    return jmesh;
}
#endif
