#include "MixedIntegerQuads.hpp"

#ifdef USE_IGL
using namespace igl::copyleft;
#endif

void JMixedIntegerQuads :: line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                                        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                                        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
{
    // Create a texture that hides the integer translation in the parametrization
    unsigned size = 128;
    unsigned size2 = size/2;
    unsigned lineWidth = 3;
    texture_R.setConstant(size, size, 255);
    for (unsigned i=0; i<size; ++i)
        for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
            texture_R(i,j) = 0;
    for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
        for (unsigned j=0; j<size; ++j)
            texture_R(i,j) = 0;

    texture_G = texture_R;
    texture_B = texture_R;
}

///////////////////////////////////////////////////////////////////////////////

void JMixedIntegerQuads :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

/*

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);

    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    // Compute face barycenters
    igl::barycenter(V, F, B);

    // Compute scale for visualizing fields
    global_scale =  .5*igl::avg_edge_length(V, F);

//  MatrixXd B1,B2,B3;
    igl::local_basis(V,F,B1,B2,B3);
*/
}

///////////////////////////////////////////////////////////////////////////////

void JMixedIntegerQuads::genRandomConstraints( int nrand)
{
    if( mesh == nullptr) return;

    b.resize(nrand);
    bc.resize(nrand,3);

    int numFaces = mesh->getSize(2);

    for( int i = 0; i < nrand; i++) {
        b[i]    =   JMath::random_value(0, numFaces-1);
        bc.coeffRef(i,0) =   drand48();
        bc.coeffRef(i,1) =   drand48();
        bc.coeffRef(i,2) =   drand48();
    }
}

///////////////////////////////////////////////////////////////////////////////

std::pair<JMeshPtr,JMeshPtr> JMixedIntegerQuads :: getCrossField()
{
    JMeshPtr minK = JMesh::newObject();
    addEdges( minK, B -global_scale*X1, B + global_scale*X1);

    JMeshPtr maxK = JMesh::newObject();
    addEdges( maxK, B- global_scale*X2, B + global_scale*X2);

    std::pair<JMeshPtr, JMeshPtr> result;
    result.first  = minK;
    result.second = maxK;
    return result;
}

///////////////////////////////////////////////////////////////////////////////

std::pair<JMeshPtr,JMeshPtr> JMixedIntegerQuads :: getBisectorField()
{
    JMeshPtr minK = JMesh::newObject();
    addEdges( minK, B - global_scale*BIS1, B + global_scale*BIS1);

    JMeshPtr maxK = JMesh::newObject();
    addEdges( maxK, B + global_scale*BIS2, B + global_scale*BIS2);

    std::pair<JMeshPtr, JMeshPtr> result;
    result.first  = minK;
    result.second = maxK;
    return result;
}
///////////////////////////////////////////////////////////////////////////////

std::pair<JMeshPtr,JMeshPtr> JMixedIntegerQuads :: getBisectorCombinedField()
{
    JMeshPtr minK = JMesh::newObject();
    addEdges( minK, B, global_scale*BIS1_combed);

    JMeshPtr maxK = JMesh::newObject();
    addEdges( maxK, B, global_scale*BIS2_combed);

    std::pair<JMeshPtr, JMeshPtr> result;
    result.first  = minK;
    result.second = maxK;
    return result;
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMixedIntegerQuads :: getSeams()
{
    int l_count = Seams.sum();
    Eigen::MatrixXd P1(l_count,3);
    Eigen::MatrixXd P2(l_count,3);

    for (unsigned i=0; i<Seams.rows(); ++i)
    {
        for (unsigned j=0; j<Seams.cols(); ++j)
        {
            if (Seams(i,j) != 0)
            {
                P1.row(l_count-1) = V.row(F(i,j));
                P2.row(l_count-1) = V.row(F(i,(j+1)%3));
                l_count--;
            }
        }
    }
    JMeshPtr seam = JMesh::newObject();
    addEdges( seam, P1,P2);
    return seam;
}


///////////////////////////////////////////////////////////////////////////////

void JMixedIntegerQuads :: getVecField( int  key)
{

    if (key == 4)
    {
        /*
                // Plot the singularities as colored dots (red for negative, blue for positive)
                for (unsigned i=0; i<singularityIndex.size(); ++i)
                {
                    if (singularityIndex(i) < 2 && singularityIndex(i) > 0)
                        viewer.data.add_points(V.row(i),Eigen::RowVector3d(1,0,0));
                    else if (singularityIndex(i) > 2)
                        viewer.data.add_points(V.row(i),Eigen::RowVector3d(0,1,0));
                }
        */

    }

    if (key == 5)
    {
        // Singularities and cuts, original field
        // Singularities and cuts
        /*

                // Plot the singularities as colored dots (red for negative, blue for positive)
                for (unsigned i=0; i<singularityIndex.size(); ++i)
                {
                    if (singularityIndex(i) < 2 && singularityIndex(i) > 0)
                        viewer.data.add_points(V.row(i),Eigen::RowVector3d(1,0,0));
                    else if (singularityIndex(i) > 2)
                        viewer.data.add_points(V.row(i),Eigen::RowVector3d(0,1,0));
                }
        */
    }

    /*
        if (key == '6')
        {
            // Global parametrization UV
            viewer.data.set_mesh(UV, FUV);
            viewer.data.set_uv(UV);
            viewer.core.show_lines = true;
        }

        if (key == '7')
        {
            // Global parametrization in 3D
            viewer.data.set_mesh(V, F);
            viewer.data.set_uv(UV,FUV);
            viewer.core.show_texture = true;
        }

        if (key == '8')
        {
            // Global parametrization in 3D with seams
            viewer.data.set_mesh(V, F);
            viewer.data.set_uv(UV_seams,FUV_seams);
            viewer.core.show_texture = true;
        }

        viewer.data.set_colors(Eigen::RowVector3d(1,1,1));

        // Replace the standard texture with an integer shift invariant texture
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;
        line_texture(texture_R, texture_G, texture_B);
        viewer.data.set_texture(texture_R, texture_B, texture_G);

        viewer.core.align_camera_center(viewer.data.V,viewer.data.F);
        return false;
    */
}

//////////////////////////////////////////////////////////////////////////////////////////
JMeshPtr JMixedIntegerQuads :: extractQuadMesh()
{
    /*

        qex_TriMesh triMesh;
        qex_QuadMesh quadMesh;
        qex_Point3* vertex;

        int numnodes = mesh->getSize(0);
        int numTris  = mesh->getSize(2);

        triMesh.vertex_count = numnodes;
        triMesh.tri_count    = numTris;

        triMesh.vertices = (qex_Point3*)malloc(sizeof(qex_Point3) * numnodes);
        triMesh.tris = (qex_Tri*)malloc(sizeof(qex_Tri) * numTris);
        triMesh.uvTris = (qex_UVTri*)malloc(sizeof(qex_UVTri) * numTris);

        qex_Point3  qp;
        for( int i = 0; i < numnodes; i++) {
             const JNodePtr &v  = mesh->getNodeAt(i);
             const Point3D  &p  = v->getXYZCoords();
             qp.x[0]  = p[0];
             qp.x[1]  = p[1];
             qp.x[2]  = p[2];
             triMesh.vertices[i] = qp;
        }

        qex_Tri conn;
        for( int i = 0; i < numnodes; i++) {
             const JFacePtr &f  = mesh->getFaceAt(i);
             conn.indices[0] = f->getNodeAt(0)->getID();
             conn.indices[1] = f->getNodeAt(1)->getID();
             conn.indices[2] = f->getNodeAt(3)->getID();
            triMesh.tris[i]  = conn;
        }


        qex_UVTri  uvt;
        qex_Point2 uv;
        for( int i = 0; i < numTris; i++) {
             for( int j = 0; j < 3; j++) {
                  int vid  = FUV.coeff(i,j);
                  uv.x[0] = UV.coeff(vid,0);
                  uv.x[1] = UV.coeff(vid,1);
                  uvt.uvs[j] = uv;
              }
              triMesh.uvTris[i] = uvt;
         }
        qex_extractQuadMesh(&triMesh, NULL, &quadMesh);
    */
    JMeshOFFExporter mexp;
    mexp.writeFile(mesh, "xyz.off");

    JMeshEigenMatrix mat;
    JMeshPtr uvmesh = mat.getMesh(UV, FUV);

    mexp.writeFile(uvmesh, "uv.off");

    int err = system("quadextract  xyz.off uv.off uvquads.off");
    if( err < 0) return nullptr;

    JMeshOFFImporter mimp;
    JMeshPtr uvQuads = mimp.readFile("uvquads.off");
    return uvQuads;

}

JMeshPtr JMixedIntegerQuads :: getUVMesh()
{
#ifdef CSV
    if( b.size() < 1) genRandomConstraints(1);

    // Create a smooth 4-RoSy field
    VectorXd S;
    comiso::nrosy(V,F,b,bc,VectorXi(),VectorXd(),MatrixXd(),4,0.5,X1,S);

    // Find the the orthogonal vector
    X2 = igl::rotate_vectors(X1, VectorXd::Constant(1,M_PI/2), B1, B2);

    // Always work on the bisectors, it is more general
    igl::compute_frame_field_bisectors(V, F, X1, X2, BIS1, BIS2);

    // Comb the field, implicitly defining the seams
    igl::comb_cross_field(V, F, BIS1, BIS2, BIS1_combed, BIS2_combed);

    // Find the integer mismatches
    igl::cross_field_missmatch(V, F, BIS1_combed, BIS2_combed, true, MMatch);

    // Find the singularities
    igl::find_cross_field_singularities(V, F, MMatch, isSingularity, singularityIndex);

    // Cut the mesh, duplicating all vertices on the seams
    igl::cut_mesh_from_singularities(V, F, MMatch, Seams);

    // Comb the frame-field accordingly
    igl::comb_frame_field(V, F, X1, X2, BIS1_combed, BIS2_combed, X1_combed, X2_combed);

    double  gradient_size = 50;
    double  iter = 0;
    double  stiffness = 5.0;
    bool    direct_round = 0;
    // Global parametrization
    comiso::miq(V, F, X1_combed, X2_combed, MMatch, isSingularity,
                Seams, UV, FUV, gradient_size, stiffness, direct_round,
                iter, 5, true);

    JMeshPtr uvQuads = extractQuadMesh();
    return uvQuads;
#endif

    /*
    // Global parametrization (with seams, only for demonstration)
        igl::comiso::miq(V,
                         F,
                         X1_combed,
                         X2_combed,
                         MMatch,
                         isSingularity,
                         Seams,
                         UV_seams,
                         FUV_seams,
                         gradient_size,
                         stiffness,
                         direct_round,
                         iter,
                         5,
                         false);
    */
}
