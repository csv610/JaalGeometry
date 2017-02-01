#include "IntegrableField.hpp"

void JIntegrableField :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

#ifdef USE_IGL
    JMeshEigenMatrix mat;
    mat.setMesh(mesh);

    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    V_border = igl::is_border_vertex(V,F);
    igl::adjacency_list(F, VV);
    igl::vertex_triangle_adjacency(V,F,VF,VFi);
    igl::triangle_triangle_adjacency(F,TT,TTi);
    igl::edge_topology(V,F,E,F2E,E2F);

    // Generate "subdivided" mesh for visualization of curl terms
    igl::false_barycentric_subdivision(V, F, Vbs, Fbs);

    // Compute scale for visualizing fields
    global_scale =  0.5*igl::avg_edge_length(V, F);

    //Compute scale for visualizing texture
    uv_scale = 0.6/igl::avg_edge_length(V, F);

    // Compute face barycenters
    igl::barycenter(V, F, B);

    // Compute local basis for faces
    igl::local_basis(V,F,B1,B2,B3);
#endif
}

//////////////////////////////////////////////////////////////////////////
void JIntegrableField :: genRandomConstraints(int nRandom)
{
    /*
        b.resize(nRandom);

        int numFaces = mesh->getSize(2);

        for( int i = 0; i < nRandom; i++)
            b[i]    =   JMath::random_value(0, numFaces-1);

        bc.resize(nRandom,3*numVecPerFace);
        for(unsigned i=0; i<b.size(); ++i)
        {
            VectorXd t = random_constraints(B1.row(b(i)),B2.row(b(i)),numVecPerFace);
            bc.row(i) = t;
        }
    */
}

//////////////////////////////////////////////////////////////////////////

void JIntegrableField :: getCuts(const Eigen::MatrixXi &cuts)
{
    /*
      int maxCutNum = cuts.sum();
      Eigen::MatrixXd start(maxCutNum,3);
      Eigen::MatrixXd end(maxCutNum,3);
      int ind = 0;
      for (unsigned int i=0;i<F.rows();i++)
        for (int j=0;j<3;j++)
          if (cuts(i,j))
          {
            start.row(ind) = V.row(F(i,j));
            end.row(ind) = V.row(F(i,(j+1)%3));
            ind++;
          }
      viewer.data.add_edges(start, end , Eigen::RowVector3d(1.,0,1.));
    */
}

void JIntegrableField :: getField( const Eigen::MatrixXd &field,
                                   const Eigen::RowVector3d &color)
{
    /*
      for (int n=0; n<2; ++n)
      {
        Eigen::MatrixXd VF = field.block(0,n*3,F.rows(),3);
        Eigen::VectorXd c = VF.rowwise().norm();
        viewer.data.add_edges(B - global_scale*VF, B + global_scale*VF , color);
      }
    */
}

/*
void JIntegrableField :: getConstraints()
{
  for (int n=0; n<2; ++n)
  {
    Eigen::MatrixXd Bc = igl::slice(B, b, 1);
    Eigen::MatrixXd color;
    color.setZero(b.rows(),3);
    color.col(2).setOnes();
    for (int i =0; i<b.rows(); ++i)
      if (blevel[i] ==1 && n>0)
        color.row(i)<<0.7,0.7,0.7;
    // Eigen::RowVector3d color; color<<0.5,0.5,0.5;
    viewer.data.add_edges(Bc - global_scale*bc.block(0,n*3,bc.rows(),3), Bc + global_scale*bc.block(0,n*3,bc.rows(),3) , color);
  }
}
*/


void JIntegrableField :: colorEdgeMeshFaces(const Eigen::VectorXd &values,
        const double &minimum, const double &maximum, Eigen::MatrixXd &C)
{
    /*
      C.setConstant(Fbs.rows(),3,1);

      Eigen::MatrixXd colors;
      igl::jet(values, minimum, maximum, colors);

      for (int ei = 0; ei<E.rows(); ++ei)
      {
        const Eigen::RowVector3d &this_color = colors.row(ei);
        int f0 = E2F(ei,0);
        int f1 = E2F(ei,1);
        if(f0 != -1)
        {
          int i0 = -1;
          for (int k = 0; k<3; k++)
            if (F2E(f0,k)== ei)
            {
              i0 = k;
              break;
            }
          C.row(3*f0+i0) = this_color;
        }
        if(f1 != -1)
        {
          int i1 = -1;
          for (int k = 0; k<3; k++)
            if (F2E(f1,k)== ei)
            {
              i1 = k;
              break;
            }
          C.row(3*f1+i1) = this_color;
        }
      }
    */

}

/*
void update_display()
{

  if (display_mode == 1)
  {
    cerr<< "Displaying original field, its singularities and its cuts"  <<endl;

    viewer.data.set_mesh(V, F);

    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;
    viewer.data.set_colors(C);

    //Draw constraints
    drawConstraints(viewer);

    // Draw Field
    Eigen::RowVector3d color; color<<0,0,1;
    drawField(viewer,two_pv_ori,color);

    // Draw Cuts
    drawCuts(viewer,cuts_ori);

    //Draw Singularities
    Eigen::MatrixXd singular_points = igl::slice(V, singularities_ori, 1);
    viewer.data.add_points(singular_points,Eigen::RowVector3d(239./255.,205./255.,57./255.));

  }

  if (display_mode == 2)
  {
    cerr<< "Displaying current field, its singularities and its cuts"  <<endl;

    viewer.data.set_mesh(V, F);

    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;
    viewer.data.set_colors(C);

    //Draw constraints
    drawConstraints(viewer);

    // Draw Field
    Eigen::RowVector3d color; color<<0,0,1;
    drawField(viewer,two_pv,color);

    // Draw Cuts
    drawCuts(viewer,cuts);

    //Draw Singularities
    Eigen::MatrixXd singular_points = igl::slice(V, singularities, 1);
    viewer.data.add_points(singular_points,Eigen::RowVector3d(239./255.,205./255.,57./255.));
  }

  if (display_mode == 3)
  {
    cerr<< "Displaying original field and its curl"  <<endl;

    viewer.data.set_mesh(Vbs, Fbs);
    Eigen::MatrixXd C;
    colorEdgeMeshFaces(curl_ori, 0, 0.2, C);
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<1,1,1;
    drawField(viewer,two_pv_ori,color);

  }

  if (display_mode == 4)
  {
    cerr<< "Displaying current field and its curl"  <<endl;

    viewer.data.set_mesh(Vbs, Fbs);
    Eigen::MatrixXd C;
    colorEdgeMeshFaces(curl, 0, 0.2, C);
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<1,1,1;
    drawField(viewer,two_pv,color);
  }

  if (display_mode == 5)
  {
    cerr<< "Displaying original poisson-integrated field and original poisson error"  <<endl;

    viewer.data.set_mesh(V, F);
    Eigen::MatrixXd C;
    igl::jet(poisson_error_ori, 0, 0.5, C);
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<1,1,1;
    drawField(viewer,two_pv_poisson_ori,color);
  }

  if (display_mode == 6)
  {
    cerr<< "Displaying current poisson-integrated field and current poisson error"  <<endl;

    viewer.data.set_mesh(V, F);
    Eigen::MatrixXd C;
    igl::jet(poisson_error, 0, 0.5, C);
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<1,1,1;
    drawField(viewer,two_pv_poisson,color);
  }

  if (display_mode == 7)
  {
    cerr<< "Displaying original texture with cuts and singularities"  <<endl;

    viewer.data.set_mesh(V, F);
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    viewer.data.set_colors(C);
    viewer.data.set_uv(uv_scale*scalars_ori, Fcut_ori);
    viewer.data.set_texture(texture_R, texture_B, texture_G);
    viewer.core.show_texture = true;

    // Draw Cuts
    drawCuts(viewer,cuts_ori);

    //Draw Singularities
    Eigen::MatrixXd singular_points = igl::slice(V, singularities_ori, 1);
    viewer.data.add_points(singular_points,Eigen::RowVector3d(239./255.,205./255.,57./255.));

  }
  if (display_mode == 8)
  {
    cerr<< "Displaying current texture with cuts and singularities"  <<endl;

    viewer.data.set_mesh(V, F);
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    viewer.data.set_colors(C);
    viewer.data.set_uv(uv_scale*scalars, Fcut);
    viewer.data.set_texture(texture_R, texture_B, texture_G);
    viewer.core.show_texture = true;

    // Draw Cuts
    drawCuts(viewer,cuts);

    //Draw Singularities
    Eigen::MatrixXd singular_points = igl::slice(V, singularities, 1);
    viewer.data.add_points(singular_points,Eigen::RowVector3d(239./255.,205./255.,57./255.));

  }

  if (display_mode == 9)
  {
    cerr<< "Displaying original field overlayed onto the current integrated field"  <<endl;

    viewer.data.set_mesh(V, F);

    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<0,0,1;
    drawField(viewer,two_pv_ori,color);

    // Draw Integrated Field
    color<<.2,.2,.2;
    drawField(viewer,two_pv_poisson_ori,color);

  }

  if (display_mode == 0)
  {
    cerr<< "Displaying current field overlayed onto the current integrated field"  <<endl;

    viewer.data.set_mesh(V, F);

    // Highlight in red the constrained faces
    MatrixXd C = MatrixXd::Constant(F.rows(),3,1);
    for (unsigned i=0; i<b.size();++i)
      C.row(b(i)) << 1, 0, 0;
    viewer.data.set_colors(C);

    // Draw Field
    Eigen::RowVector3d color; color<<0,0,1;
    drawField(viewer,two_pv,color);

    // Draw Integrated Field
    color<<.2,.2,.2;
    drawField(viewer,two_pv_poisson,color);
  }
}

*/

int JIntegrableField :: smoothField()
{
#ifdef IGL_LIB
    //do a batch of iterations
    printf("--Improving Curl--\n");
    for (int bi = 0; bi<5; ++bi)
    {
        igl::integrable_polyvector_fields_solve(ipfdata, params, two_pv, iter ==0);
        iter++;
        params.wSmooth *= params.redFactor_wsmooth;
    }
    // Post process current field
    // Compute curl_minimizing matchings and curl
    printf("--Matchings and curl--\n");
    Eigen::MatrixXi match_ab, match_ba;  // matchings across interior edges
    double avgCurl = igl::polyvector_field_matchings(two_pv, V, F, true, match_ab, match_ba, curl);
    double maxCurl = curl.maxCoeff();
    printf("curl -- max: %.5g, avg: %.5g\n", maxCurl,  avgCurl);

    // Compute singularities
    printf("--Singularities--\n");
    igl::polyvector_field_singularities_from_matchings(V, F, match_ab, match_ba, singularities);
    printf("#singularities: %ld\n", singularities.rows());
    // Get mesh cuts based on singularities
    printf("--Cuts--\n");
    igl::polyvector_field_cut_mesh_with_singularities(V, F, singularities, cuts);
    // Comb field
    printf("--Combing--\n");
    Eigen::MatrixXd combed;
    igl::polyvector_field_comb_from_matchings_and_cuts(V, F, two_pv, match_ab, match_ba, cuts, combed);
    // Reconstruct integrable vector fields from combed field
    printf("--Cut mesh--\n");
    igl::cut_mesh(V, F, cuts, Vcut, Fcut);
    printf("--Poisson--\n");
    double avgPoisson = igl::polyvector_field_poisson_reconstruction(Vcut, Fcut, combed, scalars, two_pv_poisson, poisson_error);
    double maxPoisson = poisson_error.maxCoeff();
    printf("poisson error -- max: %.5g, avg: %.5g\n", maxPoisson, avgPoisson);
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////

int JIntegrableField :: genField()
{
#ifdef IGL_LIB
    // Interpolate a 2-PolyVector field to be used as the original field
    igl::n_polyvector(V, F, b, bc, two_pv_ori);

    // Post process original field
    // Compute curl_minimizing matchings and curl
    Eigen::MatrixXi match_ab, match_ba;  // matchings across interior edges
    printf("--Matchings and curl--\n");
    double avgCurl = igl::polyvector_field_matchings(two_pv_ori, V, F, true, match_ab, match_ba, curl_ori);
    double maxCurl = curl_ori.maxCoeff();
    printf("original curl -- max: %.5g, avg: %.5g\n", maxCurl,  avgCurl);

    printf("--Singularities--\n");
    // Compute singularities
    igl::polyvector_field_singularities_from_matchings(V, F, V_border, VF, TT, E2F, F2E, match_ab, match_ba, singularities_ori);
    printf("original #singularities: %ld\n", singularities.rows());

    printf("--Cuts--\n");
    // Get mesh cuts based on singularities
    igl::polyvector_field_cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities_ori, cuts_ori);

    printf("--Combing--\n");
    // Comb field
    Eigen::MatrixXd combed;
    igl::polyvector_field_comb_from_matchings_and_cuts(V, F, TT, E2F, F2E, two_pv_ori, match_ab, match_ba, cuts_ori, combed);

    printf("--Cut mesh--\n");
    // Reconstruct integrable vector fields from combed field
    igl::cut_mesh(V, F, VF, VFi, TT, TTi, V_border, cuts_ori, Vcut_ori, Fcut_ori);

    printf("--Poisson--\n");
    double avgPoisson = igl::polyvector_field_poisson_reconstruction(Vcut_ori, Fcut_ori, combed, scalars_ori, two_pv_poisson_ori, poisson_error_ori);
    double maxPoisson = poisson_error_ori.maxCoeff();
    printf("poisson error -- max: %.5g, avg: %.5g\n", maxPoisson, avgPoisson);

    // Set the curl-free 2-PolyVector to equal the original field
    two_pv = two_pv_ori;
    singularities = singularities_ori;
    curl = curl_ori;
    cuts = cuts_ori;
    two_pv_poisson = two_pv_poisson_ori;
    poisson_error = poisson_error_ori;
    Vcut = Vcut_ori;
    Fcut = Fcut_ori;
    scalars = scalars_ori;

    printf("--Integrable - Precomputation--\n");
    // Precompute stuff for solver
    igl::integrable_polyvector_fields_precompute(V, F, b, bc, blevel, two_pv_ori, ipfdata);
#endif

    return 0;
}
