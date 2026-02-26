#pragma once

///////////////////////////////////////////////////////////////////////////////
//                      Quad-Cleanup
//
// Developed by:  Chaman Singh Verma
//                Department of Computer Sciences.
//                The University of Wisconsin, Madison
//
// Work Supported by:
//                 Dr. Tim Tautges
//                 Argonne National Lab, Chicago
//
//
// Objective:  Given a quadrilateral mesh, this class implements various strategies
// to improve the quadrilateral mesh both geometrically and topologically. The
// Laplacian ( local and global ) is used for geometric quality improvement, and for
// topological improvements various operations are used.
// The two basis operations for topological improvements are
//  1)   Face close
//  2)   doublet insertion and removal.
//

// Reference Papers:
//  1) Topological Improvement Procedures for Quadrilateral Finite Element Meshes
//     S.A. Canann,  S.N. Muthikrishnan and R.K. Phillips

//  2) Automated All Quadrilateral Mesh Adaptation Through Refinment and Coarsening
//     Bret Dallas Anderson
//     Master Thesis, Brigham Young University.
//
//  3) Non-Local Topological Clean-Up ( The idea of yring  is from this paper)
//     Guy Bunin.
//
// For suggestios, bugs and criticisms, please send e-mail to
//                      csverma@cs.wisc.edu
//
// Last Date update:  16th Feb 2010.
//
///////////////////////////////////////////////////////////////////////////////////

#include "Mesh.hpp"
#include "BinaryTreeMatch.hpp"
#include "basic_math.hpp"
#include "StopWatch.hpp"
#include "tfiblend.hpp"

#include "MeshGeodesics.hpp"
#include "QuadDefectPatch.hpp"
#include "MeshRefine.hpp"
#include "Doublet.hpp"
#include "Singlet.hpp"
#include "Diamond.hpp"

extern double area_of_poly3d(int n, double *x, double *y, double *z);

namespace Jaal {

struct IrregularNodeFilter : public JMeshFilter {
    bool pass( const JNodePtr &vertex ) const {
        if( vertex->isBoundary() ) return 1;
        if( vertex->getNumRelations(2) != 4 ) return 0;
        return 1;
    }
};

class JQuadCleanUp {
public:
    static bool isRegular( const JNodePtr &v);

    JQuadCleanUp() {
        initialize();
    }

    JQuadCleanUp(const JMeshPtr &m) {
        setMesh(m);
        initialize();
    }

    ~JQuadCleanUp() {
        finalize();
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    int getNumSingularNodes();

    // Insert new faces at the boundary edges. The boundary must be simple i.e.
    // when both the nodes of an edge are on the boundary, then the edge must
    // have only one neighbour face.
    int  insertPillows();

    vector<QDefectivePatch> search_defective_patches();
    QDefectivePatchPtr search_defective_patch();
    QDefectivePatchPtr build_defective_patch(const JNodePtr &vertex);

    // Global Cleanup methods ..
    int remesh_defective_patch( const QDefectivePatchPtr &p);
    int remesh_defective_patches();

    int collapseConcave( const JFacePtr &f);

    // Local Cleanup methods ..
    int  reduceDegree( const JNodePtr &v );
    int  vertex_degree_reduction();
    int  removeDiamonds();
    int  removeDoublets();
    int  removeSinglets();

    JNodePtr insertDoublet(const JFacePtr &face, const JNodePtr &v0, const JNodePtr &v2);
    JFaceSequence  search_flat_quads();

    int  automatic();
    void report();

#ifdef CSV
    /*
    // Query methods ...
    JNodeSequence  search_restricted_nodes();
    JFaceSequence  search_restricted_faces();
    JFaceSequence   search_diamonds(int type = 33 );
    JNodeSequence   search_singlets();
    JNodeSequence   search_doublets();
    vector<Edge>   search_tunnels();


    // int swap_concave_faces();
    // int  remove_bridges();
    // int  shift_irregular_nodes();
    // int  remove_tunnels();
    // void remove_ynodes();
    // int  clean_layer(int id);
    // void cleanup_boundary(double cutOffAngle = 100.0);
    // void advancing_front_cleanup();
    // void advancing_front_edges_swap();
    */


    // Some Feature that may be obsolete in the next version...
    JNodePtr insertDoublet(JFacePtr face);
    JNodePtr insert_boundary_doublet(JFacePtr face);

    int refine_restricted_node(Vertex *resnode, Vertex *bndnode);
    int refine_degree3_faces();
    int refine_bridges_face();
    void get_strips(Face *face, JFaceSequence &strip1, JFaceSequence strip2);
#endif

private:
    JMeshPtr mesh;
    JMeshFilterPtr filter;
    JMeshNonlinearOptimization mopt;

    JMeshGeodesics *djkpath; // Used in one defect remeshing ....
    QDefectivePatchPtr defective_patch;

    JNodeSequence  irregular_nodes;
    vector<JSinglet> vSinglets;
    vector<JDoublet> vDoublets;
    vector<JDiamond> vDiamonds; // Diamonds in the mesh;

    vector<QDefectivePatchPtr>  vDefectPatches;

    void initialize();
    void finalize();

//   int  has_interior_nodes_degree_345();
    vector<JDiamond> search_diamonds_in_layer(int l);
    int clean_layer_once(int id);
    int face_close(const JFacePtr &face, const JNodePtr &v0, const JNodePtr &v2);
    int diamond_collapse(JFaceClose &d);
    int remove_interior_doublet(JDoublet &d);

    int remove_bridges_in_layer( int l);
    int remove_bridges_once();
    void remove_diamonds_once();
    int remove_diamonds_in_layer( int l);
    int advance_front_edges_swap_once(int layerid);

    int apply_advance_front_bridge_rule( const JNodePtr &v0, const JNodePtr &v1);
    int apply_advance_front_excess_rule( const JNodePtr &v);
    int apply_advance_front_triplet_rule( const JNodePtr &v);
    int apply_advance_front_singlet_rule( const JNodePtr &v);

    int remove_doublets_once();
    int remove_interior_doublets_once();

    int boundary_vertex_degree_reduction_once();
    int internal_vertex_degree_reduction_once();

    int  reduce_internal_vertex_degree(const JNodePtr &v);
    int  reduce_boundary_vertex_degree(const JNodePtr &v);

    /*
    //   void cleanup_internal_boundary_face();
         // May become obsolere
         int refine_3434_pattern( Face *face, int pos);
         int refine_3454_pattern( Face *face, int pos);
         int refine_3444_pattern( Face *face, int pos);
         int apply_shift_node3_rule( Vertex *vertex);
    */
};

////////////////////////////////////////////////////////////////////////////////

inline bool
JQuadCleanUp::isRegular (const JNodePtr &v)
{
    assert(v);
    // Any interior vertex having four nodes( or faces ) is a regular node.
    if (!v->isBoundary () && (v->getNumRelations(2) == 4)) return 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Jaal


///////////////////////////////////////////////////////////////////////////////
