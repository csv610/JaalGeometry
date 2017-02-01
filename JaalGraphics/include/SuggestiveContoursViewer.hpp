#pragma once

#include <QDialog>
#include "MeshViewer.hpp"

#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "apparentridge.h"
#include <vector>

class JSuggestiveContoursViewer : public JViewComponent {
public:
    explicit JSuggestiveContoursViewer( JaalViewer *p);

    ~JSuggestiveContoursViewer();

    void draw();

    void loadNewMesh( const string &s);

    void smoothNormals()   {
        filter_normals();
    }
    void smoothCurvature() {
        filter_curv();
    }
    void smoothCurvatureDeriv() {
        filter_dcurv();
    }
    void subdivisionSurface() {
        subdivide_mesh();
    }

// Set to false for hardware that has problems with display lists
    const bool use_dlists = true;
// Set to false for hardware that has problems with supplying 3D texture coords
    const bool use_3dtexc = false;

    JaalViewer *viewManager;

// Globals: mesh...
    trimesh::TriMesh *themesh;
    trimesh:: point viewpos;    // Current view position

// Two cameras: the primary one, and an alternate one to fix the lines
// and see them from a different direction
    int dual_vpmode = false, mouse_moves_alt = false;
//  xform xf, xf_alt;

// Toggles for drawing various lines
    int draw_extsil = 0, draw_c = 1, draw_sc = 1;
    int draw_sh = 0, draw_phridges = 0, draw_phvalleys = 0;
    int draw_ridges = 0, draw_valleys = 0, draw_apparent = 0;
    int draw_K = 0, draw_H = 0, draw_DwKr = 0;
    int draw_bdy = 0, draw_isoph = 0, draw_topo = 0;
    int niso = 20, ntopo = 20;
    float topo_offset = 0.0f;

// Toggles for tests we perform
    int draw_hidden = 0;
    int test_c = 1, test_sc = 1, test_sh = 1, test_ph = 1, test_rv = 1, test_ar = 1;
    float sug_thresh = 0.01, sh_thresh = 0.02, ph_thresh = 0.04;
    float rv_thresh = 0.1, ar_thresh = 0.1;

// Toggles for style
    int use_texture = 0;
    int draw_faded = 1;
    int draw_colors = 0;
    int use_hermite = 0;

// Mesh colorization
    enum { COLOR_WHITE, COLOR_GRAY, COLOR_CURV, COLOR_GCURV, COLOR_MESH };
    const int ncolor_styles = 5;
    int color_style = COLOR_WHITE;
    vector<trimesh::Color> curv_colors, gcurv_colors;
    int draw_edges = false;

// Lighting
    enum { LIGHTING_NONE, LIGHTING_LAMBERTIAN, LIGHTING_LAMBERTIAN2,
           LIGHTING_HEMISPHERE, LIGHTING_TOON, LIGHTING_TOONBW, LIGHTING_GOOCH
         };
    const int nlighting_styles = 7;
    GLuint texture_contexts[7];
    int lighting_style = LIGHTING_NONE;
    float lightdir_matrix[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    int light_wrt_camera = true;

// Per-vertex vectors
    int draw_norm = 0, draw_curv1 = 0, draw_curv2 = 0, draw_asymp = 0;
    int draw_w = 0, draw_wperp = 0;

// Other miscellaneous variables
    float feature_size;     // Used to make thresholds dimensionless
    float currsmooth;       // Used in smoothing
    trimesh::vec currcolor;          // Current line color

    void init();
    void draw_tstrips();
    void make_texture(float width);

// Draw contours and suggestive contours using texture mapping
    void draw_c_sc_texture(const vector<float> &ndotv,
                           const vector<float> &kr,
                           const vector<float> &sctest_num,
                           const vector<float> &sctest_den);

// Color the mesh by curvatures
    void compute_curv_colors();

// Similar, but grayscale mapping of mean curvature H
    void compute_gcurv_colors();

// Set up textures to be used for the lighting.
// These are indexed by (n dot l), though they are actually 2D textures
// with a height of 1 because some hardware (cough, cough, ATI) is
// thoroughly broken for 1D textures...
    void make_light_textures(GLuint *texture_contexts);


// Draw the basic mesh, which we'll overlay with lines
    void draw_base_mesh();


// Compute per-vertex n dot l, n dot v, radial curvature, and
// derivative of curvature for the current view
    void compute_perview(vector<float> &ndotv, vector<float> &kr,
                         vector<float> &sctest_num, vector<float> &sctest_den,
                         vector<float> &shtest_num, vector<float> &q1,
                         vector<vec2> &t1, vector<float> &Dt1q1,
                         bool extra_sin2theta = false);

// Compute gradient of (kr * sin^2 theta) at vertex i
    vec gradkr(int i);

    float find_zero_linear(float val0, float val1);

    float find_zero_hermite(int v0, int v1, float val0, float val1,
                            const vec &grad0, const vec &grad1);

    void draw_face_isoline2(int v0, int v1, int v2,
                            const vector<float> &val,
                            const vector<float> &test_num,
                            const vector<float> &test_den,
                            bool do_hermite, bool do_test, float fade);


    void draw_face_isoline(int v0, int v1, int v2,
                           const vector<float> &val,
                           const vector<float> &test_num,
                           const vector<float> &test_den,
                           const vector<float> &ndotv,
                           bool do_bfcull, bool do_hermite,
                           bool do_test, float fade);


    void draw_isolines(const vector<float> &val,
                       const vector<float> &test_num,
                       const vector<float> &test_den,
                       const vector<float> &ndotv,
                       bool do_bfcull, bool do_hermite,
                       bool do_test, float fade);

    void draw_segment_ridge(int v0, int v1, int v2,
                            float emax0, float emax1, float emax2,
                            float kmax0, float kmax1, float kmax2,
                            float thresh, bool to_center);

    void draw_face_ridges(int v0, int v1, int v2,
                          bool do_ridge,
                          const vector<float> &ndotv,
                          bool do_bfcull, bool do_test, float thresh);

    void draw_face_ph(int v0, int v1, int v2, bool do_ridge,
                      const vector<float> &ndotv, bool do_bfcull,
                      bool do_test, float thresh);

    void draw_mesh_ridges(bool do_ridge, const vector<float> &ndotv,
                          bool do_bfcull, bool do_test, float thresh);

    void draw_mesh_ph(bool do_ridge, const vector<float> &ndotv, bool do_bfcull,
                      bool do_test, float thresh);

    void draw_silhouette(const vector<float> &ndotv);

    void draw_boundaries(bool do_hidden);

    void draw_topolines(const vector<float> &ndotv);

    void draw_misc(const vector<float> &ndotv, const vector<float> &DwKr,
                   bool do_hidden);
    void draw_mesh();

    void cls();
    void draw_isophotes(const vector<float> &ndotv);
    void filter_mesh(int dummy = 0);

    void filter_normals(int dummy = 0);
    void filter_curv(int dummy = 0);
    void filter_dcurv(int dummy = 0);
    void subdivide_mesh(int dummy = 0);

    void compute_feature_size();

    void largest_eig_2x2(float m1, float m12, float m2, vec2 &e1, float &l1);

    void draw_segment_app_ridge(int v0, int v1, int v2,
                                float emax0, float emax1, float emax2,
                                float kmax0, float kmax1, float kmax2,
                                const vec &tmax0, const vec &tmax1, const vec &tmax2,
                                float thresh, bool to_center, bool do_test);

// Compute principal view dependent curvatures and directions at vertex i.
    void compute_viewdep_curv(const TriMesh *mesh, int i, float ndotv,
                              float u2, float uv, float v2,
                              float &q1, vec2 &t1);

// Compute D_{t_1} q_1 - the derivative of max view-dependent curvature
// in the principal max view-dependent curvature direction.
    void compute_Dt1q1(const TriMesh *mesh, int i, float ndotv,
                       const vector<float> &q1, const vector<vec2> &t1,
                       float &Dt1q1);

// Draw apparent ridges of the mesh
    void draw_mesh_app_ridges(const vector<float> &ndotv, const vector<float> &q1,
                              const vector<vec2> &t1, const vector<float> &Dt1q1,
                              bool do_bfcull, bool do_test, float thresh);
    void draw_face_app_ridges(int v0, int v1, int v2,
                              const vector<float> &ndotv, const vector<float> &q1,
                              const vector<vec2> &t1, const vector<float> &Dt1q1,
                              bool do_bfcull, bool do_test, float thresh);

};
