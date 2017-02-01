#include "meshfiles.h"

#include "MeshImpl.hpp"
#include "MsqTimer.hpp"
#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "InstructionQueue.hpp"
#include "PatchData.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"

// algorythms
#include "Randomize.hpp"
#include "ConditionNumberQualityMetric.hpp"
#include "UntangleBetaQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "LInfTemplate.hpp"
#include "SteepestDescent.hpp"
#include "ConjugateGradient.hpp"
#include "PlanarDomain.hpp"

#include <iostream>
using std::cout;
using std::endl;
#include <cstdlib>

using namespace Mesquite;


int main()
{
    Mesquite::MeshImpl mesh;
    MsqPrintError err(cout);
    mesh.read_vtk(MESH_FILES_DIR "2D/VTK/tangled_quad.vtk", err);
    if (err) return 1;

    // Set Domain Constraint
    Vector3D pnt(0,0,0);
    Vector3D s_norm(0,0,1);
    PlanarDomain msq_geom(s_norm, pnt);

    // creates an intruction queue
    InstructionQueue queue1;

    // creates a mean ratio quality metric ...
    ConditionNumberQualityMetric shape_metric;
    UntangleBetaQualityMetric untangle(2);
    Randomize pass0(.05);
    // ... and builds an objective function with it
    //LInfTemplate* obj_func = new LInfTemplate(shape_metric);
    LInfTemplate obj_func(&untangle);
    LPtoPTemplate obj_func2(&shape_metric, 2, err);
    if (err) return 1;
    // creates the steepest descent optimization procedures
    ConjugateGradient pass1( &obj_func, err );
    if (err) return 1;

    //SteepestDescent* pass2 = new SteepestDescent( obj_func2 );
    ConjugateGradient pass2( &obj_func2, err );
    if (err) return 1;
    pass2.use_element_on_vertex_patch();
    if (err) return 1;
    pass2.use_global_patch();
    if (err) return 1;
    QualityAssessor stop_qa=QualityAssessor(&shape_metric);
    QualityAssessor stop_qa2=QualityAssessor(&shape_metric);

    stop_qa.add_quality_assessment(&untangle);
    // **************Set stopping criterion**************
    //untangle beta should be 0 when untangled
    TerminationCriterion sc1;
    sc1.add_relative_quality_improvement( 0.000001 );
    TerminationCriterion sc3;
    sc3.add_iteration_limit( 10 );
    TerminationCriterion sc_rand;
    sc_rand.add_iteration_limit( 1 );

    //StoppingCriterion sc1(&stop_qa,-1.0,.0000001);
    //StoppingCriterion sc3(&stop_qa2,.9,1.00000001);
    //StoppingCriterion sc2(StoppingCriterion::NUMBER_OF_PASSES,10);
    //StoppingCriterion sc_rand(StoppingCriterion::NUMBER_OF_PASSES,1);
    //either until untangled or 10 iterations
    pass0.set_outer_termination_criterion(&sc_rand);
    pass1.set_outer_termination_criterion(&sc1);
    pass2.set_inner_termination_criterion(&sc3);

    // adds 1 pass of pass1 to mesh_set1
    queue1.add_quality_assessor(&stop_qa,err);
    if (err) return 1;
    //queue1.add_preconditioner(pass0,err);MSQ_CHKERR(err);
    //queue1.add_preconditioner(pass1,err);MSQ_CHKERR(err);
    //queue1.set_master_quality_improver(pass2, err); MSQ_CHKERR(err);
    queue1.set_master_quality_improver(&pass1, err);
    if (err) return 1;
    queue1.add_quality_assessor(&stop_qa2,err);
    if (err) return 1;
    mesh.write_vtk("original_mesh.vtk", err);
    if (err) return 1;

    // launches optimization on mesh_set1
    queue1.run_instructions(&mesh, &msq_geom, err);
    if (err) return 1;

    mesh.write_vtk("smoothed_mesh.vtk", err);
    if (err) return 1;

    print_timing_diagnostics(cout);
    return 0;
}
