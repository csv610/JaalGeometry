#pragma once

#include "Ui_MeshToolsDialog.hpp"

#include <QDialog>
#include <QColorDialog>
#include <QMessageBox>

#include "MeshAffineTransformsDialog.hpp"
#include "FontsDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshNodesDialog.hpp"
#include "MeshEdgesDialog.hpp"
#include "MeshFacesDialog.hpp"
#include "MeshCellsDialog.hpp"
#include "ObjectsListDialog.hpp"
#include "MeshEntityPickerDialog.hpp"
#include "MeshDualGrapherDialog.hpp"
#include "MeshSlicerDialog.hpp"
#include "ContourEditingDialog.hpp"

#include "MeshRelationsTableDialog.hpp"
#include "MeshTopologyQueryDialog.hpp"
#include "LightsDialog.hpp"
#include "SuggestiveContoursDialog.hpp"
#include "MagnifyLensDialog.hpp"

#include "TriMesherDialog.hpp"
#include "QuadMesherDialog.hpp"
#include "TetMesherDialog.hpp"
#include "HexMesherDialog.hpp"
#include "GenSimpleShapeDialog.hpp"
#include "MeshBooleanDialog.hpp"
#include "GlobalSettingsDialog.hpp"
#include "MeshPartitionDialog.hpp"
#include "MeshSegmentationDialog.hpp"
#include "MeshGeometricQualityDialog.hpp"
#include "MeshTopologyQualityDialog.hpp"
#include "MeshOptDialog.hpp"
#include "MeshContoursDialog.hpp"
#include "StructuredMeshDialog.hpp"
#include "TriMeshCleanupDialog.hpp"
#include "TriSimplificationDialog.hpp"
#include "MeshSpacePartitionsDialog.hpp"
#include "MeshDeformationDialog.hpp"
#include "MeshDualsDialog.hpp"
#include "MeshSkeletonDialog.hpp"
#include "MeshFeaturesDialog.hpp"
#include "MeshUtilsDialog.hpp"
#include "QuadCleanUp.hpp"
#include "SurfaceParameterizationDialog.hpp"
#include "MeshSurfaceVectorFieldDialog.hpp"
#include "MeshInterpolationDialog.hpp"
#include "MeshSpectrumDialog.hpp"
#include "MeshRenderDialog.hpp"
#include "MeshRefine2DDialog.hpp"
#include "MeshRefine3DDialog.hpp"
#include "PatchQuadmeshingDialog.hpp"
#include "MorseAnalysisDialog.hpp"
#include "MeshMeanCurvatureFlowDialog.hpp"

/*
#include "CongruentTessellationDialog.hpp"
#include "InteractiveMeshingDialog.hpp"
#include "HeatConductionDialog.hpp"
#include "MeshStackDialog.hpp"
#include "MeshTangleDialog.hpp"
#include "ObjectsListDialog.hpp"
#include "SaveAnimationDialog.hpp"
#include "TangleFEMTestsDialog.hpp"
*/

#include <igl/embree/ambient_occlusion.h>
#include <igl/per_vertex_normals.h>

class JMeshToolsDialog : public QDialog, public Ui::MeshToolsDialog {
    Q_OBJECT

public:
    JMeshToolsDialog( QWidget *parent = 0);
    ~JMeshToolsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

protected:
    void  showEvent(QShowEvent *e);
    void  keyPressEvent( QKeyEvent *e);

private slots:
    void  openAffineDialog();
    void  openRelationsDialog();
    void  openTopoQueryDialog();
    void  openPickEntityDialog();
    void  openMagnifyingLensDialog();
    void  openMeshDualGrapherDialog();
    void  openMeshSlicerDialog();

    void  ambientOcclusion();
    void  openSuggestiveContoursDialog();

    void openTriMesherDialog();
    void openQuadMesherDialog();
    void openTetMesherDialog();
    void openHexMesherDialog();
    void openMeshSkeletonDialog();
    void openMeshRenderDialog();
    void openMeshOptDialog();
    void openGenSimpleShapeDialog();
    void openMeshBooleanDialog();
    void openMeshDualDialog();
    void openMeshFeaturesDialog();
    void openMeshGeomQualityDialog();
    void openMeshTopoQualityDialog();
    void openMeshSegmentationDialog();
    void openMeshRefine2DDialog();
    void openMeshRefine3DDialog();
    void openMeshSpectrumDialog();
    void openSurfVecFieldDialog();
    void openMeshDeformationDialog();
    void openMeshInterpolationDialog();
    void openMeshPartitionDialog();
    void openMeshGeodesicsDialog();
    void openSurfaceParameterizationDialog();
    void openMeshMeanCurvatureFlowDialog();
    void openContourEditingDialog();

/*
    void openCongruentTessellationDialog();
    void openHeatConductionDialog();
    void openInteractiveMeshingDialog();

    void openMeshTangleDialog();
    void openMeshUtilsDialog();
    void openMorseAnalysisDialog();
    void openObjectsListDialog();
    void openPatchRemeshDialog();
    void openQuadPartitionDialog();
    void openSaveAnimationDialog();
    void openTrimeshCleanupDialog();
    void openTriSimplifyMeshDialog();
    void openTangleFEMTestsDialog();
    void openMeshContoursDialog();
    void sphericalMap();
*/

    void  setMesh();
    void  selectMesh(QModelIndex index);
    void  resetCamera();
    void  closeDialog();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;

    boost::scoped_ptr<JObjectsListDialog> meshlistDialog;
    boost::scoped_ptr<JMeshEntityPickerDialog> pickEntityDialog;
    boost::scoped_ptr<JMeshAffineTransformsDialog>    affineDialog;
    boost::scoped_ptr<JMeshRelationsTableDialog> relationsDialog;
    boost::scoped_ptr<JMeshTopologyQueryDialog>  topoQueryDialog;
    boost::scoped_ptr<JSuggestiveContoursDialog>  suggestiveContoursDialog;
    boost::scoped_ptr<JMeshDualGrapherDialog>  dualMeshDialog;
    boost::scoped_ptr<JMeshSlicerDialog>  meshSlicerDialog;
    boost::scoped_ptr<JMagnifyingLensDialog> magnifyingLensDialog;

    boost::scoped_ptr<JTriMesherDialog>   triMesherDialog;
    boost::scoped_ptr<JQuadMesherDialog>  quadMesherDialog;
    boost::scoped_ptr<JTetMesherDialog>   tetMesherDialog;
    boost::scoped_ptr<JHexMesherDialog>   hexMesherDialog;
    boost::scoped_ptr<JMeshSkeletonDialog>  meshSkeletonDialog;
    boost::scoped_ptr<JMeshRenderDialog> meshRenderDialog;
    boost::scoped_ptr<JMeshOptDialog>  meshOptDialog;
    boost::scoped_ptr<JGenSimpleShapeDialog> genSimpleShapeDialog;
    boost::scoped_ptr<JMeshRefine2DDialog>     refine2dDialog;
    boost::scoped_ptr<JMeshRefine3DDialog>     refine3dDialog;
    boost::scoped_ptr<JMeshTopologyQualityDialog>  meshTopoQualityDialog;
    boost::scoped_ptr<JMeshSegmentationDialog> segmentationDialog;
    boost::scoped_ptr<JMeshBooleanDialog>     meshBooleanDialog;
    boost::scoped_ptr<JMeshDualsDialog>        meshdualsDialog;
    boost::scoped_ptr<JMeshFeaturesDialog>     meshFeaturesDialog;
    boost::scoped_ptr<JMeshPartitionDialog>    meshpartitionDialog;
    boost::scoped_ptr<JMeshSurfaceVectorFieldDialog> surfVecFieldDialog;
    boost::scoped_ptr<JMeshSpectrumDialog>     meshSpectrumDialog;
    boost::scoped_ptr<JMeshGeometricQualityDialog> meshGeomQualityDialog;
    boost::scoped_ptr<JMeshInterpolationDialog> meshInterpolationDialog;
    boost::scoped_ptr<JMeshDeformationDialog>  meshDeformationDialog;
    boost::scoped_ptr<JMeshPartitionDialog>  meshPartitionDialog;
    boost::scoped_ptr<JMeshGeodesicsDialog>    meshGeodesicsDialog;
    boost::scoped_ptr<JSurfaceParameterizationDialog> surfaceParameterizationDialog;
    boost::scoped_ptr<JMeshMeanCurvatureFlowDialog> meanCurvatureFlowDialog;
    boost::scoped_ptr<JContourEditingDialog> contourEditingDialog;
/*

    boost::scoped_ptr<JCongruentTessellationDialog> congruentTessellationDialog;
    boost::shared_ptr<JMeshDualViewer>  dualViewer;
    boost::scoped_ptr<JMeshTangleDialog>       meshTangleDialog;
    boost::scoped_ptr<JMeshUtilsDialog>        meshUtilsDialog;
    boost::scoped_ptr<JMorseAnalysisDialog>    morseAnalysisDialog;
    boost::scoped_ptr<JTrimeshCleanupDialog>   trimeshCleanupDialog;
    boost::scoped_ptr<JObjectsListDialog>      objectsListDialog;
    boost::scoped_ptr<JSaveAnimationDialog>    saveAnimationDialog;
    boost::scoped_ptr<JTriSimplificationDialog>  triSimplifyDialog;
    boost::scoped_ptr<JHeatConductionDialog>  heatConductionDialog;
    boost::scoped_ptr<JTangleFEMTestsDialog>     tanglefemtestsDialog;
    boost::scoped_ptr<JInteractiveMeshingDialog>  interactiveMeshingDialog;
    boost::scoped_ptr<JMeshContoursDialog>     meshContoursDialog;
*/


    void init();
    void makeConnections();
    void setListView();
    void warnMessage();
};

