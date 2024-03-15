set(classes
  vtkComputeEigenVectorField
  vtkExtractCenterLine
  vtkExtractCenterLineFromRegion
  vtkExtractMidsurface
  vtkExtractMidsurface3D
  vtkGeneratePointCloudFromSegmentationMask
  vtkSignedDistanceField
  vtkConnectedCommponentsBinaryImage
  vtkConnectedCommponentsBinaryImage3D
)

set(sources
  #signed_distance.cxx
)

vtk_module_add_module(MidsurfaceExtractor::vtkMidsurfaceExtractor
  FORCE_STATIC
  CLASSES ${classes}
  SOURCES ${sources}
  HEADERS ${headers}
)

paraview_add_server_manager_xmls(
  XMLS  MidsurfaceExtractor.xml)
