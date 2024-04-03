set(classes
  vtkComputeEigenVectorField
  vtkExtractCenterLine
  vtkExtractCenterLineFromRegion
  vtkExtractMidsurface
  vtkGeneratePointCloudFromSegmentationMask
  vtkSignedDistanceField
  vtkConnectedCommponentsBinaryImage
  vtkZipperTriangulation
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
  XMLS  MidSurfer.xml)
