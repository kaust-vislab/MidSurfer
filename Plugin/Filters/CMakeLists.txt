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

vtk_module_add_module(MidSurfer::vtkMidSurfer
  FORCE_STATIC
  CLASSES ${classes}
  SOURCES ${sources}
  HEADERS ${headers}
)

paraview_add_server_manager_xmls(
  XMLS  MidSurfer.xml)
