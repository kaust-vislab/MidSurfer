paraview_add_plugin(MidsurfaceExtractor
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  VERSION "0.1"
  MODULES MidsurfaceExtractor::vtkMidsurfaceExtractor
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filters/vtk.module")

if (BUILD_TESTING)
  add_subdirectory(Testing)
endif ()
