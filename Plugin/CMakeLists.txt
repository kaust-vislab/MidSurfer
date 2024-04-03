paraview_add_plugin(MidSurfer
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  VERSION "0.1"
  MODULES MidSurfer::vtkMidSurfer
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filters/vtk.module")

if (BUILD_TESTING)
  add_subdirectory(Testing)
endif ()
