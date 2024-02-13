#message(STATUS "Fetch content for svOneDSolver")

#find_package(OneDSolver_lib HINTS "/home/a11bmafr/software/SimVascular/svOneDSolver/cmake-build-debug" REQUIRED)

#find_package(OneDSolver_lib HINTS

#message(${BACI_OneDSolver_ROOT})

#if(NOT OneDSolver_lib_found)
#    message(
#            FATAL_ERROR
#            "OneDSolverv2 could not be found. "
#    )
#endif()
#add_subdirectory("/home/a11bmafr/software/SimVascular/install/")
add_library(svOneDSolver UNKNOWN IMPORTED)

# Set the path to your library and include directories
set(OneDSolver_INCLUDE_DIR "/home/a11bmafr/software/SimVascular/install/include/onedsolver")
set(OneDSolver_LIBRARY
    "/home/a11bmafr/software/SimVascular/install/usr/local/lib/libOneDSolver_core.a"
    )
set_target_properties(
  svOneDSolver
  PROPERTIES IMPORTED_LOCATION "${OneDSolver_LIBRARY}"
             INTERFACE_INCLUDE_DIRECTORIES "${OneDSolver_INCLUDE_DIR}"
  )

#if(NOT TARGET OneDSolver_lib)
#fetchcontent_declare(
#        OneDSolver
#        GIT_REPOSITORY git@github.com:knarfnitram/svOneDSolver.git
#        GIT_TAG d30d825f5259a72782de9fd0318cef65fa8a2955
#)
#fetchcontent_makeavailable(OneDSolver)
#baci_add_dependency(baci_all_enabled_external_dependencies OneDSolver)
#endif ()
#find_package(Qhull REQUIRED)
#baci_add_dependency(baci_all_enabled_external_dependencies OneDSolver_lib)

message(STATUS "onedsolver include directory: ${OneDSolver_INCLUDE_DIR}")
message(STATUS "onedsolver library directory: ${OneDSolver_LIBRARY}")
if(svOneDSolver)
  message(STATUS "onedsolver found")
endif()
target_link_libraries(baci_all_enabled_external_dependencies INTERFACE svOneDSolver)
#baci_add_dependency(baci_all_enabled_external_dependencies onedsolver)
